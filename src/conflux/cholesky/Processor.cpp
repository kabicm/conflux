/**
 * Copyright (C) 2020, ETH Zurich
 *
 * This product includes software developed at the Scalable Parallel Computing
 * Lab (SPCL) at ETH Zurich, headed by Prof. Torsten Hoefler. It was developped
 * as a part of the "Design of Parallel- and High-Performance Computing"
 * lecture. For more information, visit http://spcl.inf.ethz.ch/. Unless you 
 * have an agreement with ETH Zurich for a separate license, the following
 * terms and conditions apply:
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along 
 * with this program. If not, see http://www.gnu.org/licenses.

 * @file Processor.cpp
 * 
 * @brief implementation of the Processor class
 * 
 * @authors Jens Eirik Saethre, Andre Gaillard, Pascal Müller, Marc Styger
 * Contact: {saethrej, andrega, stygerma, pamuelle}@student.ethz.ch
 * 
 * @date 23.11.2020
 */

#include <sstream>
#include <mpi.h>
#include <math.h>

#include "Processor.h"

/**
 * @brief creates a Processor object and allocates all the buffers
 * 
 * This constructor creates all data structures that a processor needs through-
 * out the entirety of the algorithm. This includes basic information like
 * rank, maximum indices, but also buffers. The latter are all allocated on the
 * heap during the call of the constructor.
 * 
 * @pre The MPI environment was initialized. Otherwise, an exception is thrown.
 * @post all processor properties are computed and all buffers are allocated
 * 
 * @param prop pointer to properties of the Cholesky algorithm
 */
conflux::Processor::Processor(CholeskyProperties *prop)
{
    // check if MPI environment was already initialized, throw exception if not
    int init;
    MPI_Initialized(&init);
    if (!init) {
        throw CholeskyException(CholeskyException::errorCode::FailedMPIInit);
    }

    // get rank and grid position
    int procRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    this->rank = static_cast<ProcRank>(procRank);
    this->grid = prop->globalToGrid(this->rank);
    this->px = this->grid.px;
    this->py = this->grid.py;
    this->pz = this->grid.pz;

    // compute maximal local tile indices for this processor
    this->maxIndexA10 = (prop->Kappa - 1) / prop->P
                        + ((prop->Kappa - 1) % prop->P > this->rank ? 1 : 0);
    this->maxIndexA11i = (prop->Kappa - 1) / prop->PX
                         + ((prop->Kappa - 1) % prop->PX > this->px ? 1 : 0);
    this->maxIndexA11j = (prop->Kappa - 1) / prop->PY 
                         + ((prop->Kappa - 1) % prop->PY > this->py ? 1 : 0);

    // allocate permanent buffers for this processor
    this->A00 = new double[prop->v * prop->v];
    this->A10 = new TileMatrix(MatrixType::VECTOR, prop->v, prop->v, this->maxIndexA10);
    this->A11 = new TileMatrix(MatrixType::MATRIX, prop->v, prop->v, this->maxIndexA11i,
                               this->maxIndexA11j);

    // allocate temporary receive buffers
    this->A10rcv = new TileMatrix(MatrixType::VECTOR, prop->v, prop->l, this->maxIndexA11i);
    this->A01rcv = new TileMatrix(MatrixType::VECTOR, prop->v, prop->l, this->maxIndexA11j);

    // set the request counters to zero
    this->cntUpdateA10 = 0;
    this->cntScatterA11 = 0;

    // compute the upper bounds for the number of requests on sub-tiles
    this->sndBound = this->maxIndexA10 * (prop->PZ * (prop->PX + prop->PY));
    this->rcvBound = this->maxIndexA11i + this->maxIndexA11j;

    // reserve memory for the MPI request vectors. Note that all these sizes are
    // upper bounds, and not exact by any means. It is only important to 
    reqUpdateA10.resize(this->rcvBound);
    reqScatterA11.resize(this->maxIndexA10 + 1);
    
    // create a new communicator for all processors along the same z-axis as
    // the current processor, i.e. processors that share (px,py) coordinates
    // we define color = px * PY + py, i.e. rank on XY-plane in row-major order
    // and rank as the pz coordinate.
    MPI_Comm_split(MPI_COMM_WORLD, this->px * prop->PY + this->py, this->pz, &this->zAxisComm);

    // set the private pointer to the properties object
    m_prop = prop;
}

/**
 * @brief destroys a processor object and frees all buffers
 */
conflux::Processor::~Processor()
{
    // delete permanent buffers
    delete[] A00;
    delete A10;
    delete A11;

    // delete temporary receive buffers
    delete A10rcv;
    delete A01rcv;

    // the processor with rank 0 within its local z-axis communicator (i.e. 
    // the one with pz=0) has to free the communicator at the end of the
    // execution.
    if (this->pz == 0) {
        MPI_Comm_free(&this->zAxisComm);
    }
}

/**
 * @brief updates the broadcast communicator
 * @param remTiles the remaining number of tiles in A10
 */
void conflux::Processor::updateBcastComm(uint32_t remTiles)
{
    // TODO: implement
}

/**
 * @brief initializes the broadcast communicators
 */
void conflux::Processor::initializeBroadcastComms()
{
    // if we have 8 or less processors, it's not worth it to generate new communicators
    if (m_prop->P <= 8) {
        // inser the global communicator into the lists
        m_bcastComms.push_back(MPI_COMM_WORLD);
        m_bcastSizes.push_back(m_prop->P);
        m_membershipFlags.push_back(true);
        m_curIdx = 0;
        std::set<ProcRank> tO;
        for (size_t i = 0; i < m_prop->P; ++i) {
            tO.insert((ProcRank) i);
        }
        m_tileOwners.push_back(tO);

        // update the public variables
        bcastComm = m_bcastComms[0];
        inBcastComm = true;

        return;
    }

    // there are PX * PZ processors potentially owning tiles that become A00
    uint32_t numBcastProcessors = m_prop->PX * m_prop->PZ;
    uint32_t maxBcastSize = std::min(m_prop->P, m_prop->Kappa); // this is a power of two in our case
    uint8_t numComms = 0;

    int idx = 0;
    // if the maxBcast size is P (i.e. num processors) then use MPI_COMM_WORLD
    if (maxBcastSize == m_prop->P) {
        
    }













}

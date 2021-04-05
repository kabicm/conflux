/**
 * Copyright (C) 2020, ETH Zurich
 *
 * This product includes software developed at the Scalable Parallel Computing
 * Lab (SPCL) at ETH Zurich, headed by Prof. Torsten Hoefler. For more information
 * visit http://spcl.inf.ethz.ch/. Unless you have an agreement with ETH Zurich
 * for a separate license, the following terms and conditions apply:
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

 * @file TileMatrix.h
 * 
 * @brief declaration & implementation of class for matrices of tiles, where there is no need
 * for accessing individual elements of a single tile, but only the tile as
 * a whole, i.e. the pointer to it. More over, tiles are stored contiguously.
 * 
 * @authors Jens Eirik Saethre, André Gaillard
 * Contact: {saethrej, andrega, stygerma, pamuelle}@student.ethz.ch
 * 
 * @date 19.12.2020
 */

#ifndef TILEMATRIX_H
#define TILEMATRIX_H

#include <stdint.h>
#include <exception>

namespace conflux {

/**
 * @brief defines a type for a TileMatrix object, if it's a tile-column vector
 * only (e.g for A10snd and A01snd), which is 1-dimensional, or a tile matrix 
 * (e.g. for A11), which is 2-dimensional.
 */
enum class MatrixType {VECTOR, MATRIX};

/**
 * @brief class that represents a matrix of tiles. This is a contiguous buffer
 * that consists of multiple tiles that are stored either in a 1d vector or a
 * 2d matrix. Each tile can be accessed individually via the get() function.
 * 
 * @tparam T the template parameter used for this class. Will most likely be double
 * for the matrix elements and TileReady for the dgemm interleaving data structure.
 */
template <class T>
class TileMatrix
{
public:
    // constructors and destructors
    TileMatrix(MatrixType type, uint32_t rowsPerTile, uint32_t colsPerTile,
               uint32_t tileRows, uint32_t tileCols = 1);
    ~TileMatrix();

    // methods to access a tile in the tile matrix
    T* get(uint32_t i);
    T* get(uint32_t i, uint32_t j);

private:
    MatrixType type;        //!< type of the matrix
    T *data;                //!< pointer to the start of the tile matrix
    uint32_t rowsPerTile;   //!< number of rows per tile
    uint32_t colsPerTile;   //!< number of cols per tile
    uint32_t tileRows;      //!< number of rows in this tile matrix
    uint32_t tileCols;      //!< number of cols in this tile matrix
}; 

} // namespace conflux

/******************************* IMPLEMENTATION *******************************
 *                 (since this is a class with template params) 
 */

/**
 * @brief constructor of a tile matrix
 * 
 * constructs a new tile matrix, which is either a tile-column vector in 1
 * dimension or a tile matrix in 2 dimensions. The tiles will be stored 
 * contiguously in memory, which makes it useful to perform larger GEMM 
 * operations with LAPACK, which should lead to performance increases.
 * 
 * @param type the type of the tile matrix, either VECTOR or MATRIX
 * @param rowsPerTile the number of rows per tile
 * @param colsPerTile the number of cols per tile
 * @param tileRows the number of tile rows in this matrix
 * @param tileCols the number of tile cols in this matrix, defaults to 1
 */
template <class T>
conflux::TileMatrix<T>::TileMatrix(MatrixType type, uint32_t rowsPerTile, uint32_t colsPerTile,
                       uint32_t tileRows, uint32_t tileCols)
{
    // asign parameters to local variables
    this->type = type;
    this->rowsPerTile = rowsPerTile;
    this->colsPerTile = colsPerTile;
    this->tileRows = tileRows;
    this->tileCols = tileCols;

    // create the data buffer and zero the entries
    data = new T[tileRows * tileCols * rowsPerTile * colsPerTile]();
}

/** 
 * @brief destructor of a tile matrix instance, deletes the buffer
 */
template <class T>
conflux::TileMatrix<T>::~TileMatrix()
{
    // delete the data buffer
    delete[] data;
}

/**
 * @brief tile access operator if type is VECTOR
 * 
 * @param i the tile row index
 * @throws CholeskyException if type is MATRIX
 * @returns pointer to the start of the tile
 */
template <class T>
T* conflux::TileMatrix<T>::get(uint32_t i)
{
    // throw CholeskyException if type is not vector, since one parameter
    // does not suffice in this case
    if (type != MatrixType::VECTOR) {
        throw CholeskyException(CholeskyException::errorCode::MatrixInsteadOfVector);
    }

    return data + i * (rowsPerTile * colsPerTile);
}

/**
 * @brief tile accessor operator if type is MATRIX
 *  
 * @param i the tile row index
 * @param j the tile column index
 * @throws CholeskyException if type is VECTOR
 * @returns pointer to start of the tile
 */
template <class T>
T* conflux::TileMatrix<T>::get(uint32_t i, uint32_t j)
{
    // throw CholeskyException if type is not vector, since the second parameter
    // would bear no meaning and this is clearly unintended behaviour
    if (this->type != MatrixType::MATRIX) {
        throw CholeskyException(CholeskyException::errorCode::VectorInsteadOfMatrix);
    }

    return data + i * (tileCols * rowsPerTile * colsPerTile)
                      + j * (rowsPerTile * colsPerTile);
}

#endif // TILEMATRIX_H
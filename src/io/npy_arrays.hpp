/*
 * =====================================================================================
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * -------------------------------------------------------------------------------------
 *
 *       Filename:  npy_arrays.hpp
 *
 *        Version:  1.0
 *        Created:  6/12/13 10:44:29
 *       Revision:  none
 *       Compiler:  gcc/intel compiler
 *
 *         Author:  Changxi Zheng (cz), cxz@cs.columbia.edu
 *                  Columbia University
 *
 * =====================================================================================
 */
 
#ifndef NPY_ARRAYS_HPP
#   define NPY_ARRAYS_HPP

#include <boost/multi_array.hpp>
#include "npy.h"
#include "multiarray_io.h"

template <size_t D>
static void npy_write_array_double(const boost::multi_array<double, D>& array, const char* file)
{
    int shape[D];
    const size_t* s = array.shape();
    for(std::size_t i = 0;i < D;++ i) shape[i] = static_cast<int>(s[i]);
    npy_save_double(file, 0, D, shape, array.data());
}

template <size_t D>
inline int mad_save_array_double(const boost::multi_array<double, D>& array, const char* file)
{
    int shape[D];
    const size_t* s = array.shape();
    for(std::size_t i = 0;i < D;++ i) shape[i] = static_cast<int>(s[i]);
    return save_multiarray_double(file, D, shape, array.data());
}

template <size_t D>
inline int mad_load_array_double(boost::multi_array<double, D>& array, const char* file)
{
    fprintf(stderr, "ERROR: Do not support multi_array dimension for %d\n", D);
    return -100;
}

template <>
inline int mad_load_array_double<2>(boost::multi_array<double, 2>& array, const char* file)
{
    int cnt;
    int ndims;
    FILE* fp = fopen(file, "rb");
    if ( fp == NULL ) return -1;

    char bytes[4];
    if ( fread(bytes, sizeof(char), 4, fp) != 4 ) goto ERR_READ;
    if ( bytes[0] != MAGIC_DOUBLE[0] || bytes[1] != MAGIC_DOUBLE[1] ||
         bytes[2] != MAGIC_DOUBLE[2] )
        goto ERR_FMT;

    if ( fread(&ndims, sizeof(int), 1, fp) != 1 ) goto ERR_READ;
    if ( (is_bigendian() && bytes[3] != '\x21') || 
         (!is_bigendian() && bytes[3] != '\x12') )
    {
        fprintf(stderr, "ERROR: Not support converting endianness now\n");
        fclose(fp);
        return -4;
    }
    if ( ndims != 2 ) goto ERR_FMT;
    int dims[2];
    if ( fread(dims, sizeof(int), 2, fp) != 2 ) goto ERR_READ;
    array.resize(boost::extents[dims[0]][dims[1]]);
    cnt = dims[0] * dims[1];
    if ( (int)fread(array.data(), sizeof(double), cnt, fp) != cnt ) goto ERR_READ;

    fclose(fp);
    return 0;

ERR_READ:
    fclose(fp);
    return -2;
ERR_FMT:
    fclose(fp);
    return -3;
}

#endif


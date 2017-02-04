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
 *       Filename:  npy.h
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
 
/*
Functions to create native numpy data files (.npy).
*/
#ifndef NPY_LIB
#   define NPY_LIB

#include<stdio.h>
#include<string.h>
#include<stdlib.h>

#ifdef __cplusplus
#include<complex>

extern "C" {
#else
#include<complex.h>
#endif


static const char MAGIC[] = "\x93NUMPY";
static const int  MAJOR = 1;
static const int  MINOR = 0;
static const int  MAX_HDR_LEN = 256 * 256;
static const int  MAX_INT_LEN = 32;
static const int  PREAMBLE_LEN = 6 + 1 + 1 + 2;

#if __BYTE_ORDER == __LITTLE_ENDIAN
static const char ENDIAN_CHAR = '<';
#else
static const char ENDIAN_CHAR = '>';
#endif

int create_metadata(char preamble[PREAMBLE_LEN], char header[MAX_HDR_LEN],
                    char* descr, int fortran_order, 
		    int ndims, int* shape);

void npy_save(const char* fname, char* descr, int fortran_order,
              int ndims, int* shape, size_t sz, const void* data);

void npy_save_double(const char* fname, int fortran_order, 
                     int ndims, int* shape, const double* data);

void npy_save_float(const char* fname, int fortran_order,
                    int ndims, int* shape, const float* data);

void npy_save_int(const char* fname, int fortran_order,
                  int ndims, int* shape, const int* data);

#ifndef __cplusplus
void npy_save_float_complex(const char* fname, int fortran_order,
                            int ndims, int* shape, 
                            const float complex *data);

void npy_save_double_complex(const char* fname, int fortran_order,
                            int ndims, int* shape,
                            const double complex *data);
#endif

#ifdef __cplusplus
}
#endif

#endif


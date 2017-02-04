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
 *       Filename:  multiarray_io.h
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
 
#ifndef MULTIARRAY_IO_H
#   define MULTIARRAY_IO_H

#undef __BEGIN_DECLS
#undef __END_DECLS

#ifdef	__cplusplus
#   define __BEGIN_DECLS    extern "C" {
#   define __END_DECLS      }
#else
#   define __BEGIN_DECLS
#   define __END_DECLS
#endif

#undef __p
#if defined (__STDC__) || defined (_AIX) \
    || (defined (__mips) && defined (_SYSTYPE_SVR4)) \
    || defined(WIN32) || defined(__cplusplus)
#   define __p(protos) protos
#else
#   define __p(protos) ()
#endif

__BEGIN_DECLS

static const int ONE = 1;
#define is_bigendian() ( (*(char*)&ONE) == 0 )

static const char MAGIC_DOUBLE[] = "MA\x9";

/*!
 * Format: <MAGIC> <type[1byte]> <1byte> <uint32_t> <uint32_t> ... [data]
 *                 <double/float/...> <endianness> <ndims> <shape...> <padding> [data]
 *
 * type: 0: int8_t
 * type: 1: uint8_t
 * type: 2: int16_t
 *       3: uint16_t
 *       4: int32_t
 *       5: uint32_t
 *       6: int64_t
 *       7: uint64_t
 *       8: float
 *       9: double
 *
 * return -1   Cannot open file
 */
int save_multiarray_double(const char* file, int ndims, const int* shape, const double* data);
//int load_multiarray_double();

__END_DECLS
#endif


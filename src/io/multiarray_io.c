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
 *       Filename:  multiarray_io.c
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
 
#include "multiarray_io.h"
#include "IOEndian.h"
#include <stdio.h>

int save_multiarray_double(const char* file, 
        int ndims, const int* shape, const double* data)
{
    int i, nbytes;
    char endian = is_bigendian() ? '\x21' : '\x12';
    FILE* fp = fopen(file, "wb");

    if ( fp == NULL ) return -1;
    /*** write header ***/
    /*  magic string  */
    if ( fwrite(MAGIC_DOUBLE, sizeof(char), 3, fp) != 3 ) goto ERR_RET;
    // write uint32_t for number of dimention
    if ( fwrite(&endian, sizeof(char), 1, fp) != 1 ) goto ERR_RET;
    if ( fwrite(&ndims, sizeof(int), 1, fp) != 1 ) goto ERR_RET;
    if ( fwrite(shape, sizeof(int), ndims, fp) != ndims ) goto ERR_RET;
    nbytes = (8+4*ndims) % 16;
    if ( nbytes )
    {
        char pad = '\0';
        for(;nbytes < 16;++ nbytes) 
            if ( fwrite(&pad, sizeof(char), 1, fp) != 1 ) goto ERR_RET;
    }
    nbytes = 1;
    for(i = 0;i < ndims;++ i) nbytes *= shape[i];
    if ( fwrite(data, sizeof(double), nbytes, fp) != nbytes ) goto ERR_RET;

    fflush(fp);
    fclose(fp);
    return 0;

ERR_RET:
    fclose(fp);
    return -2;
}



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
 *       Filename:  npy.c
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
 
/* Copyright 2009 William McLean.  */


#include "npy.h"
#include "IOEndian.h"

int create_metadata(char preamble[PREAMBLE_LEN], char header[MAX_HDR_LEN],
                    char * descr, int fortran_order,                     
		    int ndims, int* shape)
{
    unsigned char byte;
    uint16_t      hdrlen;
    int n, m, l, topad;
    /* 
    See numpy/lib/format.py for details of the .npy file format.
    */ 
    strcpy(header, "{'descr': '");
    strcat(header, descr);
    strcat(header, "', 'fortran_order': ");
    if ( fortran_order )
        strcat(header, "True, ");
    else
        strcat(header, "False, ");
    strcat(header, "'shape': (");
    for ( m=0; m<ndims; m++ ) {
        l = strlen(header);
	if ( shape[m] < 0 ) {
	    printf("shape[%d] = %d is negative!\n", m, shape[m]);
	    abort();
        }
	if ( l + MAX_INT_LEN + 4 >= MAX_HDR_LEN ) {
	    printf("header too long\n");
	    abort();
        }
	sprintf(header+l, "%d,", shape[m]);
    }
    l = strlen(header);
    if ( ndims > 1 ) header[l-1] = '\0'; // remove comma
    strcat(header, "), }");

    l = strlen(header);
    topad = 16 - ( PREAMBLE_LEN + l + 1 ) % 16;
    if ( l + topad + 1 > MAX_HDR_LEN ) {
        printf("header too long\n");
	abort();
    }
    for ( m=0; m<topad; m++ ) header[l+m] = ' ';
    l += topad;
    header[l] = '\n';
    header[++l] = '\0';

    strcpy(preamble, MAGIC);
    n = strlen(preamble);
    byte = MAJOR;
    preamble[n++] = byte;
    byte = MINOR;
    preamble[n++] = byte;
    hdrlen = (uint16_t)htole16(l);
    memcpy((void*)(preamble+n), (void*)&hdrlen, 2);
    return l;
}

void npy_save(const char* fname, char * descr, int fortran_order, 
              int ndims, int* shape, size_t sz, const void* data)
{
    char preamble[PREAMBLE_LEN], header[MAX_HDR_LEN];
    FILE *fp;
    int l, m, mtd_len;

    l = create_metadata(preamble, header, descr, fortran_order, ndims, shape);
    mtd_len = PREAMBLE_LEN + l;
    if ( mtd_len % 16 != 0 ) {
        printf("formatting error: metadata length %d not divisible by 16\n", 
	       mtd_len);
        abort();
    }
    fp = fopen(fname, "w");
    fwrite(preamble, sizeof(char), PREAMBLE_LEN, fp);
    fwrite(header,   sizeof(char), l, fp);
    l = 1;
    for ( m=0; m<ndims; m++ ) l *= shape[m];
    fwrite(data, sz, l, fp);
    fclose(fp);
}

void npy_save_double(const char* fname, int fortran_order,
                     int ndims, int* shape, const double* data)
{
    char descr[5];
    descr[0] = ENDIAN_CHAR;
    descr[1] = 'f';
    sprintf(descr+2, "%d", (int) sizeof(double));
    npy_save(fname, descr, fortran_order, ndims, shape, sizeof(double), data);
}

void npy_save_float(const char* fname, int fortran_order, 
                    int ndims, int* shape, const float* data)
{
    char descr[5];
    descr[0] = ENDIAN_CHAR;
    descr[1] = 'f';
    sprintf(descr+2, "%d", (int) sizeof(float));
    npy_save(fname, descr, fortran_order, ndims, shape, sizeof(float), data);
}

void npy_save_int(const char* fname, int fortran_order, 
                  int ndims, int* shape, const int* data)
{
    char descr[5];
    descr[0] = ENDIAN_CHAR;
    descr[1] = 'i';
    sprintf(descr+2, "%d", (int) sizeof(int));
    npy_save(fname, descr, fortran_order, ndims, shape, sizeof(int), data);
}

#ifndef __cplusplus
void npy_save_float_complex(const char* fname, int fortran_order, 
                            int ndims, int* shape, 
                            const float complex *data)
{
    char descr[5];
    size_t sz;

    descr[0] = ENDIAN_CHAR;
    descr[1] = 'c';
    sz = sizeof(float complex);
    sprintf(descr+2, "%d", (int) sz);
    npy_save(fname, descr, fortran_order, ndims, shape, sz, data);
}

void npy_save_double_complex(const char* fname, int fortran_order, 
                            int ndims, int* shape, 
                            const double complex *data)
{
    char descr[5];
    size_t sz;

    descr[0] = ENDIAN_CHAR;
    descr[1] = 'c';
    sz = sizeof(double complex);
    sprintf(descr+2, "%d", (int) sz);
    npy_save(fname, descr, fortran_order, ndims, shape, sz, data);
}
#endif


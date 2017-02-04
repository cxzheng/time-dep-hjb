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
 *       Filename:  InterpErr.cpp
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
 
#include "InterpErr.h"
#include <boost/filesystem.hpp>
#include <cstdlib>
#include "io/npy_arrays.hpp"

bool InterpErr::operator() (const boost::multi_array<double, 2>& u,
                            double ts, double& L1, double& Linf)
{
    namespace fs = boost::filesystem;

    char fname[128];
    sprintf(fname, fileFmt_.c_str(), ts);
    // check whether or not the file exists
    if ( !fs::exists(fs::path(fname)) ) return false;
    // load the data
    if ( mad_load_array_double(accU_, fname) != 0 )
    {
        fprintf(stderr, "Fail to load .mad file\n");
        exit(1);
    }
    int strX = (accU_.shape()[1]-1) / nx_;
    int strY = (accU_.shape()[0]-1) / ny_;
    if ( strX * nx_ != accU_.shape()[1]-1 || strY * ny_ != accU_.shape()[0]-1 )
    {
        fprintf(stderr, "Incorrect 2D array size (%d,%d) vs (%d,%d)\n",
                (int)accU_.shape()[1], (int)accU_.shape()[0], nx_, ny_);
        exit(1);
    }

    double sum = 0; Linf = 0;
#ifdef USE_OPENMP
    #pragma omp parallel default(none) shared(strX, strY, u, sum, Linf)
    {
    double mySum = 0, myMaxErr = 0;
    #pragma omp for nowait
    for(int iy = 0;iy <= ny_;++ iy)
    for(int ix = 0;ix <= nx_;++ ix)
    {
        double ev = fabs(accU_[iy*strY][ix*strX] - u[iy][ix]);
        mySum += ev;
        myMaxErr = fmax(ev, myMaxErr);
    }

    #pragma omp critical (UPDATE_ERR)
    {
        sum += mySum;
        Linf = fmax(Linf, myMaxErr);
    }
    }
#else
    for(int iy = 0;iy <= ny_;++ iy)
    for(int ix = 0;ix <= nx_;++ ix)
    {
        double ev = fabs(accU_[iy*strY][ix*strX] - u[iy][ix]);
        sum += ev;
        Linf = fmax(Linf, ev);
    }
#endif
    L1 = sum / double((nx_+1)*(ny_+1));
    return true;
}



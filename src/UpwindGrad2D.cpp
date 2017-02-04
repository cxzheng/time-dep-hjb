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
 *       Filename:  UpwindGrad2D.cpp
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
 
#include "UpwindGrad2D.h"
#include "utils/macros.h"

double UpwindGrad2D::dx(int ix, int iy, const boost::multi_array<double, 2>& u) const
{
    MSG_ASSERT(ix >= 0 && ix <= RES_.x && iy >= 0 && iy <= RES_.y, 
            "Given coordinate is out of domain");
    /* If node is on boundary, use one-side difference in any case */
    if ( ix == 0 ) return fmin(0, (u[iy][1] - u[iy][0]) / H_.x);
    if ( ix == RES_.x ) return fmax(0, (u[iy][ix] - u[iy][ix-1]) / H_.x);

    double ret = 0;
    if ( u[iy][ix] > u[iy][ix-1] ) ret = (u[iy][ix] - u[iy][ix-1])/H_.x;
    if ( u[iy][ix] > u[iy][ix+1] && u[iy][ix+1] < u[iy][ix-1] ) 
        ret = (u[iy][ix+1] - u[iy][ix]) / H_.x;
    return ret;
}

double UpwindGrad2D::dy(int ix, int iy, const boost::multi_array<double, 2>& u) const
{
    MSG_ASSERT(ix >= 0 && ix <= RES_.x && iy >= 0 && iy <= RES_.y, 
            "Given coordinate is out of domain");

    if ( iy == 0 ) return fmin(0, (u[1][ix] - u[0][ix]) / H_.y);
    if ( iy == RES_.y ) return fmax(0, (u[iy][ix] - u[iy-1][ix]) / H_.y);

    double ret = 0;
    if ( u[iy][ix] > u[iy-1][ix] ) ret = (u[iy][ix] - u[iy-1][ix])/H_.y;
    if ( u[iy][ix] > u[iy+1][ix] && u[iy+1][ix] < u[iy-1][ix] )
        ret = (u[iy+1][ix] - u[iy][ix]) / H_.y;
    return ret;
}



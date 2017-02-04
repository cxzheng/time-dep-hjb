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
 *       Filename:  UpwindGrad2D.h
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
 
#ifndef UPWIND_GRADIENT_2D_H
#   define UPWIND_GRADIENT_2D_H

#include "utils/arrays.hpp"
#include "vector2.hpp"

/*
 * Compute the first order upwind derivative
 */
class UpwindGrad2D
{
    public:
        UpwindGrad2D(double lx, double ly, int nx, int ny):
                RES_(nx, ny), H_(lx/nx, ly/ny) { }

        double dx(int ix, int iy, const boost::multi_array<double, 2>& u) const;
        double dy(int ix, int iy, const boost::multi_array<double, 2>& u) const;

        inline double dx(const vector2i& pos, const boost::multi_array<double, 2>& u) const
        {   return dx(pos.x, pos.y, u);  }
        inline double dy(const vector2i& pos, const boost::multi_array<double, 2>& u) const
        {   return dy(pos.x, pos.y, u);  }

    private:
        const vector2i      RES_;       // lattice resolution
        const vector2d      H_;
};

#endif


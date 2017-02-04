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
 *       Filename:  LatticeData2D.h
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
 
#ifndef LATTICE_DATA_2D_H
#   define LATTICE_DATA_2D_H

#include "vector2.hpp"
#include "utils/arrays.hpp"

//! Basis data structure to maintain a 2d lattice
class LatticeData2D
{
    public:
        typedef boost::multi_array<double, 2>   TValField;

        LatticeData2D(double lx, double ly, int nx, int ny):
            SZ_(lx, ly), RES_(nx, ny), H_(lx/nx, ly/ny), curUId_(0)
        {
            u_[0].resize(boost::extents[ny+1][nx+1]);
            u_[1].resize(boost::extents[ny+1][nx+1]);
        }

        int num_nodes() const
        {   return (RES_.x+1)*(RES_.y+1); }

        inline bool in_bound(int ix, int iy) const
        {   return ix >= 0 && ix <= RES_.x && iy >= 0 && iy <= RES_.y; }

        const boost::multi_array<double, 2>& current_u() const
        {   return u_[curUId_]; }

        inline void get_node_pos(int ix, int iy, vector2d& pos) const
        {
            assert(in_bound(ix, iy));
            pos.set(ix*H_.x, iy*H_.y);
        }

    protected:
        const vector2d                  SZ_;
        const vector2i                  RES_;       // lattice resolution
        const vector2d                  H_;
        int                             curUId_;
        boost::multi_array<double, 2>   u_[2];      // U value in both current time step and 

        static const int NEIGH_DIRS[4][2];
        static const int SOLVE_DIR[4][2];
};

const int LatticeData2D::NEIGH_DIRS[4][2] = { {1, 0}, {0, 1}, {-1, 0}, {0, -1} };
const int LatticeData2D::SOLVE_DIR[4][2] = { {1, -1}, {0, -1}, {1, 1}, {0, 1} };

#endif


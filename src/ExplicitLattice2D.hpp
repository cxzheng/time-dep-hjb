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
 *       Filename:  ExplicitLattice2D.hpp
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
 
#ifndef EXPLICIT_LATTICE_HPP
#   define EXPLICIT_LATTICE_HPP

#include <stdint.h>
#include "LatticeData2D.h"

#include "utils/nano_timer.h"   // $test

/*!
 * The equation to solve is the time dependent HJ equation in 2D
 * \f[
 * u_t + F(x,t)|\nabla u| = g(x,t)
 * \f]
 * This class is to solve this equation using explicit time step method.
 *
 * - _F:     the F function in the equation
 * - _G:     the G function in the equation
 * - _TGrad: the class to compute numeric gradient
 * - _TBC:   the boundary condition
 */
template <typename _F, typename _G, class _TGrad, class _TBC>
class ExplicitLattice2D : public LatticeData2D
{
    public:
        ExplicitLattice2D(double lx, double ly, int nx, int ny, _TGrad* gd, _TBC* bc);

        void init(double ts);

        void advance(double dt);

        double time() const
        {   return ts_; }

    private:
        double      ts_;
        boost::multi_array<uint8_t, 2>  labels_;

        _TGrad*     gd_;
        _TBC*       bc_;
        _F          f_;
        _G          g_;
};

/////////////////////////////////////////////////////////////////////////////////////////

template <typename _F, typename _G, class _TGrad, class _TBC>
ExplicitLattice2D<_F,_G,_TGrad,_TBC>::ExplicitLattice2D(
        double lx, double ly, int nx, int ny, _TGrad* gd, _TBC* bc):
    LatticeData2D(lx, ly, nx, ny), labels_(boost::extents[ny+1][nx+1]), 
    gd_(gd), bc_(bc)
{
    assert(gd);
    assert(bc);
}

template <typename _F, typename _G, class _TGrad, class _TBC>
void ExplicitLattice2D<_F,_G,_TGrad,_TBC>::init(double ts)
{
    ts_ = ts;
#ifdef USE_OPENMP
    #pragma omp parallel for default(none) shared(ts)
#endif
    for(int iy = 0;iy <= RES_.y;++ iy)
    for(int ix = 0;ix <= RES_.x;++ ix)
        u_[curUId_][iy][ix] = bc_->boundary_val(ix, iy, ts);
}

template <typename _F, typename _G, class _TGrad, class _TBC>
void ExplicitLattice2D<_F,_G,_TGrad,_TBC>::advance(double dt)
{
    const int lastUId =  curUId_;
    curUId_ = 1 - curUId_;

    zero_multi_array(labels_);
    //// setup boundary values
    const int BN = bc_->num_boundary_nodes();
    const vector2i* pos = bc_->boundary_pos();
    const double*   bv  = bc_->boundary_val(ts_+dt);
#ifdef USE_OPENMP
    #pragma omp parallel for default(none) shared(bv, pos)
#endif
    for(int i = 0;i < BN;++ i)
    {
        labels_[pos[i].y][pos[i].x] = 1;
        u_[curUId_][pos[i].y][pos[i].x] = bv[i];
    }

#ifdef USE_OPENMP
    #pragma omp parallel for default(none) shared(dt)
#endif
    for(int iy = 0;iy <= RES_.y;++ iy)
    for(int ix = 0;ix <= RES_.x;++ ix)
        if ( !labels_[iy][ix] )
        {
            //// compute |\nabla u| using previous time step
            double dx = gd_->dx(ix, iy, u_[lastUId]);
            double dy = gd_->dy(ix, iy, u_[lastUId]);
            double gdu = sqrt(dx*dx + dy*dy);

            vector2d cpos(ix*H_.x, iy*H_.y);
            double fval = f_(cpos, ts_);    // f value
            double gval = g_(cpos, ts_);    // g value
            u_[curUId_][iy][ix] = (gval - fval*gdu)*dt + u_[lastUId][iy][ix];
        }

    ts_ += dt;
}

#endif



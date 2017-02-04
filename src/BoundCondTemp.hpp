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
 *       Filename:  BoundCondTemp.hpp
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
 
#ifndef BOUND_COND_TEMP_HPP
#   define BOUND_COND_TEMP_HPP

#include <vector>
#include <string.h>

/*!
 * Define the boundary condition as
 * \f[ u(x,y,t) = y + K(t + y) \f]
 * boundary is just one side, specified only on [0..1]x{0}
 */
template <class _K>
class BoundCond2D_1Side
{
    public:
        BoundCond2D_1Side(double lx, double ly, int nx, int ny);

        int num_boundary_nodes() const
        {   return totN_; }

        const vector2i* boundary_pos() const        // all the boundary positions
        {   return &bpos_[0]; }

        const double* boundary_val(double);         // boundary values

        /*!
         * time dependent boundary value
         * return u(ix*h, iy*h, ts)
         */
        double boundary_val(int , int iy, double ts) const
        {
            double dy = iy*sy_;
            return dy + kfunc_(ts + dy);
        }

    private:
        int     nx_, ny_;
        int     totN_;
        double  lx_, ly_;
        double  sx_, sy_;
        std::vector<vector2i>   bpos_;
        std::vector<double>     bv_;
        _K      kfunc_;
};

template <class _K>
BoundCond2D_1Side<_K>::BoundCond2D_1Side(double lx, double ly, int nx, int ny):
        nx_(nx), ny_(ny), totN_(nx_+1), lx_(lx), ly_(ly), sx_(lx/nx), sy_(ly/ny)
{
    bpos_.resize(totN_);
    bv_.resize(totN_);
    memset(&bv_[0], 0, sizeof(double)*totN_);

    for(int i = 0;i <= nx;++ i) bpos_[i].set(i, 0);
}

template <class _K>
const double* BoundCond2D_1Side<_K>::boundary_val(double ts)
{
#ifdef USE_OPENMP
    #pragma omp parallel for default(none) shared(ts)
#endif
    for(int i = 0;i < totN_;++ i)
        bv_[i] = kfunc_(ts);
    return &bv_[0];
}

///////////////////////////////////////////////////////////////////////////////////////////////////

/*!
 * Define the boundary condition as
 * \f[ u(x,y,t) = y + y^2 + K(t + y) \f]
 * boundary is just one side, specified only on [0..1]x{0}
 * u is the solution of equation:
 * \f[ u_t -|\nabla u| = -(2y + 1)\f]
 */
template <class _K>
class BoundCond2D_1Side_V2
{
    public:
        BoundCond2D_1Side_V2(double lx, double ly, int nx, int ny);

        int num_boundary_nodes() const
        {   return totN_; }

        const vector2i* boundary_pos() const        // all the boundary positions
        {   return &bpos_[0]; }

        const double* boundary_val(double);         // boundary values

        /*!
         * time dependent boundary value
         * return u(ix*h, iy*h, ts)
         */
        double boundary_val(int, int iy, double ts) const
        {
            double dy = iy*sy_;
            return dy + dy*dy + kfunc_(ts + dy);
        }

    private:
        int     nx_, ny_;
        int     totN_;
        double  lx_, ly_;
        double  sx_, sy_;
        std::vector<vector2i>   bpos_;
        std::vector<double>     bv_;
        _K      kfunc_;
};

template <class _K>
BoundCond2D_1Side_V2<_K>::BoundCond2D_1Side_V2(double lx, double ly, int nx, int ny):
        nx_(nx), ny_(ny), totN_(nx_+1), lx_(lx), ly_(ly), sx_(lx/nx), sy_(ly/ny)
{
    bpos_.resize(totN_);
    bv_.resize(totN_);
    memset(&bv_[0], 0, sizeof(double)*totN_);

    for(int i = 0;i <= nx;++ i) bpos_[i].set(i, 0);
}

template <class _K>
const double* BoundCond2D_1Side_V2<_K>::boundary_val(double ts)
{
    const double v = kfunc_(ts);
#ifdef USE_OPENMP
    #pragma omp parallel for default(none) shared(ts)
#endif
    for(int i = 0;i < totN_;++ i)
        bv_[i] = v;
    return &bv_[0];
}

///////////////////////////////////////////////////////////////////////////////////////////////////

/*!
 * Define the boundary condition as
 * \f[ u(x,y,t) = y + y^2 + K(t + y^2 + y) \f]
 * boundary is just one side, specified only on [0..1]x{0}
 * u is the solution of equation:
 * \f[ u_t - \frac{|\nabla u|}{2y+1} = -1 \f]
 */
template <class _K>
class BoundCond2D_1Side_V3
{
    public:
        BoundCond2D_1Side_V3(double lx, double ly, int nx, int ny);

        int num_boundary_nodes() const
        {   return totN_; }

        const vector2i* boundary_pos() const        // all the boundary positions
        {   return &bpos_[0]; }

        const double* boundary_val(double);         // boundary values

        /*!
         * time dependent boundary value
         * return u(ix*h, iy*h, ts)
         */
        double boundary_val(int , int iy, double ts) const
        {
            double dy = iy*sy_;
            double yval = dy + dy*dy;
            return yval + kfunc_(ts + yval);
        }

    private:
        int     nx_, ny_;
        int     totN_;
        double  lx_, ly_;
        double  sx_, sy_;
        std::vector<vector2i>   bpos_;
        std::vector<double>     bv_;
        _K      kfunc_;
};

template <class _K>
BoundCond2D_1Side_V3<_K>::BoundCond2D_1Side_V3(double lx, double ly, int nx, int ny):
        nx_(nx), ny_(ny), totN_(nx_+1), lx_(lx), ly_(ly), sx_(lx/nx), sy_(ly/ny)
{
    bpos_.resize(totN_);
    bv_.resize(totN_);
    memset(&bv_[0], 0, sizeof(double)*totN_);

    for(int i = 0;i <= nx;++ i) bpos_[i].set(i, 0);
}

template <class _K>
const double* BoundCond2D_1Side_V3<_K>::boundary_val(double ts)
{
#ifdef USE_OPENMP
    #pragma omp parallel for default(none) shared(ts)
#endif
    for(int i = 0;i < totN_;++ i)
        bv_[i] = kfunc_(ts);
    return &bv_[0];
}
#endif



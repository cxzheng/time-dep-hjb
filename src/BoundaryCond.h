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
 *       Filename:  BoundaryCond.h
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
 
#ifndef BOUNDARY_COND_H
#   define BOUNDARY_COND_H

#include <vector>
#include <math.h>
#include "vector2.hpp"

class SimpleBoundaryCond2D          // TESTED
{
    public:
        SimpleBoundaryCond2D(int nx, int ny);

        int num_boundary_nodes() const
        {   return totN_; }

        const vector2i* boundary_pos() const    // all the boundary positions
        {   return &bpos_[0]; }

        const double* boundary_val(double) const    // boundary values
        {   return &bv_[0]; }

        /*!
         * Boundary condition on the top slide
         * time dependent boundary value
         * return u(ix*h, iy*h, ts)
         */
        double boundary_val(int, int, double) const
        {   return 0.; }

    private:
        int     nx_, ny_;
        int     totN_;
        std::vector<vector2i>   bpos_;
        std::vector<double>     bv_;
};

class SimpleBoundaryCond2D_01          // TESTED
{
    public:
        SimpleBoundaryCond2D_01(int nx, int ny);

        int num_boundary_nodes() const
        {   return totN_; }

        const vector2i* boundary_pos() const    // all the boundary positions
        {   return &bpos_[0]; }

        const double* boundary_val(double) const    // boundary values
        {   return &bv_[0]; }

        /*!
         * Boundary condition on the top slide
         * time dependent boundary value
         * return u(ix*h, iy*h, ts)
         */
        double boundary_val(int, int, double) const
        {   return 0.; }

    private:
        int     nx_, ny_;
        int     totN_;
        std::vector<vector2i>   bpos_;
        std::vector<double>     bv_;
};

// ----------------------------------------------------------------------------

class SimpleBoundaryCond2D_V1
{
    public:
        SimpleBoundaryCond2D_V1(double lx, double ly, int nx, int ny);

        int num_boundary_nodes() const
        {   return totN_; }

        const vector2i* boundary_pos() const    // all the boundary positions
        {   return &bpos_[0]; }

        const double* boundary_val(double) const    // boundary values
        {   return &bv_[0]; }

        /*!
         * time dependent boundary value
         * return u(ix*h, iy*h, ts)
         */
        double boundary_val(int ix, int iy, double) const
        {
            double dx = ix*sx_, dy = iy*sy_;
            return fmin(dx, fmin(dy, fmin(lx_-dx, ly_-dy)));
        }

    private:
        int     nx_, ny_;
        int     totN_;
        double  lx_, ly_;
        double  sx_, sy_;
        std::vector<vector2i>   bpos_;
        std::vector<double>     bv_;
};

class SimpleBoundaryCond2D_V2
{
    public:
        SimpleBoundaryCond2D_V2(int nx, int ny);

        int num_boundary_nodes() const
        {   return totN_; }

        const vector2i* boundary_pos() const    // all the boundary positions
        {   return &bpos_[0]; }

        const double* boundary_val(double);     // boundary values

        /*!
         * time dependent boundary value
         * return u(ix*h, iy*h, ts)
         */
        double boundary_val(int, int, double ts) const
        {   return ts; }

    private:
        int     nx_, ny_;
        int     totN_;
        std::vector<vector2i>   bpos_;
        std::vector<double>     bv_;
};

/*!
 * Define the boundary condition as
 * \f[ u(x,y,t) = (x+1)(y+1) + t \f]
 */
class SimpleBoundaryCond2D_V3
{
    public:
        SimpleBoundaryCond2D_V3(double lx, double ly, int nx, int ny);

        int num_boundary_nodes() const
        {   return totN_; }

        const vector2i* boundary_pos() const        // all the boundary positions
        {   return &bpos_[0]; }

        const double* boundary_val(double);         // boundary values

        /*!
         * time dependent boundary value
         * return u(ix*h, iy*h, ts)
         */
        double boundary_val(int ix, int iy, double ts) const
        {
            double dx = ix*sx_, dy = iy*sy_;
            return (dx+1.)*(dy+1.) + ts;
        }

    private:
        int     nx_, ny_;
        int     totN_;
        double  lx_, ly_;
        double  sx_, sy_;
        std::vector<vector2i>   bpos_;
        std::vector<double>     bv_;
};

/*!
 */
class SimpleBoundaryCond2D_V4
{
    public:
        SimpleBoundaryCond2D_V4(int nx, int ny);

        int num_boundary_nodes() const
        {   return 1; }

        const vector2i* boundary_pos() const        // all the boundary positions
        {   return &bpos_; }

        const double* boundary_val(double) const         // boundary values
        {   return &bv_; }

        /*!
         * time dependent boundary value
         * return u(ix*h, iy*h, ts)
         */
        double boundary_val(int, int, double) const
        {   return 0.; }

    private:
        int     nx_, ny_;
        vector2i   bpos_;
        double     bv_;
};

class SimpleBoundaryCond2D_V42
{
    public:
        SimpleBoundaryCond2D_V42(int nx, int ny):nx_(nx), ny_(ny)
        {   bpos_.set(nx_/2, ny_/2); }

        int num_boundary_nodes() const
        {   return 1; }

        const vector2i* boundary_pos() const        // all the boundary positions
        {   return &bpos_; }

        const double* boundary_val(double t)        // boundary values
        {   
            //bv_ = cos(4*M_PI*t+M_PI)+1.;
            bv_ = 0;
            return &bv_; 
        }

        /*!
         * time dependent boundary value
         * return u(ix*h, iy*h, ts)
         */
        double boundary_val(int, int, double) const
        {   return 0; }

    private:
        int         nx_, ny_;
        vector2i    bpos_;
        double      bv_;
};

/*!
 * Define the boundary condition as
 * \f[ u(x,y,t) = (x+1)^2(y+1)^2 + t \f]
 */
class SimpleBoundaryCond2D_V5
{
    public:
        SimpleBoundaryCond2D_V5(double lx, double ly, int nx, int ny);

        int num_boundary_nodes() const
        {   return totN_; }

        const vector2i* boundary_pos() const        // all the boundary positions
        {   return &bpos_[0]; }

        const double* boundary_val(double);         // boundary values

        /*!
         * time dependent boundary value
         * return u(ix*h, iy*h, ts)
         */
        double boundary_val(int ix, int iy, double ts) const
        {
            double dx = ix*sx_+K, dy = iy*sy_+K;
            return dx*dx*dy*dy + ts;
        }

    private:
        static const double K;

        int     nx_, ny_;
        int     totN_;
        double  lx_, ly_;
        double  sx_, sy_;
        std::vector<vector2i>   bpos_;
        std::vector<double>     bv_;
};

/*!
 * Same as BoundaryCond2D_V5 but use zero boundary condition on 
 * top slide
 */
class SimpleBoundaryCond2D_V50
{
    public:
        SimpleBoundaryCond2D_V50(double lx, double ly, int nx, int ny);

        int num_boundary_nodes() const
        {   return totN_; }

        const vector2i* boundary_pos() const        // all the boundary positions
        {   return &bpos_[0]; }

        const double* boundary_val(double);         // boundary values

        /*!
         * time dependent boundary value
         * return u(ix*h, iy*h, ts)
         */
        double boundary_val(int ix, int iy, double ts) const
        {   return 28.; }

    private:
        static const double K;

        int     nx_, ny_;
        int     totN_;
        double  lx_, ly_;
        double  sx_, sy_;
        std::vector<vector2i>   bpos_;
        std::vector<double>     bv_;
};

/*!
 * Define the boundary condition as
 * \f[ u(x,y,t) = (x+1)^2(y+1)^2 e^t \f]
 */
class SimpleBoundaryCond2D_V6           // TESTED
{
    public:
        SimpleBoundaryCond2D_V6(double beta, double lx, double ly, int nx, int ny);

        int num_boundary_nodes() const
        {   return totN_; }

        const vector2i* boundary_pos() const        // all the boundary positions
        {   return &bpos_[0]; }

        const double* boundary_val(double);         // boundary values

        /*!
         * time dependent boundary value
         * return u(ix*h, iy*h, ts)
         */
        double boundary_val(int ix, int iy, double ts) const
        {
            double dx = ix*sx_+1., dy = iy*sy_+1.;
            return dx*dx*dy*dy*exp(beta_*ts);
        }

    private:
        const double beta_;
        int     nx_, ny_;
        int     totN_;
        double  lx_, ly_;
        double  sx_, sy_;
        std::vector<vector2i>   bpos_;
        std::vector<double>     bv_;
};

/*!
 * Define the boundary condition as
 * \f[ u(x,y,t) = (x+1)^2(y+1)^2 e^t \f]
 * boundary is specified only on [0..1]x{0} and {0}x[0..1]
 */
class SimpleBoundaryCond2D_V7   // TESTED
{
    public:
        SimpleBoundaryCond2D_V7(double lx, double ly, int nx, int ny);

        int num_boundary_nodes() const
        {   return totN_; }

        const vector2i* boundary_pos() const        // all the boundary positions
        {   return &bpos_[0]; }

        const double* boundary_val(double);         // boundary values

        /*!
         * time dependent boundary value
         * return u(ix*h, iy*h, ts)
         */
        double boundary_val(int ix, int iy, double ts) const
        {
            double dx = ix*sx_+1., dy = iy*sy_+1.;
            return dx*dx*dy*dy*exp(ts);
        }

    private:
        int     nx_, ny_;
        int     totN_;
        double  lx_, ly_;
        double  sx_, sy_;
        std::vector<vector2i>   bpos_;
        std::vector<double>     bv_;
};

/*!
 * Define the boundary condition as
 * \f[ u(x,y,t) = (x+1)^2(y+1)^2 e^{t_max-t} \f]
 */
class SimpleBoundaryCond2D_V8
{
    public:
        SimpleBoundaryCond2D_V8(double tmax, double lx, double ly, int nx, int ny);

        int num_boundary_nodes() const
        {   return totN_; }

        const vector2i* boundary_pos() const        // all the boundary positions
        {   return &bpos_[0]; }

        const double* boundary_val(double);         // boundary values

        /*!
         * time dependent boundary value
         * return u(ix*h, iy*h, ts)
         */
        double boundary_val(int ix, int iy, double ts) const
        {
            double dx = ix*sx_+1., dy = iy*sy_+1.;
            return dx*dx*dy*dy*exp(tmax_ - ts);
        }

    private:
        const double tmax_;
        int     nx_, ny_;
        int     totN_;
        double  lx_, ly_;
        double  sx_, sy_;
        std::vector<vector2i>   bpos_;
        std::vector<double>     bv_;
};

class SimpleBoundaryCond2D_V9
{
    public:
        SimpleBoundaryCond2D_V9(double lx, double ly, int nx, int ny);

        int num_boundary_nodes() const
        {   return totN_; }

        const vector2i* boundary_pos() const
        {   return &bpos_[0]; }

        const double* boundary_val(double);

        /*! time dependent boundary value */
        inline double boundary_val(int ix, int iy, double ts) const
        {
            double tx = ix*sx_, ty = iy*sy_;
            double xy = tx + ty;
            return sin(xy*ts) + xy;
        }

    private:
        int     nx_, ny_;
        int     totN_;
        double  lx_, ly_;
        double  sx_, sy_;
        std::vector<vector2i>   bpos_;
        std::vector<double>     bv_;
};

class SimpleBoundaryCond2D_V10
{
    public:
        SimpleBoundaryCond2D_V10(double lx, double ly, int nx, int ny);

        int num_boundary_nodes() const
        {   return totN_; }

        const vector2i* boundary_pos() const        // all the boundary positions
        {   return &bpos_[0]; }

        const double* boundary_val(double);         // boundary values

        /*!
         * time dependent boundary value
         * return u(ix*h, iy*h, ts)
         */
        double boundary_val(int ix, int, double ts) const
        {
            double dx = ix*sx_+1.;
            return dx*dx + ts;
        }

    private:
        int     nx_, ny_;
        int     totN_;
        double  lx_, ly_;
        double  sx_, sy_;
        std::vector<vector2i>   bpos_;
        std::vector<double>     bv_;
};

class SimpleBoundaryCond2D_V11
{
    public:
        SimpleBoundaryCond2D_V11(double lx, double ly, int nx, int ny);

        int num_boundary_nodes() const
        {   return totN_; }

        const vector2i* boundary_pos() const        // all the boundary positions
        {   return &bpos_[0]; }

        const double* boundary_val(double);         // boundary values

        /*!
         * time dependent boundary value
         * return u(ix*h, iy*h, ts)
         */
        double boundary_val(int ix, int iy, double ts) const
        {
            const double K = 0.1;
            double vx = ix*sx_ + K;
            double vy = iy*sy_ + K;
            return vx*vx*vx*vx + vy + ts;
        }

    private:
        int     nx_, ny_;
        int     totN_;
        double  lx_, ly_;
        double  sx_, sy_;
        std::vector<vector2i>   bpos_;
        std::vector<double>     bv_;
};
#endif



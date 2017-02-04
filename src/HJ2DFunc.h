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
 *       Filename:  HJ2DFunc.h
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
 
#ifndef HJ_2D_FUNC_H
#   define HJ_2D_FUNC_H

#include <math.h>
#include "vector2.hpp"

struct func_F1
{
    inline double operator() (const vector2d&, double) const
    {   return -1; }

    inline double operator() (const vector2d&) const
    {   return -1; }
};

struct func_G1
{
    inline double operator() (const vector2d&, double) const
    {   return -1; }
};

//! Vanishing g function, to model minimum-time to the boundary
struct func_G0
{
    inline double operator() (const vector2d&, double) const
    {   return 0; }
};

/*!
 * Function \f$F\f$ is defined as
 * \f[
 * F(x,y) = \frac{1}{\sqrt{x^2+y^2}}
 * \f]
 */
struct func_F_V3
{
    inline double operator() (const vector2d& pos, double) const
    {   return -1./sqrt((pos.x+1.)*(pos.x+1.) + (pos.y+1.)*(pos.y+1.)); }
};


template <int _C>
struct func_G_const
{
    inline double operator() (const vector2d&, double) const
    {   return (double)_C; }
};

/*!
 * Oscillatory speed profile
 */
struct func_F_V4
{
    inline double operator() (const vector2d& pos, double) const
    {
        return -1 - 0.8 * sin(8.*M_PI*pos.x) * sin(8.*M_PI*pos.y);
    }

    inline double operator() (const vector2d& pos) const
    {
        return -1 - 0.8 * sin(8.*M_PI*pos.x) * sin(8.*M_PI*pos.y);
    }
};

struct func_F_V5
{
    inline double operator() (const vector2d& pos, double t) const
    {
        //return -1 - 0.8 * sin(8.*M_PI*pos.x) * sin(8.*M_PI*pos.y) * sin(8.*M_PI*t);
        return -1 - 0.8 * sin(8.*M_PI*pos.x) * sin(8.*M_PI*pos.y) * sin(M_PI*t);
        //return -1 - 0.8 * sin(8.*M_PI*pos.x) * sin(8.*M_PI*pos.y) * t * 0.5;
    }
};

struct func_F_V51
{
    inline double operator() (const vector2d& pos, double t) const
    {
        //return -1 - 0.8 * sin(8.*M_PI*pos.x) * sin(M_PI*t);
        return -1 - 0.8 * sin(8.*M_PI*pos.x); // * sin(M_PI*t);
    }
};

/*!
 * suppose u = (x+1)^2(y+1)^2 + t
 * then |\nabla u| = 2(x+1)(y+1)\sqrt((x+1)^2+(y+1)^2)
 *      f = 1 / |\nabla u|
 */
struct func_F_V6
{
    inline double operator() (const vector2d& pos, double) const
    {
        static const double K = 1.;
        const vector2d vv = pos + K;
        return -1./(2.*vv.x*vv.y*vv.length());
    }
};

/*!
 * suppose u = (x+1)^2(y+1)^2 + t
 * then |\nabla u| = 2(x+1)(y+1)\sqrt((x+1)^2+(y+1)^2)
 *      f = 2 / |\nabla u|
 */
struct func_F_V7
{
    inline double operator() (const vector2d& pos, double) const
    {
        const vector2d vv = pos + 1.;
        return -1./(vv.x*vv.y*vv.length());
    }
};

/*!
 * suppose u = (x+1)^2(y+1)^2 e^t
 * then |\nabla u| = 2(x+1)(y+1)e^t\sqrt((x+1)^2+(y+1)^2)
 *      f = -(x+1)(y+1) / (2*\sqrt{(x+1)^2+(y+1)^2})
 */
struct func_F_V8
{
    inline double operator() (const vector2d& pos, double) const
    {
        const vector2d vv = pos + 1.;
        return -vv.x*vv.y*BETA/(2.*vv.length());
    }

    static double BETA;
};

/*!
 * suppose u = (x+1)^2(y+1)^2 e^t
 * then |\nabla u| = 2(x+1)(y+1)e^t\sqrt((x+1)^2+(y+1)^2)
 *      f = -(x+1)(y+1) / (\sqrt{(x+1)^2+(y+1)^2})
 */
struct func_F_V9
{
    inline double operator() (const vector2d& pos, double) const
    {
        const vector2d vv = pos + 1.;
        return -vv.x*vv.y/vv.length();
    }
};

/*
 * f = -1/(2y+1)
 * the PDE is  u_t - f|\nabla u| = -1
 */
struct func_F_V10
{
    inline double operator() (const vector2d& pos, double) const
    {   return -1. / (2*pos.y + 1.); }

    inline double operator() (const vector2d& pos) const
    {   return -1. / (2*pos.y + 1.); }
};

struct func_F_V11
{
    inline double operator() (const vector2d& pos, double ts) const
    {
        static const double S = sqrt(2);
        return -(pos.x + pos.y) / (S*ts);
    }
};

struct func_F_V12
{
    inline double operator() (const vector2d& pos, double) const
    {
        return -1. / (2*(pos.x + 1));
    }
};

struct func_F_V13
{
    inline double operator() (const vector2d& pos, double) const
    {
        const double K = 0.1;
        const vector2d vv = pos + K;
        const double ux = 3.*vv.x*vv.x*vv.x;
        return -sqrt(ux*ux + 1.);
    }
};

/*!
 * suppose u = (x+1)^2(y+1)^2 e^t
 * then |\nabla u| = 2(x+1)(y+1)e^t\sqrt((x+1)^2+(y+1)^2)
 *      f = -(x+1)(y+1) / (\sqrt{(x+1)^2+(y+1)^2})
 * u_t - f|\nabla u| = -u
 */
struct func_G_V1
{
    inline double operator() (const vector2d& pos, double ts) const
    {   
        const vector2d vv = pos + 1.;
        return -vv.x*vv.x*vv.y*vv.y*exp(ts);
    }
};

/*!
 * suppose u = (x+1)^2(y+1)^2 e^(T-t)
 * then |\nabla u| = 2(x+1)(y+1)e^t\sqrt((x+1)^2+(y+1)^2)
 *      f = -(x+1)(y+1) / (\sqrt{(x+1)^2+(y+1)^2})
 * u_t - f|\nabla u| = -u
 */
struct func_G_V2
{
    inline double operator() (const vector2d& pos, double ts) const
    {   
        const vector2d vv = pos + 1.;
        return -2.*vv.x*vv.x*vv.y*vv.y*exp(T_MAX - ts);
    }

    static double T_MAX;
};

/* f = exp(\beta * s) */
struct func_exp
{
    inline double operator() (double s) const
    {   return exp(BETA*s); }

    static double BETA;
};

/*!
 * g = -(2y + 1)
 */
struct func_G_V3
{
    inline double operator() (const vector2d& pos, double) const
    {   return -(2.*pos.y + 1); }
};

struct func_G_V4
{
    inline double operator() (const vector2d& pos, double t) const
    {   return - (pos.x + pos.y) / t; }
};

#endif



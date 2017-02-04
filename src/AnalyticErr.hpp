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
 *       Filename:  AnalyticErr.hpp
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
 
#ifndef ANALYTIC_ERR_HPP
#   define ANALYTIC_ERR_HPP

#include <cmath>
#include "vector2.hpp"
#include "utils/arrays.hpp"

/*
 * The analytic solution:
 *     u = (x+1)^2(y+1)^2*exp(t)
 */
struct AnalySol_V1
{
    inline double operator() (const vector2d& pos, double ts) const
    {
        const vector2d vv = pos + 1.;
        const double aaa = vv.x * vv.y;     // (pos.x+1)*(pos.y+1)
        return aaa*aaa*exp(ts);             // (pos.x+1)^2 (pos.y+1)^2 * exp(ts)
    }
};

/*
 * The analytic solution:
 *     u = (x+1)^2(y+1)^2 + t
 */
struct AnalySol_V2
{
    inline double operator() (const vector2d& pos, double ts) const
    {
        static const double K = 1.;
        const vector2d vv = pos + K;    // (x+1), (y+1)
        const double aaa = vv.x * vv.y; // (x+1)(y+1)
        return aaa*aaa + ts;            // (x+1)^2(y+1)^2 + t
    }
};

/*
 * The analytic solution:
 *     u = (x+1)^2(y+1)^2*exp(T-t)
 */
struct AnalySol_V3
{
    inline double operator() (const vector2d& pos, double ts) const
    {
        const vector2d vv = pos + 1.;
        const double aaa = vv.x * vv.y;
        return aaa*aaa*exp(ts);
    }
    static double T_MAX;
};

template <class _K>
struct AnalySol_V4
{
    inline double operator() (const vector2d& pos, double ts) const
    {   return pos.y + kfunc_(pos.y + ts); }

    _K      kfunc_;
};

/*!
 * solution of the equation: 
 *   u_t - |\nabla u| = -(2y+1)
 * the solution is 
 *   u = K(t+y) + y^2 + y
 */
template <class _K>
struct AnalySol_V5
{
    inline double operator() (const vector2d& pos, double ts) const
    {   return pos.y + pos.y*pos.y + kfunc_(pos.y + ts); }

    _K      kfunc_;
};

/*!
 * solution of the equation: 
 *   u_t - \frac{|\nabla u|}{2y+1} = -1
 * the solution is 
 *   u = K(t+y+y^2) + y^2 + y
 */
template <class _K>
struct AnalySol_V6
{
    inline double operator() (const vector2d& pos, double ts) const
    {   
        double yval = pos.y + pos.y*pos.y;
        return yval + kfunc_(yval + ts); 
    }

    _K      kfunc_;
};

struct AnalySol_V7
{
    inline double operator() (const vector2d& pos, double ts) const
    {
        return fmin(T_MAX-ts, fmin(pos.x, fmin(pos.y, fmin(1.-pos.x, 1.-pos.y))));
    }
    static double T_MAX;
};

struct AnalySol_V8
{
    inline double operator() (const vector2d& pos, double ts) const
    {
        const double xy = pos.x + pos.y;
        return sin(xy*ts) + xy;
    }
};

/*
 * The analytic solution:
 *     u = (x+1)^2(y+1)^2 + t
 */
struct AnalySol_V9
{
    inline double operator() (const vector2d& pos, double ts) const
    {
        return (pos.x+1.)*(pos.x+1.) + ts;
    }
};

/*
 * u = (x+k)^4 + (y+k) + t
 */
struct AnalySol_V10
{
    inline double operator() (const vector2d& pos, double t) const
    {
        const double K = 0.1;
        const vector2d vv = pos + K;
        return vv.x*vv.x*vv.x*vv.x + vv.y + t;
    }
};

struct AnalySol_V11
{
    inline double operator() (const vector2d& pos, double) const
    {
        return fmin(pos.x, fmin(pos.y, fmin(1.-pos.x, 1.-pos.y)));
    }
};

//=======================================================================================
template <typename _Sol>
class AnalyticErr
{
    public:
        AnalyticErr(double lx, double ly, int nx, int ny):
                nx_(nx), ny_(ny), sx_(lx/nx), sy_(ly/ny)
        { }

        inline bool operator()(const boost::multi_array<double, 2>& u, double ts,
                               double& L1, double& Linf) const
        {
            double sum = 0; Linf = 0;
#ifdef USE_OPENMP
            #pragma omp parallel default(none) shared(ts,u,sum,Linf)
            {
            double mySum = 0, myMaxErr = 0;
            #pragma omp for nowait
            for(int iy = 0;iy <= ny_;++ iy)
            for(int ix = 0;ix <= nx_;++ ix)
            {
                vector2d pos(ix*sx_, iy*sy_);
                double ev = fabs(sol_(pos, ts) - u[iy][ix]);
                mySum += ev;
                myMaxErr = fmax(myMaxErr, ev);
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
                vector2d pos(ix*sx_, iy*sy_);
                double ev = fabs(sol_(pos, ts) - u[iy][ix]);
                sum += ev;
                Linf = fmax(Linf, ev);
            }
#endif

            L1 = sum / double((nx_+1)*(ny_+1));
            return true;
        }

    private:
        int     nx_, ny_;
        double  sx_, sy_;   // grid size
        _Sol    sol_;
};

#endif


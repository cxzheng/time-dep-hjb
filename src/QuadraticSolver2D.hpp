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
 *       Filename:  QuadraticSolver2D.hpp
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
 
#ifndef QUADRATIC_SOLVER_2D_HPP
#   define QUADRATIC_SOLVER_2D_HPP

#include <math.h>
#include <stdlib.h>
#include <limits>
#include "QuadraParam.h"
#include "utils/macros.h"
#include "utils/math.hpp"

//! Parameters for time-dependent problems
struct TDQuadraParam2D : public QuadraParam2D
{
    double          ts;
    vector2<int8_t> deriv;
    double          dt;
    double          lastU;
};

// ============================================================================
/*!
 * Quadratic solver for time-dependent HJ PDE, used by Lattice2D
 *
 * The Eikonal equ. we need to solve here is
 * \f[
 * |\nabla u|\tilde{F} = g + \frac{u^n - u^{n+1}}{\Delta t}
 * \f]
 *
 * - _F is the function \f$\tilde{F}(x,t)\f$
 * - _G is the function \f$g(x,t)\f$
 */
template <typename _F, typename _G>
class TDQuadraSolver2D
{
    public:
        /*! Solve the given quadratic equation */
        double two_sided_solve(const TDQuadraParam2D* param);

        /*! Solve for the 1-sided case (linear equation)
         * the upwind value is only given for y-dir
         */
        double one_sided_solve_y(const TDQuadraParam2D* param);
        /*! Solve for the 1-sided case (linear equation)
         * the upwind value is only given for x-dir
         */
        double one_sided_solve_x(const TDQuadraParam2D* param);

    private:
        inline double one_sided_solve_y(const TDQuadraParam2D* param, 
                double fval, double gval, double ubd);
        inline double one_sided_solve_x(const TDQuadraParam2D* param, 
                double fval, double gval, double ubd);

    private:
        _F      f_;
        _G      g_;
};

///////////////////////////////////////////////////////////////////////////////

template <typename _F, typename _G>
double TDQuadraSolver2D<_F,_G>::one_sided_solve_x(const TDQuadraParam2D* param)
{
    const double F_VAL = f_(param->pos, param->ts);
    const double G_VAL = g_(param->pos, param->ts);
    const double U_BOUND = param->lastU - param->dt * G_VAL;    // because g is on the right hand side here

    double A = F_VAL * param->dt - param->H.x;
    if ( M_IS_ZERO(A, 1E-11) ) return U_BOUND;
    double B = G_VAL*param->dt*param->H.x - param->lastU*param->H.x +
               param->upU.x*F_VAL*param->dt;
    double ret = B/A;
    /* Note the solution ret >= upU.x should be satisfied */
    return ret >= param->upU.x ? ret : U_BOUND;
}

template <typename _F, typename _G>
double TDQuadraSolver2D<_F,_G>::one_sided_solve_x(const TDQuadraParam2D* param,
        double F_VAL, double G_VAL, double U_BOUND)
{
    double A = F_VAL * param->dt - param->H.x;
    if ( M_IS_ZERO(A, 1E-11) ) return U_BOUND;
    double B = G_VAL*param->dt*param->H.x - param->lastU*param->H.x +
               param->upU.x*F_VAL*param->dt;
    double ret = B/A;
    return ret >= param->upU.x ? ret : U_BOUND;
}

template <typename _F, typename _G>
double TDQuadraSolver2D<_F,_G>::one_sided_solve_y(const TDQuadraParam2D* param)
{
    const double F_VAL = f_(param->pos, param->ts);
    const double G_VAL = g_(param->pos, param->ts);
    const double U_BOUND = param->lastU - param->dt * G_VAL;    // because g is on the right hand side here

    double A = F_VAL * param->dt - param->H.x;
    if ( M_IS_ZERO(A, 1E-11) ) return U_BOUND;
    double B = G_VAL*param->dt*param->H.x - param->lastU*param->H.x +
               param->upU.y*F_VAL*param->dt;
    double ret = B/A;
    return ret >= param->upU.y ? ret : U_BOUND;
}

template <typename _F, typename _G>
double TDQuadraSolver2D<_F,_G>::one_sided_solve_y(const TDQuadraParam2D* param,
        double F_VAL, double G_VAL, double U_BOUND)
{
    double A = F_VAL * param->dt - param->H.x;
    if ( M_IS_ZERO(A, 1E-11) ) return U_BOUND;
    double B = G_VAL*param->dt*param->H.x - param->lastU*param->H.x +
               param->upU.y*F_VAL*param->dt;
    double ret = B/A;
    return ret >= param->upU.y ? ret : U_BOUND;
}

template <typename _F, typename _G>
double TDQuadraSolver2D<_F,_G>::two_sided_solve(const TDQuadraParam2D* param)
{
    const double F_VAL = f_(param->pos, param->ts);
    const double G_VAL = g_(param->pos, param->ts);
    double t1   = F_VAL*param->dt;
    double fdt2 = 1. / (t1*t1);                                 // positive value
    t1 = G_VAL*param->dt - param->lastU;                        // g*dt - U^{n+1}

    /* this is a upper bound of the U value at current node U(i,j) */
    const double U_BOUND = param->lastU - param->dt * G_VAL;          // because g is on the right hand side here
    /* the A,B,C coefficients for quadratic equation */
    double A = - fdt2;              // negative value
    double B = - 2.*t1*fdt2;
    double C = - t1*t1*fdt2;
    double hinv = 1. / (param->H.x*param->H.x);
    A += hinv*2.;
    B -= 2.*hinv*(param->upU.x + param->upU.y);
    C += (param->upU.x*param->upU.x + param->upU.y*param->upU.y)*hinv;

    if ( A < 0. ) { A = -A; B *= -1.; C *= -1.; }
//    if ( A > 1E-10 )
//    {
        double m = B*B - 4*A*C;
        //MSG_ASSERT(m >= 0., "The quadratic equation should have real solution(m=%.4g)\n", m);
        if ( m >= 0. )
        {
            m = sqrt(m);
            double tu = (-B - m) / (2.*A);
            if ( (param->deriv.x != 0 && tu < param->upU.x) ||
                 (param->deriv.y != 0 && tu < param->upU.y) )
            {
                tu = (-B + m) / (2. * A);
                if ( (param->deriv.x == 0 || tu >= param->upU.x) &&
                     (param->deriv.y == 0 || tu >= param->upU.y) )
                {
                    return tu;
                }
            }
            else
            {
                return tu;
            }
        }
//    }
//    else
//    {
//        fprintf(stderr, "ERROR Happened #2\n"); exit(1);
//    }

    if ( param->upU.x >= param->upU.y )
        return one_sided_solve_y(param, F_VAL, G_VAL, U_BOUND);
    else
        return one_sided_solve_x(param, F_VAL, G_VAL, U_BOUND);
}

#endif



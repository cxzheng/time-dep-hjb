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
 *       Filename:  HJBQuadraSolver2D.hpp
 *
 *    Description:  Quadratic Solver for HJB equation. Same as TDQuadraSolver2D but with
 *                  precomputed f/g value
 *
 *        Version:  1.0
 *        Created:  10/14/11 12:14:23
 *       Revision:  none
 *       Compiler:  gcc/intel compiler
 *
 *         Author:  Changxi Zheng (cz), cxz@cs.columbia.edu
 *                  Columbia University
 *
 * =====================================================================================
 */
#ifndef HJB_QUADRA_SOLVER_2D_INC
#   define HJB_QUADRA_SOLVER_2D_INC

#include <math.h>
#include <limits>
#include "QuadraParam.h"
#include "utils/math.hpp"

struct HJBQuadraSolver2D : public QuadraParam2D
{
    double  dt;
    double  lastU;
    double  fval;
    double  gval;
    double  ubound;

    inline void update_u_bound()
    {   ubound = lastU - dt*gval; }

    inline double one_sided_solve_x() const
    {   
        const double A = fval*dt - H.x;
        if ( M_IS_ZERO(A, 1E-11) ) return ubound;

        const double B = gval*dt*H.x - lastU*H.x + upU.x*fval*dt;
        const double ret = B / A;
        return ret >= upU.x ? ret : ubound;
    }

    inline double one_sided_solve_y() const
    {
        const double A = fval*dt - H.y;
        if ( M_IS_ZERO(A, 1E-11) ) return ubound;

        const double B = gval*dt*H.y - lastU*H.y + upU.y*fval*dt;
        const double ret = B / A;
        return ret >= upU.y ? ret : ubound;
    }

    inline double two_sided_solve() const
    {
        double t1   = fval*dt;
        double fdt2 = 1. / (t1*t1);
        t1 = gval*dt - lastU;                   // g*dt - U^{n+1}

        /* the A,B,C coefficients for quadratic equation */
        double A = - fdt2;              // negative value
        double B = - 2.*t1*fdt2;
        double C = - t1*t1*fdt2;
        double hinv = 1. / (H.x*H.y);
        A += hinv*2.;
        B -= 2.*hinv*(upU.x + upU.y);
        C += (upU.x*upU.x + upU.y*upU.y)*hinv;

        if ( A < 0. ) { A = -A; B = -B; C = -C; }

        double m = B*B - 4*A*C;
        //MSG_ASSERT(m >= 0., "The quadratic equation should have real solution(m=%.4g)\n", m);
        if ( m >= 0. )
        {
            m = sqrt(m);
            double tu = (-B - m) / (2.*A);
            if ( tu >= upU.x && tu >= upU.y ) return tu;

            tu = (-B + m) / (2. * A);
            if ( tu >= upU.x && tu >= upU.y ) return tu;
        }

        return upU.x >= upU.y ? one_sided_solve_y() : one_sided_solve_x();
    }
};


#endif


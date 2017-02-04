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
 *       Filename:  TestCases.h
 *
 *        Version:  1.0
 *        Created:  10/19/11 10:44:29
 *       Revision:  none
 *       Compiler:  gcc/intel compiler
 *
 *         Author:  Changxi Zheng (cz), cxz@cs.columbia.edu
 *                  Columbia University
 *
 * =====================================================================================
 */
#ifndef TEST_CASES_INC
#   define TEST_CASES_INC

#include <math.h>
#include <vector>
#include "vector2.hpp"

namespace Test07
{
    static const double ALPHA = 1;

    struct F
    {
        inline double operator() (const vector2d& pos, double t) const
        {
            const double tx = sin(8.*M_PI*pos.x);
            const double ty = sin(8.*M_PI*pos.y);
            const double tx2 = tx*tx;
            const double tx4 = tx2 * tx2;
            const double tx8 = tx4 * tx4;
            const double ty2 = ty*ty;
            const double ty4 = ty2 * ty2;
            const double ty8 = ty4 * ty4;
            const double tv  = sin(M_PI * t);
            const double t2  = tv * tv;
            //return -ALPHA - 0.8*ALPHA * tx2*tx2*tx * ty2*ty2*ty * sin(M_PI*t);
            //return -0.1*ALPHA - 2.4*ALPHA * tx4*tx4 * ty4*ty4 * t2;
            //return -0.1*ALPHA - 2.4*ALPHA * tx8*tx8 * ty8*ty8 * t2;
            return -0.1*ALPHA - 4.9*ALPHA * tx8*tx8 * ty8*ty8 * t2;
        }
    };
}

namespace Test40
{
    static const double LAMBDA = 8.;
    //static const double S_SCALE = 2048.;    // 2^11
    static const double S_SCALE = 32.;    // 2^5
    static const double EV = exp(LAMBDA);
    static const double EVD = 1. / (EV - 1.);

    struct F
    {
        inline double operator() (const vector2d& pos, double) const
        {
            const double d = fmin(pos.x, fmin(pos.y, fmin(1.-pos.x, 1.-pos.y)));
            const double s = 1. + 2.*d;
            //// ================== exponent = 11 ============
            //const double s2 = s*s;
            //const double s4 = s2*s2;
            //const double s8 = s4*s4;
            ////return -s2*s2*s; // s^5
            //return -(s8*s2*s)*(1./S_SCALE); // s^11/Scale: rescale the speed profile
            // ================== exponent = 5 ============
            const double s2 = s*s;
            const double s4 = s2*s2;
            return -(s4*s)*(1./S_SCALE); // s^11/Scale: rescale the speed profile
        }

        inline double operator() (const vector2d& pos) const
        {
            const double d = fmin(pos.x, fmin(pos.y, fmin(1.-pos.x, 1.-pos.y)));
            const double s = 1. + 2.*d;
            //// ================== exponent = 11 ============
            //const double s2 = s*s;
            //const double s4 = s2*s2;
            //const double s8 = s4*s4;
            ////return -s2*s2*s; // s^5
            //return -(s8*s2*s)*(1./S_SCALE); // s^11/Scale: rescale the speed profile
            // ================== exponent = 5 ============
            const double s2 = s*s;
            const double s4 = s2*s2;
            return -(s4*s)*(1./S_SCALE);
        }
    };

    class BoundaryCond2D
    {
        public:
            BoundaryCond2D(double lx, double ly, int nx, int ny);

            int num_boundary_nodes() const
            {   return totN_; }

            const vector2i* boundary_pos() const
            {   return &bpos_[0]; }

            const double* boundary_val(double);

            /*!
             * time dependent boundary value
             * return u(ix*h, iy*h, ts)
             */
            double boundary_val(int ix, int iy, double ts) const
            {
                const double x = ix*sx_;
                const double y = iy*sy_;
                const double d = fmin(x, fmin(y, fmin(1.-x, 1.-y)));
                
                const double s = (1. + 2.*d);
                // ================== exponent = 11 ============
                //const double s2 = s*s;
                //const double s4 = s2*s2;
                //const double s8 = s4*s4;
                //double dt = (1. - 1. / (s8*s2)) * S_SCALE * 0.05;
                //return dt + EV*EVD*(1. - exp(-LAMBDA*(ts +dt)));
                // double dt = 0.125 * (1. - 1./(s2*s2));
                // ================== exponent = 5 ============
                const double s2 = s*s;
                const double s4 = s2*s2;
                double dt = (1. - 1. / s4) * S_SCALE * 0.125;
                return dt + EV*EVD*(1. - exp(-LAMBDA*(ts +dt)));

                //return dt;
                //return dt + (ts + dt);
                //return dt + exp(LAMBDA * (ts + dt));
                //return dt + exp(LAMBDA * (ts + dt - 1.));
                /*
                 * q(t) = exp(lambda*T) - exp(lambda*(T-t))
                 *        ---------------------------------
                 *             exp(lambda*T) - 1
                 * q(t) is in [0, 1]  
                 */
                //return dt + EV*EVD*(1. - exp(-LAMBDA*(ts +dt)));
            }

        private:
            int     nx_, ny_;
            int     totN_;
            double  sx_, sy_;

            std::vector<vector2i>   bpos_;
            std::vector<double>     bv_;
    };
}

namespace Test41
{
    static const double LAMBDA = 1.;
    static const double EXP_LAMB = exp(LAMBDA);

    /*
     * Solution:
     *
     *                e^lamba - e^lambda(1-t-y)
     * u = y + y^2 + ---------------------------
     *                      e^lambda - 1
     */
    struct G
    {
        inline double operator() (const vector2d& pos, double) const
        {   return -(2.*pos.y + 1); }
    };

    class BoundaryCond2D
    {
        public:
            BoundaryCond2D(double lx, double ly, int nx, int ny);

            int num_boundary_nodes() const
            {   return totN_; }

            const vector2i* boundary_pos() const
            {   return &bpos_[0]; }

            const double* boundary_val(double);         // boundary values

            /*!
             * time dependent boundary value
             * return u(ix*h, iy*h, ts)
             */
            double boundary_val(int, int iy, double ts) const
            {
                const double y = iy*sy_;
                return y + y*y + (EXP_LAMB*(1.-exp(-LAMBDA*(ts + y))))/(EXP_LAMB-1.);
            }

        private:
            int     nx_, ny_;
            int     totN_;
            double  sx_, sy_;

            std::vector<vector2i>   bpos_;
            std::vector<double>     bv_;
    };
}

namespace Test30
{
    static const double ALPHA = M_PI;
    static const double BETA  = 2. * M_PI;
    static const double CVAL  = 1.0;

    struct F
    {
        inline double operator() (const vector2d& pos, double) const
        {   
            double xv = pos.x + 1;
            double yv = pos.y + 2;
            return -BETA / (ALPHA * sqrt(xv*xv + yv*yv));
        }

        inline double operator() (const vector2d& pos) const
        {   
            double xv = pos.x + 1;
            double yv = pos.y + 2;
            return -BETA / (ALPHA * sqrt(xv*xv + yv*yv));
        }
    };

    class BoundaryCond2D
    {
        public:
            BoundaryCond2D(double lx, double ly, int nx, int ny);

            int num_boundary_nodes() const
            {   return totN_; }

            const vector2i* boundary_pos() const
            {   return &bpos_[0]; }

            const double* boundary_val(double);

            /*!
             * time dependent boundary value
             * return u(ix*h, iy*h, ts)
             */
            double boundary_val(int ix, int iy, double ts) const
            {
                double xv = ix*sx_ + 1.;    // x+1
                double yv = iy*sy_ + 2.;    // y+2
                //return sin(ALPHA*xv*yv + BETA*ts) + CVAL;

                double dd = ALPHA*xv*yv + BETA*ts;
                return sin(dd) + CVAL + 2*dd;
            }

        private:
            int     nx_, ny_;
            int     totN_;
            double  sx_, sy_;

            std::vector<vector2i>   bpos_;
            std::vector<double>     bv_;
    };
}

namespace Test31
{
    static const double ALPHA = M_PI;
    static const double BETA  = 2. * M_PI;
    static const double AVAL  = 1.1;
    static const double GAMMA = 0.1;

    struct G
    {
        inline double operator() (const vector2d&, double ts) const
        {   return AVAL*GAMMA*exp(GAMMA*ts); }
    };

    struct F
    {
        inline double operator() (const vector2d& pos, double) const
        {   
            double xv = pos.x + 1;
            double yv = pos.y + 2;
            return -BETA / (ALPHA * sqrt(xv*xv + yv*yv));
        }

        inline double operator() (const vector2d& pos) const
        {   
            double xv = pos.x + 1;
            double yv = pos.y + 2;
            return -BETA / (ALPHA * sqrt(xv*xv + yv*yv));
        }
    };

    class BoundaryCond2D
    {
        public:
            BoundaryCond2D(double lx, double ly, int nx, int ny);

            int num_boundary_nodes() const
            {   return totN_; }

            const vector2i* boundary_pos() const
            {   return &bpos_[0]; }

            const double* boundary_val(double);

            /*!
             * time dependent boundary value
             * return u(ix*h, iy*h, ts)
             */
            double boundary_val(int ix, int iy, double ts) const
            {
                double xv = ix*sx_ + 1.;    // x+1
                double yv = iy*sy_ + 2.;    // y+2
                return sin(ALPHA*xv*yv + BETA*ts) + AVAL*exp(GAMMA*ts);
            }

        private:
            int     nx_, ny_;
            int     totN_;
            double  sx_, sy_;

            std::vector<vector2i>   bpos_;
            std::vector<double>     bv_;
    };
}
#endif

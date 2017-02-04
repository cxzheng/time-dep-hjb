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
 *       Filename:  TestCases.cpp
 *
 *        Version:  1.0
 *        Created:  10/19/11 10:49:24
 *       Revision:  none
 *       Compiler:  gcc/intel compiler
 *
 *         Author:  Changxi Zheng (cz), cxz@cs.columbia.edu
 *                  Columbia University
 *
 * =====================================================================================
 */
#include "TestCases.h"
#include <string.h>

namespace Test30
{
    BoundaryCond2D::BoundaryCond2D(double lx, double ly, int nx, int ny):
        nx_(nx), ny_(ny), totN_((nx+ny)*2), sx_(lx/nx), sy_(ly/ny)
    {
        bpos_.resize(totN_);
        bv_.resize(totN_);
        memset(&bv_[0], 0, sizeof(double)*totN_);

        int cnt = 0;
        for(int i = 0;i <= nx;++ i)
        {
            bpos_[cnt ++].set(i, 0);
            bpos_[cnt ++].set(i, ny);
        }
        for(int i = 1;i < ny;++ i)
        {
            bpos_[cnt ++].set(0, i);
            bpos_[cnt ++].set(nx, i);
        }
    }

    const double* BoundaryCond2D::boundary_val(double ts)
    {
        for(int i = 0;i < totN_;++ i)
            bv_[i] = boundary_val(bpos_[i].x, bpos_[i].y, ts);
        return &bv_[0];
    }
}

namespace Test31
{
    BoundaryCond2D::BoundaryCond2D(double lx, double ly, int nx, int ny):
        nx_(nx), ny_(ny), totN_((nx+ny)*2), sx_(lx/nx), sy_(ly/ny)
    {
        bpos_.resize(totN_);
        bv_.resize(totN_);
        memset(&bv_[0], 0, sizeof(double)*totN_);

        int cnt = 0;
        for(int i = 0;i <= nx;++ i)
        {
            bpos_[cnt ++].set(i, 0);
            bpos_[cnt ++].set(i, ny);
        }
        for(int i = 1;i < ny;++ i)
        {
            bpos_[cnt ++].set(0, i);
            bpos_[cnt ++].set(nx, i);
        }
    }

    const double* BoundaryCond2D::boundary_val(double ts)
    {
        for(int i = 0;i < totN_;++ i)
            bv_[i] = boundary_val(bpos_[i].x, bpos_[i].y, ts);
        return &bv_[0];
    }
}

namespace Test41
{
    BoundaryCond2D::BoundaryCond2D(double lx, double ly, int nx, int ny):
        nx_(nx), ny_(ny), totN_(nx+1), sx_(lx/nx), sy_(ly/ny)
    {
        bpos_.resize(totN_);
        bv_.resize(totN_);
        memset(&bv_[0], 0, sizeof(double)*totN_);

        for(int i = 0;i <= nx;++ i) bpos_[i].set(i, 0);
    }

    const double* BoundaryCond2D::boundary_val(double ts)
    {
        const double v = (EXP_LAMB * (1. - exp(-LAMBDA*ts))) / (EXP_LAMB - 1.);
        for(int i = 0;i < totN_;++i)
            bv_[i] = v;
        return &bv_[0];
    }
}

namespace Test40
{
    BoundaryCond2D::BoundaryCond2D(double lx, double ly, int nx, int ny):
        nx_(nx), ny_(ny), totN_((nx+ny)*2), sx_(lx/nx), sy_(ly/ny)
    {
        bpos_.resize(totN_);
        bv_.resize(totN_);
        memset(&bv_[0], 0, sizeof(double)*totN_);
        int cnt = 0;
        for(int i = 0;i <= nx;++ i)
        {
            bpos_[cnt ++].set(i, 0);
            bpos_[cnt ++].set(i, ny);
        }
        for(int i = 1;i < ny;++ i)
        {
            bpos_[cnt ++].set(0, i);
            bpos_[cnt ++].set(nx, i);
        }
    }

    const double* BoundaryCond2D::boundary_val(double ts)
    {
        const double v = EV*EVD*(1.-exp(-LAMBDA*ts)); 
        for(int i = 0;i < totN_;++ i) bv_[i] = v;
        return &bv_[0];
    }
}

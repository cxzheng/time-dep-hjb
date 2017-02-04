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
 *       Filename:  BoundaryCond.cpp
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
 
#include "BoundaryCond.h"
#include <string.h>

// --------------------- SimpleBoundaryCond2D -------------------------
SimpleBoundaryCond2D::SimpleBoundaryCond2D(int nx, int ny):
        nx_(nx), ny_(ny), totN_((nx+ny)*2)
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

// --------------------- SimpleBoundaryCond2D_01 -------------------------
SimpleBoundaryCond2D_01::SimpleBoundaryCond2D_01(int nx, int ny):
        nx_(nx), ny_(ny), totN_((nx+ny)*2)
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

    for(int i = 0;i < cnt;++ i)
        if ( bpos_[i].x > 0 ) bv_[i] = 1E+4;
}

// --------------------- SimpleBoundaryCond2D_V1 -------------------------
SimpleBoundaryCond2D_V1::SimpleBoundaryCond2D_V1(double lx, double ly, int nx, int ny):
        nx_(nx), ny_(ny),totN_((nx_+1)*2 + (ny_-1)*2), lx_(lx), ly_(ly), sx_(lx/nx), sy_(ly/ny)
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

// --------------------- SimpleBoundaryCond2D_V2 -------------------------
SimpleBoundaryCond2D_V2::SimpleBoundaryCond2D_V2(int nx, int ny):
        nx_(nx), ny_(ny),totN_((nx_+1)*2 + (ny_-1)*2)
{
    bpos_.resize(totN_);
    bv_.resize(totN_);

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

const double* SimpleBoundaryCond2D_V2::boundary_val(double ts)
{
#ifdef USE_OPENMP
    #pragma omp parallel for default(none) shared(ts)
#endif
    for(int i = 0;i < totN_;++ i)
        bv_[i] = ts;
    return &bv_[0];
}

// --------------------- SimpleBoundaryCond2D_V3 -------------------------
SimpleBoundaryCond2D_V3::SimpleBoundaryCond2D_V3(double lx, double ly, int nx, int ny):
        nx_(nx), ny_(ny),totN_((nx_+1)*2 + (ny_-1)*2), lx_(lx), ly_(ly), sx_(lx/nx), sy_(ly/ny)
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

const double* SimpleBoundaryCond2D_V3::boundary_val(double ts)
{
#ifdef USE_OPENMP
    #pragma omp parallel for default(none) shared(ts)
#endif
    for(int i = 0;i < totN_;++ i)
        bv_[i] = (bpos_[i].x*sx_+1.)*(bpos_[i].y*sy_+1.) + ts;
    return &bv_[0];
}


// --------------------- SimpleBoundaryCond2D_V4 -------------------------
SimpleBoundaryCond2D_V4::SimpleBoundaryCond2D_V4(int nx, int ny):nx_(nx), ny_(ny)
{
    bpos_.set(nx_/2, ny_/2);
    bv_ = 0.;
}

// --------------------- SimpleBoundaryCond2D_V5 -------------------------
const double SimpleBoundaryCond2D_V5::K = 1.;

SimpleBoundaryCond2D_V5::SimpleBoundaryCond2D_V5(double lx, double ly, int nx, int ny):
        nx_(nx), ny_(ny),totN_((nx+ny)*2), lx_(lx), ly_(ly), sx_(lx/nx), sy_(ly/ny)
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

const double* SimpleBoundaryCond2D_V5::boundary_val(double ts)
{
#ifdef USE_OPENMP
    #pragma omp parallel for default(none) shared(ts)
#endif
    for(int i = 0;i < totN_;++ i)
    {
        const double tx = bpos_[i].x*sx_+K;
        const double ty = bpos_[i].y*sy_+K;
        bv_[i] = tx*tx*ty*ty + ts;
    }
    return &bv_[0];
}

// --------------------- SimpleBoundaryCond2D_V50 -------------------------
const double SimpleBoundaryCond2D_V50::K = 1.;

SimpleBoundaryCond2D_V50::SimpleBoundaryCond2D_V50(double lx, double ly, int nx, int ny):
        nx_(nx), ny_(ny),totN_((nx_+1)*2 + (ny_-1)*2), lx_(lx), ly_(ly), sx_(lx/nx), sy_(ly/ny)
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

const double* SimpleBoundaryCond2D_V50::boundary_val(double ts)
{
#ifdef USE_OPENMP
    #pragma omp parallel for default(none) shared(ts)
#endif
    for(int i = 0;i < totN_;++ i)
    {
        const double tx = bpos_[i].x*sx_+K;
        const double ty = bpos_[i].y*sy_+K;
        bv_[i] = tx*tx*ty*ty + ts;
    }
    return &bv_[0];
}

// --------------------- SimpleBoundaryCond2D_V6 -------------------------
SimpleBoundaryCond2D_V6::SimpleBoundaryCond2D_V6(
        double beta, double lx, double ly, int nx, int ny):
        beta_(beta), nx_(nx), ny_(ny), totN_(nx_*2 + ny_*2), 
        lx_(lx), ly_(ly), sx_(lx/nx), sy_(ly/ny)
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

const double* SimpleBoundaryCond2D_V6::boundary_val(double ts)
{
#ifdef USE_OPENMP
    #pragma omp parallel for default(none) shared(ts)
#endif
    for(int i = 0;i < totN_;++ i)
    {
        const double tx = bpos_[i].x*sx_+1.;
        const double ty = bpos_[i].y*sy_+1.;
        bv_[i] = tx*tx*ty*ty*exp(beta_*ts);
    }
    return &bv_[0];
}

// --------------------- SimpleBoundaryCond2D_V7 -------------------------
SimpleBoundaryCond2D_V7::SimpleBoundaryCond2D_V7(
        double lx, double ly, int nx, int ny): nx_(nx), ny_(ny),
        totN_(nx_+ny_+1), lx_(lx), ly_(ly), sx_(lx/nx), sy_(ly/ny)
{
    bpos_.resize(totN_);
    bv_.resize(totN_);
    memset(&bv_[0], 0, sizeof(double)*totN_);

    int cnt = 0;
    for(int i = 0;i <= nx;++ i) bpos_[cnt ++].set(i, 0);
    for(int i = 1;i <= ny;++ i) bpos_[cnt ++].set(0, i);
}

const double* SimpleBoundaryCond2D_V7::boundary_val(double ts)
{
#ifdef USE_OPENMP
    #pragma omp parallel for default(none) shared(ts)
#endif
    for(int i = 0;i < totN_;++ i)
    {
        const double tx = bpos_[i].x*sx_+1.;
        const double ty = bpos_[i].y*sy_+1.;
        bv_[i] = tx*tx*ty*ty*exp(ts);
    }
    return &bv_[0];
}

// --------------------- SimpleBoundaryCond2D_V8 -------------------------
SimpleBoundaryCond2D_V8::SimpleBoundaryCond2D_V8(
        double tmax, double lx, double ly, int nx, int ny):
        tmax_(tmax), nx_(nx), ny_(ny), totN_(nx_*2 + ny_*2), 
        lx_(lx), ly_(ly), sx_(lx/nx), sy_(ly/ny)
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

const double* SimpleBoundaryCond2D_V8::boundary_val(double ts)
{
#ifdef USE_OPENMP
    #pragma omp parallel for default(none) shared(ts)
#endif
    for(int i = 0;i < totN_;++ i)
    {
        const double tx = bpos_[i].x*sx_+1.;
        const double ty = bpos_[i].y*sy_+1.;
        bv_[i] = tx*tx*ty*ty*exp(tmax_ - ts);
    }
    return &bv_[0];
}

// --------------------- SimpleBoundaryCond2D_V9 -------------------------
SimpleBoundaryCond2D_V9::SimpleBoundaryCond2D_V9(double lx, double ly, int nx, int ny):
        nx_(nx), ny_(ny),totN_((nx_+1)*2 + (ny_-1)*2), lx_(lx), ly_(ly), sx_(lx/nx), sy_(ly/ny)
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

const double* SimpleBoundaryCond2D_V9::boundary_val(double ts)
{
#ifdef USE_OPENMP
    #pragma omp parallel for default(none) shared(ts)
#endif
    for(int i = 0;i < totN_;++ i)
        bv_[i] = boundary_val(bpos_[i].x, bpos_[i].y, ts);

    return &bv_[0];
}

// --------------------- SimpleBoundaryCond2D_V10 -------------------------
SimpleBoundaryCond2D_V10::SimpleBoundaryCond2D_V10(double lx, double ly, int nx, int ny):
        nx_(nx), ny_(ny), totN_((ny_+1)*2), lx_(lx), ly_(ly), sx_(lx/nx), sy_(ly/ny)
{
    bpos_.resize(totN_);
    bv_.resize(totN_);
    memset(&bv_[0], 0, sizeof(double)*totN_);

    int cnt = 0;
    for(int i = 0;i <= ny;++ i)
    {
        bpos_[cnt ++].set(0, i);
        bpos_[cnt ++].set(nx, i);
    }
}

const double* SimpleBoundaryCond2D_V10::boundary_val(double ts)
{
#ifdef USE_OPENMP
    #pragma omp parallel for default(none) shared(ts)
#endif
    for(int i = 0;i < totN_;++ i)
        bv_[i] = boundary_val(bpos_[i].x, bpos_[i].y, ts);
    return &bv_[0];
}

// --------------------- SimpleBoundaryCond2D_V11 -------------------------
SimpleBoundaryCond2D_V11::SimpleBoundaryCond2D_V11(double lx, double ly, int nx, int ny):
        nx_(nx), ny_(ny), totN_(ny+nx+1), lx_(lx), ly_(ly), sx_(lx/nx), sy_(ly/ny)
{
    bpos_.resize(totN_);
    bv_.resize(totN_);
    memset(&bv_[0], 0, sizeof(double)*totN_);

    int cnt = 0;
    for(int i = 0;i <= ny;++ i)
        bpos_[cnt ++].set(0, i);
    for(int i = 1;i <= nx;++ i)
        bpos_[cnt ++].set(i, 0);
}

const double* SimpleBoundaryCond2D_V11::boundary_val(double ts)
{
#ifdef USE_OPENMP
    #pragma omp parallel for default(none) shared(ts)
#endif
    for(int i = 0;i < totN_;++ i)
        bv_[i] = boundary_val(bpos_[i].x, bpos_[i].y, ts);
    return &bv_[0];
}





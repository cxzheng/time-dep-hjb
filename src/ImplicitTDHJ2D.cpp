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
 *       Filename:  ImplicitTDHJ2D.cpp
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
 
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include "config.h"
#include "HJ2DFunc.h"
#include "BoundaryCond.h"
#include "BoundCondTemp.hpp"
#include "TimeDepHJ2DEqu.hpp"
#include "io/npy_arrays.hpp"
#include "CachedImpLattice2D.hpp"
#include "TestCases.h"

#include "utils/nano_timer.h"

using namespace std;

#ifdef USE_TEST_00
#warning ====== Use Test Case #00 ======
// --------00-- Simple F=1, G=1, time-dependent HJ equ ---------------
/* ** exp-1 **
 * Equation: u_t - |\nabla u| = -1
 * BC: u(x,t) = 0 where x \in \partial\Omega or t = T
 */
typedef SimpleBoundaryCond2D                        TBC;
typedef CachedImpLattice2D<func_F1, func_G1, TBC>   TLattice;

// ------------ Simple F=1, G=1, time-dependent HJ equ ---------------
/*
 * Equation: u_t - |\nabla u| = -1
 * BC: u(x,t) = 0 where x \in \partial\Omega
 *     u(x,T) = \phi(x) where \phi is the distance value to the boundary
 */
//typedef SimpleBoundaryCond2D_V1               TBC;
//typedef TDQuadraSolver2D<func_F1, func_G1>    TSolver;

// ------------ Simple F=1, G=0, time-dependent HJ equ ---------------
/*
 * Equation: u_t - |\nabla u| = 0
 * BC:  u(x,t) = t
 */
//typedef SimpleBoundaryCond2D_V2               TBC;
//typedef TDQuadraSolver2D<func_F1, func_G0>    TSolver;

// -------03-- Simple F=-1/|x+(1,1)|, G=0, time-dependent HJ equ ------------
/*
 * Equation: u_t - |\nabla u| /|x+(1,1)| = 0
 * the analytic solution: u(x,t) = (x+1)(y+1) + t
 */
//typedef SimpleBoundaryCond2D_V3                 TBC;
//typedef TDQuadraSolver2D<func_F_V3, func_G0>    TSolver;

// -------04----------------------------------------------------------
/*
 * Equation: u_t - |\nabla u|F = -1
 *           where F = 1 + 0.8*sin(8\pi x)sin(8\pi y)
 * BC: u(x,t) = 0 where x \in \partial\Omega or t = T
 */
//typedef SimpleBoundaryCond2D               TBC;
//typedef TDQuadraSolver2D< func_F_V4, func_G1 >  TSolver;

#elif defined(USE_TEST_05)
#warning ====== Use Test Case #05 ======
// -------05----------------------------------------------------------
/*
 * Equation: u_t - |\nabla u|F = -1
 *           where F = 1 + 0.8*sin(8\pi x)sin(8\pi y)
 * BC: u(x,t) = 0 where t = T
 *     u(x,t) = 0 where x = (lx/2, ly/2)
 */
typedef SimpleBoundaryCond2D_V4                     TBC;
typedef CachedImpLattice2D<func_F_V4, func_G1, TBC> TLattice;

#elif defined(USE_TEST_07)
#warning ====== Use Test Case #07 ======
// -------06,07,08----------------------------------------------------
/*
 * Equation: u_t - |\nabla u|F = -1
 *           where F = 1 + 0.8*sin(8\pi x)sin(8\pi y)*sin(8\pi t)
 * BC: u(x,t) = 0 where x \in \partial\Omega or t = T
 */
typedef SimpleBoundaryCond2D                        TBC;
typedef CachedImpLattice2D<Test07::F, func_G1, TBC> TLattice;

// -------09----------------------------------------------------------
/*
 * u_t - |\nabla u|F = 0
 *       where F = 1/2(x+1)(y+1)\sqrt{(x+1)^2+(y+1)^2}
 * the analytic solution: u(x,t) = (x+1)^2(y+1)^2 + t
 */
//typedef SimpleBoundaryCond2D_V5              TBC;
//typedef TDQuadraSolver2D<func_F_V6, func_G0> TSolver;

// -------10----------------------------------------------------------
/*
 * u_t - |\nabla u|F = -1
 *       where F = 1/(x+1)(y+1)\sqrt{(x+1)^2+(y+1)^2}
 * the analytic solution: u(x,t) = (x+1)^2(y+1)^2 + t
 */
//typedef SimpleBoundaryCond2D_V5              TBC;
//typedef TDQuadraSolver2D<func_F_V7, func_G1> TSolver;

// -------11----------------------------------------------------------
/*
 * u_t - |\nabla u|F = 0
 *       where F = (x+1)(y+1)/(2*\sqrt{(x+1)^2+(y+1)^2})
 * the analytic solution: u(x,t) = (x+1)^2(y+1)^2 exp(t)
 */
//typedef SimpleBoundaryCond2D_V6              TBC;
//typedef SimpleBoundaryCond2D_V7              TBC;
//typedef TDQuadraSolver2D<func_F_V8, func_G0> TSolver;

// -------12----------------------------------------------------------
/*
 * u_t - |\nabla u|F = -(x+1)^2(y+1)^2 exp(t) 
 *       where F = (x+1)(y+1)/(\sqrt{(x+1)^2+(y+1)^2})
 * the analytic solution: u(x,t) = (x+1)^2(y+1)^2 exp(t)
 */
//typedef SimpleBoundaryCond2D_V6                TBC;
//typedef TDQuadraSolver2D<func_F_V9, func_G_V1> TSolver;

// -------13----------------------------------------------------------
/*
 * u_t - |\nabla u|F = -2(x+1)^2(y+1)^2 exp(T_max-t) 
 *       where F = (x+1)(y+1)/(2.\sqrt{(x+1)^2+(y+1)^2})
 * the analytic solution: u(x,t) = (x+1)^2(y+1)^2 exp(T_max-t)
 */
//typedef SimpleBoundaryCond2D_V8                TBC;
//typedef TDQuadraSolver2D<func_F_V8, func_G_V2> TSolver;

#elif defined(USE_TEST_20)
#warning ====== Use Test Case #20 ======
// -------20----------------------------------------------------------
/*
 * u_t - |\nabla u| = -1 
 * the analytic solution: u(x,y,t) = y + exp(beta*(t + y))
 */
typedef BoundCond2D_1Side<func_exp>             TBC;
typedef TDQuadraSolver2D<func_F1, func_G1>      TSolver;

#elif defined(USE_TEST_21)
#warning ====== Use Test Case #21 ======
// -------21----------------------------------------------------------
/*
 * u_t - |\nabla u| = -(2y+1)
 * the analytic solution: u(x,y,t) = y + y^2 + exp(beta*(t + y))
 * boundary is defined as [0,1]X{0}
 */
typedef BoundCond2D_1Side_V2<func_exp>              TBC;
typedef CachedImpLattice2D<func_F1, func_G_V3, TBC> TLattice;

#elif defined(USE_TEST_22)
#warning ====== Use Test Case #22 ======
// -------22----------------------------------------------------------
/*
 * u_t - \frac{|\nabla u|}{2y+1} = -1 
 * the analytic solution: u(x,y,t) = y + y^2 + exp(beta*(t + y + y^2))
 * boundary is defined as [0,1]X{0}
 */
typedef BoundCond2D_1Side_V3<func_exp>                  TBC;
typedef CachedImpLattice2D<func_F_V10, func_G1, TBC>    TLattice;

#elif defined(USE_TEST_30)
#warning ====== Use Test Case #30 ======
typedef Test30::BoundaryCond2D                      TBC;
typedef CachedImpLattice2D<Test30::F, func_G0, TBC> TLattice;
//typedef Test30_1::BoundaryCond2D                      TBC;
//typedef CachedImpLattice2D<Test30_1::F, func_G0, TBC> TLattice;

#elif defined(USE_TEST_31)
#warning ====== Use Test Case #31 ======
typedef Test31::BoundaryCond2D                          TBC;
typedef CachedImpLattice2D<Test31::F, Test31::G, TBC>   TLattice;

#elif defined(USE_TEST_40)
#warning ====== Use Test Case #40 ======
typedef Test40::BoundaryCond2D                          TBC;
typedef CachedImpLattice2D<Test40::F, func_G1, TBC>     TLattice;

#elif defined(USE_TEST_41)
#warning ====== Use Test Case #41 ======
typedef Test41::BoundaryCond2D                          TBC;
typedef CachedImpLattice2D<func_F1, Test41::G, TBC>     TLattice;

#else
#warning ====== Use Test Case #00 ======
// --------00-- Simple F=1, G=1, time-dependent HJ equ ---------------
/* ** exp-1 **
 * Equation: u_t - |\nabla u| = -1
 * BC: u(x,t) = 0 where x \in \partial\Omega or t = T
 */
typedef SimpleBoundaryCond2D                        TBC;
typedef CachedImpLattice2D<func_F1, func_G1, TBC>   TLattice;
#endif

// -------------------------------------------------------------------
typedef TimeDepHJ2DEqu<TLattice>              TTimeDepHJEqu;

static TBC*             pbc = NULL;
static TLattice*        plattice = NULL;
static TTimeDepHJEqu*   equ = NULL;

static double lx = 1, ly = 1;
static double T = 2., dt = -0.05;
static double endT = 0;
static int nx = 100, ny = 100;
static string outFile = "out.dat";

static void clean()
{
    delete equ;
    delete plattice;
    delete pbc;
}

static void usage(const char* cmd)
{
    printf("Usage: %s [-h] [-W lx] [-H ly] [-X nx] [-Y ny] [-T tmax] [-E tend] [-d dt] [-o file pattern]\n", cmd);
    printf("options:\n");
    printf("    -W lx       domain size along x-axis\n");
    printf("    -H ly       domain size along y-axis\n");
    printf("    -X nx       resolution along x-axis\n");
    printf("    -Y nx       resolution along y-axis\n");
    printf("    -d dt       time step size\n");
    printf("    -T tmax     solve equation in time domain [0, T]\n");
    printf("    -E tend     solve the equation in time domain [tend, tmax]\n");
    printf("    -o outfile  file name pattern for output files\n");
    printf("\n");
    printf("%s  Copyright (C) 2013  Changxi Zheng\n", cmd);
    printf("This program comes with ABSOLUTELY NO WARRANTY;\n");
    printf("This is free software, and you are welcome to redistribute it\n");
    printf("under certain conditions; see http://www.gnu.org/copyleft/gpl.html for details.\n");
    printf("\n");
}

static void parse_cmd(int argc, char* argv[])
{
    int opt;
    while ( (opt = getopt(argc, argv, "hW:H:X:Y:T:d:o:E:")) != -1 )
    {
        switch (opt)
        {
            case 'h':
                usage(argv[0]);
                exit(0);
            case 'W':
                lx = atof(optarg);
                break;
            case 'H':
                ly = atof(optarg);
                break;
            case 'X':
                nx = atoi(optarg);
                break;
            case 'Y':
                ny = atoi(optarg);
                break;
            case 'T':
                T  = atof(optarg);
                break;
            case 'E':
                endT = atof(optarg);
                break;
            case 'd':
                dt = -atof(optarg);
                break;
            case 'o':
                outFile = optarg;
                break;
        }
    }
}

static void write_u(int id)
{
    char file[256];
    sprintf(file, outFile.c_str(), id);

    const TLattice::TValField& val = plattice->current_u();
    npy_write_array_double<2>(val, file);
}

int main(int argc, char* argv[])
{
    parse_cmd(argc, argv);
    cout << "==============================================" << endl;

#if defined(USE_TEST_05)
    cout << "  R U N  T E S T *05* " << endl;
    pbc = new SimpleBoundaryCond2D_V4(nx, ny);            // 5
#elif defined(USE_TEST_07)
    cout << "  R U N  T E S T *07* " << endl;
    pbc = new TBC(nx, ny);
#elif defined(USE_TEST_21)
    cout << "  R U N  T E S T *21* " << endl;
    pbc = new TBC(lx, ly, nx, ny);      // 20,21,22
#elif defined(USE_TEST_22)
    cout << "  R U N  T E S T *22* " << endl;
    pbc = new TBC(lx, ly, nx, ny);      // 20,21,22
#elif defined(USE_TEST_30)
    cout << "  R U N  T E S T *30* " << endl;
    pbc = new TBC(lx, ly, nx, ny);      // 30
#elif defined(USE_TEST_31)
    cout << "  R U N  T E S T *31* " << endl;
    pbc = new TBC(lx, ly, nx, ny);      // 31
#elif defined(USE_TEST_40)
    cout << "  R U N  T E S T *40* " << endl;
    pbc = new TBC(lx, ly, nx, ny);      // 40
#elif defined(USE_TEST_41)
    cout << "  R U N  T E S T *41* " << endl;
    pbc = new TBC(lx, ly, nx, ny);      // 41
#else                                   // by default using Test00
    cout << "  R U N  T E S T *00* " << endl;
    pbc = new SimpleBoundaryCond2D(nx, ny);                 // 0,4,6,7,8
#endif
    cout << "----------------------" << endl;
    cout << "  WxH        = " << lx << " x " << ly << endl;
    cout << "  Nx X Ny    = " << nx << " x " << ny << endl;
    cout << "  Tmax       = " << T << endl;
    cout << "  dt         = " <<-dt << endl;
    cout << "  outfile    = " << outFile << endl; 
    cout << "==============================================" << endl;

    plattice = new TLattice(lx, ly, nx, ny, pbc);
    equ      = new TTimeDepHJEqu(plattice, T);
    printf("Run simulation ...\n");
    double remainT = equ->time() - endT;

    double st = GetMilliTimed();
    while ( remainT > 0. )
    {
        if ( remainT <= -dt )
        {
            equ->advance(-remainT);
            printf("\rtime = %-10g", equ->time()); fflush(stdout);
            break;
        }
        equ->advance(dt);

#ifdef COUNT_UPDATE_NUM
        cerr << equ->time() << ' ' << plattice->average_update_cnt() << endl;
#endif
        printf("\rtime = %-10g", equ->time()); fflush(stdout);
        remainT = equ->time() - endT;
    }
    fprintf(stderr, "%lf\n", GetMilliTimed() - st);
    printf("\nDone!\n");
    write_u(0);

    clean();
    return 0;
}



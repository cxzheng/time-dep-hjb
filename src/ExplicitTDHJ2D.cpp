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
 *       Filename:  ExplicitTDHJ2D.cpp
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
#include <iomanip>
#include <fstream>
#include <string>
#include <boost/program_options.hpp>
#include "config.h"
#include "BoundaryCond.h"
#include "BoundCondTemp.hpp"
#include "HJ2DFunc.h"
#include "ExplicitLattice2D.hpp"
#include "UpwindGrad2D.h"
#include "io/npy_arrays.hpp"
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
typedef SimpleBoundaryCond2D                TBC;
typedef ExplicitLattice2D<func_F1, func_G1, UpwindGrad2D, TBC>  TLattice;

// --------01-- Simple F=1, G=1, time-dependent HJ equ ---------------
/*
 * Equation: u_t - |\nabla u| = -1
 * BC: u(x,t) = 0 where x \in \partial\Omega
 *     u(x,T) = \phi(x) where \phi is the distance value to the boundary
 */
//typedef SimpleBoundaryCond2D_V1                                 TBC;
//typedef ExplicitLattice2D<func_F1, func_G1, UpwindGrad2D, TBC>  TLattice;

// --------02-- Simple F=1, G=0, time-dependent HJ equ ---------------
/*
 * Equation: u_t - |\nabla u| = 0
 * BC:  u(x,t) = t
 */
//typedef SimpleBoundaryCond2D_V2               TBC;
//typedef ExplicitLattice2D<func_F1, func_G0, UpwindGrad2D, TBC>  TLattice;

// -------03-- Simple F=-1/|x+(1,1)|, G=0, time-dependent HJ equ ------------
/*
 * Equation: u_t - |\nabla u| /|x+(1,1)| = 0
 * the analytic solution: u(x,t) = (x+1)(y+1) + t
 */
//typedef SimpleBoundaryCond2D_V3               TBC;
//typedef ExplicitLattice2D<func_F_V3, func_G0, UpwindGrad2D, TBC>  TLattice;

// -------04----------------------------------------------------------
/*
 * Equation: u_t - |\nabla u|F = -1
 *           where F = 1 + 0.8*sin(8\pi x)sin(8\pi y)
 * BC: u(x,t) = 0 where x \in \partial\Omega or t = T
 */
//typedef SimpleBoundaryCond2D               TBC;
//typedef ExplicitLattice2D<func_F_V4, func_G1, UpwindGrad2D, TBC>  TLattice;

#elif defined(USE_TEST_05)
#warning ====== Use Test Case #05 ======
// -------05----------------------------------------------------------
/*
 * Equation: u_t - |\nabla u|F = -1
 *           where F = 1 + 0.8*sin(8\pi x)sin(8\pi y)
 * BC: u(x,t) = 0 where t = T
 *     u(x,t) = 0 where x = (lx/2, ly/2)
 */
typedef SimpleBoundaryCond2D_V4               TBC;
typedef ExplicitLattice2D<func_F_V4, func_G1, UpwindGrad2D, TBC>  TLattice;

#elif defined(USE_TEST_07)
#warning ====== Use Test Case #07 ======
// -------06,07,08----------------------------------------------------
/*
 * Equation: u_t - |\nabla u|F = -1
 *           where F = 1 + 0.8*sin(8\pi x)sin(8\pi y)*sin(8\pi t)
 * BC: u(x,t) = 0 where x \in \partial\Omega or t = T
 */
typedef SimpleBoundaryCond2D                  TBC;
typedef ExplicitLattice2D<Test07::F, func_G1, UpwindGrad2D, TBC>  TLattice;

// -------09----------------------------------------------------------
/*
 * u_t - |\nabla u|F = 0
 *       where F = 1/2(x+1)(y+1)\sqrt{(x+1)^2+(y+1)^2}
 * the analytic solution: u(x,t) = (x+1)^2(y+1)^2 + t
 */
//typedef SimpleBoundaryCond2D_V5              TBC;
//typedef ExplicitLattice2D<func_F_V6, func_G0, UpwindGrad2D, TBC>  TLattice;

// -------10----------------------------------------------------------
/*
 * u_t - |\nabla u|F = -1
 *       where F = 1/(x+1)(y+1)\sqrt{(x+1)^2+(y+1)^2}
 * the analytic solution: u(x,t) = (x+1)^2(y+1)^2 + t
 */
//typedef SimpleBoundaryCond2D_V5              TBC;
//typedef ExplicitLattice2D<func_F_V7, func_G1, UpwindGrad2D, TBC>  TLattice;

// -------11----------------------------------------------------------
/*
 * u_t - |\nabla u|F = 0
 *       where F = (x+1)(y+1)/(2*\sqrt{(x+1)^2+(y+1)^2})
 * the analytic solution: u(x,t) = (x+1)^2(y+1)^2 exp(t)
 */
//typedef SimpleBoundaryCond2D_V6              TBC;
//typedef SimpleBoundaryCond2D_V7              TBC;
//typedef ExplicitLattice2D<func_F_V8, func_G0, UpwindGrad2D, TBC>  TLattice;

// -------12----------------------------------------------------------
/*
 * u_t - |\nabla u|F = -(x+1)^2(y+1)^2 exp(t) 
 *       where F = (x+1)(y+1)/(\sqrt{(x+1)^2+(y+1)^2})
 * the analytic solution: u(x,t) = (x+1)^2(y+1)^2 exp(t)
 */
//typedef SimpleBoundaryCond2D_V6              TBC;
//typedef ExplicitLattice2D<func_F_V9, func_G_V1, UpwindGrad2D, TBC>  TLattice;

// -------13----------------------------------------------------------
/*
 * u_t - |\nabla u|F = -2(x+1)^2(y+1)^2 exp(T_max-t) 
 *       where F = (x+1)(y+1)/(2.\sqrt{(x+1)^2+(y+1)^2})
 * the analytic solution: u(x,t) = (x+1)^2(y+1)^2 exp(T_max-t)
 */
//typedef SimpleBoundaryCond2D_V8                TBC;
//typedef ExplicitLattice2D<func_F_V8, func_G_V2, UpwindGrad2D, TBC>  TLattice;

// -------20----------------------------------------------------------
/*
 * u_t - |\nabla u| = -1 
 * the analytic solution: u(x,y,t) = y + exp(beta*(t + y))
 */
//typedef BoundCond2D_1Side<func_exp>             TBC;
//typedef ExplicitLattice2D<func_F1, func_G1, UpwindGrad2D, TBC> TLattice;

#elif defined(USE_TEST_21)
#warning ====== Use Test Case #21 ======
// -------21----------------------------------------------------------
/*
 * u_t - |\nabla u| = -(2y+1)
 * the analytic solution: u(x,y,t) = y + y^2 + exp(beta*(t + y))
 * boundary is defined as [0,1]X{0}
 */
typedef BoundCond2D_1Side_V2<func_exp>          TBC;
typedef ExplicitLattice2D<func_F1, func_G_V3, UpwindGrad2D, TBC> TLattice;

#elif defined(USE_TEST_22)
#warning ====== Use Test Case #22 ======
// -------22----------------------------------------------------------
/*
 * u_t - \frac{|\nabla u|}{2y+1} = -1 
 * the analytic solution: u(x,y,t) = y + y^2 + exp(beta*(t + y + y^2))
 * boundary is defined as [0,1]X{0}
 */
typedef BoundCond2D_1Side_V3<func_exp>         TBC;
typedef ExplicitLattice2D<func_F_V10, func_G1, UpwindGrad2D, TBC> TLattice;

#elif defined(USE_TEST_23)
#warning ====== Use Test Case #23 ======
// -------23----------------------------------------------------------
/*
 * Equation: u_t - |\nabla u|F = -1
 *           where F = 1 + 0.8*sin(8\pi x)sin(8\pi y)
 * BC: u(x,t) = 0 where t = T
 *     u(x,t) = sin(8*pi*t) where x = (lx/2, ly/2)
 */
typedef SimpleBoundaryCond2D_V42               TBC;
typedef ExplicitLattice2D<func_F_V4, func_G1, UpwindGrad2D, TBC>  TLattice;

#elif defined(USE_TEST_30)
#warning ====== Use Test Case #30 ======
typedef Test30::BoundaryCond2D                                      TBC;
typedef ExplicitLattice2D<Test30::F, func_G0, UpwindGrad2D, TBC>    TLattice;

#elif defined(USE_TEST_40)
#warning ====== Use Test Case #40 ======
typedef Test40::BoundaryCond2D                                    TBC;
typedef ExplicitLattice2D<Test40::F, func_G1, UpwindGrad2D, TBC>  TLattice;

#elif defined(USE_TEST_41)
#warning ====== Use Test Case #41 ======
typedef Test41::BoundaryCond2D                                    TBC;
typedef ExplicitLattice2D<func_F1, Test41::G, UpwindGrad2D, TBC>  TLattice;

#else
#warning ====== Use Test Case #00 ======
// --------00-- Simple F=1, G=1, time-dependent HJ equ ---------------
/* ** exp-1 **
 * Equation: u_t - |\nabla u| = -1
 * BC: u(x,t) = 0 where x \in \partial\Omega or t = T
 */
typedef SimpleBoundaryCond2D                TBC;
typedef ExplicitLattice2D<func_F1, func_G1, UpwindGrad2D, TBC>  TLattice;
#endif

// -------------------------------------------------------------------
static TBC*             pbc = NULL;
static UpwindGrad2D*    pgd = NULL;
static TLattice*        plattice = NULL;

static double lx = 1, ly = 1;
static int nx = 100, ny = 100;
static double T = 2., dt = 0.05;
static double endT = 0.;
static string outFile;

static int dumpF = 0;

static void parse_cmd(int argc, char* argv[])
{
    namespace po = boost::program_options;
    po::options_description genericOpt("Generic options");
    genericOpt.add_options()
            ("help,h", "display help information");
    po::options_description configOpt("Configuration");
    configOpt.add_options()
            ("width,W", po::value<double>(&lx), "domain size along x-axis")
            ("height,H", po::value<double>(&ly), "domain size along y-axis")
            ("nx,X", po::value<int>(&nx), "resolution along x-axis")
            ("ny,Y", po::value<int>(&ny), "resolution along y-axis")
            ("dt,d", po::value<double>(&dt), "time step size")
            ("tmax,T", po::value<double>(&T), "solve equation in time domain [0, T]")
            ("tend,E", po::value<double>(&endT), "solver ending time")
            ("out,o", po::value<string>(&outFile)->default_value("out-%d.txt"),
                    "file name patter for output data dumps")
            ("dump,D", po::value<int>(&dumpF), 
                    "dump the data at each time step (This is for error estimation)");
    // use configure file to specify the option
    po::options_description cfileOpt("Configure file");
    cfileOpt.add_options()
            ("cfg-file", po::value<string>(), "configuration file");

    po::options_description cmdOpts;
    cmdOpts.add(genericOpt).add(configOpt).add(cfileOpt);

    po::variables_map vm;
    store(po::parse_command_line(argc, argv, cmdOpts), vm);
    if ( vm.count("cfg-file") )
    {
        ifstream ifs(vm["cfg-file"].as<string>().c_str());
        store(parse_config_file(ifs, configOpt), vm);
    }
    po::notify(vm);

    if ( vm.count("help") )
    {
        printf("Usage: %s [options] \n", argv[0]);
        cout << cmdOpts << endl << endl;
        printf("%s  Copyright (C) 2013  Changxi Zheng\n", argv[0]);
        printf("This program comes with ABSOLUTELY NO WARRANTY;\n");
        printf("This is free software, and you are welcome to redistribute it\n");
        printf("under certain conditions; see http://www.gnu.org/copyleft/gpl.html for details.\n");
        exit(0);
    }
}

static void clean()
{
    delete pbc;
    delete pgd;
    delete plattice;
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
    dt *= -1;

    cout << "==============================================" << endl;
    //pbc = new SimpleBoundaryCond2D_V1(lx, ly, nx, ny);
#ifdef USE_TEST_00
    cout << "  R U N  T E S T *00* " << endl;
    pbc = new SimpleBoundaryCond2D(nx, ny);     // 0,4
#elif defined(USE_TEST_05)
    cout << "  R U N  T E S T *05* " << endl;
    pbc = new SimpleBoundaryCond2D_V4(nx, ny);  // 5
#elif defined(USE_TEST_07)
    cout << "  R U N  T E S T *07* " << endl;
    pbc = new TBC(nx, ny);  // 07
#elif defined(USE_TEST_23)
    cout << "  R U N  T E S T *23* " << endl;
    pbc = new TBC(nx, ny);  // 23
#elif defined(USE_TEST_21)
    cout << "  R U N  T E S T *21* " << endl;
    pbc = new TBC(lx, ly, nx, ny);      // 20,21,22
#elif defined(USE_TEST_22)
    cout << "  R U N  T E S T *22* " << endl;
    pbc = new TBC(lx, ly, nx, ny);      // 20,21,22
#elif defined(USE_TEST_30)
    cout << "  R U N  T E S T *30* " << endl;
    pbc = new TBC(lx, ly, nx, ny);      // 20,21,22
#elif defined(USE_TEST_40)
    cout << "  R U N  T E S T *40* " << endl;
    pbc = new TBC(lx, ly, nx, ny);      // 40
#elif defined(USE_TEST_41)
    cout << "  R U N  T E S T *41* " << endl;
    pbc = new TBC(lx, ly, nx, ny);      // 41
#else
    cout << "  R U N  T E S T *00* " << endl;
    pbc = new SimpleBoundaryCond2D(nx, ny);     // 0,4
#endif
    cout << "----------------------" << endl;
    cout << "  WxH        = " << lx << " x " << ly << endl;
    cout << "  Nx X Ny    = " << nx << " x " << ny << endl;
    cout << "  Tmax       = " << T << endl;
    cout << "  dt         = " <<-dt << endl;
    cout << "  outfile    = " << outFile << endl; 
    if ( dumpF > 0 ) cout  << "  Dump DATA = " << dumpF << endl;
#ifdef USE_OPENMP
    cout << "   OPENMP:      enabled" << endl;
#endif
    cout << "==============================================" << endl;

    pgd = new UpwindGrad2D(lx, ly, nx, ny);
    plattice = new TLattice(lx, ly, nx, ny, pgd, pbc);
    plattice->init(T);

    printf("Run simulation ...\n");
    double remainT = plattice->time() - endT;

    double st = GetMilliTimed();
    while ( remainT > 0. )
    {
        if ( remainT <= -dt )
        {
            plattice->advance(-remainT);
            printf("\rtime = %-10g", plattice->time()); fflush(stdout);
            break;
        }
        plattice->advance(dt);

        printf("\rtime = %-10g", plattice->time()); fflush(stdout);
        remainT = plattice->time() - endT;
    }
    fprintf(stderr, "\nTotal time: %lf (sec)\n", GetMilliTimed() - st);
    printf("\nDone!\n");
    write_u(0);

    clean();
    return 0;
}



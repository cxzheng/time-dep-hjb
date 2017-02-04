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
 *       Filename:  MixedTDHJ2D.cpp
 *
 *        Version:  1.0
 *        Created:  10/15/11 09:18:00
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
#include "io/npy_arrays.hpp"
#include "HJ2DFunc.h"
#include "BoundaryCond.h"
#include "BoundCondTemp.hpp"
#include "TimeIndepSpeedMixedLattice2D.hpp"
#include "MixedLattice2D.hpp"
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
typedef TimeIndepSpeedMixedLattice2D<func_F1, func_G1, TBC>  TLattice;

#elif defined(USE_TEST_05)
#warning ====== Use Test Case #05 ======
typedef SimpleBoundaryCond2D_V4                                 TBC;
typedef TimeIndepSpeedMixedLattice2D<func_F_V4, func_G1, TBC>   TLattice;

#elif defined(USE_TEST_07)
#warning ====== Use Test Case #07 ======
typedef SimpleBoundaryCond2D                    TBC;
typedef MixedLattice2D<Test07::F, func_G1, TBC> TLattice;

#elif defined(USE_TEST_21)
#warning ====== Use Test Case #21 ======
typedef BoundCond2D_1Side_V2<func_exp>                          TBC;
typedef TimeIndepSpeedMixedLattice2D<func_F1, func_G_V3, TBC>   TLattice;

#elif defined(USE_TEST_22)
#warning ====== Use Test Case #22 ======
typedef BoundCond2D_1Side_V3<func_exp>                          TBC;
typedef TimeIndepSpeedMixedLattice2D<func_F_V10, func_G1, TBC>  TLattice;

#elif defined(USE_TEST_23)
#warning ====== Use Test Case #23 ======
typedef SimpleBoundaryCond2D_V42                                TBC;
typedef TimeIndepSpeedMixedLattice2D<func_F_V4, func_G1, TBC>   TLattice;

#elif defined(USE_TEST_30)
#warning ====== Use Test Case #30 ======
typedef Test30::BoundaryCond2D                                  TBC;
typedef TimeIndepSpeedMixedLattice2D<Test30::F, func_G0, TBC>   TLattice;

#elif defined(USE_TEST_40)
#warning ====== Use Test Case #40 ======
typedef Test40::BoundaryCond2D                                  TBC;
typedef TimeIndepSpeedMixedLattice2D<Test40::F, func_G1, TBC>   TLattice;

#else
#warning ====== Use Test Case #00 ======
typedef SimpleBoundaryCond2D                TBC;
typedef TimeIndepSpeedMixedLattice2D<func_F1, func_G1, TBC>  TLattice;
#endif

// -------------------------------------------------------------------

static TBC*         pbc = NULL;
static TLattice*    plattice = NULL;

static double lx = 1, ly = 1;
static int nx = 100, ny = 100;
static double T = 2., dt = 0.05;
static double endT = 0.;
static string outFile = "out.npy";

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
            ("out,o", po::value<string>(&outFile)->default_value("out.npy"),
                    "file name patter for output data dumps");
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

    cout << "==============================================" << endl;
#if defined(USE_TEST_05)
    cout << "  R U N  T E S T *05* " << endl;
    pbc = new TBC(nx, ny);
#elif defined(USE_TEST_07)
    cout << "  R U N  T E S T *07* " << endl;
    pbc = new TBC(nx, ny);
#elif defined(USE_TEST_21)
    cout << "  R U N  T E S T *21* " << endl;
    pbc = new TBC(lx, ly, nx, ny);      // 20,21,22
#elif defined(USE_TEST_22)
    cout << "  R U N  T E S T *22* " << endl;
    pbc = new TBC(lx, ly, nx, ny);      // 20,21,22
#elif defined(USE_TEST_23)
    cout << "  R U N  T E S T *23* " << endl;
    pbc = new TBC(nx, ny);
#elif defined(USE_TEST_30)
    cout << "  R U N  T E S T *30* " << endl;
    pbc = new TBC(lx, ly, nx, ny);
#elif defined(USE_TEST_40)
    cout << "  R U N  T E S T *40* " << endl;
    pbc = new TBC(lx, ly, nx, ny);
#else
    cout << "  R U N  T E S T *00* " << endl;
    pbc = new SimpleBoundaryCond2D(nx, ny);     // 0,4
#endif
    cout << "----------------------" << endl;
    cout << "  WxH        = " << lx << " x " << ly << endl;
    cout << "  Nx X Ny    = " << nx << " x " << ny << endl;
    cout << "  Tmax       = " << T << endl;
    cout << "  dt         = " << dt << endl;
    cout << "  outfile    = " << outFile << endl; 
#ifdef USE_OPENMP
    cout << "  OPENMP:      enabled" << endl;
#endif
    cout << "==============================================" << endl;

    plattice = new TLattice(lx, ly, nx, ny, pbc, dt);
    plattice->init(T);

    dt = -dt;
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

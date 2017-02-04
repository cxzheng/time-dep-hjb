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
 *       Filename:  InterpErr.h
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
 
#ifndef INTERP_ERR_H
#   define INTERP_ERR_H

#include <string>
#include "utils/arrays.hpp"

class InterpErr
{
    public:
        InterpErr(double lx, double ly, int nx, int ny, 
                  const std::string& fileptn):
                nx_(nx), ny_(ny), sx_(lx/nx), sy_(ly/ny),
                fileFmt_(fileptn)
        { }

        bool operator() (const boost::multi_array<double, 2>& u,
                         double ts, double& L1, double& Linf);
    private:
        int     nx_, ny_;
        double  sx_, sy_;   // grid size
        const std::string   fileFmt_;
        boost::multi_array<double, 2>   accU_;
};

#endif


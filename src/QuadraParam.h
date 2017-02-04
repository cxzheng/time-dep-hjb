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
 *       Filename:  QuadraParam.h
 *
 *    Description:  Parameters for quadratic equation
 *
 *        Version:  1.0
 *        Created:  10/14/11 12:17:33
 *       Revision:  none
 *       Compiler:  gcc/intel compiler
 *
 *         Author:  Changxi Zheng (cz), cxz@cs.columbia.edu
 *                  Columbia University
 *
 * =====================================================================================
 */
#ifndef QUADRA_PARAM_INC
#   define QUADRA_PARAM_INC

#include <stdint.h>
#include "vector2.hpp"

struct QuadraParam2D
{
    vector2d            pos;
    vector2d            upU;
    vector2d            H;
};

#endif

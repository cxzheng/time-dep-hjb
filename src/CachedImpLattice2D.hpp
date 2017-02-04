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
 *       Filename:  CachedImpLattice2D.hpp
 *
 *    Description:  The same as ImplicitLattice2D but using cached f/g values
 *
 *        Version:  1.0
 *        Created:  10/14/11 12:07:52
 *       Revision:  none
 *       Compiler:  gcc/intel compiler
 *
 *         Author:  Changxi Zheng (cz), cxz@cs.columbia.edu
 *                  Columbia University
 *
 * =====================================================================================
 */
#ifndef CACHED_IMP_LATTICE_2D_INC
#   define CACHED_IMP_LATTICE_2D_INC

#include <stdint.h>
#include <limits>
#include "HJBQuadraSolver2D.h"
#include "LatticeData2D.h"
#include "generic/PriorityQueue.hpp"

template<typename _F, typename _G, class _TBC>
class CachedImpLattice2D : public LatticeData2D
{
    public:
        struct TNodeRec
        {
            int         qIdx;
            vector2i    ipos;
            double      u;

            TNodeRec() { }
            TNodeRec(int ix, int iy, double u):ipos(ix, iy), u(u) { }

            bool operator < (const TNodeRec& b) const
            {   return b.u > u; }
        };

        struct TNode
        {
            TNodeRec            rec;
            HJBQuadraSolver2D   solver;
        };

        typedef PriorityQueue<TNodeRec>     TPQ;

    public:
        CachedImpLattice2D(double lx, double ly, int nx, int ny, _TBC* bc);

        //! Initiate with time boundary condition
        void init(double ts);

        void fm_begin_next_step(double ts, double dt, TPQ& que);
        //! Solve and add new neighbors to the queue 
        void fm_add_neighbors(const TNodeRec* curnode, TPQ& que);

    private:
        double solve_node(int tx, int ty, int cdir);

    private:
        boost::multi_array<uint8_t, 2>      labels_;
        boost::multi_array<TNode, 2>        nodes_;

        double  ts_;        // current time
        double  dt_;
        _TBC*   bc_;        // define the boundary condition
        _F      f_;
        _G      g_;
};

// --------------------------------------------------------------------------------------

template<typename _F, typename _G, class _TBC>
CachedImpLattice2D<_F, _G, _TBC>::CachedImpLattice2D(double lx, double ly, 
        int nx, int ny, _TBC* bc):LatticeData2D(lx, ly, nx, ny),
        labels_(boost::extents[ny+1][nx+1]),
        nodes_(boost::extents[ny+1][nx+1]),
        bc_(bc)
{
    assert(bc);

    for(int iy = 0;iy <= RES_.y;++ iy)
    for(int ix = 0;ix <= RES_.x;++ ix)
    {
        TNode& node = nodes_[iy][ix];
        node.rec.ipos.set(ix, iy);
        node.solver.pos.set(ix*H_.x, iy*H_.y);
        node.solver.H = H_;
    }
}

template<typename _F, typename _G, class _TBC>
void CachedImpLattice2D<_F, _G, _TBC>::init(double ts)
{
    for(int iy = 0;iy <= RES_.y;++ iy)
    for(int ix = 0;ix <= RES_.x;++ ix)
        u_[curUId_][iy][ix] = bc_->boundary_val(ix, iy, ts);
}

/*!
 * Fast Marching method:
 */
template<typename _F, typename _G, class _TBC>
void CachedImpLattice2D<_F, _G, _TBC>::fm_begin_next_step(double ts, double dt, TPQ& que)
{
    curUId_ = 1 - curUId_;
    zero_multi_array(labels_);  // initialize labels
    ts_ = ts;
    dt_ = fabs(dt);

    // set boundary values
    const int BN = bc_->num_boundary_nodes();
    const vector2i* pos = bc_->boundary_pos();
    const double*   bv  = bc_->boundary_val(ts);

    for(int i = 0;i < BN;++ i)
    {
        const int ty = pos[i].y;
        const int tx = pos[i].x;
        labels_[ty][tx] = 2;
        u_[curUId_][ty][tx] = nodes_[ty][tx].rec.u = bv[i];
        que.push(&(nodes_[ty][tx].rec));
    }
}

template<typename _F, typename _G, class _TBC>
void CachedImpLattice2D<_F, _G, _TBC>::fm_add_neighbors(const TNodeRec* curnode, TPQ& que)
{
    labels_[curnode->ipos.y][curnode->ipos.x] = 2;

    for(int i = 0;i < 4;++ i)
    {
        const int tx = curnode->ipos.x + NEIGH_DIRS[i][0];
        const int ty = curnode->ipos.y + NEIGH_DIRS[i][1];

        if ( in_bound(tx, ty) && labels_[ty][tx] < 2 ) 
        {
            if ( labels_[ty][tx] == 0 )
            {
                // update f, g and u_bound
                HJBQuadraSolver2D& solver = nodes_[ty][tx].solver;
                TNodeRec&             rec = nodes_[ty][tx].rec;
                solver.fval  = f_(solver.pos, ts_);
                solver.gval  = g_(solver.pos, ts_);
                solver.dt    = dt_;
                solver.lastU = u_[1 - curUId_][ty][tx];
                solver.update_u_bound();

                double umin = solve_node(tx, ty, i);
                labels_[ty][tx] = 1;
                u_[curUId_][ty][tx] = rec.u = umin;
                que.push(&rec);
            }
            else if ( u_[curUId_][curnode->ipos.y][curnode->ipos.x] < u_[curUId_][ty][tx] )
            {
                double umin = solve_node(tx, ty, i);
                if ( umin < u_[curUId_][ty][tx] )    // labels_[ty][tx] == 1
                {
                    TNodeRec& rec = nodes_[ty][tx].rec;
                    u_[curUId_][ty][tx] = rec.u = umin;
                    que.update_node(&rec);
                }
            }
        }
    }   // end for
}

template<typename _F, typename _G, class _TBC>
double CachedImpLattice2D<_F, _G, _TBC>::solve_node(int tx, int ty, int cdir)
{
    double ret = std::numeric_limits<double>::infinity();
    HJBQuadraSolver2D& solver = nodes_[ty][tx].solver;

    if ( SOLVE_DIR[cdir][0] )   // y dir
    {
        solver.upU.x = u_[curUId_][ty][tx+SOLVE_DIR[cdir][1]];

        int ttyMinus = ty - 1;
        int ttyPlus  = ty + 1;
        if ( ttyMinus >= 0 && labels_[ttyMinus][tx] > 1 )
        {
            if ( ttyPlus <= RES_.y && labels_[ttyPlus][tx] > 1 && 
                 u_[curUId_][ttyPlus][tx] < u_[curUId_][ttyMinus][tx] )
            {   // use ty+1
                solver.upU.y = u_[curUId_][ttyPlus][tx];
                ret = fmin(ret, solver.two_sided_solve());
            }
            else
            {   // use ty-1
                solver.upU.y = u_[curUId_][ttyMinus][tx];
                ret = fmin(ret, solver.two_sided_solve());
            }
        }
        else if ( ttyPlus <= RES_.y && labels_[ttyPlus][tx] > 1 )
        {
            solver.upU.y = u_[curUId_][ttyPlus][tx];
            ret = fmin(ret, solver.two_sided_solve());
        }
        else
        {
            ret = fmin(ret, solver.one_sided_solve_x());
        }
    }
    else                        // x dir
    {
        solver.upU.y = u_[curUId_][ty+SOLVE_DIR[cdir][1]][tx];

        int ttxMinus = tx - 1;
        int ttxPlus  = tx + 1;
        if ( ttxMinus >= 0 && labels_[ty][ttxMinus] > 1 )
        {
            if ( ttxPlus <= RES_.x && labels_[ty][ttxPlus] > 1 &&
                 u_[curUId_][ty][ttxPlus] < u_[curUId_][ty][ttxMinus] )
            {
                solver.upU.x = u_[curUId_][ty][ttxPlus];
                ret = fmin(ret, solver.two_sided_solve());
            }
            else
            {
                solver.upU.x = u_[curUId_][ty][ttxMinus];
                ret = fmin(ret, solver.two_sided_solve());
            }
        }
        else if ( ttxPlus <= RES_.x && labels_[ty][ttxPlus] > 1 )
        {
            solver.upU.x = u_[curUId_][ty][ttxPlus];
            ret = fmin(ret, solver.two_sided_solve());
        }
        else
        {
            ret = fmin(ret, solver.one_sided_solve_y());
        }
    }

    return ret;
}

#endif

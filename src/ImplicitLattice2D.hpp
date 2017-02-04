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
 *       Filename:  ImplicitLattice2D.hpp
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
 
#ifndef IMPLICIT_LATTICE_2D_HPP
#   define IMPLICIT_LATTICE_2D_HPP

#include <limits>
#include "LatticeData2D.h"
#include "BoundaryCond.h"
#include "generic/PriorityQueue.hpp"

//! 2D Lattice for time dependent HJ solver
template <typename _TSolver2D, class _TBC>
class ImplicitLattice2D : public LatticeData2D
{
    public:
        //! Record structure for using in Priority queue
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

        typedef PriorityQueue<TNodeRec>         TPQ;

        // ------------------------------------------------
        ImplicitLattice2D(double lx, double ly, int nx, int ny, _TBC* bc); 

        //! Initiate with time boundary condition
        void init(double ts);

        // ------------------------------------------------
        void fm_begin_next_step(double ts, double dt, TPQ& que);
        //! Solve and add new neighbors to the queue 
        void fm_add_neighbors(const TNodeRec* curnode, TPQ& que);
        //! Move the given node to the alive set
        void fm_node_alive(const TNodeRec* curnode)
        {   labels_[curnode->ipos.y][curnode->ipos.x] = 2;  }

#ifdef COUNT_UPDATE_NUM
        double  average_update_cnt() const;
#endif

    private:
        double solve_node(int tx, int ty, int cdir);

    private:
        boost::multi_array<uint8_t, 2>  labels_;
        boost::multi_array<TNodeRec,2>  nodeRec_;
#ifdef COUNT_UPDATE_NUM
        boost::multi_array<uint8_t, 2>  updateCnt_;
#endif

        double                          ts_;        // current time
        _TBC*                           bc_;        // define the boundary condition
        _TSolver2D                      solver_;
        TDQuadraParam2D                 param_;
};

///////////////////////////////////////////////////////////////////////////////
template <typename _TSolver2D, class _TBC>
ImplicitLattice2D<_TSolver2D,_TBC>::ImplicitLattice2D(
        double lx, double ly, int nx, int ny, _TBC* bc):
    LatticeData2D(lx, ly, nx, ny),
    labels_(boost::extents[ny+1][nx+1]),
    nodeRec_(boost::extents[ny+1][nx+1]), 
#ifdef COUNT_UPDATE_NUM
    updateCnt_(boost::extents[ny+1][nx+1]),
#endif
    bc_(bc)
{
    assert(bc);
    param_.H = H_;

#ifdef USE_OPENMP
    #pragma omp parallel for default(none) 
#endif
    for(int iy = 0;iy <= RES_.y;++ iy)
    for(int ix = 0;ix <= RES_.x;++ ix)
        nodeRec_[iy][ix].ipos.set(ix, iy);
}

template <typename _TSolver2D, class _TBC>
void ImplicitLattice2D<_TSolver2D,_TBC>::init(double ts)
{
#ifdef USE_OPENMP
    #pragma omp parallel for default(none) shared(ts)
#endif
    for(int iy = 0;iy <= RES_.y;++ iy)
    for(int ix = 0;ix <= RES_.x;++ ix)
        u_[curUId_][iy][ix] = bc_->boundary_val(ix, iy, ts);
}

/*!
 * Fast Marching method:
 * - put all the boundary points into the priority queue
 * - set the label into 2 (accepted, so their value can't change anymore)
 */
template <typename _TSolver2D, class _TBC>
void ImplicitLattice2D<_TSolver2D,_TBC>::fm_begin_next_step(double ts, double dt, TPQ& que)
{
    curUId_ = 1 - curUId_;
    zero_multi_array(labels_);  // initialize labels
#ifdef COUNT_UPDATE_NUM
    zero_multi_array(updateCnt_);
#endif
    ts_ = ts;
    // set time parameter
    param_.ts = ts;
    param_.dt = fabs(dt);
    // set boundary values
    const int BN = bc_->num_boundary_nodes();
    const vector2i* pos = bc_->boundary_pos();
    const double*   bv  = bc_->boundary_val(ts);
#ifdef USE_OPENMP
    #pragma omp parallel for default(none) shared(bv, pos)
#endif
    for(int i = 0;i < BN;++ i)
    {
        labels_[pos[i].y][pos[i].x] = 2;
        u_[curUId_][pos[i].y][pos[i].x] = nodeRec_[pos[i].y][pos[i].x].u = bv[i];
        que.push(&nodeRec_[pos[i].y][pos[i].x]);
    }
}

template <typename _TSolver2D, class _TBC>
double ImplicitLattice2D<_TSolver2D,_TBC>::solve_node(int tx, int ty, int cdir)
{
    double ret = std::numeric_limits<double>::infinity();

    get_node_pos(tx, ty, param_.pos);
    param_.lastU = u_[1 - curUId_][ty][tx];
    
    if ( SOLVE_DIR[cdir][0] )   // y dir
    {
        param_.deriv.x = (int8_t)SOLVE_DIR[cdir][1];    // deriv.x != 0
        param_.upU.x   = u_[curUId_][ty][tx+param_.deriv.x];

        int ttyMinus = ty - 1;
        int ttyPlus  = ty + 1;
        if ( ttyMinus >= 0 && labels_[ttyMinus][tx] > 1 )
        {
            if ( ttyPlus <= RES_.y && labels_[ttyPlus][tx] > 1 && 
                 u_[curUId_][ttyPlus][tx] < u_[curUId_][ttyMinus][tx] )
            {   // use ty+1
                param_.deriv.y = 1;
                param_.upU.y   = u_[curUId_][ttyPlus][tx];
                ret = fmin(ret, solver_.two_sided_solve(&param_));      // two-sided solve
#ifdef COUNT_UPDATE_NUM
                ++ updateCnt_[ty][tx];
#endif
            }
            else
            {   // use ty-1
                param_.deriv.y = -1;
                param_.upU.y   = u_[curUId_][ttyMinus][tx];
                ret = fmin(ret, solver_.two_sided_solve(&param_));      // two-sided solve
#ifdef COUNT_UPDATE_NUM
                ++ updateCnt_[ty][tx];
#endif
            }
        }
        else if ( ttyPlus <= RES_.y && labels_[ttyPlus][tx] > 1 )
        {   // ty-1 is not available, ty+1 is available
            param_.deriv.y = 1;
            param_.upU.y   = u_[curUId_][ttyPlus][tx];
            ret = fmin(ret, solver_.two_sided_solve(&param_));          // two-sided solve
#ifdef COUNT_UPDATE_NUM
            ++ updateCnt_[ty][tx];
#endif
        }
        else
        {
            param_.deriv.y = 0;
            ret = fmin(ret, solver_.one_sided_solve_x(&param_));          // one-sided solve
#ifdef COUNT_UPDATE_NUM
            ++ updateCnt_[ty][tx];
#endif
        }
    }
    else                        // x dir
    {
        param_.deriv.y = (int8_t)SOLVE_DIR[cdir][1];
        param_.upU.y   = u_[curUId_][ty+param_.deriv.y][tx];

        int ttxMinus = tx - 1;
        int ttxPlus  = tx + 1;
        if ( ttxMinus >= 0 && labels_[ty][ttxMinus] > 1 )
        {
            if ( ttxPlus <= RES_.x && labels_[ty][ttxPlus] > 1 &&
                 u_[curUId_][ty][ttxPlus] < u_[curUId_][ty][ttxMinus] )
            {
                param_.deriv.x = 1;
                param_.upU.x   = u_[curUId_][ty][ttxPlus];
                ret = fmin(ret, solver_.two_sided_solve(&param_));
#ifdef COUNT_UPDATE_NUM
                ++ updateCnt_[ty][tx];
#endif
            }
            else
            {
                param_.deriv.x = -1;
                param_.upU.x   = u_[curUId_][ty][ttxMinus];
                ret = fmin(ret, solver_.two_sided_solve(&param_));
#ifdef COUNT_UPDATE_NUM
                ++ updateCnt_[ty][tx];
#endif
            }
        }
        else if ( ttxPlus <= RES_.x && labels_[ty][ttxPlus] > 1 )
        {
            param_.deriv.x = 1;
            param_.upU.x   = u_[curUId_][ty][ttxPlus];
            ret = fmin(ret, solver_.two_sided_solve(&param_));
#ifdef COUNT_UPDATE_NUM
            ++ updateCnt_[ty][tx];
#endif
        }
        else
        {
            param_.deriv.x = 0;
            ret = fmin(ret, solver_.one_sided_solve_y(&param_));
#ifdef COUNT_UPDATE_NUM
            ++ updateCnt_[ty][tx];
#endif
        }
    }

    return ret;
}

template <typename _TSolver2D, class _TBC>
void ImplicitLattice2D<_TSolver2D,_TBC>::fm_add_neighbors(const TNodeRec* curnode, TPQ& que)
{
    labels_[curnode->ipos.y][curnode->ipos.x] = 2;

    for(int i = 0;i < 4;++ i)
    {   // iterate 4 neighbors
        const int tx = curnode->ipos.x + NEIGH_DIRS[i][0];
        const int ty = curnode->ipos.y + NEIGH_DIRS[i][1];

        if ( in_bound(tx, ty) && labels_[ty][tx] < 2 ) 
        {
            if ( labels_[ty][tx] == 0 )
            {
                double umin = solve_node(tx, ty, i);
                labels_[ty][tx]     = 1;
                u_[curUId_][ty][tx] = nodeRec_[ty][tx].u = umin;
                que.push(&nodeRec_[ty][tx]);
            }
            else if ( u_[curUId_][curnode->ipos.y][curnode->ipos.x] < u_[curUId_][ty][tx] )
            {   // labels_[ty][tx] == 1
                double umin = solve_node(tx, ty, i);
                if ( umin < u_[curUId_][ty][tx] )    // labels_[ty][tx] == 1
                {
                    u_[curUId_][ty][tx] = nodeRec_[ty][tx].u = umin;
                    que.update_node(&nodeRec_[ty][tx]);
                }
            }
        }
    } // end for
}

#ifdef COUNT_UPDATE_NUM
template <typename _TSolver2D, class _TBC>
double ImplicitLattice2D<_TSolver2D,_TBC>::average_update_cnt() const
{
    int sum = 0;
#ifdef USE_OPENMP
    #pragma omp parallel for default(none) reduction(+:sum)
#endif
    for(int iy = 0;iy <= RES_.y;++ iy)
    for(int ix = 0;ix <= RES_.x;++ ix)
        sum += updateCnt_[iy][ix];
    return (double)sum/(double)((RES_.x+1)*(RES_.y+1));
}
#endif

#endif



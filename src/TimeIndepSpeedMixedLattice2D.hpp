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
 *       Filename:  TimeIndepSpeedMixedLattice2D.hpp
 *
 *    Description:  2D lattice for mixed solver where the speed function (F) is assumed
 *                  to be constant
 *
 *        Version:  1.0
 *        Created:  10/14/11 11:44:02
 *       Revision:  none
 *       Compiler:  gcc/intel compiler
 *
 *         Author:  Changxi Zheng (cz), cxz@cs.columbia.edu
 *                  Columbia University
 *
 * =====================================================================================
 */
#ifndef TIME_INDEP_SPEED_MIXED_LATTICE_2D_INC
#   define TIME_INDEP_SPEED_MIXED_LATTICE_2D_INC

#include <vector>
#include "LatticeData2D.h"
#include "UpwindGrad2D.h"
#include "HJBQuadraSolver2D.h"
#include "utils/macros.h"
#include "generic/PriorityQueue.hpp"

template <typename _F, typename _G, class _TBC>
class TimeIndepSpeedMixedLattice2D : public LatticeData2D
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
            double              fval;
            TNodeRec            rec;
            HJBQuadraSolver2D   solver;
        };

        typedef PriorityQueue<TNodeRec>     TPQ;

    public:
        TimeIndepSpeedMixedLattice2D(double lx, double ly, int nx, int ny, 
                                     _TBC* bc, double dt);

        //! Initiate with time boundary condition
        void init(double T);

        void advance(double dt);

        double time() const
        {   return ts_; }

    private:
        void begin_next_step(double ts);
        void explicit_advance(double dt);
        void begin_implicit_step();
        //! Solve and add new neighbors to the queue 
        void add_neighbors(TNodeRec* curnode);
        double solve_node(int tx, int ty, int cdir);

    private:
        TPQ             que_;
        UpwindGrad2D    gd_;

        boost::multi_array<uint8_t, 2>      labels_;
        boost::multi_array<TNode, 2>        nodes_;

        // --- list of points using implicit/explicit method ---
        std::vector<TNode*>                 expNodes_;
        int                                 nImpNodes_;

        double      ts_;
        double      dt_;
        _TBC*       bc_;
        _F          f_;
        _G          g_;
};

// --------------------------------------------------------------------------------------

template <typename _F, typename _G, class _TBC>
TimeIndepSpeedMixedLattice2D<_F, _G, _TBC>::TimeIndepSpeedMixedLattice2D(
        double lx, double ly, int nx, int ny, _TBC* bc, double dt):
        LatticeData2D(lx, ly, nx, ny),
        gd_(lx, ly, nx, ny),
        labels_(boost::extents[ny+1][nx+1]),
        nodes_(boost::extents[ny+1][nx+1]),
        nImpNodes_(0), bc_(bc)
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
    que_.resize(num_nodes());

    // initialize the implicit/explicit list
    zero_multi_array(labels_);
    const int BN = bc_->num_boundary_nodes();
    const vector2i* pos = bc_->boundary_pos();
    for(int i = 0;i < BN;++ i)
    {
        const int ty = pos[i].y;
        const int tx = pos[i].x;
        labels_[ty][tx] = 1;
    }
    for(int iy = 0;iy <= RES_.y;++ iy)
    for(int ix = 0;ix <= RES_.x;++ ix)
    {
        if ( labels_[iy][ix] ) continue;    // ignore boundary node
        TNode& node = nodes_[iy][ix];
        node.fval = f_(node.solver.pos);
        if ( -node.fval * dt * (M_SQRT2+1e-12) < H_.x )        // NOTE: ASSUME H_.x == H_.y
            expNodes_.push_back(&node);
        else
            ++ nImpNodes_;
    }
    printf("# of explicit nodes: %d\n", (int)expNodes_.size());
    printf("# of implicit nodes: %d\n", nImpNodes_);
}

template <typename _F, typename _G, class _TBC>
void TimeIndepSpeedMixedLattice2D<_F, _G, _TBC>::init(double T)
{
    ts_ = T;
    for(int iy = 0;iy <= RES_.y;++ iy)
    for(int ix = 0;ix <= RES_.x;++ ix)
        u_[curUId_][iy][ix] = bc_->boundary_val(ix, iy, T);
}

template <typename _F, typename _G, class _TBC>
void TimeIndepSpeedMixedLattice2D<_F, _G, _TBC>::advance(double dt)
{
    curUId_ = 1 - curUId_;

    begin_next_step(ts_ + dt);
    // ----- explict part ------
    explicit_advance(dt);
    ts_ += dt;

    if ( unlikely(nImpNodes_==0) ) return;

    // ------ implicit part ------
    dt_ = fabs(dt);
    begin_implicit_step();
    TNodeRec* curNode;
    while ( !que_.empty() ) 
    {
        curNode = que_.pop();
        add_neighbors(curNode);
    }
}

//! update boundary nodes
template <typename _F, typename _G, class _TBC>
void TimeIndepSpeedMixedLattice2D<_F, _G, _TBC>::begin_next_step(double ts)
{
    const int BN = bc_->num_boundary_nodes();
    const vector2i* pos = bc_->boundary_pos();
    const double*   bv  = bc_->boundary_val(ts);
    for(int i = 0;i < BN;++ i)
        u_[curUId_][pos[i].y][pos[i].x] = bv[i];
}

template <typename _F, typename _G, class _TBC>
void TimeIndepSpeedMixedLattice2D<_F, _G, _TBC>::explicit_advance(double dt)
{
    int lastUId = 1 - curUId_;
    for(size_t i = 0;i < expNodes_.size();++ i)
    {
        const vector2i& ipos = expNodes_[i]->rec.ipos;
        double dx = gd_.dx(ipos, u_[lastUId]); 
        double dy = gd_.dy(ipos, u_[lastUId]);
        double gdu = sqrt(dx*dx + dy*dy);
        double gval = g_(expNodes_[i]->solver.pos, ts_);
        u_[curUId_][ipos.y][ipos.x] = (gval - expNodes_[i]->fval*gdu)*dt + 
                                      u_[lastUId][ipos.y][ipos.x];
    }
}

template <typename _F, typename _G, class _TBC>
void TimeIndepSpeedMixedLattice2D<_F, _G, _TBC>::begin_implicit_step()
{
    zero_multi_array(labels_);

    // ------ push Boundary into queue ------
    const int BN = bc_->num_boundary_nodes();
    const vector2i* pos = bc_->boundary_pos();
    for(int i = 0;i < BN;++ i)
    {
        const int ty = pos[i].y;
        const int tx = pos[i].x;
        labels_[ty][tx] = 2;
        nodes_[ty][tx].rec.u = u_[curUId_][ty][tx];
        que_.push(&(nodes_[ty][tx].rec));
    }

    // ------ push Explicit nodes ------
    for(size_t i = 0;i < expNodes_.size();++ i)
    {
        const vector2i& ipos = expNodes_[i]->rec.ipos;
        labels_[ipos.y][ipos.x] = 2;
        expNodes_[i]->rec.u = u_[curUId_][ipos.y][ipos.x];
        que_.push(&(expNodes_[i]->rec));
    }
}

template <typename _F, typename _G, class _TBC>
void TimeIndepSpeedMixedLattice2D<_F, _G, _TBC>::add_neighbors(TNodeRec* curnode)
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
                HJBQuadraSolver2D& solver = nodes_[ty][tx].solver;
                TNodeRec&             rec = nodes_[ty][tx].rec;
                solver.fval  = nodes_[ty][tx].fval;
                solver.gval  = g_(solver.pos, ts_);
                solver.dt    = dt_;
                solver.lastU = u_[1 - curUId_][ty][tx];
                solver.update_u_bound();

                double umin = solve_node(tx, ty, i);
                labels_[ty][tx] = 1;
                u_[curUId_][ty][tx] = rec.u = umin;
                que_.push(&rec);
            }
            else if ( u_[curUId_][curnode->ipos.y][curnode->ipos.x] < u_[curUId_][ty][tx] )
            {
                double umin = solve_node(tx, ty, i);
                if ( umin < u_[curUId_][ty][tx] )    // labels_[ty][tx] == 1
                {
                    TNodeRec& rec = nodes_[ty][tx].rec;
                    u_[curUId_][ty][tx] = rec.u = umin;
                    que_.update_node(&rec);
                }
            }
        }
    }
}

template <typename _F, typename _G, class _TBC>
double TimeIndepSpeedMixedLattice2D<_F, _G, _TBC>::solve_node(int tx, int ty, int cdir)
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

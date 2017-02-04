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
 *       Filename:  MixedLattice2D.hpp
 *
 *    Description:  2D lattice for mixed solver
 *
 *        Version:  1.0
 *        Created:  10/14/11 11:11:39
 *       Revision:  none
 *       Compiler:  gcc/intel compiler
 *
 *         Author:  Changxi Zheng (cz), cxz@cs.columbia.edu
 *                  Columbia University
 *
 * =====================================================================================
 */
#ifndef MIXED_LATTICE_2D_INC
#   define MIXED_LATTICE_2D_INC

#include "LatticeData2D.h"

template <typename _F, typename _G, class _TBC>
class MixedLattice2D : public LatticeData2D
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
        MixedLattice2D(double lx, double ly, int nx, int ny, 
                       _TBC* bc, double dt);
        
        void init(double T);

        void advance(double dt);

        double time() const
        {   return ts_; }

    private:
        void begin_next_step(double ts);
        void explicit_advance();
        void add_neighbors(TNodeRec* curnode);  // Solve and add new neighbors to the queue 
        void begin_implicit_step();
        double solve_node(int tx, int ty, int cdir);

    private:
        double                  ts_;
        double                  dt_;
        _TBC*                   bc_;
        _F                      f_;
        _G                      g_;
        int                     lastUId_;

        TPQ                     que_;
        UpwindGrad2D            gd_;

        boost::multi_array<uint8_t, 2> labels_;
        boost::multi_array<TNode, 2>   nodes_;

        std::vector<TNode*>     expNodes_[2];
        int                     nExpNodes_[2];
        int                     nLastExpNodes_[2];

        int                     nFreeNodes_;
};

// --------------------------------------------------------------------------------------

template <typename _F, typename _G, class _TBC>
MixedLattice2D<_F, _G, _TBC>::MixedLattice2D(
        double lx, double ly, int nx, int ny, _TBC* bc, double dt): // dt > 0
        LatticeData2D(lx, ly, nx, ny),
        dt_(dt), bc_(bc), gd_(lx, ly, nx, ny),
        labels_(boost::extents[ny+1][nx+1]),
        nodes_(boost::extents[ny+1][nx+1])
{
    assert(bc);

    for(int iy = 0;iy <= RES_.y;++ iy)
    for(int ix = 0;ix <= RES_.x;++ ix)
    {
        TNode&  node = nodes_[iy][ix];
        node.rec.ipos.set(ix, iy);
        node.solver.pos.set(ix*H_.x, iy*H_.y);
        node.solver.H = H_;
    }

    que_.resize(num_nodes());
    expNodes_[0].resize(num_nodes());
    expNodes_[1].resize(num_nodes());
    nLastExpNodes_[0] = nLastExpNodes_[1] = 0;
}

template <typename _F, typename _G, class _TBC>
void MixedLattice2D<_F, _G, _TBC>::init(double T)
{
    ts_ = T;

    zero_multi_array(labels_);
    const int BN = bc_->num_boundary_nodes();
    nFreeNodes_ = num_nodes() - BN;
    const vector2i* pos = bc_->boundary_pos();

    for(int i = 0;i < BN;++ i)
    {
        const int ty = pos[i].y;
        const int tx = pos[i].x;
        labels_[ty][tx] = 1;
    }

    nExpNodes_[curUId_] = 0;
    for(int iy = 0;iy <= RES_.y;++ iy)
    for(int ix = 0;ix <= RES_.x;++ ix)
    {
        u_[curUId_][iy][ix] = bc_->boundary_val(ix, iy, T);

        // initialize the explicit/impliict list
        if ( unlikely(labels_[iy][ix]) ) continue;
        TNode& node = nodes_[iy][ix];
        node.solver.fval = f_(node.solver.pos, T);
        if ( -node.solver.fval * dt_ * (M_SQRT2+1e-12) < H_.x )  // dt > 0
            expNodes_[curUId_][nExpNodes_[curUId_] ++] = &node;
    }
    printf("# of explicit nodes: %d\n", nExpNodes_[curUId_]);
    printf("# of implicit nodes: %d\n", nFreeNodes_ - nExpNodes_[curUId_]);
}

template <typename _F, typename _G, class _TBC>
void MixedLattice2D<_F, _G, _TBC>::advance(double dt)
{
    lastUId_ = curUId_;
    curUId_ = 1 - curUId_;
    dt_ = fabs(dt);

    begin_next_step(ts_ + dt);
    // ------ explict part -------
    explicit_advance();
    ts_ += dt;
    if ( unlikely(nExpNodes_[curUId_] == nFreeNodes_) ) return;

    // ------ implicit part ------
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
void MixedLattice2D<_F, _G, _TBC>::begin_next_step(double ts)
{
    const int BN = bc_->num_boundary_nodes();
    const vector2i* pos = bc_->boundary_pos();
    const double*   bv  = bc_->boundary_val(ts);
    for(int i = 0;i < BN;++ i)
        u_[curUId_][pos[i].y][pos[i].x] = bv[i];
}

template <typename _F, typename _G, class _TBC>
void MixedLattice2D<_F, _G, _TBC>::explicit_advance()
{
    nExpNodes_[curUId_] = 0;
    for(int i = 0;i < nLastExpNodes_[lastUId_];++ i)
    {
        TNode* node = expNodes_[lastUId_][i];
        double fval = f_(node->solver.pos, ts_);
        if ( -fval * dt_ * (M_SQRT2+1e-12) < H_.x ) // explicit
        {
            const vector2i& ipos = node->rec.ipos;
            double dx = gd_.dx(ipos, u_[lastUId_]);
            double dy = gd_.dy(ipos, u_[lastUId_]);
            double gdu = sqrt(dx*dx + dy*dy);
            double gval = g_(node->solver.pos, ts_);

            u_[curUId_][ipos.y][ipos.x] = (fval*gdu - gval)*dt_ + u_[lastUId_][ipos.y][ipos.x];
            expNodes_[curUId_][nExpNodes_[curUId_]++] = node;
        }
    }

    for(int i = nLastExpNodes_[lastUId_];i < nExpNodes_[lastUId_];++ i)
    {
        TNode* node = expNodes_[lastUId_][i];
        const vector2i& ipos = node->rec.ipos;
        double dx = gd_.dx(ipos, u_[lastUId_]);
        double dy = gd_.dy(ipos, u_[lastUId_]);
        double gdu = sqrt(dx*dx + dy*dy);
        double gval = g_(node->solver.pos, ts_);
        u_[curUId_][ipos.y][ipos.x] = (node->solver.fval*gdu - gval)*dt_ + 
                                      u_[lastUId_][ipos.y][ipos.x];
        expNodes_[curUId_][nExpNodes_[curUId_]++] = node;
    }
    nLastExpNodes_[curUId_] = nExpNodes_[curUId_];
}

template <typename _F, typename _G, class _TBC>
void MixedLattice2D<_F, _G, _TBC>::begin_implicit_step()
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
    for(int i = 0;i < nLastExpNodes_[curUId_];++ i)
    {
        TNodeRec& rec = expNodes_[curUId_][i]->rec;
        const vector2i& ipos = rec.ipos;
        labels_[ipos.y][ipos.x] = 2;
        rec.u = u_[curUId_][ipos.y][ipos.x];
        que_.push(&rec);
    }
}

template <typename _F, typename _G, class _TBC>
void MixedLattice2D<_F, _G, _TBC>::add_neighbors(TNodeRec* curnode)
{
    labels_[curnode->ipos.y][curnode->ipos.x] = 2;
    
    for(int i = 0;i < 4;++ i)
    {
        const int tx = curnode->ipos.x + NEIGH_DIRS[i][0];
        const int ty = curnode->ipos.y + NEIGH_DIRS[i][1];

        if ( in_bound(tx, ty) && labels_[ty][tx] < 2 ) 
        {
            if ( labels_[ty][tx] == 0 ) 
            {   // first time for this node
                HJBQuadraSolver2D& solver = nodes_[ty][tx].solver;
                TNodeRec&             rec = nodes_[ty][tx].rec;
                solver.fval  = f_(solver.pos, ts_);
                solver.gval  = g_(solver.pos, ts_);
                solver.dt    = dt_;
                solver.lastU = u_[lastUId_][ty][tx];
                solver.update_u_bound();
                if ( -solver.fval * dt_ * (M_SQRT2+1e-12) < H_.x )  // use explicit method for the next timestep
                    expNodes_[curUId_][nExpNodes_[curUId_]++] = &(nodes_[ty][tx]);

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
double MixedLattice2D<_F, _G, _TBC>::solve_node(int tx, int ty, int cdir)
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


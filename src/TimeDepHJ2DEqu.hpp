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
 *       Filename:  TimeDepHJ2DEqu.hpp
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
 
#ifndef TIME_DEP_HJ_2D_EQU_HPP
#   define TIME_DEP_HJ_2D_EQU_HPP

/*
 * Time dependent HJ PDE solver in 2D
 */
template <class _TLattice>
class TimeDepHJ2DEqu
{
    public:
        TimeDepHJ2DEqu(_TLattice* lattice, double T);

        void advance(double dt);            // advance to the next time step

        double time() const
        {   return ts_; }

    private:
        _TLattice*                  plattice_;
        typename _TLattice::TPQ     pq_;    // priority queue
        double                      ts_;
};

///////////////////////////////////////////////////////////////////////////////

template <class _TLattice>
TimeDepHJ2DEqu<_TLattice>::TimeDepHJ2DEqu(_TLattice* lattice, double T):
        plattice_(lattice), ts_(T)
{
    plattice_->init(T);
    pq_.resize(plattice_->num_nodes());
}

template <class _TLattice>
void TimeDepHJ2DEqu<_TLattice>::advance(double dt)
{
    typename _TLattice::TNodeRec* curptr;
    ts_ += dt;   // advance time step
    plattice_->fm_begin_next_step(ts_, dt, pq_);

    while ( !pq_.empty() )
    {
        curptr = pq_.pop();
        plattice_->fm_add_neighbors(curptr, pq_);
    }
}

#endif


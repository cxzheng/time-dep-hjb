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
 *       Filename:  PriorityQueue.hpp
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
 
#ifndef PRIORITY_QUEUE
#   define PRIORITY_QUEUE

#include <functional>
#include <memory>
#include <assert.h>
#include "utils/macros.h"

#ifdef USE_OPENMP_PQ
#include <omp.h>
#endif

/*!
 * Data should have a field Data.qIdx
 * Compare returns true if the first argument is less than the second one
 */
template<
        typename Data, 
        typename Compare = std::less<Data>,
        typename Alloc = std::allocator<Data*>
        >
class PriorityQueue
{
    public:
        // =========== Constructors ==========
        PriorityQueue():mptr_queue(NULL), m_used(0) 
        { 
#ifdef USE_OPENMP_PQ
            omp_init_lock(&m_lock);
#endif
        }
        PriorityQueue(size_t s):m_size(s), m_used(0)
        {
            mptr_queue = _alloc.allocate(s);
#ifdef USE_OPENMP_PQ
            omp_init_lock(&m_lock);
#endif
        }
        ~PriorityQueue()
        {
            if ( m_size > 0 ) 
                _alloc.deallocate(mptr_queue, m_size);
#ifdef USE_OPENMP_PQ
            omp_destroy_lock(&m_lock);
#endif
        }

        void  resize(size_t size);
        void  push(Data* ptr);
        void  update_node(Data* val);
        Data* pop();
        bool  empty() const
        {
            return m_used == 0; 
        }
        void clear() 
        {
            m_used = 0;
        }
        size_t size() const { return m_used; }

    private:
        //! this one is the same with update_node
        //  but is not thread-safe
        void  private_update_node(Data* val);

        Data**          mptr_queue;
        size_t          m_size;
        size_t          m_used;

#ifdef USE_OPENMP_PQ
        omp_lock_t      m_lock;
#endif
        Compare         _comp;
        Alloc           _alloc;
};

// ==================== implementation ===================
template<typename Data, typename Compare, typename Alloc>
void PriorityQueue<Data, Compare, Alloc>::resize(size_t size)
{
#ifdef USE_OPENMP_PQ
    omp_set_lock(&m_lock);
#endif
    if ( m_size > 0 ) _alloc.deallocate(mptr_queue, m_size);
    m_size = size;
    m_used = 0;
    mptr_queue = _alloc.allocate(m_size);
#ifdef USE_OPENMP_PQ
    omp_unset_lock(&m_lock);
#endif
}

template<typename Data, typename Compare, typename Alloc>
void PriorityQueue<Data, Compare, Alloc>::push(Data* ptr)
{
#ifdef USE_OPENMP_PQ
    omp_set_lock(&m_lock);
#endif
    assert(m_used < m_size && ptr != NULL);
    mptr_queue[m_used] = ptr;
    ptr->qIdx = m_used ++;
    private_update_node(ptr);
#ifdef USE_OPENMP_PQ
    omp_unset_lock(&m_lock);
#endif
}

template<typename Data, typename Compare, typename Alloc>
void PriorityQueue<Data, Compare, Alloc>::private_update_node(Data* ptr)
{
    int parent = (ptr->qIdx - 1) / 2;
    while ( parent >= 0 && _comp(*ptr, *mptr_queue[parent]) )
    {
        mptr_queue[parent]->qIdx = ptr->qIdx;
        std::swap(mptr_queue[parent], mptr_queue[ptr->qIdx]);
        ptr->qIdx = parent;
        parent = (ptr->qIdx - 1) / 2;
    }
}

template<typename Data, typename Compare, typename Alloc>
void PriorityQueue<Data, Compare, Alloc>::update_node(Data* ptr)
{
#ifdef USE_OPENMP_PQ
    omp_set_lock(&m_lock);
#endif
    int parent = (ptr->qIdx - 1) / 2;
    while ( parent >= 0 && _comp(*ptr, *mptr_queue[parent]) )
    {
        mptr_queue[parent]->qIdx = ptr->qIdx;
        std::swap(mptr_queue[parent], mptr_queue[ptr->qIdx]);
        ptr->qIdx = parent;
        parent = (ptr->qIdx - 1) / 2;
    }
#ifdef USE_OPENMP_PQ
    omp_unset_lock(&m_lock);
#endif
}

template<typename Data, typename Compare, typename Alloc>
Data* PriorityQueue<Data, Compare, Alloc>::pop()
{
    if ( !m_used ) return NULL;

#ifdef USE_OPENMP_PQ
    omp_set_lock(&m_lock);
#endif
    Data* ret = mptr_queue[0];
    mptr_queue[0] = mptr_queue[-- m_used];
    mptr_queue[0]->qIdx = 0;

    Data* cur = mptr_queue[0];
    size_t child = 1;
    child += (child+1 < m_used && 
            _comp(*mptr_queue[child+1], *mptr_queue[child]) );
    while ( child < m_used && _comp(*mptr_queue[child], *cur) )
    {
        mptr_queue[child]->qIdx = (child - 1) / 2;
        std::swap(mptr_queue[child], mptr_queue[cur->qIdx]);
        cur->qIdx = child;
        child = child * 2 + 1;
        child += (child+1 < m_used && 
                _comp(*mptr_queue[child+1], *mptr_queue[child]) );
    }
#ifdef USE_OPENMP_PQ
    omp_unset_lock(&m_lock);
#endif
    return ret;
}

#endif


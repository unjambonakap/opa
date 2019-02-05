/*
    Copyright (c) 2011-2014 University of Zurich
    
    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
    THE SOFTWARE.
*/

#ifndef PLLL_INCLUDE_GUARD__LLL2_DEEPINSERTIONS_CPP
#define PLLL_INCLUDE_GUARD__LLL2_DEEPINSERTIONS_CPP

namespace plll
{
    template<class RealTypeContext, class IntTypeContext>
    class DoClassicDeepInsertion
    {
    private:
        LatticeReduction::DIMode d_mode;
        LatticeReduction::DIChoice d_choice;
        unsigned d_bs;
        LatticeReduction::Statistics & d_stats;
    
    public:
        DoClassicDeepInsertion(LatticeReduction::DIMode mode, LatticeReduction::DIChoice choice, unsigned bs, LatticeReduction::Statistics & stats)
            : d_mode(mode), d_choice(choice), d_bs(bs), d_stats(stats)
        {
        }
    
        template<class Enumerator, class CallbackFunction>
        bool operator() (Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction> & workspace,
                         Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned & k, signed modified_up_to,
                         bool before_sizereduce, bool size_reduce_changed)
        {
            // Check if something should be done?
            switch (d_mode)
            {
            case LatticeReduction::DIM_AfterSR: if (before_sizereduce) return false; break;
            case LatticeReduction::DIM_BeforeSR: if (!before_sizereduce) return false; break;
            case LatticeReduction::DIM_Both: if (!before_sizereduce && !size_reduce_changed) return false; break;
            }
        
            // Do deep insertions
            if (modified_up_to < (signed)k)
                modified_up_to = (signed)k;
            unsigned new_k = k;
            bool modified = false;
            for (unsigned begin_k = k; begin_k <= (unsigned)modified_up_to; ++begin_k)
            {
                k = begin_k;
            
                unsigned dest = lattice.range().begin();
                typename RealTypeContext::Real res(lattice.rc());
                for (; dest < k; ++dest)
                {
                    if (dest >= lattice.range().begin() + d_bs)
                    {
                        if (d_choice == LatticeReduction::DIC_First)
                        {
                            dest = k;
                            break;
                        }
                        if ((d_choice != LatticeReduction::DIC_All) && (dest + d_bs < k))
                            dest = k - d_bs;
                    }
                    else
                        if ((d_choice == LatticeReduction::DIC_Block) && (dest + d_bs < k))
                            dest = k - d_bs;
                
                    lattice.computeProjectionLengthBV(res, dest, k);
                    if (workspace.verbose().yieldsOutput(LatticeReduction::VL_Chatter))
                        workspace.verbose()(LatticeReduction::VL_Chatter) << "[" << dest << ":" << lattice.getNormSq(dest) << ", " << k << ":" << res << "]";
                    if (lattice.getLLLalpha() * lattice.getNormSq(dest) > res)
                        break;
                }
            
                if (dest < k)
                {
                    if (k - dest > 1)
                    {
                        if (workspace.verbose().yieldsOutput(LatticeReduction::VL_Chatter))
                            workspace.verbose()(LatticeReduction::VL_Chatter) << "Moving vector " << k << " to " << dest;
                        ++d_stats.deepinsertions;
                    }
                    while (k > dest)
                    {
                        lattice.swap(k, k - 1);
                        --k;
                    }
                    modified = true;
                    if (k < new_k)
                        new_k = k;
                }
            }
            k = new_k;
            return modified;
        }
    };
}

#include "ensurehugeexponent.hpp"

namespace plll
{
    template<class RealTypeContext, class IntTypeContext>
    class DoMinPotDeepInsertion
    {
    private:
        unsigned d_rank;
        LatticeReduction::DIMode d_mode;
        LatticeReduction::DIChoice d_choice;
        unsigned d_bs;
        LatticeReduction::Statistics & d_stats;
        bool d_find_smallest;
    
        typedef typename arithmetic::HugeExponent<RealTypeContext>::Type Huge;
//    typedef typename RealTypeContext::Real Huge;
    
    public:
        DoMinPotDeepInsertion(unsigned rank, LatticeReduction::DIMode mode, LatticeReduction::DIChoice choice, unsigned bs,
                              LatticeReduction::Statistics & stats, bool findSmallest = true)
            : d_rank(rank), d_mode(mode), d_choice(choice), d_bs(bs), d_stats(stats), d_find_smallest(findSmallest)
        {
        }
    
        template<class Enumerator, class CallbackFunction>
        bool operator() (Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction> & workspace,
                         Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned & k, signed modified_up_to,
                         bool before_sizereduce, bool size_reduce_changed)
        /*
          The potential is defined as   D := \prod_{i=begin}^end \|b_i^*\|_2^{2(A.rows()-i)}.
      
          In case   begin <= k < end,  swapping b_k with b_{k+1} changes D by a factor of
          \|\pi_k(b_{k+1})\|_2^2 / \|b_k^*\|_2^2.
      
          PROOF: D_new / D_old
          = \|b_{k,new}\|_2^{2(A.rows()-k)}*\|b_{k+1,new}\|_2^{2(A.rows()-k-1)}
          / \|b_{k,old}\|_2^{2(A.rows()-k)}*\|b_{k+1,old}\|_2^{2(A.rows()-k-1)}
          = \|proj(b_{k+1,old})\|_2^{2(A.rows()-k)}*(\|b_{k,old}\|_2*\|b_{k+1,old}\|_2\|_2/\|proj(b_{k+1,old})\|_2)^{2(A.rows()-k-1)}
          / \|b_{k,old}\|_2^{2(A.rows()-k)}*\|b_{k+1,old}\|_2^{2(A.rows()-k-1)}
          = \|proj(b_{k+1,old})\|_2^{2(A.rows()-k)-2(A.rows()-k-1)}
          / \|b_{k,old}\|_2^{2(A.rows()-k)-2(A.rows()-k-1)}*\|b_{k+1,old}\|_2^{2(A.rows()-k-1)-2(A.rows()-k-1)}
          = \|proj(b_{k+1,old})\|_2^2 / \|b_{k,old}\|_2^2
      
          Therefore, if   begin <= k <= ell <= end,   and we move the ell-th basis vector to the k-th
          position, then the overall potential change is
          \prod_{i=k}^{ell-1} (\|\pi_i(b_{ell})\|_2^2 / \|b_i^*\|_2^2).
      
          We compute this factor iteratively for k = ell-1, ell-2, ..., until the desired maximal
          insertion depth, and keep track where it reaches a minimum.
        */
        {
            // Only do potential deep insertions if the number of basis vectors equals the number of basis vectors at the beginning
            if (lattice.dimension() != d_rank)
                return false;
        
            // Check if something should be done?
            switch (d_mode)
            {
            case LatticeReduction::DIM_AfterSR: if (before_sizereduce) return false; break;
            case LatticeReduction::DIM_BeforeSR: if (!before_sizereduce) return false; break;
            case LatticeReduction::DIM_Both: if (!before_sizereduce && !size_reduce_changed) return false; break;
            }
        
            // Do deep insertions
            if (modified_up_to < (signed)k)
                modified_up_to = (signed)k;
            unsigned new_k = k;
            bool modified = false;
            Huge minpot(lattice.rc()), pot(lattice.rc());
            typename RealTypeContext::Real t(lattice.rc());
            for (unsigned begin_k = k; begin_k <= (unsigned)modified_up_to; ++begin_k)
            {
                k = begin_k;
            
                unsigned blockbegin = lattice.range().begin();
                if ((d_choice == LatticeReduction::DIC_Block) && (blockbegin + d_bs < k))
                    blockbegin = k - d_bs;
                if (blockbegin == k)
                    continue;
            
                // Find position where potential is minimal
                unsigned minimal_index = k;
                minpot = Huge(lattice.getLLLalpha());
                setOne(pot);
                lattice.update(k + 1);
                for (; blockbegin < k; --k)
                {
                    // Update potential
                    pot /= Huge(lattice.getNormSq(k - 1));
                    lattice.computeProjectionLengthBV(t, k - 1, begin_k);
                    pot *= Huge(t);
                
                    // Check whether we want to actually compare the potential for this position
                    bool doComp = true;
                    if (k > lattice.range().begin() + d_bs)
                    {
                        if (d_choice == LatticeReduction::DIC_First)
                            doComp = false;
                        if ((d_choice != LatticeReduction::DIC_All) && (k + d_bs <= begin_k))
                            doComp = false;
                    }
                    else
                        if ((d_choice == LatticeReduction::DIC_Block) && (k + d_bs <= begin_k))
                            doComp = false;
                
                    // Compare potential
                    if (doComp)
                        if (pot < minpot)
                        {
                            if (d_find_smallest)
                                minpot = pot;
                            minimal_index = k - 1;
                        }
                }
            
                // Move vector to right position
                if (minimal_index < begin_k)
                {
                    if (begin_k - minimal_index > 1)
                    {
                        ++d_stats.deepinsertions;
                        if (workspace.verbose().yieldsOutput(LatticeReduction::VL_Chatter))
                            workspace.verbose()(LatticeReduction::VL_Chatter) << "Moving vector " << begin_k << " to " << minimal_index;
                    }
                    k = begin_k;
                    while (k > minimal_index)
                    {
                        lattice.swap(k, k - 1);
                        --k;
                    }
                    // Note that we did something, and store smallest destination
                    modified = true;
                    if (k < new_k)
                        new_k = k;
                }
            }
            k = new_k;
            lattice.update(k + 1);
            return modified;
        }
    };
}

#endif

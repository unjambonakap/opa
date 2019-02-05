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

#ifndef PLLL_INCLUDE_GUARD__BKZ_IMPL_CPP
#define PLLL_INCLUDE_GUARD__BKZ_IMPL_CPP

#include <algorithm>

namespace plll
{
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    template<class Reorderer, class Annealer>
    void Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    partialBKZ_Classical(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize, Reorderer & reorder, Annealer & anneal)
    // Applies Classical (Schnorr-Euchner) BKZ to a subset of the vectors.
    {
        typename LatticeHelper<RealTypeContext, IntTypeContext>::RangePopper rangepopper(lattice);
        
        // Apply LLL
        lattice.range().dupRange();
        partialLLL(lattice, lattice.range().begin() + 1, reorder);
        
        // Prepare GS
        lattice.update(lattice.range().end() + 1);
        bool avoidInsertions = lattice.avoidInsertions();
        TransformDataChangeNotifyer<IntTypeContext> TT;
        typename LatticeHelper<RealTypeContext, IntTypeContext>::AddNotifier addnotifier(lattice, &TT);
        
        // Main loop
        linalg::math_rowvector<typename IntTypeContext::Integer> result;
        
        unsigned z = 0, j = lattice.range().end();
        while (z < lattice.range().dimension() - 1)
        {
            // Callback function call
            callback(lattice, j);
            
            // Do reorder (Deep Insertion)
            if (reorder(*this, lattice, j, j, true, false))
            {
                // Apply LLL
                lattice.range().addRange(lattice.range().begin(), j);
                partialLLL(lattice, j, reorder);
                
                // Set z to 0
                z = 0;
                continue;
            }
            
            // Choose next j
            if (++j > lattice.range().end())
                j = lattice.range().begin();
            unsigned k = std::min(lattice.range().end(), j + blocksize - 1);
            
            // Find shortest vector
            lattice.range().addRange(j, k);
            solveSVP(lattice, result, reorder, false);
            
            bool doreplace;
            if (result.size() > 0)
            {
                // Compute length of projection
                typename RealTypeContext::Real l(lattice.rc());
                lattice.computeProjectionLength(l, j, j, result);
                
                // Check Schnorr-Euchner BKZ condition
                doreplace = (lattice.getLLLalpha() * lattice.getNormSq(j) > l);
            }
            else
                // Didn't find anything
                doreplace = false;
            
            // Do annealing
            if (!doreplace)
                doreplace = anneal(lattice, k, blocksize, result);
            
            if (doreplace)
            {
                // Yes? Include new vector
                if (avoidInsertions)
                {
                    lattice.range().addRange(j, std::min(k + 1, lattice.range().end()));
                    rearrange(lattice, result);
                    lattice.range().addRange(lattice.range().begin(), std::min(k + 1, lattice.range().end()));
                    partialLLL(lattice, j, reorder);
                }
                else
                {
                    unsigned jj = lattice.range().begin();
                    lattice.range().addRange(j, std::min(k + 1, lattice.range().end()));
                    rearrangeLLL(lattice, jj, result, reorder);
                }
                
                // Set z to 0
                z = 0;
                TT.resetChangeNotifyer();
            }
            else
            {
                // No new shortest vector to be inserted?
                if (TT.wasChanged())
                {
                    // If something was changed, reset z. This is necessary since solveSVP() might
                    // modify the local basis before enumeration, and thereby change the property that
                    // the BKZ property holds (or even size reduction).
                    z = 0;
                }
                else
                    // If nothing was changed, increase z
                    ++z;
                
                // Apply LLL
                lattice.range().addRange(lattice.range().begin(), std::min(k + 1, lattice.range().end()));
                unsigned ms = partialLLL(lattice, k, reorder);
                
                // Check if something before k was modified
                if ((unsigned)ms < k)
                {
                    // If yes: we cannot assume anything anymore ==> set z = 0
                    z = 0;
                }
                TT.resetChangeNotifyer();
            }
        }
        // Now the basis is BKZ reduced.
    }
    
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    template<class Reorderer>
    void Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    partialBKZ_Classical(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize, Reorderer & reorder)
    {
        NoAnnealBKZ<RealTypeContext, IntTypeContext> anneal;
        partialBKZ_Classical(lattice, blocksize, reorder, anneal);
    }
    
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    template<class Reorderer, class Annealer>
    void Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    partialBKZ_Simplified(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize, Reorderer & reorder, Annealer & anneal)
    // Applies Simplified (Hanrot-Pujol-Stehle) BKZ to a subset of the vectors.
    {
        // Prepare GS
        lattice.update(lattice.range().end() + 1);
        TransformDataChangeNotifyer<IntTypeContext> TT;
        typename LatticeHelper<RealTypeContext, IntTypeContext>::AddNotifier addnotifier(lattice, &TT);
        typename LatticeHelper<RealTypeContext, IntTypeContext>::RangePopper rangepopper(lattice);
        bool avoidInsertions = lattice.avoidInsertions();
        
        linalg::math_rowvector<typename IntTypeContext::Integer> result;
        do
        {
            TT.resetChangeNotifyer();
            for (unsigned k = lattice.range().begin(); k < lattice.range().end(); ++k)
            {
                callback(lattice, k);
                
                // Do reorder (Deep Insertion)
                reorder(*this, lattice, k, k, true, false);
                
                // Find shortest vector
                unsigned kk = std::min(lattice.range().end(), k + blocksize - 1);
                lattice.range().addRange(k, kk);
                solveSVP(lattice, result, reorder);
                
                bool doreplace;
                if (result.size() > 0)
                {
                    // Compute length of projection
                    typename RealTypeContext::Real l(lattice.rc());
                    lattice.computeProjectionLength(l, k, k, result);
                    
                    // Check Schnorr-Euchner BKZ condition
                    doreplace = (lattice.getLLLalpha() * lattice.getNormSq(k) > l);
                }
                else
                    // Didn't find anything
                    doreplace = false;
                
                // Do annealing
                if (!doreplace)
                    doreplace = anneal(lattice, k, blocksize, result);
                
                if (doreplace)
                {
                    // Yes? Include new vector
                    if (avoidInsertions)
                    {
                        lattice.range().addRange(k, std::min(kk + 1, lattice.range().end()));
                        rearrange(lattice, result);
                        lattice.range().addRange(lattice.range().begin(), std::min(kk + 1, lattice.range().end()));
                        unsigned ms = partialLLL(lattice, k, reorder);
                        // Go back to last change
                        if ((unsigned)ms < k)
                            k = ms - 1;
                    }
                    else
                    {
                        unsigned j = lattice.range().begin();
                        lattice.range().addRange(k, std::min(kk + 1, lattice.range().end()));
                        unsigned ms = rearrangeLLL(lattice, j, result, reorder);
                        // Go back to last change
                        if ((unsigned)ms < k)
                            k = ms - 1;
                    }
                }
                else
                {
                    // Apply LLL
                    lattice.range().addRange(lattice.range().begin(), std::min(kk + 1, lattice.range().end()));
                    unsigned ms = partialLLL(lattice, kk, reorder); 
                    // Go back to last change
                    if ((unsigned)ms <= k)
                        k = ms - 1;
                }
            }
        } while (TT.wasChanged());
    }
    
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    template<class Reorderer>
    void Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    partialBKZ_Simplified(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize, Reorderer & reorder)
    {
        NoAnnealBKZ<RealTypeContext, IntTypeContext> noanneal;
        partialBKZ_Simplified(lattice, blocksize, reorder, noanneal);
    }
    
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    template<class Reorderer>
    void Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    partialBKZ_Terminating(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize,
                           bool doHKZ, const arithmetic::Integer & tours, Reorderer & reorder)
    // doHKZ: true = apply HKZ to local bases, false = apply SVP to local bases
    // tours: number of tours (i.e. runs)
    //
    // Applies Terminating BKZ (Hanrot-Pujol-Stehle) to a subset of the vectors.
    {
        // Prepare GS
        lattice.update(lattice.range().end() + 1);
        TransformDataChangeNotifyer<IntTypeContext> TT;
        typename LatticeHelper<RealTypeContext, IntTypeContext>::AddNotifier addnotifier(lattice, &TT);
        typename LatticeHelper<RealTypeContext, IntTypeContext>::RangePopper rangepopper(lattice);
        
        // Main loop
        arithmetic::Integer tour;
        do
        {
            if (lattice.range().dimension() > 40)
                verbose()(LatticeReduction::VL_Information) << tour << "/" << tours << " tours";
            TT.resetChangeNotifyer();
            unsigned end = lattice.range().end();
            if (doHKZ)
            {
                if (end < blocksize - 1)
                    end = lattice.range().begin();
                else
                {
                    end = end + 1 - blocksize;
                    if (end < lattice.range().begin())
                        end = lattice.range().begin();
                }
            }
            for (unsigned k = lattice.range().begin(); k <= end; ++k)
            {
                callback(lattice, k);
                
                // Compute SVP/HKZ basis
                if (doHKZ)
                {
                    lattice.range().addRange(k, std::min(k + blocksize - 1, lattice.range().end()));
                    computeHKZbasis(lattice);
                }
                else
                {
                    lattice.range().addRange(k, std::min(k + blocksize - 1, lattice.range().end()));
                    computeSVPbasis(lattice, true);
                }
                
                // Do reorder (Deep Insertion)
                unsigned kk = k, ek = std::min(k + blocksize - 1, lattice.range().end());
                if (reorder(*this, lattice, k, ek, true, false))
                {
                    lattice.update(kk + 1);
                    for (unsigned j = k; j <= ek; ++j)
                        lattice.sizereduce(j);
                    k = kk;
                }
                
                // Size reduction is already done except last vector (if in range):
                if (k + blocksize <= lattice.range().end())
                {
                    lattice.update(k + blocksize + 1);
                    lattice.sizereduce(k + blocksize);
                }
            }
            
            // Check if enough tours were done
            ++tour;
            if (tour >= tours)
                break;
        } while (TT.wasChanged());
    }
    
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    void Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    computeProjectedDeterminant(typename RealTypeContext::Real & result,
                                const Lattice<RealTypeContext, IntTypeContext> & lattice,
                                const std::pair<unsigned, unsigned> & block)
    {
        if (block.second < block.first)
        {
            setOne(result);
        }
        else
        {
            result = lattice.getNormSq(block.first);
            for (unsigned i = block.first + 1; i <= block.second; ++i)
                result *= lattice.getNormSq(i);
        }
    }
    
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    template<class Reorderer>
    void Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    partialBKZ_Terminating(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize, bool doHKZ, Reorderer & reorder)
    // doHKZ: true = apply HKZ to local bases, false = apply SVP to local bases
    //
    // Applies Terminating BKZ (Hanrot-Pujol-Stehle) to a subset of the vectors.
    {
        typename LatticeHelper<RealTypeContext, IntTypeContext>::RangePopper rangepopper(lattice);
        // Prepare GS
        lattice.update(lattice.range().end() + 1);
        
        // Compute approximation of squared determinant
        typename RealTypeContext::Real det(lattice.rc());
        computeProjectedDeterminant(det, lattice, std::make_pair(lattice.range().begin(), lattice.range().end()));
        if (isZero(det))
            throw reduction_error("(Numerical GS) Determinant of lattice is zero, probably due to linearly dependent vectors or precision problems. Run LLL first!");
        
        // Compute maximum of base vectors square norms
        arithmetic::RealContext irc;
        arithmetic::Real max(irc);
        if (lattice.range().begin() == 0)
        {
            arithmetic::Real t(irc);
            arithmetic::convert(max, lattice.getNormSqUP(lattice.range().begin()), irc);
            for (unsigned i = lattice.range().begin() + 1; i <= lattice.range().end(); ++i)
            {
                arithmetic::convert(t, lattice.getNormSqUP(i), irc);
                if (t > max)
                    max = t;
            }
        }
        else
        {
            typename RealTypeContext::Real t(lattice.rc()), max_(lattice.rc());
            max_ = lattice.getNormSq(lattice.range().begin());
            for (unsigned i = lattice.range().begin() + 1; i <= lattice.range().end(); ++i)
            {
                lattice.computeProjectionLengthBV(t, lattice.range().begin(), i);
                if (t > max_)
                    max_ = t;
            }
            arithmetic::convert(max, max_, irc);
        }
        
        // Compute maximal number of tours
        arithmetic::IntegerContext iic;
        arithmetic::Integer tours;
        arithmetic::Real C = arithmetic::convert(1.0, irc); // ??? !!! ... FIND OUT WHAT C SHOULD BE !!!!!!!!!
        unsigned dim = lattice.range().end() - lattice.range().begin() + 1;
        max = arithmetic::sqrt(max);
        max *= arithmetic::power(arithmetic::convert(det, irc), arithmetic::convert(-1, irc) / arithmetic::convert(2*dim, irc));
        arithmetic::Real tau = C * arithmetic::convert(dim * dim, irc) * arithmetic::convert(dim, irc)
            / arithmetic::convert(blocksize * blocksize, irc)
            * (arithmetic::log(arithmetic::convert(dim, irc)) + arithmetic::log(arithmetic::log(max)));
        convert_ceil(tours, tau, iic);
        ++tours; // for rounding problems...
        verbose()(LatticeReduction::VL_Information) << "Terminating BKZ: doing " << tours << " tour" << (isOne(tours) ? "" : "s") << "...";
        
        // Apply BKZ
        lattice.range().dupRange();
        partialBKZ_Terminating(lattice, blocksize, doHKZ, tours, reorder);
    }
    
    template<class RealTypeContext, class IntTypeContext>
    std::pair<unsigned, unsigned> getBlock(unsigned block, Lattice<RealTypeContext, IntTypeContext> & lattice,
                                           unsigned blocksize, unsigned slide = 0, unsigned enlarge = 0)
    {
        unsigned
            ib = std::min(lattice.range().begin() + block * blocksize + slide, lattice.range().end()),
            ie = std::min(lattice.range().begin() + (block + 1) * blocksize - 1 + slide + enlarge, lattice.range().end());
        return std::make_pair(ib, ie);
    }
    
    template<class RealTypeContext, class IntTypeContext>
    std::pair<unsigned, unsigned> getBlock2(unsigned block, Lattice<RealTypeContext, IntTypeContext> & lattice,
                                            unsigned blocksize, unsigned slide = 0, unsigned enlarge = 0)
    // Get indices for blocks [block, block+1]
    {
        unsigned
            ib = std::min(lattice.range().begin() + block * blocksize + slide, lattice.range().end()),
            ie = std::min(lattice.range().begin() + (block + 2) * blocksize - 1 + slide + enlarge, lattice.range().end());
        return std::make_pair(ib, ie);
    }
    
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    template<class Reorderer, class Annealer>
    void Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    partialBKZ_SemiBlock2k(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize, Reorderer & reorder, Annealer & anneal)
    /*
      Applies Schnorr's Semi Block 2k Reduction to a subset of the vectors.
      
      This implementation is a mix of four descriptions of this algorithm:
      * [S] Schnorr: "A Hierarchy of Polynomial Time Lattice Basis Reduction" (1987)
      * [SB] Schnorr: "Blockwise Lattice Basis Reduction Revisited" (2006)
      (http://www.mi.informatik.uni-frankfurt.de/research/papers/Pridu1.pdf)
      * [SL] Schnorr: "Gitter und Kryptographie"; lecture notes
      (http://www.mi.informatik.uni-frankfurt.de/teaching/lecture_notes/gitterAKTUELL.pdf)
      * [SP] Schnorr: "Progress on LLL and Lattice Reduction"
      
      In [S], first all blocks are HKZ reduced. Then, one selects the smallest i such that the
      determinant or the Siegel condition are violated. In the latter case, an LLL step is applied to
      the two vectors, and both adjacent blocks are HKZ reduced. In the former case, one applies HKZ to
      the combined block (of double size). This is repeated until the two conditions are always
      satisfies.
      
      In [SB], one first computes a HKZ basis of the first block, then applies LLL to the whole
      lattice. Then, for each ell, one first HKZ-reduces block ell+1. Then one tests whether swapping
      the adjacent vectors from the two blocks shortens the GS norm of the last vector of the ell-th
      block by at least a factor of alpha. Finally, the double block is always HKZ reduced, though the
      changes are not directly applied. (It is not totally clear to me if the double block HKZ reduction
      is always applied, or only if the swap occured.) Then, the new determinant for block ell (after
      HKZ-reduction of double block) is computed and compared to the old; if it decreased by a factor of
      at least sqrt(alpha), the changes from the double block HKZ reduction are applied and ell is
      decreased. Otherwise, the double-block HKZ-reduction changes are ignored and ell is increased.
      
      In [SL], one first computes an LLL basis of the whole lattice, and then a HKZ basis of the first
      block. Then, for each ell, one first HKZ-reduces block ell+1. Then, one applies HKZ reduction to
      the double block ell, ell+1, but does not directly apply the changes. As in [SB], they are only
      applied if the determinant of block ell decreases by a certain factor (which can be different from
      the LLL alpha). At the very end of the algorithm, size reduction is applied to the basis.
      
      In [SP], one first computes an LLL basis of the whole lattice, and then a HKZ basis of the first
      block. Then, for each ell (starting with ell=2), one first HKZ-reduces block ell+1. Then, one
      applies LLL reduction to the double block ell, ell+1, but does not directly apply the changes. If
      an LLL swap connects the two blocks, the LLL reduction is applied, ell decreased, and one restarts
      with this ell. In case no connecting swap appears, one throws away the changes and instead
      HKZ-reduces the double block ell,ell+1. As in [SB] and [SL], one only applies the changes from the
      HKZ reduction in case the determinant decreases by a certain factor; in that case, ell is
      decreased, otherwise increased. In case ell is 1, one restarts the whole algorithm (beginning with
      LLL-reducing the basis); otherwise one continues with the new value of ell.
    */
    {
        typename LatticeHelper<RealTypeContext, IntTypeContext>::RangePopper rangepopper(lattice);
        if (lattice.range().dimension() <= 1)
            return;
        if (lattice.range().dimension() % blocksize)
            verbose()(LatticeReduction::VL_Warning) << "WARNING: Dimension (" << lattice.range().dimension() << ") not divisible by blocksize (" << blocksize << ") in Semi Block 2k Reduction!";
        
        unsigned blockno = (lattice.range().dimension() + blocksize - 1) / blocksize;
        
        // Here, we first run LLL and then HKZ, as the LLL run could destroy the HKZ property of the
        // first block. [SL,SP]
        
        // LLL-reduce everything [SL,SP]
        lattice.range().dupRange();
        partialLLL(lattice, lattice.range().begin() + 1, reorder);
        // HKZ-reduce first block [S,SL,SP]
        lattice.range().addRange(getBlock(0, lattice, blocksize));
        computeHKZbasis(lattice);
        
        unsigned ell = 0;
        while (ell + 1 < blockno)
        {
            // HKZ-reduce block ell+1 [SB,SL,SP]
            lattice.range().addRange(getBlock(ell + 1, lattice, blocksize));
            computeHKZbasis(lattice);
            
            // // Swap boundary between blocks ell and ell+1? (Test for Lovász condition.) [SB]
            // if (LovaszCondition(lattice, lattice.range().begin() + (ell + 1) * blocksize))
            //     lattice.swap(lattice.range().begin() + (ell + 1) * blocksize - 1, lattice.range().begin() + (ell + 1) * blocksize);
            //
            // This swap sometimes destroys the HKZ property of the blocks. Therefore, we won't do it.
            
            // Create backup
            typename Lattice<RealTypeContext, IntTypeContext>::DuplicateStorage backup(lattice);
            
            // Compute determinant (to the fourth power)
            typename RealTypeContext::Real Dell(lattice.rc());
            std::pair<unsigned, unsigned> b = getBlock(ell, lattice, blocksize);
            lattice.update(b.second + 1);
            computeProjectedDeterminant(Dell, lattice, b);
            square(Dell, Dell);
            
            // HKZ reduce larger block
            lattice.range().addRange(getBlock2(ell, lattice, blocksize));
            computeHKZbasis(lattice);
            
            // Compute determinant (to the fourth power)
            typename RealTypeContext::Real Dellnew(lattice.rc());
            b = getBlock(ell, lattice, blocksize);
            lattice.update(b.second + 1);
            computeProjectedDeterminant(Dellnew, lattice, b);
            square(Dellnew, Dellnew);
            
            // Check whether new determinant is small enough
            if (Dellnew <= lattice.getLLLalpha() * Dell)
            {
                // Short enough: keep result and go back
                if (ell > 0)
                    --ell;
            }
            else
            {
                // Not short enough: restore backup and go to next block
                lattice.copyFrom(backup);
                ++ell;
            }
        }
    }
    
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    template<class Reorderer>
    void Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    partialBKZ_SemiBlock2k(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize, Reorderer & reorder)
    {
        NoAnnealBKZ<RealTypeContext, IntTypeContext> noanneal;
        partialBKZ_SemiBlock2k(lattice, blocksize, reorder, noanneal);
    }
    
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    template<class Reorderer, class Annealer>
    void Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    partialBKZ_PrimalDual(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize, Reorderer & reorder, Annealer & anneal)
    /*
      Applies Koy's Primal Dual BKZ Reduction to a subset of the vectors.
      
      This implementation is a mix of three descriptions of this algorithm:
      * [K] The description by Koy in the slides from 2004
      (http://www.mi.informatik.uni-frankfurt.de/research/papers/primdual.ps)
      * [SB] Schnorr: "Blockwise Lattice Basis Reduction Revisited" (2006)
      (http://www.mi.informatik.uni-frankfurt.de/research/papers/Pridu1.pdf)
      See also Schnorr: "Progress on LLL and Lattice Reduction"
      * [SL] Schnorr: "Gitter und Kryptographie"; lecture notes
      (http://www.mi.informatik.uni-frankfurt.de/teaching/lecture_notes/gitterAKTUELL.pdf)
      
      The idea to first run LLL on the whole basis and then HKZ on the first block is from [SB]. That
      each round starts by HKZ-reducing the block ell+1 is shared by [SB] and [SL], while [K] only
      computes an SVP basis.
      
      As the next step, [K] computes a DSVP basis for block ell, while [SL] computes a dual HKZ basis
      and [SB] computes a dual HKZ basis but only applies the transformation conditionally. In the three
      sources, the aim is to maximize the norm of the last GS vector in the ell-th block.
      
      In [SB] and [SL], a LLL step is applied to both blocks ell and ell+1; in case a swap appeared
      between the two blocks, the changes are taken. (In [SB], the HKZ reduction of the dual of block
      ell is only applied to the basis in this case. In [SL], the HKZ reduction of the dual is always
      applied.) In [K], one checks the Lovász condition for the two vectors where the blocks meet, and
      if it is not satisfied LLL reduction is applied to both blocks.
      
      Finally, in [SB], no size reduction is applied. In [SL], size reduction is applied only at the
      very end of the algorithm (and the algorithm never increases ell). In [K], size reduction is not
      mentioned.
      
      This implementation is close to [K], with an additional swap for the adjacent vectors before
      running LLL. Additionally, we run size reduction at the end as in [SL].
    */
    {
        typename LatticeHelper<RealTypeContext, IntTypeContext>::RangePopper rangepopper(lattice); 
        if (lattice.range().dimension() <= 1)
            return;
        if (lattice.range().dimension() % blocksize)
            verbose()(LatticeReduction::VL_Warning) << "WARNING: Dimension (" << lattice.range().dimension() << ") not divisible by blocksize (" << blocksize << ") in Primal Dual BKZ Reduction!";
        
        unsigned blockno = (lattice.range().dimension() + blocksize - 1) / blocksize;
        
        // LLL-reduce everything
        lattice.range().dupRange();
        partialLLL(lattice, lattice.range().begin() + 1, reorder);
        // HKZ-reduce first block
        if (verbose().yieldsOutput(LatticeReduction::VL_Chatter))
            verbose()(LatticeReduction::VL_Chatter) << "makeHKZ(" << 0 << ")";
        lattice.range().addRange(getBlock(0, lattice, blocksize));
        computeHKZbasis(lattice);
        
        unsigned ell = 0;
        while (ell + 1 < blockno)
        {
            // HKZ-reduce block ell+1
            if (verbose().yieldsOutput(LatticeReduction::VL_Chatter))
                verbose()(LatticeReduction::VL_Chatter) << "makeHKZ(" << ell + 1 << ")";
            lattice.range().addRange(getBlock(ell + 1, lattice, blocksize));
            computeHKZbasis(lattice);
            
            // Create backup (proceed as in [SL])
            typename Lattice<RealTypeContext, IntTypeContext>::DuplicateStorage backup(lattice);
            
            // Compute dual SVP basis for block ell
            if (verbose().yieldsOutput(LatticeReduction::VL_Chatter))
                verbose()(LatticeReduction::VL_Chatter) << "makeDSVP(" << ell << ")\n";
            std::pair<unsigned, unsigned> b = getBlock(ell, lattice, blocksize);
            lattice.range().addRange(b);
            computeDSVPbasis(lattice);
            
            // Check Lovász condition for adjacent block connecting vectors
            lattice.update(b.first + blocksize + 1);
            if (verbose().yieldsOutput(LatticeReduction::VL_Chatter))
                verbose()(LatticeReduction::VL_Chatter) << lattice.getNormSq(b.first + blocksize - 1) << " "
                                                        << lattice.getNormSq(b.first + blocksize) + square(lattice.getCoeff(b.first + blocksize, b.first + blocksize - 1)) * lattice.getNormSq(b.first + blocksize - 1)
                                                        << " (" << lattice.getCoeff(b.first + blocksize, b.first + blocksize - 1) << ", " << lattice.getNormSq(b.first + blocksize) << ")";
            if (LovaszCondition(lattice, b.first + blocksize))
            {
                // Swap vectors
                if (verbose().yieldsOutput(LatticeReduction::VL_Chatter))
                    verbose()(LatticeReduction::VL_Chatter) << "SwapAdjacent(" << ell << ", " << ell + 1 << ")";
                lattice.swap(b.first + blocksize - 1, b.first + blocksize);
                
                // LLL-reduce blocks ell and ell+1
                b = getBlock2(ell, lattice, blocksize);
                lattice.range().addRange(b);
                partialLLL(lattice, lattice.range().begin() + 1, reorder);
                
                // Keep changes and go back one block
                if (ell > 0)
                    --ell;
                else
                {
                    // HKZ-reduce first block
                    if (verbose().yieldsOutput(LatticeReduction::VL_Chatter))
                        verbose()(LatticeReduction::VL_Chatter) << "makeHKZ(" << 0 << ")";
                    lattice.range().addRange(getBlock(0, lattice, blocksize));
                    computeHKZbasis(lattice);
                }
            }
            else
            {
                // Restore from backup and go to next block
                lattice.copyFrom(backup);
                ++ell;
            }
        }
        // Do size reduction
        sizereduction(lattice);
    }
    
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    template<class Reorderer>
    void Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    partialBKZ_PrimalDual(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize, Reorderer & reorder)
    {
        NoAnnealBKZ<RealTypeContext, IntTypeContext> noanneal;
        partialBKZ_PrimalDual(lattice, blocksize, reorder, noanneal);
    }
    
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    template<class Reorderer>
    void Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    partialBKZ_PrimalDualBHO(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize, Reorderer & reorder)
    {
        typename LatticeHelper<RealTypeContext, IntTypeContext>::RangePopper rangepopper(lattice); 
        NoAnnealBKZ<RealTypeContext, IntTypeContext> noanneal;
        // Count how many times we have to double the block size
        unsigned c = 0;
        while ((blocksize & 1) == 0)
        {
            ++c;
            blocksize >>= 1;
        }
        // Compute
        for (unsigned i = 0; i <= c; ++i)
        {
            verbose()(LatticeReduction::VL_Information) << "Blocksize " << blocksize << "\n";
            lattice.range().dupRange();
            partialBKZ_PrimalDual(lattice, blocksize, reorder, noanneal);
            blocksize <<= 1;
        }
    }
    
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    template<class Reorderer, class Annealer>
    void Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    partialBKZ_Slide(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize, Reorderer & reorder, Annealer & anneal)
    // Applies Gama-Nguyen Slide Reduction to a subset of the vectors.
    {
        typename LatticeHelper<RealTypeContext, IntTypeContext>::RangePopper rangepopper(lattice); 
        if (lattice.range().dimension() <= 1)
            return;
        if (lattice.range().dimension() % blocksize)
            verbose()(LatticeReduction::VL_Warning) << "WARNING: Dimension (" << lattice.range().dimension() << ") not divisible by blocksize (" << blocksize << ") in Slide Reduction!";
        
        // Prepare watchdogs
        TransformDataChangeNotifyer<IntTypeContext> TT1, TT2;
        typename LatticeHelper<RealTypeContext, IntTypeContext>::AddNotifier addnotifier1(lattice, &TT1);
        typename LatticeHelper<RealTypeContext, IntTypeContext>::AddNotifier addnotifier2(lattice, &TT2);
        
        unsigned blockno = (lattice.range().dimension() + blocksize - 1) / blocksize;
        do
        {
            // While basis is not slide reduced
            TT2.resetChangeNotifyer();
            do
            {
                // While one of the primal conditions does not hold
                TT1.resetChangeNotifyer();
                // LLL-reduce
                lattice.range().dupRange();
                if (verbose().yieldsOutput(LatticeReduction::VL_Chatter))
                    verbose()(LatticeReduction::VL_Chatter) << "LLL" << lattice.range() << "\n";
                partialLLL(lattice, lattice.range().begin() + 1, reorder);
                // HKZ-reduce blocks
                for (unsigned b = blockno; b > 0; --b)
                {
                    // Compute index range to work on
                    std::pair<unsigned, unsigned> block = getBlock(b - 1, lattice, blocksize);
                    if (verbose().yieldsOutput(LatticeReduction::VL_Chatter))
                        verbose()(LatticeReduction::VL_Chatter) << "HKZ[" << block.first << "," << block.second << "]\n";
                    // Make sure callback is called from time to time
                    callback(lattice, block.first);
                    // Compute HKZ basis
                    lattice.range().addRange(block);
                    computeHKZbasis(lattice);
                }
            } while (TT1.wasChanged());
            // Ensure dual conditions
            for (unsigned b = blockno - 1; b > 0; --b)
            {
                // Compute index range to work on
                std::pair<unsigned, unsigned> block = getBlock(b - 1, lattice, blocksize, 1);
                if (block.first == lattice.range().end())
                    continue;
                if (verbose().yieldsOutput(LatticeReduction::VL_Chatter))
                    verbose()(LatticeReduction::VL_Chatter) << "DSVP[" << block.first << "," << block.second << "]\n";
                // Make sure callback is called from time to time
                callback(lattice, block.first);
                // Compute DSVP basis
                lattice.range().addRange(block);
                computeDSVPbasis(lattice);
            }
        } while (TT2.wasChanged());
    }
    
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    template<class Reorderer>
    void Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    partialBKZ_Slide(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize, Reorderer & reorder)
    {
        NoAnnealBKZ<RealTypeContext, IntTypeContext> noanneal;
        partialBKZ_Slide(lattice, blocksize, reorder, noanneal);
    }
    
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    typename Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::ImprovedSlideReturn
    Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    partialBKZ_ImprovedSlide_Original(Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction> & workspace,
                                      unsigned ell, Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize)
    {
        // Compute SVP bases for blocks ell and ell+1. Note that this does not changes D[ell]
        // and D[ell+1]!
        lattice.range().addRange(getBlock(ell, lattice, blocksize));
        workspace.computeSVPbasis(lattice, true);
        lattice.range().addRange(getBlock(ell + 1, lattice, blocksize));
        workspace.computeSVPbasis(lattice, true);
        
        // Now check the dual SVP condition.
        lattice.range().addRange(getBlock(ell, lattice, blocksize, 1));
        bool r = workspace.computeDSVPbasis(lattice);
        return r ? ISR_Continue : ISR_Stop;
    }
    
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    typename Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::ImprovedSlideReturn
    Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    partialBKZ_ImprovedSlide_LargerDSVP(Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction> & workspace,
                                        unsigned ell, Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize)
    {
        // Compute SVP bases for block ell+1. Note that this does not change D[ell]
        // and D[ell+1]!
        lattice.range().addRange(getBlock(ell + 1, lattice, blocksize));
        workspace.computeSVPbasis(lattice, true);
        
        // Now check the dual SVP condition for the extended block ell.
        std::pair<unsigned, unsigned> b = getBlock(ell, lattice, blocksize, 1);
        --b.first;
        lattice.range().addRange(b);
        bool r = workspace.computeDSVPbasis(lattice);
        return r ? ISR_Continue : ISR_Stop;
    }
    
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    typename Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::ImprovedSlideReturn
    Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    partialBKZ_ImprovedSlide_LargerSVP(Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction> & workspace,
                                       unsigned ell, Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize)
    {
        // Now check the dual SVP condition for the non-slided block ell.
        lattice.range().addRange(getBlock(ell, lattice, blocksize));
        bool r = workspace.computeDSVPbasis(lattice);
        
        // Compute SVP bases for block ell+1.
        std::pair<unsigned, unsigned> b = getBlock(ell + 1, lattice, blocksize, 1);
        --b.first;
        lattice.range().addRange(b);
        workspace.computeSVPbasis(lattice, true);
        
        return r ? ISR_Continue : ISR_Stop;
    }
    
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    template<class Reduction, class Reorderer, class Annealer>
    void Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    partialBKZ_ImprovedSlide(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize,
                             Reduction reduction, Reorderer & reorder, Annealer & anneal)
    // Applies Schnorr's Improved Slide Reduction to a subset of the vectors.
    {
        typename LatticeHelper<RealTypeContext, IntTypeContext>::RangePopper rangepopper(lattice); 
        if (lattice.range().dimension() <= 1)
            return;
        if (lattice.range().dimension() % blocksize)
            verbose()(LatticeReduction::VL_Warning) << "WARNING: Dimension (" << lattice.range().dimension() << ") not divisible by blocksize (" << blocksize << ") in Improved Slide Reduction!";
        
        // First do LLL reduction
        lattice.range().dupRange();
        partialLLL(lattice, lattice.range().begin() + 1, reorder);
        
        // Compute number of blocks
        unsigned blockno = (lattice.range().dimension() + blocksize - 1) / blocksize;
        if (blockno >= 2)
        {
            // In case there is only one block, directly proceed to compute HKZ basis of it
            
            // Compute determinants
            lattice.update(lattice.range().end() + 1);
            linalg::math_rowvector<typename RealTypeContext::Real> D(blockno);
            for (unsigned i = 0; i < blockno; ++i)
            {
                D[i].setContext(lattice.rc());
                computeProjectedDeterminant(D[i], lattice, getBlock(i, lattice, blocksize));
            }
            
            // Now do slide reduction rounds
            while (true)
            {
                // Find maximizing index
                unsigned ell = 0;
                typename RealTypeContext::Real curr(lattice.rc()), t(lattice.rc());
                curr = D[0] / D[1];
                for (unsigned i = 2; i < blockno; ++i)
                {
                    t = D[i - 1] / D[i];
                    if (t > curr)
                    {
                        ell = i - 1;
                        curr = t;
                    }
                }
                
                ImprovedSlideReturn ret = reduction(*this, ell, lattice, blocksize);
                if (ret == ISR_Stop)
                    // Stop?
                    break;
                
                // Polish blocks using LLL
                lattice.range().addRange(getBlock(ell, lattice, blocksize, 1));
                partialLLL(lattice, lattice.range().begin() + 1, reorder);
                lattice.range().addRange(getBlock(ell + 1, lattice, blocksize, 1));
                partialLLL(lattice, lattice.range().begin() + 1, reorder);
                
                // Update determinants
                lattice.update(getBlock(ell + 1, lattice, blocksize).second + 1);
                computeProjectedDeterminant(D[ell], lattice, getBlock(ell, lattice, blocksize));
                computeProjectedDeterminant(D[ell + 1], lattice, getBlock(ell + 1, lattice, blocksize));
            }
        }
        
        // Compute SVP basis for first block
        lattice.range().addRange(getBlock(0, lattice, blocksize));
        computeSVPbasis(lattice, true);
        
        // Do size reduction
        sizereduction(lattice);
    }
    
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    template<class Reduction, class Reorderer>
    void Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    partialBKZ_ImprovedSlide(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize,
                             Reduction reduction, Reorderer & reorder)
    {
        NoAnnealBKZ<RealTypeContext, IntTypeContext> noanneal;
        return partialBKZ_ImprovedSlide(lattice, blocksize, reduction, reorder, noanneal);
    }
    
    namespace
    {
        template<class T>
        class descending_sorter
        {
        public:
            bool operator() (const T & a, const T & b) const
            {
                return a > b;
            }
        };
    }
    
    template<class RealTypeContext>
    long computeLog2SuccessProbability(std::vector<typename RealTypeContext::Real> & lv, unsigned umax, typename RealTypeContext::Real & scaledfirstnorm, unsigned dim, RealTypeContext & rc)
    {
        // Sort first dim - umax + 1 lengths
        std::sort(lv.begin(), lv.begin() + dim - umax, descending_sorter<typename RealTypeContext::Real>());
        // Compute the success probability for k = 0, ..., dim - umax - 1 and takes the maximum
        long max = 0, cur;
        // There is no need to initialize max, but we want to avoid compiler warnings...
        typename RealTypeContext::Real tmp(rc), q(rc);
        for (unsigned k = 0; k < dim - umax; ++k)
        {
            // Compute return value of LogSuccessProbBound(lv, k, umax, gamma)
            
            // We want to solve the equation:
            //     1/12 * \sum_{i=0}^{k-1} q^{k-i} lv[i]
            //   + 1/12 * \sum_{i=k}^{lv.size() - umax - 2} lv[i])
            //   + 1/3 * \sum_{i=lv.size() - umax - 1}^{lv.size() - 2} lv[i]
            //   + lv[lv.size() - 1]
            //   = scaledfirstnorm.
            // This is equivalent to solve:
            //     \sum_{i=0}^{k-1} q^{k-i} lv[i]
            //   = 12 * scaledfirstnorm
            //   - 12 * lv[lv.size() - 1]
            //   - 4 * \sum_{i=lv.size() - umax - 1}^{lv.size() - 2} lv[i]
            //   - \sum_{i=k}^{lv.size() - umax - 2} lv[i]).
            // We first compute the right-hand side; then solving this is easier.
            typename RealTypeContext::Real rhs_plus(rc), rhs_minus(rc);
            arithmetic::convert(rhs_plus, 12, rc);
            rhs_plus *= scaledfirstnorm;
            // First take lv[lv.size() - 1], the only (negative) term with factor 12
            rhs_minus = lv[lv.size() - 1];
            arithmetic::convert(tmp, 3, rc);
            rhs_minus *= tmp; // multiply by 3 = 12/4
            // Take care of terms with factor 4
            for (unsigned i = lv.size() - umax - 1; i < lv.size() - 1; ++i)
                rhs_minus += lv[i];
            arithmetic::convert(tmp, 4, rc);
            rhs_minus *= tmp; // multiply by 4
            // Now take care of the terms with factor 1
            for (unsigned i = k; i < lv.size() - umax - 1; ++i)
                rhs_minus += lv[i];
            // Now our equation is:   \sum_{i=0}^{k-1} q^{k-i} lv[i]  +  rhs_minus == rhs_plus
            
            // First check the equation for q == 1. If q == 1, then \sum_{i=0}^{k-1} q^{k-i}
            // lv[i] == \sum_{i=0}^{k-1} lv[i].
            setZero(tmp);
            for (unsigned i = 0; i < k; ++i)
                tmp += lv[i];
            tmp += rhs_minus;
            if (tmp <= rhs_plus)
                cur = -1;
            else
            {
                // Next check the equation for q == 0. If q == 0, then \sum_{i=0}^{k-1} q^{k-i}
                // lv[i] == 0.
                if (rhs_minus >= rhs_plus)
                {
                    cur = -(long)umax - 1; // should be -infinity, but this is enough :)
                }
                else
                {
                    // Now find q such that
                    //     \sum_{i=0}^{k-1} q^{k-i} lv[i]  +  rhs_minus == rhs_plus
                    // using regula falsi. Define
                    //     f(q) = \sum_{i=0}^{k-1} q^{k-i} lv[i] + rhs_minus - rhs_plus.
                    // Then we know that f(1) > 0 and f(0) < 0.
                    
                    typename RealTypeContext::Real lower(rc), upper(rc), middle(rc);
                    setZero(lower);
                    setOne(upper);
                    for (unsigned it = 0; it < 60; ++it)
                        // Do 60 iterations
                    {
                        // Compute middle position
                        middle = lower + upper;
                        middle >>= 1;
                        // Evaluate at middle position
                        setZero(tmp);
                        setOne(q);
                        for (unsigned i = 0; i < k; ++i)
                        {
                            tmp += lv[i];
                            tmp *= middle;
                        }
                        tmp += rhs_minus;
                        // Now we have to compare tmp to rhs_plus
                        int c = compare(tmp, rhs_plus);
                        if (c < 0)
                            lower = middle;
                        else if (c > 0)
                            upper = middle;
                        else
                        {
                            // We hit the target! Unlikely, but can happen...
                            lower = upper = middle;
                            break;
                        }
                    }
                    // Take q as the element in the middle of the interval
                    q = lower + upper;
                    q >>= 1;
                    // Compute log_2 of success probability
                    log2(q, q);
                    arithmetic::convert(tmp, k * (k + 1), rc);
                    tmp >>= 2;
                    q *= tmp;
                    --q;
                    cur = arithmetic::convert_floor<long>(q);
                }
            }
            // Compare
            if ((k == 0) || (cur > max))
                max = cur;
        }
        return max;
    }
    
    class IsFP { };
    class IsNotFP { };
    
    template<class RealTypeContext, class IntTypeContext>
    long computeLog2SuccessProbability(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned umax,
                                       typename RealTypeContext::Real & scaledfirstnorm, IsFP)
    {
        std::vector<typename RealTypeContext::Real> lv;
        lv.resize(lattice.range().dimension());
        for (unsigned i = 0; i < lv.size(); ++i)
        {
            lv[i].setContext(lattice.rc());
            lv[i] = lattice.getNormSq(lattice.range().begin() + i);
        }
        return computeLog2SuccessProbability<RealTypeContext>(lv, umax, scaledfirstnorm, lattice.range().dimension(), lattice.rc());
    }
    
    template<class RealTypeContext, class IntTypeContext>
    long computeLog2SuccessProbability(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned umax,
                                       typename RealTypeContext::Real & scaledfirstnorm, IsNotFP)
    {
        arithmetic::RealContext rc;
        std::vector<arithmetic::Real> lv;
        arithmetic::Real scaledfirstnorm_(rc);
        scaledfirstnorm_ = arithmetic::convert(scaledfirstnorm, rc);
        lv.resize(lattice.range().dimension());
        for (unsigned i = 0; i < lv.size(); ++i)
        {
            lv[i].setContext(rc);
            lv[i] = arithmetic::convert(lattice.getNormSq(lattice.range().begin() + i), rc);
        }
        return computeLog2SuccessProbability<arithmetic::RealContext>(lv, umax, scaledfirstnorm_, lattice.range().dimension(), rc);
    }
    
    template<class RealTypeContext, class IntTypeContext>
    long computeLog2SuccessProbability(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned umax,
                                       typename RealTypeContext::Real & scaledfirstnorm)
    {
        return computeLog2SuccessProbability(lattice, umax, scaledfirstnorm,
                                             typename helper::SelectFirstType<RealTypeContext::has_special_fns, IsFP, IsNotFP>::result());
    }
    
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    template<class Reorderer, class Annealer>
    void Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    partialBKZ_SamplingReduction(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize, Reorderer & reorder, Annealer & anneal)
    // Applies Buchmann and Ludwig's Sampling Reduction to a subset of the vectors. For details, see
    // Buchmann, Ludwig: "Practical Lattice Basis Sampling Reduction".
    {
        typename LatticeHelper<RealTypeContext, IntTypeContext>::RangePopper rangepopper(lattice); 
        if (lattice.range().dimension() <= 1)
            return;
        
        typename RealTypeContext::Real gamma(lattice.rc());
        arithmetic::convert(gamma, 0.99, lattice.rc());
        unsigned umax = 30;
        if (umax > lattice.range().dimension())
            umax = lattice.range().dimension();
        
        while (true)
        {
            verbose()(LatticeReduction::VL_Information) << "Starting SamplingReduction iteration";
            // Run BKZ
            lattice.range().dupRange();
            partialBKZ_Classical(lattice, blocksize, reorder, anneal);
            
            // Compute GS coefficients
            lattice.update(lattice.range().end() + 1);
            
            // Compute scaled first norm
            typename RealTypeContext::Real scaledfirstnorm(lattice.rc());
            scaledfirstnorm = lattice.getNormSq(lattice.range().begin());
            scaledfirstnorm *= gamma;
            
            // Check best bound condition
            long max = computeLog2SuccessProbability(lattice, umax, scaledfirstnorm);
            if (-max > (long)umax)
            {
                verbose()(LatticeReduction::VL_Information) << "Sampling Reduction termination reason: success probability too small.";
                break;
            }
            
            // Do sampling
            verbose()(LatticeReduction::VL_Information) << "  Do sampling (floor(log_2(success probability)) = " << max << ")...";
            unsigned count = 1 << umax;
            bool inserted_vector = false;
            typename IntTypeContext::Integer x;
            typename RealTypeContext::Real t(lattice.rc());
            linalg::math_rowvector<typename IntTypeContext::Integer> v;
            linalg::math_rowvector<typename RealTypeContext::Real> w;
            v.resize(lattice.range().dimension());
            w.resize(lattice.range().dimension());
            for (unsigned i = 0; i < w.size(); ++i)
                w[i].setContext(lattice.rc());
            for (unsigned ell = 1; ell <= count; ++ell)
            {
                // Callback function call
                callback(lattice, lattice.range().begin());
                // Sample vector
                for (unsigned i = 0; i + 1 < w.size(); ++i)
                {
                    w[i] = lattice.getCoeff(lattice.range().end(), lattice.range().begin() + i);
                    setZero(v[i]);
                }
                setOne(w[w.size() - 1]);
                setOne(v[v.size() - 1]);
                unsigned l = ell;
                for (unsigned i = w.size() - 1; i > 0; --i)
                {
                    // Use vector lattice.range().begin() + i - 1
                    bool rounded_up;
                    arithmetic::convert_round(x, w[i - 1], rounded_up, lattice.ic());
                    // Is bit 0 set in l?
                    if (l & 1)
                    {
                        // Yes: modify x by +-1
                        if (rounded_up)
                            --x;
                        else
                            ++x;
                    }
                    // Go to next bit
                    l >>= 1;
                    // Subtract x times the (lattice.range().begin()+i-1)-th basis vector
                    v[i - 1] = -x;
                    if (isPMOne(x))
                    {
                        if (sign(x) > 0)
                            for (unsigned j = 0; j + 1 < i; ++j)
                                w[j] -= lattice.getCoeff(lattice.range().begin() + i - 1, lattice.range().begin() + j);
                        else
                            for (unsigned j = 0; j + 1 < i; ++j)
                                w[j] += lattice.getCoeff(lattice.range().begin() + i - 1, lattice.range().begin() + j);
                    }
                    else
                    {
                        arithmetic::convert(t, x, lattice.rc());
                        for (unsigned j = 0; j + 1 < i; ++j)
                            w[j] -= t * lattice.getCoeff(lattice.range().begin() + i - 1, lattice.range().begin() + j);
                        // We can ignore w[i-1] since this won't be needed again
                    }
                }
                // Check if its squared norm is shorter than the squared norm of the first basis vector
                // times gamma.
                lattice.computeProjectionLength(t, lattice.range().begin(), lattice.range().begin(), v);
                if (t < scaledfirstnorm)
                {
                    lattice.insertVectorLC(lattice.range().begin(), v);
                    inserted_vector = true;
                    break;
                }
            }
            
            // Did we insert something?
            if (!inserted_vector)
            {
                // No: stop
                verbose()(LatticeReduction::VL_Information) << "Sampling Reduction termination reason: search space exhausted.";
                break;
            }
        }
    }
    
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    template<class Reorderer>
    void Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    partialBKZ_SamplingReduction(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize, Reorderer & reorder)
    {
        NoAnnealBKZ<RealTypeContext, IntTypeContext> noanneal;
        partialBKZ_SamplingReduction(lattice, blocksize, reorder, noanneal);
    }
    
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    template<class Reorderer, class Annealer>
    void Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    partialBKZ_Experimental(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize,
                            Reorderer & reorder, Annealer & anneal, unsigned id)
    {
        typename LatticeHelper<RealTypeContext, IntTypeContext>::RangePopper rangepopper(lattice); 
        switch(id)
        {
        case 0:
        // Some experimental block reduction code. This is reversed simplified BKZ with SVP steps.
        {
            // First do LLL reduction
            lattice.range().dupRange();
            partialLLL(lattice, lattice.range().begin() + 1, reorder);
            
            TransformDataChangeNotifyer<IntTypeContext> TT;
            typename LatticeHelper<RealTypeContext, IntTypeContext>::AddNotifier addnotifier(lattice, &TT);
            bool avoidInsertions = lattice.avoidInsertions();
            linalg::math_rowvector<typename IntTypeContext::Integer> result;
            do
            {
                TT.resetChangeNotifyer();
                for (signed k = lattice.range().end(); k >= (signed)lattice.range().begin(); --k)
                {
                    callback(lattice, k);
                    
                    // Do reorder (Deep Insertion)
                    unsigned kkk = k;
                    reorder(*this, lattice, kkk, kkk, true, false);
                    k = kkk;
                    
                    // Find shortest vector
                    unsigned kk = std::min(lattice.range().end(), k + blocksize - 1);
                    lattice.range().addRange(k, kk);
                    solveSVP(lattice, result, reorder);
                    
                    bool doreplace;
                    if (result.size() > 0)
                    {
                        // Compute length of projection
                        typename RealTypeContext::Real l(lattice.rc());
                        lattice.computeProjectionLength(l, k, k, result);
                        
                        // Check Schnorr-Euchner BKZ condition
                        doreplace = (lattice.getLLLalpha() * lattice.getNormSq(k) > l);
                    }
                    else
                        // Didn't find anything
                        doreplace = false;
                    
                    // Do annealing
                    if (!doreplace)
                        doreplace = anneal(lattice, k, blocksize, result);
                    
                    if (doreplace)
                    {
                        // Yes? Include new vector
                        if (avoidInsertions)
                        {
                            lattice.range().addRange(k, std::min(kk + 1, lattice.range().end()));
                            rearrange(lattice, result);
                            lattice.range().addRange(lattice.range().begin(), std::min(kk + 1, lattice.range().end()));
                            unsigned ms = partialLLL(lattice, k, reorder);
                            // Go back to last change
                            if ((signed)ms < k)
                                k = ms - 1;
                        }
                        else
                        {
                            unsigned j = lattice.range().begin();
                            lattice.range().addRange(k, std::min(kk + 1, lattice.range().end()));
                            unsigned ms = rearrangeLLL(lattice, j, result, reorder);
                            // Go back to last change
                            if ((signed)ms < k)
                                k = ms - 1;
                        }
                    }
                    else
                    {
                        // Apply LLL
                        lattice.range().addRange(lattice.range().begin(), std::min(kk + 1, lattice.range().end()));
                        unsigned ms = partialLLL(lattice, kk, reorder);
                        // Go back to last change
                        if ((signed)ms <= k)
                            k = ms - 1;
                    }
                }
            } while (TT.wasChanged());
        }
        break;
        case 1:
        // Some experimental block reduction code. This is reversed simplified BKZ with SVP steps and dual enums to.
        {
            // First do LLL reduction
            lattice.range().dupRange();
            partialLLL(lattice, lattice.range().begin() + 1, reorder);
            
            TransformDataChangeNotifyer<IntTypeContext> TT;
            typename LatticeHelper<RealTypeContext, IntTypeContext>::AddNotifier addnotifier(lattice, &TT);
            bool avoidInsertions = lattice.avoidInsertions();
            
            linalg::math_rowvector<typename IntTypeContext::Integer> result;
            do
            {
                TT.resetChangeNotifyer();
//            for (signed k = lattice.range().end(); k >= (signed)lattice.range().begin(); --k)
                for (signed k = lattice.range().begin(); k <= (signed)lattice.range().end(); ++k)
                {
                    callback(lattice, k);
                    
                    // Do reorder (Deep Insertion)
                    unsigned ku = k;
                    reorder(*this, lattice, ku, ku, true, false);
                    k = ku;
                    
                    unsigned
                        kk = std::min(lattice.range().end(), k + blocksize - 1),
                        l = kk,
                        ll = std::min(lattice.range().end(), l + blocksize - 1);
                    
                    if (ll > kk)
                    {
                        // compute DHKZ basis for block [l, ll]
                        lattice.range().addRange(l, ll);
                        computeDHKZbasis(lattice);
                    }
                    
                    // Find shortest vector
                    lattice.range().addRange(k, kk);
                    solveSVP(lattice, result, reorder);
                    
                    bool doreplace;
                    if (result.size() > 0)
                    {
                        // Compute length of projection
                        typename RealTypeContext::Real l(lattice.rc());
                        lattice.computeProjectionLength(l, k, k, result);
                        
                        // Check Schnorr-Euchner BKZ condition
                        doreplace = (lattice.getLLLalpha() * lattice.getNormSq(k) > l);
                    }
                    else
                        // Didn't find anything
                        doreplace = false;
                    
                    // Do annealing
                    if (!doreplace)
                        doreplace = anneal(lattice, k, blocksize, result);
                    
                    if (doreplace)
                    {
                        // Yes? Include new vector
                        if (avoidInsertions)
                        {
                            lattice.range().addRange(k, std::min(kk + 1, lattice.range().end()));
                            rearrange(lattice, result);
                            lattice.range().addRange(lattice.range().begin(), std::min(kk + 1, lattice.range().end()));
                            unsigned ms = partialLLL(lattice, k, reorder);
                            // Go back to last change
                            if ((signed)ms < k)
                                k = ms - 1;
                        }
                        else
                        {
                            unsigned j = lattice.range().begin();
                            lattice.range().addRange(k, std::min(kk + 1, lattice.range().end()));
                            unsigned ms = rearrangeLLL(lattice, j, result, reorder);
                            // Go back to last change
                            if ((signed)ms < k)
                                k = ms - 1;
                        }
                    }
                    else
                    {
                        // Apply LLL
                        lattice.range().addRange(lattice.range().begin(), std::min(kk + 1, lattice.range().end()));
                        unsigned ms = partialLLL(lattice, kk, reorder);
                        // Go back to last change
                        if ((signed)ms <= k)
                            k = ms - 1;
                    }
                }
            } while (TT.wasChanged());
        }
        break;
        default:
            return;
        }
    }
    
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    template<class Reorderer>
    void Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    partialBKZ_Experimental(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize, Reorderer & reorder, unsigned id)
    {
        NoAnnealBKZ<RealTypeContext, IntTypeContext> noanneal;
        partialBKZ_Experimental(lattice, blocksize, reorder, noanneal, id);
    }
    
    template<class AnnealFunction, class IntTypeContext>
    class AnnealerBKZ
    {
    private:
        AnnealFunction d_af;
        arithmetic::RealContext & d_rc;
        arithmetic::RandomNumberGenerator & d_rng;
        arithmetic::Real & d_T;
        
    public:
        AnnealerBKZ(AnnealFunction af, arithmetic::RealContext & rc, arithmetic::RandomNumberGenerator & rng, arithmetic::Real & T)
            : d_af(af), d_rc(rc), d_rng(rng), d_T(T)
        {
        }
        
        template<class Lattice>
        bool operator() (Lattice & lattice, int k, int windowsize, linalg::math_rowvector<typename IntTypeContext::Integer> & lincomb)
        {
            std::pair<linalg::math_matrix<typename IntTypeContext::Integer> *, unsigned> ll = lattice.getMatrixIndex(k);
            return d_af(lattice.ic(), *ll.first, ll.second, windowsize, lincomb, d_rc, d_rng, d_T);
        }
    };
    
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    template<class AnnealFunction, class Reorderer>
    void Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    partialBKZ_Anneal(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize,
                      LatticeReduction::AnnealCallbackFunction ACF, AnnealFunction AF,
                      Reorderer & reorder, bool simplified)
    // Applies BKZ to a subset of the vectors.
    {
        typename LatticeHelper<RealTypeContext, IntTypeContext>::RangePopper rangepopper(lattice); 
        if (ACF == NULL)
            ACF = &DefaultAnnealCallbackFunction;
        
        arithmetic::RealContext arc;
        arithmetic::Real t(arc);
        setOne(t);
        t = -t;
        arithmetic::RandomNumberGenerator rng;
        rng.randomizeSeed();
        
        lattice.range().dupRange();
        partialLLL(lattice, lattice.range().begin() + 1, reorder);
        assert(lattice.canGetIntegerVector());
        while (true)
        {
            {
                MatrixConversionConst2<IntTypeContext> m(lattice.ic(), *lattice.getMatrix());
                if (!ACF(m.matrix(), arc, t))
                    break;
            }
            verbose()(LatticeReduction::VL_Information) << "Temperature " << t;
            arithmetic::Real tt(arc);
            tt = t;
            AnnealerBKZ<AnnealFunction, IntTypeContext> annealer(AF, arc, rng, tt);
            lattice.range().dupRange();
            if (simplified)
                partialBKZ_Simplified(lattice, blocksize, reorder, annealer);
            else
                partialBKZ_Classical(lattice, blocksize, reorder, annealer);
        }
    }
}

#endif

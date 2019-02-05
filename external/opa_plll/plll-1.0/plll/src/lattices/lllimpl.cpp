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

#ifndef PLLL_INCLUDE_GUARD__LLL_IMPL_CPP
#define PLLL_INCLUDE_GUARD__LLL_IMPL_CPP

namespace plll
{
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    void Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::sizereduction(Lattice<RealTypeContext, IntTypeContext> & lattice)
    // Applies size reduction: A.row(begin) to A.row(end) will be size-reduced.
    {
        for (unsigned dest = lattice.range().begin(); dest <= lattice.range().end(); ++dest)
        {
            lattice.update(dest + 1);
            lattice.sizereduce(dest);
        }
    }
    
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    bool Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    LovaszCondition(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned k)
    {
        return lattice.getLLLalpha() * lattice.getNormSq(k - 1) > lattice.getNormSq(k) + square(lattice.getCoeff(k, k - 1)) * lattice.getNormSq(k - 1);
    }
    
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    bool Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    UnprojectedLovaszCondition(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned k)
    {
        return lattice.getNormSqUP(k - 1) > lattice.getNormSqUP(k);
    }
    
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    bool Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    SiegelCondition(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned k)
    {
        return convert(3, lattice.rc()) * lattice.getLLLalpha() * lattice.getNormSq(k - 1) > convert(4, lattice.rc()) * lattice.getNormSq(k);
    }
    
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    template<class Reorderer, class Annealer, class LCondition>
    unsigned Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    partialLLL(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned beginstage,
               Reorderer & reorder, Annealer & anneal, LCondition & lcondition)
    // Applies LLL to a subset of the vectors.
    {
        typename LatticeHelper<RealTypeContext, IntTypeContext>::RangePopper rangepopper(lattice);
        if (lattice.range().dimension() == 0)
            return lattice.range().begin() + 1;
        // Initialize variables
        unsigned k = beginstage; // stage: b_begin, ..., b_{k-1} are reduced
        unsigned minstage = k;
        // Main loop
        bool jump = false;
        if (k > 0)
        {
            --k;
            jump = true;
        }
        // Remove zero vectors from the beginning: these will never be found by size reduction
        while (lattice.range().dimension() > 0)
        {
            if (lattice.isRowZero(lattice.range().begin()))
                lattice.removeZeroVector(lattice.range().begin());
            else
                break;
        }
        while (k <= lattice.range().end())
        {
            // Call callback function
            callback(lattice, k);
            // Make sure enough GS coefficients are computed
            lattice.update(k + 1);
            // Do reorder (Deep Insertion)
            if (reorder(*this, lattice, k, k, true, false))
            {
                if (k < lattice.range().begin() + 1)
                {
                    k = lattice.range().begin();
                    jump = true;
                }
                continue;
            }
            // Do size reduction (if possible)
            bool src = lattice.sizereduce(k);
            if (lattice.isRowZero(k))
            {
                // Remove zero vector
                lattice.removeZeroVector(k);
                continue;
            }
            // Do reorder again (Deep Insertion)
            if (reorder(*this, lattice, k, k, false, src))
            {
                if (k < lattice.range().begin() + 1)
                {
                    k = lattice.range().begin();
                    jump = true;
                }
                continue;
            }
            // Do jump to valid stage (if requested)
            if (jump)
            {
                ++k;
                jump = false;
                continue;
            }
            // Update minimum stage
            if (minstage > k)
                minstage = k;
            // Ensure Lovász condition
            bool doswap = lcondition(lattice, k);
            // Do annealing
            if (!doswap)
                doswap = anneal(lattice, k);
            // Do swap?
            if (doswap)
            {
                // Swap k and k-1
                lattice.swap(k, k - 1);
                // Go back
                if (k > lattice.range().begin() + 1)
                    --k;
            }
            else
            {
                // Condition is satisfied
                ++k;
            }
        }
        return minstage;
    }
    
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    template<class Reorderer>
    unsigned Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    partialLLL(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned beginstage, Reorderer & reorder)
    {
        NoAnnealLLL<RealTypeContext, IntTypeContext> anneal;
        return partialLLL(lattice, beginstage, reorder, anneal,
                          Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::LovaszCondition);
    }
    
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    template<class Reorderer>
    unsigned Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    partialLLL2(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned beginstage, Reorderer & reorder)
    {
        NoAnnealLLL<RealTypeContext, IntTypeContext> anneal;
        return partialLLL(lattice, beginstage, reorder, anneal,
                          Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::UnprojectedLovaszCondition);
    }
    
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    template<class Reorderer>
    unsigned Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    partialLLL3(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned beginstage, Reorderer & reorder)
    {
        NoAnnealLLL<RealTypeContext, IntTypeContext> anneal;
        return partialLLL(lattice, beginstage, reorder, anneal,
                          Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::SiegelCondition);
    }
    
    template<class AnnealFunction, class IntTypeContext>
    class AnnealerLLL
    {
    private:
        AnnealFunction d_af;
        arithmetic::RealContext & d_rc;
        arithmetic::RandomNumberGenerator & d_rng;
        arithmetic::Real & d_T;
    
    public:
        AnnealerLLL(AnnealFunction af, arithmetic::RealContext & rc, arithmetic::RandomNumberGenerator & rng, arithmetic::Real & T)
            : d_af(af), d_rc(rc), d_rng(rng), d_T(T)
        {
        }
    
        template<class Lattice>
        bool operator() (Lattice & lattice, int k)
        {
            std::pair<linalg::math_matrix<typename IntTypeContext::Integer> *, unsigned> ll = lattice.getMatrixIndex(k);
            return d_af(lattice.ic(), *ll.first, ll.second, d_rc, d_rng, d_T);
        }
    };

    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    template<class AnnealFunction, class Reorderer>
    void Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    partialLLL_Anneal(Lattice<RealTypeContext, IntTypeContext> & lattice, LatticeReduction::AnnealCallbackFunction ACF, AnnealFunction AF, Reorderer & reorder)
    // Applies LLL to a subset of the vectors.
    {
        typename LatticeHelper<RealTypeContext, IntTypeContext>::RangePopper rangepopper(lattice); 
        if (ACF.empty())
            ACF = &DefaultAnnealCallbackFunction;
        arithmetic::RealContext arc;
        arithmetic::Real t(arc);
        setOne(t);
        t = -t;
        arithmetic::RandomNumberGenerator rng;
        rng.randomizeSeed();
    
        NoAnnealLLL<RealTypeContext, IntTypeContext> noanneal;
        lattice.range().dupRange();
        partialLLL(lattice, lattice.range().begin() + 1, reorder, noanneal,
                   Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::LovaszCondition);
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
            AnnealerLLL<AnnealFunction, IntTypeContext> annealer(AF, arc, rng, tt);
            lattice.range().dupRange();
            partialLLL(lattice, lattice.range().begin() + 1, reorder, annealer,
                       Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::LovaszCondition);
        }
    }
}

#endif

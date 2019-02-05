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

#ifndef PLLL_INCLUDE_GUARD__ENUMIMPL_BOUNDSEL_CPP
#define PLLL_INCLUDE_GUARD__ENUMIMPL_BOUNDSEL_CPP

namespace plll
{
    template<class RealTypeContext, class IntTypeContext>
    void StandardEnumBoundSelection(GaussianFactorComputer & gfc, Lattice<RealTypeContext, IntTypeContext> & lattice,
                                    typename RealTypeContext::Real & bound)
    {
        typename RealTypeContext::Real t(lattice.rc());
        bound = lattice.getNormSq(lattice.range().begin());
        // Compute lengths of projected basis vectors
        for (unsigned i = lattice.range().begin() + 1; i <= lattice.range().end(); ++i)
        {
            lattice.computeProjectionLengthBV(t, lattice.range().begin(), i);
            // If we find a shorter length, take that one
            if (t < bound)
                bound = t;
        }
    }

    template<class RealTypeContext, class IntTypeContext, bool hasFullPower>
    class GHEBSImpl;

    template<class RealTypeContext, class IntTypeContext>
    class GHEBSImpl<RealTypeContext, IntTypeContext, true>
    {
    public:
        void operator() (GaussianFactorComputer & gfc, Lattice<RealTypeContext, IntTypeContext> & lattice,
                         typename RealTypeContext::Real & bound, typename RealTypeContext::Real factor) const
        {
            typename RealTypeContext::Real det(lattice.rc());
            // Compute determinant (more precisely: its square)
            setOne(det);
            for (unsigned i = lattice.range().begin(); i != lattice.range().end(); ++i)
                det *= lattice.getNormSq(i);
            // If determinant is zero: fallback
            if (isZero(det))
            {
                StandardEnumBoundSelection(gfc, lattice, bound);
                return;
            }
            // Compute dim-th root of determinant and multiply with Hermite prefactor
            det = power(det, convert(1, lattice.rc()) / convert(2*lattice.range().dimension(), lattice.rc()));
            convert(bound, gfc.HPF(lattice.range().dimension()), lattice.rc());
            bound *= factor;
            bound *= det;
            // Compare with "traditional" bound:
            StandardEnumBoundSelection(gfc, lattice, det); // we ignore naming conventions and store the length in "det"
            if (bound > det)
                bound = det;
        }
    };

    template<class RealTypeContext, class IntTypeContext>
    class GHEBSImpl<RealTypeContext, IntTypeContext, false>
    {
    public:
        void operator() (GaussianFactorComputer & gfc, Lattice<RealTypeContext, IntTypeContext> & lattice,
                         typename RealTypeContext::Real & bound, typename RealTypeContext::Real factor) const
        {
            arithmetic::RealContext rrc;
            arithmetic::RealContext::Real det(rrc);
            // Compute determinant (more precisely: its square)
            setOne(det);
            for (unsigned i = lattice.range().begin(); i != lattice.range().end(); ++i)
                det *= convert(lattice.getNormSq(i), rrc);
            // If determinant is zero: fallback
            if (isZero(det))
            {
                StandardEnumBoundSelection(gfc, lattice, bound);
                return;
            }
            // Compute dim-th root of determinant and multiply with Hermite prefactor
            det = power(det, convert(1, rrc) / convert(2*lattice.range().dimension(), rrc));
            arithmetic::convert(bound, gfc.HPF(lattice.range().dimension()), lattice.rc());
            bound *= factor;
            bound *= convert(det, lattice.rc());
            // Compare with "traditional" bound:
            typename RealTypeContext::Real t(lattice.rc());
            StandardEnumBoundSelection(gfc, lattice, t);
            if (bound > t)
                bound = t;
        }
    };

    template<class RealTypeContext, class IntTypeContext>
    void GaussianHeuristicEnumBoundSelection<RealTypeContext, IntTypeContext>::operator()
        (GaussianFactorComputer & gfc, Lattice<RealTypeContext, IntTypeContext> & lattice, typename RealTypeContext::Real & bound) const
    {
        GHEBSImpl<RealTypeContext, IntTypeContext, RealTypeContext::has_full_power>()(gfc, lattice, bound, d_factor);
    }
}

#endif

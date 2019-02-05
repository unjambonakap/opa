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

#ifndef PLLL_INCLUDE_GUARD__GRAMSCHMIDT_CLASSIC_CPP
#define PLLL_INCLUDE_GUARD__GRAMSCHMIDT_CLASSIC_CPP

#include "gramschmidt-generic.cpp"

namespace plll
{
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Gram-Schmidt computers: classical Gram-Schmidt
    
    template<class RealTypeContext, class IntTypeContext>
    class GSClassic
    // classical Gram-Schmidt orthogonalization
    {
    private:
        RealTypeContext & d_rc;
        IntTypeContext & d_ic;
        linalg::math_matrix<typename IntTypeContext::Integer> & d_A;
    
        bool d_recompute;
    
        unsigned d_level, d_vectors_used;
        linalg::math_matrix<typename RealTypeContext::Real> d_coeffs;
        linalg::math_rowvector<typename RealTypeContext::Real> d_sqnorms;
        linalg::base_rowvector<bool> d_dependent;
    
        void doGS(unsigned start, unsigned stop)
        {
            // Apply (modified) Gram-Schmidt.
            typename RealTypeContext::Real t(d_rc), norm(d_rc);
            for (unsigned i = start; i <= stop; ++i)
            {
                // Start with norm (unprojected)
                convert(norm, normSq(d_A.row(i)), d_rc);
                // Compute coefficients and projected norm
                for (unsigned j = 0; j < i; ++j)
                {
                    if (d_dependent[j])
                        setZero(d_coeffs(i, j));
                    else
                    {
                        // Compute coefficient
                        convert(t, dot(d_A.row(i), d_A.row(j)), d_rc);
                        for (unsigned k = 0; k < j; ++k)
                            t -= d_coeffs(j, k) * d_coeffs(i, k) * d_sqnorms[k];
                        t /= d_sqnorms[j];
                        d_coeffs(i, j) = t;
                    
                        // Continue with computation of squared norm
                        square(t, t);
                        t *= d_sqnorms[j];
                        norm -= t;
                    }
                }
            
                // Store norm
                d_sqnorms[i] = norm;
                // Is dependent on the previous vectors?
                d_dependent[i] = isZero(norm);
            }
        }
        
        std::deque<std::pair<typename RealTypeContext::Real *, typename RealTypeContext::Real *> > d_precupdater_stack;
        
        mutable typename IntTypeContext::Integer d_tmp;
        
        Verbose * d_verbose;
        
        inline void adjustZerodotfive(typename RealTypeContext::Real & zerodotfive, const typename RealTypeContext::Real & ralpha, helper::BoolToType<true>)
        // Adjust the value of "zerodotfive", which should be 0.5 for exact arithmetic, and slightly
        // larger for floating point approximations
        {
            zerodotfive += d_rc.getEpsilon();
        }
        
        inline void adjustZerodotfive(typename RealTypeContext::Real & zerodotfive, const typename RealTypeContext::Real & ralpha, helper::BoolToType<false>)
        // Adjust the value of "zerodotfive", which should be 0.5 for exact arithmetic, and slightly
        // larger for floating point approximations
        {
        }
        
    public:
        class DuplicateStorage;
        friend class DuplicateStorage;
        
        class DuplicateStorage
        {
        private:
            // Prohibit copying
            DuplicateStorage(const DuplicateStorage &);
            DuplicateStorage & operator = (const DuplicateStorage &);
            
        public:
            inline DuplicateStorage(const GSClassic & gs, linalg::math_matrix<typename IntTypeContext::Integer> & m, RealTypeContext & rc, IntTypeContext & ic)
                : d_gs(gs, m, rc, ic)
            {
            }
        
            GSClassic d_gs;
        };
        
        void copyFrom(const GSClassic & gs)
        {
            d_level = gs.d_level;
            d_vectors_used = gs.d_vectors_used;
            d_coeffs = gs.d_coeffs;
            d_sqnorms = gs.d_sqnorms;
            d_dependent = gs.d_dependent;
            d_precupdater_stack = gs.d_precupdater_stack;
        }
        
        inline void copyFrom(const DuplicateStorage & storage)
        {
            copyFrom(storage.d_gs);
        }
        
        GSClassic(const GSClassic & source, RealTypeContext & rc, IntTypeContext & ic, linalg::math_matrix<typename IntTypeContext::Integer> & A)
            : d_rc(rc), d_ic(ic), d_A(A), d_recompute(source.d_recompute), d_level(source.d_level), d_vectors_used(source.d_vectors_used),
              d_coeffs(source.d_coeffs), d_sqnorms(source.d_sqnorms), d_dependent(source.d_dependent),
              d_precupdater_stack(source.d_precupdater_stack), d_verbose(source.d_verbose)
        {
        }
        
        GSClassic(Verbose & verbose, RealTypeContext & rc, IntTypeContext & ic,
                  linalg::math_matrix<typename IntTypeContext::Integer> & A, bool recompute = false)
            : d_rc(rc), d_ic(ic), d_A(A), d_recompute(recompute), d_level(0), d_vectors_used(d_A.rows()), d_verbose(&verbose)
        {
            d_coeffs.resize(A.rows(), A.rows(), linalg::Initialize(d_rc));
            d_sqnorms.resize(A.rows(), linalg::Initialize(d_rc));
            d_dependent.resize(A.rows());
            for (unsigned i = 0; i < A.rows(); ++i)
                d_dependent[i] = false;
        }
        
        inline linalg::math_matrix<typename RealTypeContext::Real> * getCoeffsStorage()
        {
            return &d_coeffs;
        }
        
        inline linalg::math_rowvector<typename RealTypeContext::Real> * getSqNormsStorage()
        {
            return &d_sqnorms;
        }
        
        inline void registerZDFAlpha(typename RealTypeContext::Real & zerodotfive, typename RealTypeContext::Real & ralpha)
        // Informs the Gram-Schmidt computer about the current values of zerodotfive and ralpha used. In
        // case the precision is adjusted, the RealTypeContext's of these values will be reset as well.
        {
            d_precupdater_stack.push_back(std::make_pair(&zerodotfive, &ralpha));
        }
        
        inline void unregisterZDFAlpha()
        // Removes the previously added values of zerodotfive and ralpha from the precision updating stack.
        {
            d_precupdater_stack.pop_back();
        }
        
        inline void adjustZerodotfive(typename RealTypeContext::Real & zerodotfive, const typename RealTypeContext::Real & ralpha)
        // Adjust the value of "zerodotfive", which should be 0.5 for exact arithmetic, and slightly
        // larger for floating point approximations
        {
            adjustZerodotfive(zerodotfive, ralpha, helper::BoolToType<!RealTypeContext::is_exact && RealTypeContext::has_constants>());
        }
        
        inline const typename IntTypeContext::Integer & getNormSqUP(unsigned i) const
        {
            normSq(d_tmp, d_A.row(i));
            return d_tmp;
        }
        
        inline void getNormSqUP(typename IntTypeContext::Integer & r, unsigned i) const
        {
            normSq(r, d_A.row(i));
        }
        
        inline void getDotProduct(typename IntTypeContext::Integer & r, unsigned i, unsigned j) const
        {
            dot(r, d_A.row(i), d_A.row(j));
        }
        
        void swap(unsigned i, unsigned j)
        {
            typename RealTypeContext::Real t(d_rc);
            if (i < j)
                std::swap(i, j);
            if (j + 1 < i)
            {
                // More work to do! Just recompute...
                if (d_level > j)
                    d_level = j;
                return;
            }
            if (d_level <= i)
            {
                // Early exit: just recompute the GS coefficients
                if (d_level > j)
                    d_level = j;
                return;
            }
            GSGeneric::swap(d_rc, d_ic, d_coeffs, d_sqnorms, d_dependent, d_level, i, j);
        }
        
        void add(unsigned i, const typename IntTypeContext::Integer & m, unsigned j)
        {
            if (i < j)
            {
                // Wrong order. Recompute!
                if (d_level > i)
                    d_level = i;
                return;
            }
            if (d_level <= i)
                return;
            GSGeneric::add(d_rc, d_ic, d_coeffs, d_sqnorms, d_dependent, d_level, i, m, j);
        }
        
        void flip(unsigned i)
        {
            if (d_level <= i)
                return;
            GSGeneric::flip(d_rc, d_ic, d_coeffs, d_sqnorms, d_dependent, d_level, i);
        }
        
        void trans(unsigned k, unsigned ell,
                   const typename IntTypeContext::Integer & B00, const typename IntTypeContext::Integer & B01,
                   const typename IntTypeContext::Integer & B10, const typename IntTypeContext::Integer & B11)
        {
            if (k + 1 != ell)
            {
                // Wrong order or distance: recompute
                if (d_level > std::min(k, ell))
                    d_level = std::min(k, ell);
                return;
            }
            if (d_dependent[k] || d_dependent[ell] || (d_level <= ell))
            {
                // Early exit: just recompute the GS coefficients
                if (d_level > k)
                    d_level = k;
                return;
            }
            GSGeneric::trans(d_rc, d_ic, d_coeffs, d_sqnorms, d_dependent, d_level, k, ell, B00, B01, B10, B11);
        }
        
        void update(unsigned level)
        {
            assert(level <= d_vectors_used);
            if (level <= d_level)
                return;
            unsigned start = d_level;
            if (d_recompute && !RealTypeContext::is_exact) // no need to recompute if exact arithmetic is used
                start = 0;
            
            // Apply Gram-Schmidt
            doGS(start, level - 1);
            
            d_level = level;
        }
        
        void reset()
        {
            d_level = 0;
        }
        
        template<class Lattice> // stores information on the transformations done to A
        bool sizereduce(Lattice & lattice, typename RealTypeContext::Real & zerodotfive, typename RealTypeContext::Real & ralpha, unsigned begin, unsigned end, unsigned dest)
        // zerodotfive: something in [0.5, 1)
        //
        // Applies size reduction: A.row(begin) to A.row(end) will be used to size-reduce A.row(dest).
        // Assumes begin <= end < dest.
        {
            return GSGeneric::sizereduce(d_rc, d_ic, d_coeffs, d_sqnorms, d_dependent, d_level, lattice, zerodotfive, begin, end, dest);
        }
        
        void computeDotProductProjected(typename RealTypeContext::Real & r, unsigned k, unsigned i, unsigned j) const
        {
            GSGeneric::computeDotProductProjected(d_rc, d_ic, d_coeffs, d_sqnorms, d_dependent, d_level, r, k, i, j);
        }
        
        void computeProjectionLength(typename RealTypeContext::Real & result, unsigned k, unsigned b,
                                     const linalg::math_rowvector<typename IntTypeContext::Integer> & vec) const
        // Computes the squared length of the projection of \sum_{i=0}^{m-1} vec[i] B_{b+i} onto the
        // orthogonal complement of B_0, ..., B_{k-1}. (Here, B_0, ..., B_{n-1} is the basis of the
        // lattice.)
        {
            GSGeneric::computeProjectionLength(d_rc, d_ic, d_coeffs, d_sqnorms, d_dependent, d_level, result, k, b, vec);
        }
        
        void computeProjectionLengthBV(typename RealTypeContext::Real & result, unsigned k, unsigned b) const
        // Computes the squared length of the projection of B_b onto the orthogonal complement of B_0, ...,
        // B_{k-1}. (Here, B_0, ..., B_{n-1} is the basis of the lattice.)
        {
            GSGeneric::computeProjectionLengthBV(d_rc, d_ic, d_coeffs, d_sqnorms, d_dependent, d_level, result, k, b);
        }
        
        void computeProjectionLength(typename RealTypeContext::Real & result, unsigned k,
                                     const linalg::math_rowvector<typename IntTypeContext::Integer> & vec) const
        // Computes the squared length of the projection of vec onto the orthogonal complement of B_0,
        // ..., B_{k-1}. (Here, B_0, ..., B_{n-1} is the basis of the lattice.)
        {
            GSGeneric::computeProjectionLength(d_rc, d_ic, d_coeffs, d_sqnorms, d_dependent, d_A, d_level, result, k, vec);
        }
        
        void insertVectorLC(unsigned ofs, const linalg::math_rowvector<typename IntTypeContext::Integer> & result)
        // Inserts a new vector at index ofs, which will equal \sum_{i=0}^{result.size()-1} result[i] * b[ofs+i]
        {
            GSGeneric::insertVector(d_rc, d_ic, d_coeffs, d_sqnorms, d_dependent, d_vectors_used, ofs, false);
            if (d_level > ofs)
                d_level = ofs;
        }
        
        void insertVector(unsigned ofs, const linalg::math_rowvector<typename IntTypeContext::Integer> & result)
        // Inserts a new vector <result> at index ofs
        {
            GSGeneric::insertVector(d_rc, d_ic, d_coeffs, d_sqnorms, d_dependent, d_vectors_used, ofs, false);
            if (d_level > ofs)
                d_level = ofs;
        }
        
        void removeZeroVector(unsigned ofs)
        // Removes the zero vector at position ofs
        {
            if (ofs < d_level)
            {
                GSGeneric::removeVector(d_rc, d_ic, d_coeffs, d_sqnorms, d_dependent, d_vectors_used, d_level, ofs);
                --d_level;
            }
            else
                --d_vectors_used;
        }
        
        static void compactify()
        {
            // GSGeneric::compactify(d_rc, d_ic, d_coeffs, d_sqnorms, d_dependent, d_vectors_used);
        }
        
        void changeOfPrecision()
        // Notify the Gram-Schmidt computation engine that the precision of d_rc changed
        {
            for (unsigned i = 0; i < d_coeffs.rows(); ++i)
            {
                for (unsigned j = 0; j < d_coeffs.rows(); ++j)
                    d_coeffs(i, j).setContext(d_rc);
                d_sqnorms[i].setContext(d_rc);
            }
            if (d_level > 0)
                doGS(0, d_level - 1);
            for (typename std::deque<std::pair<typename RealTypeContext::Real *, typename RealTypeContext::Real *> >::iterator i = d_precupdater_stack.begin(); i != d_precupdater_stack.end(); ++i)
            {
                i->first->setContext(d_rc);
                i->second->setContext(d_rc);
            }
        }
    };
}

#endif

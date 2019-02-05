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

#ifndef PLLL_INCLUDE_GUARD__GRAMSCHMIDT_NUMSTABLELLL_CPP
#define PLLL_INCLUDE_GUARD__GRAMSCHMIDT_NUMSTABLELLL_CPP

#include <iostream>
#include <iomanip>
#include <deque>

#include "gramschmidt-generic.cpp"

namespace plll
{
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Gram-Schmidt computers: incremental numerically stable Gram-Schmidt for LLL reductions

    template<class RealTypeContext, class IntTypeContext>
    class GSNumStableLLL
    // Incremental numerically stable Gram-Schmidt for LLL reductions as described in "Floating-Point
    // LLL Revisited" by P. Q. Nguyen and D. Stehle, 2005.
    {
    private:
        RealTypeContext & d_rc;
        IntTypeContext & d_ic;
        linalg::math_matrix<typename IntTypeContext::Integer> & d_A;
    
        unsigned d_level, d_computed, d_G_computed, d_relaxed, d_vectors_used;
        linalg::math_matrix<typename IntTypeContext::Integer> d_G;
        linalg::math_matrix<typename RealTypeContext::Real> d_r;
        linalg::math_matrix<typename RealTypeContext::Real> d_coeffs;
        linalg::math_rowvector<typename RealTypeContext::Real> d_sqnorms;
        linalg::math_matrix<typename RealTypeContext::Real> d_s;
        linalg::base_rowvector<bool> d_dependent;
    
        mutable typename IntTypeContext::Integer d_tmp; // for getNormSqUP(unsigned)
    
        void doGS(unsigned start, unsigned stop)
        {
            // Compute dot products if necessary
            if (stop >= d_G_computed)
            {
                for (unsigned ii = d_G_computed; ii <= stop; ++ii)
                    for (unsigned jj = 0; jj <= ii; ++jj)
                        dot(d_G(ii, jj), d_A.row(ii), d_A.row(jj));
                d_G_computed = stop + 1;
            }
            typename RealTypeContext::Real t(d_rc);
            for (unsigned i = start; i <= stop; ++i)
            {
                // Orthogonalize vector
                for (unsigned j = 0; j <= i; ++j)
                {
                    // Compute dot product
                    if (d_dependent[j])
                    {
                        setZero(d_r(i, j));
                        setZero(d_coeffs(i, j));
                    }
                    else
                    {
                        arithmetic::convert(d_r(i, j), d_G(i, j), d_rc);
                        for (unsigned k = 0; k < j; ++k)
                        {
                            t = d_r(i, k) * d_coeffs(j, k);
                            d_r(i, j) -= t;
                        }
                        if (j < i)
                            d_coeffs(i, j) = d_r(i, j) / d_r(j, j);
                        else
                            setOne(d_coeffs(i, i));
                    }
                }
                // Compute lengths of projections
                arithmetic::convert(d_s(i, 0), d_G(i, i), d_rc);
                for (unsigned j = 0; j < i; ++j)
                {
                    t = d_coeffs(i, j) * d_r(i, j);
                    d_s(i, j + 1) = d_s(i, j) - t;
                }
                d_sqnorms[i] = d_r(i, i) = d_s(i, i);
                // Detect dependences
                d_dependent[i] = isZero(d_sqnorms[i]);
                if (d_dependent[i])
                {
                    setZero(d_sqnorms[i]);
                    setZero(d_r(i, i));
                    setZero(d_s(i, i));
                }
            }
        
            d_computed = stop + 1;
        }

        inline void add_G(unsigned i, const typename IntTypeContext::Integer & m, unsigned j)
        {
            if (i > j)
            {
                // Usual add (happening in LLL etc.)
                if (isPMOne(m))
                {
                    if (sign(m) > 0)
                    {
                        // Update G(i, i)
                        d_G(i, i) += d_G(j, j);
                        d_G(i, i) += d_G(i, j);
                        d_G(i, i) += d_G(i, j);
                        assert(d_G(i, i) == dot(d_A.row(i), d_A.row(i)));
                        // Update G(k, i) for k > i
                        for (unsigned k = i + 1; k < d_G_computed; ++k)
                            d_G(k, i) += d_G(k, j);
                        // Update G(i, k) for k < i
                        for (unsigned k = 0; k < j; ++k)
                            d_G(i, k) += d_G(j, k);
                        for (unsigned k = j; k < i; ++k)
                            d_G(i, k) += d_G(k, j);
                    }
                    else
                    {
                        // Update G(i, i)
                        d_G(i, i) += d_G(j, j);
                        d_G(i, i) -= d_G(i, j);
                        d_G(i, i) -= d_G(i, j);
                        assert(d_G(i, i) == dot(d_A.row(i), d_A.row(i)));
                        // Update G(k, i) for k > i
                        for (unsigned k = i + 1; k < d_G_computed; ++k)
                            d_G(k, i) -= d_G(k, j);
                        // Update G(i, k) for k < i
                        for (unsigned k = 0; k < j; ++k)
                            d_G(i, k) -= d_G(j, k);
                        for (unsigned k = j; k < i; ++k)
                            d_G(i, k) -= d_G(k, j);
                    }
                }
                else
                {
                    // Update G(i, i)
                    d_G(i, i) += d_G(j, j) * square(m);
                    d_G(i, i) += d_G(i, j) * (m << 1);
                    assert(d_G(i, i) == dot(d_A.row(i), d_A.row(i)));
                    // Update G(k, i) for k > i
                    for (unsigned k = i + 1; k < d_G_computed; ++k)
                        d_G(k, i) += m * d_G(k, j);
                    // Update G(i, k) for k < i
                    for (unsigned k = 0; k < j; ++k)
                        d_G(i, k) += m * d_G(j, k);
                    for (unsigned k = j; k < i; ++k)
                        d_G(i, k) += m * d_G(k, j);
                }
            }
            else
            {
                // "Wrong" add (happening when randomizing basis f.ex.), where j > i. This one is not
                // optimized for m = +-1.
            
                // Update G(i, i)
                d_G(i, i) += d_G(j, j) * square(m);
                d_G(i, i) += d_G(j, i) * (m << 1);
                // Update G(k, i) for k > i
                for (unsigned k = i + 1; k < j; ++k)
                    d_G(k, i) += m * d_G(j, k);
                for (unsigned k = j; k < d_G_computed; ++k)
                    d_G(k, i) += m * d_G(k, j);
                // Update G(i, k) for k < i
                for (unsigned k = 0; k < i; ++k)
                    d_G(i, k) += m * d_G(j, k);
            }
        }
    
        void verifyG()
        // Debug function
        {
            bool total_ok = true;
            for (unsigned i = 0; i < d_G_computed; ++i)
                for (unsigned j = 0; j <= i; ++j)
                {
                    bool ok = dot(d_A.row(i), d_A.row(j)) == d_G(i, j);
                    if (!ok)
                    {
                        std::cout << i << " " << j << ": " << dot(d_A.row(i), d_A.row(j)) << " vs. " << d_G(i, j) << "\n";
                        total_ok = false;
                    }
                }
            assert(total_ok);
        }
    
        std::deque<std::pair<typename RealTypeContext::Real *, typename RealTypeContext::Real *> > d_precupdater_stack;
    
        Verbose * d_verbose;
    
        void do_add_row()
        {
            // Add row
            d_r.resize(d_r.rows() + 1, d_r.cols() + 1, linalg::Initialize(d_rc));
            d_coeffs.resize(d_coeffs.rows() + 1, d_coeffs.cols() + 1, linalg::Initialize(d_rc));
            d_sqnorms.resize(d_sqnorms.size() + 1, linalg::Initialize(d_rc));
            d_s.resize(d_s.rows() + 1, d_s.cols() + 1, linalg::Initialize(d_rc));
            d_dependent.resize(d_dependent.size() + 1);
        }
        
        void adjustZerodotfive(typename RealTypeContext::Real & zerodotfive, const typename RealTypeContext::Real & ralpha, helper::BoolToType<false>)
        // Adjust the value of "zerodotfive", which should be 0.5 for exact arithmetic, and slightly
        // larger for floating point approximations
        {
        }
        
        void adjustZerodotfive(typename RealTypeContext::Real & zerodotfive, const typename RealTypeContext::Real & ralpha, helper::BoolToType<true>)
        // Adjust the value of "zerodotfive", which should be 0.5 for exact arithmetic, and slightly
        // larger for floating point approximations
        {
            if (d_relaxed)
                zerodotfive += d_rc.getEpsilon() << ((signed long)d_relaxed - 1);
            while (true)
            {
                if (zerodotfive * zerodotfive >= ralpha)
                    (*d_verbose)(LatticeReduction::VL_Warning) << "WARNING: zerodotfive < sqrt(alpha) not satisfied!";
                else
                {
                    arithmetic::Real_precision_check_disable();
                    arithmetic::RealContext rc;
                    unsigned long prec = arithmetic::convert_ceil<long>(log(arithmetic::convert(square(arithmetic::convert(1, d_rc) + zerodotfive)
                                                                                                / (ralpha - square(zerodotfive)), rc))
                                                                        * arithmetic::convert(d_vectors_used + 1, rc));
                    arithmetic::Real_precision_check_enable();
                    if ((unsigned long)prec > d_rc.getRealPrecision())
                    {
                        (*d_verbose)(LatticeReduction::VL_Warning) << "WARNING: precision is " << d_rc.getRealPrecision()
                                                                   << ", but should better be at least " << prec << "";
                        if (RealTypeContext::is_variable_precision)
                        {
                            (*d_verbose)(LatticeReduction::VL_Information) << "Adjusting precision to " << prec << "...";
                            d_relaxed += prec - d_rc.getRealPrecision();
                            d_rc.setRealPrecision(prec);
                            zerodotfive.setContext(d_rc);
                            changeOfPrecision();
                            continue; // repeat loop
                        }
                    }
                }
                return;
            }
        }
    
    public:
        class DuplicateStorage;
        friend class DuplicateStorage;
    
        class DuplicateStorage
        {
            friend class GSNumStableLLL;
        
        private:
            // Prohibit copying
            DuplicateStorage(const DuplicateStorage &);
            DuplicateStorage & operator = (const DuplicateStorage &);
        
        public:
            DuplicateStorage(const GSNumStableLLL & gs, linalg::math_matrix<typename IntTypeContext::Integer> & m, RealTypeContext & rc, IntTypeContext & ic)
                : d_gs(gs, rc, ic, m)
            {
            }
        
            GSNumStableLLL d_gs;
        };
        
        void copyFrom(const GSNumStableLLL & gs)
        {
            d_level = gs.d_level;
            d_computed = gs.d_computed;
            d_G_computed = gs.d_G_computed;
            d_relaxed = gs.d_relaxed;
            d_vectors_used = gs.d_vectors_used;
            d_G = gs.d_G;
            d_r = gs.d_r;
            d_coeffs = gs.d_coeffs;
            d_sqnorms = gs.d_sqnorms;
            d_s = gs.d_s;
            d_dependent = gs.d_dependent;
            d_precupdater_stack = gs.d_precupdater_stack;
        }
        
        inline void copyFrom(const DuplicateStorage & storage)
        {
            copyFrom(storage.d_gs);
        }
        
        GSNumStableLLL(Verbose & verbose, RealTypeContext & rc, IntTypeContext & ic, linalg::math_matrix<typename IntTypeContext::Integer> & A)
            : d_rc(rc), d_ic(ic), d_A(A), d_level(0), d_computed(0), d_G_computed(0), d_relaxed(d_rc.is_exact ? 0 : 1), d_vectors_used(d_A.rows()), d_verbose(&verbose)
        {
            d_s.resize(A.rows(), A.rows(), linalg::Initialize(d_rc));
            d_coeffs.resize(A.rows(), A.rows(), linalg::Initialize(d_rc));
            d_r.resize(A.rows(), A.rows(), linalg::Initialize(d_rc));
            d_G.resize(A.rows(), A.rows(), linalg::Initialize(d_ic));
            d_sqnorms.resize(A.rows(), linalg::Initialize(d_rc));
            d_dependent.resize(A.rows());
            for (unsigned i = 0; i < A.rows(); ++i)
                d_dependent[i] = false;
        }
        
        GSNumStableLLL(const GSNumStableLLL & source, RealTypeContext & rc, IntTypeContext & ic, linalg::math_matrix<typename IntTypeContext::Integer> & A)
            : d_rc(rc), d_ic(ic), d_A(A), d_level(source.d_level), d_computed(source.d_computed),
              d_G_computed(source.d_G_computed), d_relaxed(source.d_relaxed),
              d_vectors_used(source.d_vectors_used), d_G(source.d_G),
              d_r(source.d_r), d_coeffs(source.d_coeffs), d_sqnorms(source.d_sqnorms),
              d_s(source.d_s), d_dependent(source.d_dependent),
            d_precupdater_stack(source.d_precupdater_stack), d_verbose(source.d_verbose)
        {
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
            adjustZerodotfive(zerodotfive, ralpha, helper::BoolToType<RealTypeContext::has_constants && !RealTypeContext::is_exact>());
        }
        
        inline const typename IntTypeContext::Integer & getNormSqUP(unsigned i) const
        {
            if (d_G_computed > i)
                return d_G(i, i);
            else
            {
                normSq(d_tmp, d_A.row(i));
                return d_tmp;
            }
        }
    
        inline void getNormSqUP(typename IntTypeContext::Integer & r, unsigned i) const
        {
            if (d_G_computed > i)
            {
                r = d_G(i, i);
                assert(sign(r) >= 0);
            }
            else
            {
                normSq(r, d_A.row(i));
                assert(sign(r) >= 0);
            }
        }
    
        inline void getDotProduct(typename IntTypeContext::Integer & r, unsigned i, unsigned j) const
        {
            if ((d_G_computed > i) && (d_G_computed > j))
                r = d_G(std::max(i, j), std::min(i, j));
            else
                dot(r, d_A.row(i), d_A.row(j));
        }
    
        void swap(unsigned i, unsigned j)
        {
            typename RealTypeContext::Real t(d_rc);
            if (i < j)
                std::swap(i, j);
            if (d_computed > j)
                d_computed = j;
            // Update G(., .)
            if (d_G_computed <= i)
            {
                if (d_G_computed > j)
                    d_G_computed = j;
            }
            else
            {
                arithmetic::swap(d_G(i, i), d_G(j, j));
                for (unsigned k = 0; k < j; ++k)
                {
                    arithmetic::swap(d_G(i, k), d_G(j, k));
                }
                for (unsigned k = j + 1; k < i; ++k)
                {
                    arithmetic::swap(d_G(i, k), d_G(k, j));
                }
                for (unsigned k = i + 1; k < d_G_computed; ++k)
                {
                    arithmetic::swap(d_G(k, i), d_G(k, j));
                }
            }
        }
    
        void add(unsigned i, const typename IntTypeContext::Integer & m, unsigned j)
        {
            if (d_computed > i)
                d_computed = i;
            if (d_G_computed > i)
            {
                if (d_G_computed > j)
                    add_G(i, m, j);
                else
                    d_G_computed = i;
            }
        }
    
        void flip(unsigned i)
        {
            for (unsigned k = 0; k < i; ++k)
            {
                neg(d_coeffs(i, k), d_coeffs(i, k));
                neg(d_r(i, k), d_r(i, k));
            }
            for (unsigned k = i + 1; k < d_computed; ++k)
            {
                neg(d_coeffs(k, i), d_coeffs(k, i));
                neg(d_r(k, i), d_r(k, i));
            }
            // Update G(., .)
            if (d_G_computed > i)
            {
                for (unsigned k = 0; k < i; ++k)
                {
                    neg(d_G(i, k), d_G(i, k));
                }
                for (unsigned k = i + 1; k < d_G_computed; ++k)
                {
                    neg(d_G(k, i), d_G(k, i));
                }
            }
        }
    
        void trans(unsigned k, unsigned ell,
                   const typename IntTypeContext::Integer & B00, const typename IntTypeContext::Integer & B01,
                   const typename IntTypeContext::Integer & B10, const typename IntTypeContext::Integer & B11)
        {
            if (k > ell)
            {
                trans(ell, k, B01, B00, B11, B10);
                return;
            }
            // Now we are sure that k < ell
            if (d_computed > k)
                d_computed = k;
            // Update G(., .)
            if (d_G_computed > k)
            {
                if ((d_G_computed <= ell) || (k + 1 != ell))
                    d_G_computed = k;
                else
                {
                    typename IntTypeContext::Integer x, y, t1, t2;
                    // compute new d_G(k, k)
                    square(x, B00);
                    x *= d_G(k, k);
                    square(y, B01);
                    y *= d_G(ell, ell);
                    x += y;
                    y = d_G(ell, k);
                    y <<= 1;
                    y *= B00;
                    y *= B01;
                    t1 = x + y;
                    // compute new d_G(ell, ell)
                    square(x, B10);
                    x *= d_G(k, k);
                    square(y, B11);
                    y *= d_G(ell, ell);
                    x += y;
                    y = d_G(ell, k);
                    y <<= 1;
                    y *= B10;
                    y *= B11;
                    t2 = x + y;
                    // compute new d_G(ell, k)
                    x = B00 * B10;
                    x *= d_G(k, k);
                    y = B01 * B11;
                    y *= d_G(ell, ell);
                    x += y;
                    y = B00 * B11;
                    y += B01 * B10;
                    y *= d_G(ell, k);
                    // Store results
                    d_G(ell, k) = x + y;
                    d_G(k, k) = t1;
                    d_G(ell, ell) = t2;
//                assert(d_G(k, k) == dot(d_A.row(k), d_A.row(k)));
//                assert(d_G(ell, k) == dot(d_A.row(k), d_A.row(ell)));
//                assert(d_G(ell, ell) == dot(d_A.row(ell), d_A.row(ell)));
                    // Now update other entries
                    for (unsigned i = 0; i < k; ++i)
                    {
                        x = B00 * d_G(k, i);
                        x += B01 * d_G(ell, i);
                        y = B10 * d_G(k, i);
                        y += B11 * d_G(ell, i);
                        d_G(k, i) = x;
                        d_G(ell, i) = y;
                    }
                    for (unsigned i = k + 1; i < ell; ++i)
                    {
                        x = B00 * d_G(i, k);
                        x += B01 * d_G(ell, i);
                        y = B10 * d_G(i, k);
                        y += B11 * d_G(ell, i);
                        d_G(i, k) = x;
                        d_G(ell, i) = y;
                    }
                    for (unsigned i = ell + 1; i < d_G_computed; ++i)
                    {
                        x = B00 * d_G(i, k);
                        x += B01 * d_G(i, ell);
                        y = B10 * d_G(i, k);
                        y += B11 * d_G(i, ell);
                        d_G(i, k) = x;
                        d_G(i, ell) = y;
                    }
                }
            }
        }
    
        void print()
        {
            std::cout << "Computed: " << d_computed << "\n";
            for (unsigned i = 0; i < d_computed; ++i)
            {
                std::cout << i << ":";
                for (unsigned j = 0; j < i; ++j)
                {
                    std::cout << " " << d_coeffs(i, j) << " [" << d_r(i, j) << "]";
                }
                std::cout << " -> " << d_sqnorms[i] << " [" << d_r(i, i) << "]\n";
            }
            std::cout << d_A << "\n";
        }
    
        void update(unsigned level)
        {
//        print();
        
            assert(level <= d_vectors_used);
            if ((level <= d_level) && (level <= d_computed))
                return;
            unsigned start = d_computed;
        
            // Apply Gram-Schmidt
            if (level > start)
                doGS(start, level - 1);
        
            d_level = level;
        }
    
        void reset()
        {
            d_level = 0;
            d_computed = 0;
        }
    
        template<class Lattice> // stores information on the transformations done to A
        bool sizereduce(Lattice & lattice, typename RealTypeContext::Real & zerodotfive, typename RealTypeContext::Real & ralpha, unsigned begin, unsigned end, unsigned dest)
        // zerodotfive: something in [0.5, 1)
        //
        // Applies size reduction: A.row(begin) to A.row(end) will be used to size-reduce A.row(dest).
        // Assumes begin <= end < dest.
        {
            typename RealTypeContext::Real etatilde(d_rc), t(d_rc), x(d_rc);
            bool changed = false;
//        unsigned bigiters = 0;
            while (true)
            {
//            ++bigiters;
                setOne(etatilde);
                etatilde >>= 1;
                etatilde += zerodotfive;
                etatilde >>= 1;
                typename IntTypeContext::Integer X;
                if (dest >= d_computed)
                    doGS(d_computed, dest);
                unsigned iterations = 0, extra_iterations = 0;
                bool relax = false;
//            if (bigiters > 1)
//            {
//                std::cout << "One big round of size reduction of #" << dest << ": eta = " << std::setprecision(20) << etatilde << "\n";
//                print();
//            }
                long last_exp = std::numeric_limits<long>::max();
                while (true)
                {
                    long exp = d_coeffs(dest, begin).getApproxExponent();
                    for (unsigned i = begin + 1; i <= end; ++i)
                        exp = std::max(exp, d_coeffs(dest, i).getApproxExponent());
                    // We want to gain at least 5 bits
                    if (exp > last_exp - 5)
                    {
                        if ((exp > 2) || (extra_iterations < 3))
                            ++extra_iterations;
                        else
                        {
                            relax = true;
                            break;
                        }
                    }
                    last_exp = exp;
                    ++iterations;
                    // Apply one iteration of Babai's Nearest Plane Algorithm
                    bool allzero = true;
                    for (signed i = end; i >= (signed)begin; --i)
                    {
                        if (arithmetic::abs(d_coeffs(dest, i)) >= etatilde)
                        {
                            allzero = false;
                            arithmetic::convert_round(X, d_coeffs(dest, i), d_ic);
                            neg(X, X);
                            arithmetic::convert(x, X, d_rc);
                            for (unsigned j = 0; j < (unsigned)i; ++j)
                            {
                                t = x * d_coeffs(i, j);
                                d_coeffs(dest, j) += t;
                            }
                            lattice.add(dest, X, i);
                            changed = true;
                        }
                    }
                    // In case nothing changed, we are done!
                    if (allzero)
                        break;
                    // Update GS coefficients
                    doGS(dest, dest);
//            print();
                }
                // Relax precision?
                if (!relax)
                    break;
                ++d_relaxed;
                adjustZerodotfive(zerodotfive, ralpha);
                (*d_verbose)(LatticeReduction::VL_Warning) << "WARNING: increasing relaxiation level to " << d_relaxed << " (" << std::setprecision(20) << zerodotfive << ")!";
                typename RealTypeContext::Real one(d_rc);
                setOne(one);
                if (zerodotfive >= one)
                    // Too big!
                    throw reduction_error("Relaxiation level too high, aborting!");
            }
            if (d_level > dest)
                d_level = dest + 1;
            return changed;
        }
    
        void computeDotProductProjected(typename RealTypeContext::Real & r, unsigned k, unsigned i, unsigned j) const
        {
            if ((i < k) || (j < k))
            {
                // Projection is 0
                setZero(r);
                return;
            }
            if (i < j)
                std::swap(i, j);
            typename RealTypeContext::Real t(d_rc);
            if (i == j)
            {
                r = d_sqnorms[j];
                for (unsigned ell = k; ell < j; ++ell)
                {
                    square(t, d_coeffs(i, ell));
                    r += t * d_sqnorms[ell];
                }
            }
            else
            {
                r = d_r(i, j);
                for (unsigned ell = k; ell < j; ++ell)
                    r += d_coeffs(i, ell) * d_r(j, ell);
            }
        }
    
        void computeProjectionLength(typename RealTypeContext::Real & result, unsigned k, unsigned b,
                                     const linalg::math_rowvector<typename IntTypeContext::Integer> & vec) const
        // Computes the squared length of the projection of \sum_{i=0}^{m-1} vec[i] B_{b+i} onto the
        // orthogonal complement of B_0, ..., B_{k-1}. (Here, B_0, ..., B_{n-1} is the basis of the
        // lattice.)
        {
            typename RealTypeContext::Real t(d_rc);
            setZero(result);
            unsigned beg = k < b ? b : k;
            unsigned end = b + vec.size() - 1;
            if ((end < k) || (vec.size() == 0))
                // Projection is 0
                return;
            for (unsigned i = k; i < beg; ++i)
            {
                setZero(t);
                unsigned begg = i < beg ? beg : i + 1;
                for (unsigned j = begg; j <= end; ++j)
                    t += arithmetic::convert(vec[j - b], d_rc) * d_coeffs(j, i);
                square(t, t);
                result += t * d_sqnorms[i];
            }
            for (unsigned i = beg; i <= end; ++i)
            {
                arithmetic::convert(t, vec[i - b], d_rc);
                for (unsigned j = i + 1; j <= end; ++j)
                    t += arithmetic::convert(vec[j - b], d_rc) * d_coeffs(j, i);
                square(t, t);
                result += t * d_sqnorms[i];
            }
        }
    
        void computeProjectionLengthBV(typename RealTypeContext::Real & result, unsigned k, unsigned b) const
        // Computes the squared length of the projection of B_b onto the orthogonal complement of B_0, ...,
        // B_{k-1}. (Here, B_0, ..., B_{n-1} is the basis of the lattice.)
        {
            if (b < k)
                // Projection is 0
                setZero(result);
            else
                result = d_s(b, k);
        }
    
        void computeProjectionLength(typename RealTypeContext::Real & result, unsigned k,
                                     const linalg::math_rowvector<typename IntTypeContext::Integer> & vec) const
        // Computes the squared length of the projection of vec onto the orthogonal complement of B_0,
        // ..., B_{k-1}. (Here, B_0, ..., B_{n-1} is the basis of the lattice.)
        {
            typename RealTypeContext::Real t(d_rc);
            // Compute squared norm (unprojected)
            arithmetic::convert(result, normSq(vec), d_rc);
            if (k == 0)
                // In case we just need the norm, return
                return;
            // Orthogonalize vector
            linalg::math_rowvector<typename RealTypeContext::Real> r(k);
            for (unsigned j = 0; j < k; ++j)
                if (!d_dependent[j])
                {
                    r[j].setContext(d_rc);
                    // Compute dot product <B_j, vec>
                    arithmetic::convert(r[j], dot(d_A.row(j), vec), d_rc);
                    for (unsigned k = 0; k < j; ++k)
                    {
                        t = r[k] * d_coeffs(j, k);
                        r[j] -= t;
                    }
                    // Compute length of projection iteratively
                    t = r[j] / d_r(j, j); // this would be coeffs[j] in the notation of doGS()
                    t *= r[j];
                    result = result - t;
                }
        }
    
        void insertVectorLC(unsigned ofs, const linalg::math_rowvector<typename IntTypeContext::Integer> & result)
        // Inserts a new vector at index ofs, which will equal \sum_{i=0}^{result.size()-1} result[i] * b[ofs+i]
        {
            bool add_row = d_vectors_used == d_G.rows();
            if (d_computed > ofs)
                d_computed = ofs;
            if (d_level > ofs)
                d_level = ofs;
            if (add_row)
                d_G.resize(d_G.rows() + 1, d_G.cols() + 1, linalg::Initialize(d_ic));
            if (d_G_computed > ofs)
            {
                if (d_G_computed >= ofs + result.size())
                {
                    // Insert row
                    for (unsigned i = d_G_computed; i > ofs; --i)
                        linalg::swap(d_G.row(i), d_G.row(i - 1));
                    // Insert column
                    for (unsigned i = d_G_computed; i > ofs; --i)
                        linalg::swap(d_G.col(i), d_G.col(i - 1));
                    // Add new data column
                    for (unsigned i = 0; i < ofs; ++i)
                    {
                        setZero(d_G(ofs, i));
                        for (unsigned j = 0; j < result.size(); ++j)
                            d_G(ofs, i) += result[j] * d_G(ofs + j + 1, i);
//                    assert(d_G(ofs, i) == dot(d_A.row(i), d_A.row(ofs)));
                    }
                    for (unsigned i = ofs + 1; i <= d_G_computed; ++i)
                    {
                        setZero(d_G(i, ofs));
                        for (unsigned j = 0; j < result.size(); ++j)
                            if (i <= ofs + j)
                                d_G(i, ofs) += result[j] * d_G(ofs + j + 1, i);
                            else
                                d_G(i, ofs) += result[j] * d_G(i, ofs + j + 1);
//                    assert(d_G(i, ofs) == dot(d_A.row(i), d_A.row(ofs)));
                    }
                    normSq(d_G(ofs, ofs), d_A.row(ofs));
                    ++d_G_computed;
                }
                else
                    d_G_computed = ofs;
            }
            if (add_row)
                do_add_row();
            ++d_vectors_used;
        }
    
        void insertVector(unsigned ofs, const linalg::math_rowvector<typename IntTypeContext::Integer> & result)
        // Inserts a new vector <result> at index ofs
        {
            // Here, result == d_A.row(ofs) since d_A was already updated
            bool add_row = d_vectors_used == d_G.rows();
            if (d_computed > ofs)
                d_computed = ofs;
            if (d_level > ofs)
                d_level = ofs;
            if (add_row)
                d_G.resize(d_G.rows() + 1, d_G.cols() + 1, linalg::Initialize(d_ic));
            if (d_G_computed > ofs)
            {
                if (d_G_computed >= ofs + result.size())
                {
                    // Insert row
                    for (unsigned i = d_G_computed; i > ofs; --i)
                        linalg::swap(d_G.row(i), d_G.row(i - 1));
                    // Insert column
                    for (unsigned i = d_G_computed; i > ofs; --i)
                        linalg::swap(d_G.col(i), d_G.col(i - 1));
                    // Add new data column
                    for (unsigned i = 0; i < ofs; ++i)
                        dot(d_G(ofs, i), d_A.row(i), d_A.row(ofs));
                    for (unsigned i = ofs + 1; i <= d_G_computed; ++i)
                        dot(d_G(i, ofs), d_A.row(i), d_A.row(ofs));
                    normSq(d_G(ofs, ofs), d_A.row(ofs));
                    ++d_G_computed;
                }
                else
                    d_G_computed = ofs;
            }
            if (add_row)
                do_add_row();
            ++d_vectors_used;
        }
    
        void removeZeroVector(unsigned ofs)
        // Removes the zero vector at position ofs
        {
            if (ofs < d_computed)
                assert(d_dependent[ofs]);
            if (ofs < d_G_computed)
            {
                // Remove row
                for (unsigned i = ofs + 1; i < d_G_computed; ++i)
                    linalg::swap(d_G.row(i), d_G.row(i - 1));
                // Remove column
                for (unsigned i = ofs + 1; i < d_G_computed; ++i)
                    linalg::swap(d_G.col(i), d_G.col(i - 1));
                --d_G_computed;
            }
            if (ofs < d_computed)
            {
                // Remove rows
                for (unsigned i = ofs + 1; i < d_computed; ++i)
                {
                    linalg::swap(d_r.row(i), d_r.row(i - 1));
                    linalg::swap(d_coeffs.row(i), d_coeffs.row(i - 1));
                    linalg::swap(d_s.row(i), d_s.row(i - 1));
                }
                // Remove columns
                for (unsigned i = ofs + 1; i < d_computed; ++i)
                {
                    linalg::swap(d_r.col(i), d_r.col(i - 1));
                    linalg::swap(d_coeffs.col(i), d_coeffs.col(i - 1));
                    linalg::swap(d_s.col(i), d_s.col(i - 1));
                }
                // Remove norms
                for (unsigned i = ofs + 1; i < d_computed; ++i)
                    arithmetic::swap(d_sqnorms[i], d_sqnorms[i - 1]);
                --d_computed;
            }
            if (ofs < d_level)
                --d_level;
            for (unsigned i = ofs + 1; i < d_vectors_used; ++i)
                d_dependent[i - 1] = d_dependent[i];
            --d_vectors_used;
        }
    
        static void compactify()
        {
            /*
              if (d_vectors_used < d_G.rows())
              {
              d_G.resize(d_vectors_used, d_vectors_used, linalg::Initialize(d_ic));
              d_r.resize(d_vectors_used, d_vectors_used, linalg::Initialize(d_rc));
              d_coeffs.resize(d_vectors_used, d_vectors_used, linalg::Initialize(d_rc));
              d_sqnorms.resize(d_vectors_used, linalg::Initialize(d_rc));
              d_s.resize(d_vectors_used, d_vectors_used, linalg::Initialize(d_rc));
              d_dependent.resize(d_vectors_used);
              }
            */
        }
    
        void changeOfPrecision()
        // Notify the Gram-Schmidt computation engine that the precision of d_rc changed
        {
            for (unsigned i = 0; i < d_coeffs.rows(); ++i)
            {
                for (unsigned j = 0; j < d_coeffs.rows(); ++j)
                {
                    d_r(i, j).setContext(d_rc);
                    d_s(i, j).setContext(d_rc);
                    d_coeffs(i, j).setContext(d_rc);
                }
                d_sqnorms[i].setContext(d_rc);
            }
            if (d_level > 0)
                // Redo GS
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

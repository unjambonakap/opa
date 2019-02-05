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

#ifndef PLLL_INCLUDE_GUARD__ENUMIMPL_SIEVE_CPP
#define PLLL_INCLUDE_GUARD__ENUMIMPL_SIEVE_CPP

// #define PLLL_INCLUDE_GUARD__ENUMIMPL_SIEVE_PROFILING

#include "enumimpl-sieve-impl.cpp"
#include <map>

#ifdef PLLL_INCLUDE_GUARD__ENUMIMPL_SIEVE_PROFILING
  #include "profiling.hpp"
#endif

namespace plll
{
    template<class SourceRealTypeContext,
             class RealTypeContext,
             class IntTypeContext>
    class ListSieveImpl
    /*
      Implements the ListSieve-Birthday algorithm described in X. Pujol, D. Stehle: "Solving the
      Shortest Lattice Vector Problem in Time 2^{2.465n}", and the ListSieve algorithm described in
      D. Micciancio, P. Voulgaris: "Faster exponential time algorithms for the shortest vector problem",
      2009.
    */
    {
    public:
        typedef boost::function<void(const linalg::math_matrix<typename IntTypeContext::Integer> & basis,
                                     int p,
                                     const linalg::math_rowvector<typename IntTypeContext::Integer> & vec)> CallbackFunction;
        
    private:
        RealTypeContext d_rc;
        IntTypeContext d_ic;
        CallbackFunction d_cf;
    
#ifdef PLLL_INCLUDE_GUARD__ENUMIMPL_SIEVE_PROFILING
        TimedDataCollector dc_gen, dc_sizered, dc_reduce;
#endif
    
        typedef std::map<linalg::math_rowvector<typename IntTypeContext::Integer>,
                         typename RealTypeContext::Real,
                         TypeVecPComparator<IntTypeContext> > ListContainer;
    
        linalg::math_matrix<typename RealTypeContext::Real> d_dotproducts;
    
        template<class Iterator>
        void findClosest(linalg::math_rowvector<typename IntTypeContext::Integer> & res, typename RealTypeContext::Real & currnorm,
                         const linalg::math_rowvector<typename IntTypeContext::Integer> & v, const typename RealTypeContext::Real & vsqnorm,
                         const Iterator & begin, const Iterator & end, RealTypeContext & rc, IntTypeContext & ic)
        /*
          Search the range [begin, end) for a w such that ||v - w|| is minimal. Stores the difference w
          - v in <res> and ||v - w||^2 in currnorm. In case res.size() != 0, only searches for w such
          that ||v - w||^2 < currnorm.
        */
        {
            DPLinForm<RealTypeContext, IntTypeContext> dp(rc);
            dp.reset_dotprods(v, d_dotproducts);
            typename RealTypeContext::Real tmp(rc), tmp2(rc);
            for (Iterator j = begin; j != end; ++j)
            {
                // Compute squared norm of difference
                dp(tmp, j->first);
                tmp <<= 1;
                tmp = j->second - tmp;
                tmp += vsqnorm;
                // If resulting vector is shorter, take it
                if ((res.size() == 0) || (tmp < currnorm))
                {
                    res.resize(v.size());
                    res = j->first - v;
                    currnorm = tmp;
                }
            }
        }
    
        void findClosestPair(linalg::math_rowvector<typename IntTypeContext::Integer> & res, const ListContainer & list, RealTypeContext & rc, IntTypeContext & ic)
        /*
          Searches for v, w appearing in <list> such that ||v - w|| > 0 is minimal, and returns v - w in
          <res>.
        */
        {
            res.resize(0);
            typename RealTypeContext::Real currnorm(rc);
            setZero(currnorm);
            for (typename ListContainer::const_iterator i = list.begin(); i != list.end(); ++i)
            {
                typename ListContainer::const_iterator j = i;
                ++j;
                findClosest(res, currnorm, i->first, i->second, j, list.end(), rc, ic);
            }
        }
    
        void listReduce(linalg::math_rowvector<typename IntTypeContext::Integer> & u,
                        linalg::math_rowvector<typename RealTypeContext::Real> & uprime, const ListContainer & reduction_list,
                        RealTypeContext & rc, IntTypeContext & ic, const typename RealTypeContext::Real & deltasq)
        /*
          Performs a list reduction of (u, u') using reduction_list. Here, u is a lattice vector, u+u' is
          constant during the execution of this function, and reduction_list is a list of (short) lattice
          vectors. We say that (u, u') are reduced if there exists no vector w in reduction_list such that
          ||u' - w|| < sqrt(deltasq) * ||u'||.
        */
        {
            if (reduction_list.empty())
                return;
            DPLinForm<RealTypeContext, IntTypeContext> dp(rc);
            dp.reset_dotprods(uprime, d_dotproducts);
            // While there exists a w in reduction_list with ||u' - w||^2 < deltasq * ||u'||^2, replace
            // (u, u') by (u + w, u' - w).
            bool change = true;
            typename RealTypeContext::Real n(rc), ns(rc), nn(rc), t(rc);
            // Compute norm of u'
            dp(n, uprime);
            // Multiply with delta^2
            ns = n * deltasq;
            while (change)
            {
                // Search for vectors w
                change = false;
                for (typename ListContainer::const_iterator w = reduction_list.begin(); w != reduction_list.end(); ++w)
                {
                    // Compute norm of u'-w
                    dp(nn, w->first);
                    nn <<= 1;
                    nn = w->second - nn;
                    nn += n;
                    // Is it shorter?
                    if (nn < ns)
                    {
                        // Multiply with delta^2
                        ns = nn * deltasq;
                        // Add w to u and subtract w from u'
                        u += w->first;
                        for (unsigned i = 0; i < uprime.size(); ++i)
                        {
                            arithmetic::convert(t, w->first[i], rc);
                            uprime[i] -= t;
                        }
                        // Restart
                        change = true;
                        break;
                    }
                }
            }
        }
    
        linalg::math_matrix<typename RealTypeContext::Real> d_basis;
    
        void GenReduce(linalg::math_rowvector<typename IntTypeContext::Integer> & u, unsigned dim,
                       GaussianGenerator<RealTypeContext> & grng, typename RealTypeContext::UniformRNG & rng,
                       const typename RealTypeContext::Real & xi, const typename RealTypeContext::Real & mu,
                       const linalg::math_matrix<typename RealTypeContext::Real> & inverse_basis, const ListContainer & reduction_list,
                       const typename RealTypeContext::Real & deltasq)
        {
            // Generate a random point u' in B_dim(0, xi*mu).
            linalg::math_rowvector<typename RealTypeContext::Real> uprime2;
            uprime2.resize(dim, linalg::Initialize(d_rc));
            linalg::math_rowvector<typename RealTypeContext::Real> uprime;
            uprime.resize(dim, linalg::Initialize(d_rc));
            {
#ifdef PLLL_INCLUDE_GUARD__ENUMIMPL_SIEVE_PROFILING
                Profiler<TimedDataCollector> p(dc_gen);
#endif
                // Generate a random point u' in B_dim(0, xi*mu).
                GenSphere(uprime2, grng, rng, d_rc, xi * mu);
                // Transform
                uprime = uprime2 * inverse_basis;
            }
        
            // Write u' as u+u', where u is a lattice vector and u' is "reduced" modulo the fundamental
            // parallelotope of the given basis of the lattice.
            u.resize(dim);
            {
#ifdef PLLL_INCLUDE_GUARD__ENUMIMPL_SIEVE_PROFILING
                Profiler<TimedDataCollector> p(dc_sizered);
#endif
                typename RealTypeContext::Real t(d_rc);
                for (unsigned i = 0; i < dim; ++i)
                {
                    arithmetic::convert_round(u[i], uprime[i], d_ic);
                    if (!isZero(u[i]))
                    {
                        arithmetic::convert(t, u[i], d_rc);
                        uprime[i] -= t;
                    }
                }
            }
        
            // While there exists a w in reduction_list with ||u' + w|| < (1 - 1/dim) * ||u'||, replace
            // (u, u') by (u - w, u' + w).
#ifdef PLLL_INCLUDE_GUARD__ENUMIMPL_SIEVE_PROFILING
            Profiler<TimedDataCollector> p(dc_reduce);
#endif
            listReduce(u, uprime, reduction_list, d_rc, d_ic, deltasq);
        
            // Return u
        }

        Verbose & d_verbose;
    
    public:
        ListSieveImpl(Verbose & v, GaussianFactorComputer &, unsigned enumdimension)
            : d_cf(NULL),
#ifdef PLLL_INCLUDE_GUARD__ENUMIMPL_SIEVE_PROFILING
              dc_gen("Generation"), dc_sizered("Size reduction"), dc_reduce("Reduction"),
#endif
              d_verbose(v)
        {
        }
    
        static void estimate_radius_lowerbound(typename SourceRealTypeContext::Real & lower_bound,
                                               const Lattice<SourceRealTypeContext, IntTypeContext> & lattice, unsigned begin, unsigned end,
                                               const typename SourceRealTypeContext::Real & bound)
        {
            typename SourceRealTypeContext::Real t(lattice.rc());
            lower_bound = lattice.getNormSq(begin);
            for (unsigned i = begin + 1; i <= end; ++i)
            {
                lattice.computeProjectionLengthBV(t, begin, i);
                if (lower_bound > t)
                    lower_bound = t;
            }
            if ((lower_bound > bound) && !isZero(bound))
                lower_bound = bound;
        }
    
        void enumerate_LS(const Lattice<SourceRealTypeContext, IntTypeContext> & lattice, unsigned begin, unsigned end,
                          linalg::math_rowvector<typename IntTypeContext::Integer> & result,
                          typename SourceRealTypeContext::Real & bound)
        /*
          Implements the ListSieve algorithm described in D. Micciancio, P. Voulgaris: "Faster
          exponential time algorithms for the shortest vector problem", 2009.
        */
        {
            // Find a lower bound on lambda_1(Lambda)
            typename SourceRealTypeContext::Real mu(lattice.rc());
            estimate_radius_lowerbound(mu, lattice, begin, end, bound);
        
            // Initialize random number generator
            arithmetic::RandomNumberGenerator rng_;
            rng_.randomizeSeed();
            typename RealTypeContext::UniformRNG rng(rng_);
        
            // Compute bounds
            unsigned long minprec = 0;
            if (RealTypeContext::is_variable_precision)
            {
                // Determine a precision ("just guessin'")
                minprec = 2 * arithmetic::convert_ceil<long>(log2(arithmetic::convert(mu, d_rc))) + 20;
                minprec = (minprec + 63) & ~63; // round to multiple of 64
                if (minprec == 0)
                    minprec = 64;
                if (minprec > d_rc.getRealPrecision())
                    d_rc.setRealPrecision(minprec);
            }
        
            const unsigned dim = end - begin + 1;
        
            GaussianGenerator<RealTypeContext> grng(d_rc, rng);
            typename RealTypeContext::Real deltasq(d_rc);
            arithmetic::convert(deltasq, dim - 1, d_rc);
            deltasq /= arithmetic::convert(dim, d_rc);
            square(deltasq, deltasq);
        
            typename RealTypeContext::Real xi(d_rc), c(d_rc), rmu(d_rc);
            arithmetic::convert(rmu, mu, d_rc);
            sqrt(rmu, rmu); // take square root!
            arithmetic::convert(xi, 0.685, d_rc);
            c = log(xi + sqrt(square(xi) + arithmetic::convert(1, d_rc))) + arithmetic::convert(0.401, d_rc) + arithmetic::convert(0.5, d_rc) * log(square(xi) / (square(xi) - arithmetic::convert(0.25, d_rc)));
            typename RealTypeContext::Real KK = power(arithmetic::convert(2, d_rc), c * arithmetic::convert(dim, d_rc));
            arithmetic::Integer K = arithmetic::convert_ceil<arithmetic::Integer>(KK);
            unsigned long k = arithmetic::convert<long>(K);
            {
                Verbose::VerboseStream vv = d_verbose(LatticeReduction::VL_Information);
                vv << "ListSieve [" << begin << "," << end << "]: ";
                if (RealTypeContext::is_variable_precision)
                    vv << "precision = " << minprec << ", ";
                vv << "k = " << k << " (K = " << K << ", c = " << c << ")";
            }
        
            // Compute approximation of basis and its inverse
            linalg::math_matrix<typename RealTypeContext::Real> basis, inverse_basis;
            computeBasisApproximation<RealTypeContext>(basis, inverse_basis, d_rc, begin, end, lattice);
            d_basis.resize(0, 0);
            d_basis.resize(basis.rows(), basis.cols(), linalg::Initialize(d_rc));
            d_basis = basis;
            typename SourceRealTypeContext::Real st(lattice.rc());
            d_dotproducts.resize(0, 0);
            d_dotproducts.resize(dim, dim, linalg::Initialize(d_rc));
            for (unsigned i = 0; i < dim; ++i)
                for (unsigned j = 0; j <= i; ++j)
                {
                    lattice.computeDotProductProjected(st, begin, begin + i, begin + j);
                    arithmetic::convert(d_dotproducts(i, j), st, d_rc);
                    if (j < i)
                        d_dotproducts(j, i) = d_dotproducts(i, j);
                }
        
            while (true)
            {
                d_verbose(LatticeReduction::VL_Information) << "Radius: " << rmu;
                d_verbose(LatticeReduction::VL_Information) << "Starting sieving...";
            
                ListContainer L;
                linalg::math_rowvector<typename IntTypeContext::Integer> vec;
                typename RealTypeContext::Real bound(d_rc), norm(d_rc), n(d_rc);
                DPLinForm<RealTypeContext, IntTypeContext> dp(d_rc);
                // Sieve vectors
                vec.resize(dim);
                setZero(norm);
                L.insert(std::make_pair(vec, norm)); // insert 0 vector
                bound = square(rmu);
                Percentage pp1(d_verbose, k);
                for (unsigned long i = 0; i < k; ++i)
                {
                    pp1.update(i);
                    GenReduce(vec, dim, grng, rng, xi, rmu, inverse_basis, L, deltasq);
                    if (L.find(vec) == L.end())
                    {
                        // Compute norm
                        dp.reset_dotprods(vec, d_dotproducts);
                        dp(norm, vec);
                        // Look if we find difference which is ultra-short
                        for (typename ListContainer::iterator j = L.begin(); j != L.end(); ++j)
                        {
                            // Compute squared norm of vec - j->first
                            dp(n, j->first);
                            n <<= 1;
                            n = j->second - n;
                            n += norm;
                            // Shorter?
                            if (n < bound)
                            {
                                result.resize(vec.size());
                                result = vec - j->first;
                                break;
                            }
                        }
                        if (result.size())
                            break;
                        L.insert(std::make_pair(vec, norm));
                    }
                }
                d_verbose(LatticeReduction::VL_Information) << "Size of L is " << L.size();
            
                if (result.size() > 0)
                {
                    // Found something!
                    Verbose::VerboseStream vv = d_verbose(LatticeReduction::VL_Information);
                    vv << "Found vector";
                    if (lattice.canGetIntegerVector())
                    {
                        linalg::math_rowvector<typename IntTypeContext::Integer> res;
                        lattice.getLinearCombination(res, begin, result);
                        vv << " " << res << " of squared norm " << normSq(res) << " and";
                    }
                    dp.reset_dotprods(result, d_dotproducts);
                    dp(norm, result);
                    vv << " of projected norm " << sqrt(norm) << "!";
                    break;
                }
            
                // Increase lower bound on lambda_1
                rmu *= arithmetic::convert(dim + 1, d_rc) / arithmetic::convert(dim, d_rc);
            }
        }
    
        void enumerate_LSB(const Lattice<SourceRealTypeContext, IntTypeContext> & lattice, unsigned begin, unsigned end,
                           linalg::math_rowvector<typename IntTypeContext::Integer> & result,
                           typename SourceRealTypeContext::Real & bound)
        /*
          Implements the ListSieve-Birthday algorithm described in X. Pujol, D. Stehle: "Solving the
          Shortest Lattice Vector Problem in Time 2^{2.465n}".
        */
        {
            // Find a lower bound on lambda_1(Lambda)
            typename SourceRealTypeContext::Real mu(lattice.rc());
            estimate_radius_lowerbound(mu, lattice, begin, end, bound);
        
            // Initialize random number generator
            arithmetic::RandomNumberGenerator rng_;
            rng_.randomizeSeed();
            typename RealTypeContext::UniformRNG rng(rng_);
        
            // Compute bounds
            unsigned long minprec = 0;
            if (RealTypeContext::is_variable_precision)
            {
                // Determine a precision ("just guessin'")
                minprec = 2 * arithmetic::convert_ceil<long>(log2(arithmetic::convert(mu, d_rc))) + 20;;
                minprec = (minprec + 63) & ~63; // round to multiple of 64
                if (minprec == 0)
                    minprec = 64;
                if (minprec > d_rc.getRealPrecision())
                    d_rc.setRealPrecision(minprec);
            }
        
            const unsigned dim = end - begin + 1;
        
            GaussianGenerator<RealTypeContext> grng(d_rc, rng);
            typename RealTypeContext::Real deltasq(d_rc);
            arithmetic::convert(deltasq, dim - 1, d_rc);
            deltasq /= arithmetic::convert(dim, d_rc);
            square(deltasq, deltasq);
        
            typename RealTypeContext::Real xi(d_rc), r0(d_rc), rmu(d_rc);
            arithmetic::convert(xi, 0.9476, d_rc); // value from Proof of Theorem 1
            arithmetic::convert(r0, 3.0169, d_rc); // value from Proof of Theorem 1
            arithmetic::convert(rmu, mu, d_rc);
            sqrt(rmu, rmu); // take square root!
        
            // EXPLICIT VALUES FOR o(dim) ARE MISSING!!! ??? ...
            typename RealTypeContext::Real NB = power(arithmetic::convert(2, d_rc),
                                                      (log2(arithmetic::convert(r0, d_rc)) + arithmetic::convert(0.401, d_rc))
                                                      * arithmetic::convert(dim, d_rc) + arithmetic::convert(0/*o(dim)*/, d_rc)); // formula from Lemma 3
            typename RealTypeContext::Real NT = power(arithmetic::convert(2, d_rc),
                                                      (arithmetic::convert(-0.5, d_rc) * log2(arithmetic::convert(1, d_rc) - arithmetic::convert(2, d_rc) * xi / arithmetic::convert(r0, d_rc)) + arithmetic::convert(0.401, d_rc))
                                                      * arithmetic::convert(dim, d_rc) + arithmetic::convert(0/*o(dim)*/, d_rc)); // formula from Lemma 4
            typename RealTypeContext::Real NG = power(arithmetic::convert(2, d_rc),
                                                      arithmetic::convert(-0.5, d_rc) * log2(arithmetic::convert(1, d_rc) - arithmetic::convert(1, d_rc) / (arithmetic::convert(4, d_rc) * square(xi)))
                                                      * arithmetic::convert(dim, d_rc) + arithmetic::convert(0/*o(dim)*/, d_rc)); // formula from Lemma 5
            arithmetic::Integer N1max = arithmetic::convert_ceil<arithmetic::Integer>(arithmetic::convert(4, d_rc) * NG * NT);
            arithmetic::Integer N2 = arithmetic::convert_ceil<arithmetic::Integer>(arithmetic::convert(8, d_rc) * NG)
                * arithmetic::convert_ceil<arithmetic::Integer>(sqrt(NB));
        
            arithmetic::Integer N1 = rng_.random(N1max);
            unsigned long n1 = arithmetic::convert<long>(N1), n2 = arithmetic::convert<long>(N2);
            {
                Verbose::VerboseStream vv = d_verbose(LatticeReduction::VL_Information);
                vv << "ListSieve-Birthday [" << begin << "," << end << "]: ";
                if (RealTypeContext::is_variable_precision)
                    vv << "precision = " << minprec << ", ";
                vv << "n1 = " << n1 << ", n2 = " << n2 << " (NB = " << NB << ", NT = " << NT << ", NG = " << NG << ", n1_max = " << N1max << ", N2 = " << N2 << ")";
            }
        
            // Compute approximation of basis
            linalg::math_matrix<typename RealTypeContext::Real> basis, inverse_basis;
            computeBasisApproximation<RealTypeContext>(basis, inverse_basis, d_rc, begin, end, lattice);
            typename SourceRealTypeContext::Real st(lattice.rc());
            d_dotproducts.resize(0, 0);
            d_dotproducts.resize(dim, dim, linalg::Initialize(d_rc));
            for (unsigned i = 0; i < dim; ++i)
                for (unsigned j = 0; j <= i; ++j)
                {
                    lattice.computeDotProductProjected(st, begin, begin + i, begin + j);
                    arithmetic::convert(d_dotproducts(i, j), st, d_rc);
                    if (j < i)
                        d_dotproducts(j, i) = d_dotproducts(i, j);
                }
        
            while (true)
            {
                d_verbose(LatticeReduction::VL_Information) << "Radius: " << rmu;
                d_verbose(LatticeReduction::VL_Information) << "Starting sieving...";
                // Sieve vectors
                ListContainer T, U;
                linalg::math_rowvector<typename IntTypeContext::Integer> vec;
                vec.resize(dim);
                typename RealTypeContext::Real bound(d_rc), norm(d_rc);
                DPLinForm<RealTypeContext, IntTypeContext> dp(d_rc);
                bound = square(r0 * rmu);
                Percentage pp1(d_verbose, n1);
                for (unsigned long i = 0; i < n1; ++i)
                {
                    pp1.update(i);
                    GenReduce(vec, dim, grng, rng, xi, rmu, inverse_basis, T, deltasq);
                    dp.reset_dotprods(vec, d_dotproducts);
                    dp(norm, vec);
                    if (norm >= bound)
                        T.insert(std::make_pair(vec, norm));
                }
                d_verbose(LatticeReduction::VL_Information) << "Size of T is " << T.size();
                d_verbose(LatticeReduction::VL_Information) << "Starting second sieving round...";
                Percentage pp2(d_verbose, n2);
                for (unsigned long i = 0; i < n2; ++i)
                {
                    pp2.update(i);
                    GenReduce(vec, dim, grng, rng, xi, rmu, inverse_basis, T, deltasq);
                    dp.reset_dotprods(vec, d_dotproducts);
                    dp(norm, vec);
                    U.insert(std::make_pair(vec, norm));
                }
                d_verbose(LatticeReduction::VL_Information) << "Size of U is " << U.size();
            
                // Find closest pair
                d_verbose(LatticeReduction::VL_Information) << "Finding closest pair...";
                findClosestPair(result, U, d_rc, d_ic);
            
                if (result.size() > 0)
                {
                    // Found something!
                    Verbose::VerboseStream vv = d_verbose(LatticeReduction::VL_Information);
                    vv << "Found vector";
                    if (lattice.canGetIntegerVector())
                    {
                        linalg::math_rowvector<typename IntTypeContext::Integer> res;
                        lattice.getLinearCombination(res, begin, result);
                        vv << " " << res << " of squared norm " << normSq(res) << " and";
                    }
                    dp.reset_dotprods(result, d_dotproducts);
                    dp(norm, result);
                    vv << " of projected norm " << sqrt(norm) << "!";
                    break;
                }
            
                // Increase lower bound on lambda_1
                rmu *= arithmetic::convert(dim + 1, d_rc) / arithmetic::convert(dim, d_rc);
            }
        }
    
        void setCallback(CallbackFunction cf)
        {
            d_cf = cf;
        }
    };
    
    template<class SourceRealTypeContext,
             class RealTypeContext,
             class IntTypeContext,
             class RealIntTypeContext>
    class GaussSieveImpl2
    /*
      Implements the GaussSieve algorithm described in D. Micciancio, P. Voulgaris: "Faster exponential
      time algorithms for the shortest vector problem", 2009.
    */
    {
    public:
        typedef boost::function<void(const linalg::math_matrix<typename IntTypeContext::Integer> & basis,
                                     int p,
                                     const linalg::math_rowvector<typename IntTypeContext::Integer> & vec)> CallbackFunction;
        
    private:
        RealTypeContext d_rc;
        IntTypeContext d_ic;
        RealIntTypeContext & d_ric;
        CallbackFunction d_cf;
        
#ifdef PLLL_INCLUDE_GUARD__ENUMIMPL_SIEVE_PROFILING
        TimedDataCollector dc_gen, dc_reduce;
#endif
        
        linalg::math_rowvector<typename RealTypeContext::Real> d_gaussian_sprime, d_gaussian_sprime_squared, d_gaussian_coeffs;
        linalg::math_matrix<typename RealTypeContext::Real> d_gaussian_gs;
        typename RealTypeContext::Real d_gaussian_logdim, d_minus_pi;
        
        void initGaussian(const Lattice<SourceRealTypeContext, IntTypeContext> & lattice, unsigned begin, unsigned end)
        {
            unsigned dim = end - begin + 1;
            d_gaussian_sprime.resize(dim, linalg::Initialize(d_rc));
            d_gaussian_sprime_squared.resize(dim, linalg::Initialize(d_rc));
            d_gaussian_coeffs.resize(dim, linalg::Initialize(d_rc));
            d_gaussian_gs.resize(dim, dim, linalg::Initialize(d_rc));
            for (unsigned i = 0; i < d_gaussian_sprime.size(); ++i)
            {
                d_gaussian_sprime[i].setContext(d_rc);
                d_gaussian_sprime_squared[i].setContext(d_rc);
                d_gaussian_coeffs[i].setContext(d_rc);
                for (unsigned j = 0; j < d_gaussian_sprime.size(); ++j)
                {
                    d_gaussian_gs(i, j).setContext(d_rc);
                    if (j < i)
                        arithmetic::convert(d_gaussian_gs(i, j), lattice.getCoeff(begin + i, begin + j), d_rc);
                }
            }
        
            // Compute maximum
            typename RealTypeContext::Real m(d_rc), mm(d_rc);
            typename SourceRealTypeContext::Real min(lattice.rc());
            min = lattice.getNormSq(begin);
            for (unsigned i = 1; i < dim; ++i)
                if (min < lattice.getNormSq(begin + i))
                    min = lattice.getNormSq(begin + i);
            arithmetic::convert(m, min, d_rc);
            // Multiply with sqrt(log(dim))
            d_gaussian_logdim.setContext(d_rc);
            arithmetic::convert(d_gaussian_logdim, dim, d_rc);
            d_gaussian_logdim = log(d_gaussian_logdim);
            m *= d_gaussian_logdim;
            // Fill in s' and s'^2
            for (unsigned i = 0; i < d_gaussian_sprime.size(); ++i)
            {
                arithmetic::convert(mm, lattice.getNormSq(begin + i), d_rc);
                d_gaussian_sprime_squared[i] = m / mm;
                d_gaussian_sprime[i] = sqrt(d_gaussian_sprime_squared[i]);
            }
        
            // Compute -pi
            d_minus_pi.setContext(d_rc);
            d_rc.getPi(d_minus_pi);
            neg(d_minus_pi, d_minus_pi);
        }
        
        void sampleGaussianZ(typename IntTypeContext::Integer & v, const typename RealTypeContext::Real & c,
                             const typename RealTypeContext::Real & s, const typename RealTypeContext::Real & ssq,
                             typename RealTypeContext::UniformRNG & rngu, arithmetic::RandomNumberGenerator & rng)
        {
            typename IntTypeContext::Integer min_c, max_c;
            arithmetic::convert_floor(min_c, c - s * d_gaussian_logdim, d_ic);
            arithmetic::convert_ceil(max_c, c + s * d_gaussian_logdim, d_ic);
            typename IntTypeContext::Integer int_len = max_c - min_c;
            ++int_len;
            typename RealTypeContext::Real t(d_rc);
            typename IntTypeContext::UniformRNG urng(rng);
            while (true)
            {
                // Compute random integer in interval
                urng.random(v, int_len);
                v += min_c;
                // Check if it satisfies the desired distribution
                arithmetic::convert(t, v, d_rc);
                t -= c;
                square(t, t);
                t /= ssq;
                t *= d_minus_pi;
                t = exp(t);
                if (rngu.randomUniform(d_rc) <= t)
                    break;
            }
        }
        
        void sampleGaussian(linalg::math_rowvector<typename IntTypeContext::Integer> & v,
                            typename RealTypeContext::UniformRNG & rngu, arithmetic::RandomNumberGenerator & rng)
        {
            for (unsigned i = 0; i < d_gaussian_coeffs.size(); ++i)
                setZero(d_gaussian_coeffs[i]);
            for (unsigned i = d_gaussian_coeffs.size(); i > 0; )
            {
                --i;
                // Sample
                sampleGaussianZ(v[0], d_gaussian_coeffs[i], d_gaussian_sprime[i], d_gaussian_sprime_squared[i], rngu, rng);
                arithmetic::convert(d_gaussian_coeffs[i], v[0], d_rc);
                // Add vector
                for (unsigned j = 0; j < i; ++j)
                    d_gaussian_coeffs[j] -= d_gaussian_coeffs[i] * d_gaussian_gs(i, j);
            }
            // Round
            for (unsigned i = 0; i < v.size(); ++i)
                arithmetic::convert_round(v[i], d_gaussian_coeffs[i], d_ic);
        }
        
        Verbose & d_verbose;
        
        template<class IntTypeContext_, class RealIntTypeContext_, bool bothInteger>
        class ComputeRoundedQuotient;
        
        template<class IntTypeContext_, class RealIntTypeContext_>
        class ComputeRoundedQuotient<IntTypeContext_, RealIntTypeContext_, true>
        {
        public:
            static void computeRoundedQuotient(typename IntTypeContext_::Integer & result,
                                               const typename RealIntTypeContext_::Type & a, const typename RealIntTypeContext_::Type & b,
                                               IntTypeContext_ & ic, RealIntTypeContext_ & rc)
            // compute round(ii, a / b)
            {
                arithmetic::convert_round(result, a / b, ic);
            }
        };
        
        template<class IntTypeContext_, class RealIntTypeContext_>
        class ComputeRoundedQuotient<IntTypeContext_, RealIntTypeContext_, false>
        {
        public:
            static void computeRoundedQuotient(typename IntTypeContext_::Integer & result,
                                               const typename RealIntTypeContext_::Type & a, const typename RealIntTypeContext_::Type & b,
                                               IntTypeContext_ & ic, RealIntTypeContext_ & rc)
            // compute round(ii, a / b)
            {
                roundDiv(result, a, b);
            }
        };
        
        static void computeRoundedQuotient(typename IntTypeContext::Integer & result,
                                           const typename RealIntTypeContext::Type & a, const typename RealIntTypeContext::Type & b,
                                           IntTypeContext & ic, RealIntTypeContext & rc)
        // compute round(ii, a / b)
        {
            ComputeRoundedQuotient<IntTypeContext, RealIntTypeContext, RealIntTypeContext::is_realtype>::
                computeRoundedQuotient(result, a, b, ic, rc);
        }
    
    public:
        linalg::math_matrix<typename RealIntTypeContext::Type> dot_products;
    
        GaussSieveImpl2(IntTypeContext & ic, RealIntTypeContext & c, Verbose & v, unsigned dimension, CallbackFunction cf)
            : d_ic(ic), d_ric(c), d_cf(cf),
#ifdef PLLL_INCLUDE_GUARD__ENUMIMPL_SIEVE_PROFILING
              dc_gen("Generation"), dc_reduce("Reduction"),
#endif
              d_verbose(v)
        {
            dot_products.resize(dimension, dimension, linalg::Initialize(d_ric));
        }
    
        GaussSieveImpl2(IntTypeContext & ic, Verbose & v, unsigned dimension, CallbackFunction cf)
            : d_ic(ic), d_ric(d_rc), d_cf(cf),
#ifdef PLLL_INCLUDE_GUARD__ENUMIMPL_SIEVE_PROFILING
              dc_gen("Generation"), dc_reduce("Reduction"),
#endif
              d_verbose(v)
        {
            dot_products.resize(dimension, dimension, linalg::Initialize(d_ric));
        }
    
        void enumerate_GS(Lattice<SourceRealTypeContext, IntTypeContext> & lattice, unsigned begin, unsigned end,
                          linalg::math_rowvector<typename IntTypeContext::Integer> & result,
                          typename SourceRealTypeContext::Real & bound)
        /*
          Implements the GaussSieve algorithm described in D. Micciancio, P. Voulgaris: "Faster
          exponential time algorithms for the shortest vector problem", 2009.
        */
        {
            // Initialize random number generator
            arithmetic::RandomNumberGenerator rng_;
            rng_.randomizeSeed();
            typename RealTypeContext::UniformRNG rng(rng_);
        
        
            // Compute bounds
            unsigned long minprec = 0;
            if (RealTypeContext::is_variable_precision)
            {
                // First determine maximal length of basis vectors
                typename SourceRealTypeContext::Real mu(lattice.rc());
                mu = lattice.getNormSq(begin);
                for (unsigned i = begin + 1; i <= end; ++i)
                    if (mu < lattice.getNormSq(i))
                        mu = lattice.getNormSq(i);
                // Determine a precision ("just guessin'")
                minprec = 2 * arithmetic::convert_ceil<long>(log2(arithmetic::convert(mu, d_rc))) + 20;
                minprec = (minprec + 63) & ~63; // round to multiple of 64
                if (minprec == 0)
                    minprec = 64;
                if (minprec > d_rc.getRealPrecision())
                    d_rc.setRealPrecision(minprec);
            }
        
            const unsigned dim = end - begin + 1;
        
            {
                Verbose::VerboseStream vv = d_verbose(LatticeReduction::VL_Information);
                vv << "GaussSieve [" << begin << "," << end << "]";
                if (RealTypeContext::is_variable_precision)
                    vv << ": precision = " << minprec;
            }
        
            unsigned long iterations = 0, collisions = 0;
            std::list<std::pair<linalg::math_rowvector<typename IntTypeContext::Integer>, typename RealIntTypeContext::Type> > L;
            linalg::math_rowvector<typename IntTypeContext::Integer> vec;
            vec.resize(dim);
            std::list<linalg::math_rowvector<typename IntTypeContext::Integer> > S;
            for (unsigned i = 0; i < dim; ++i)
            {
                ++vec[i];
                S.push_front(vec);
                --vec[i];
            }
            std::size_t max_size = L.size();

            typename RealIntTypeContext::Type n(d_ric), nn(d_ric), shortest_n(d_ric);
            setZero(shortest_n);
            initGaussian(lattice, begin, end);
            try
            {
                do
                {
                    ++iterations;
                    if (((iterations % 1000) == 0))
                    {
                        Verbose::VerboseStream vv = d_verbose(LatticeReduction::VL_Information);
                        vv << "Iteration " << iterations << ": list has " << L.size() << " elements, with a maximum of " << max_size << " ever attained; stack has " << S.size() << " elements, and " << collisions << " collisions occured";
                        if (!L.empty()) vv << ", minimal norm^2 so far is " << L.front().second;
//                dc_gen.output();
//                dc_reduce.output();
                    }
                    // Get next vector
                    if (!S.empty())
                    {
                        vec = S.back();
                        S.pop_back();
                    }
                    else
                    {
#ifdef PLLL_INCLUDE_GUARD__ENUMIMPL_SIEVE_PROFILING
                        Profiler<TimedDataCollector> p(dc_gen);
#endif
                        sampleGaussian(vec, rng, rng_);
                    }
                    
                    // Do Gauss reduction
                    {
                        // Set up linear form and compute norm
                        DPLinForm<RealIntTypeContext, IntTypeContext> dp(d_ric);
                        dp.reset_dotprods(vec, dot_products);
                        dp(n, vec);
                        // Reduce
#ifdef PLLL_INCLUDE_GUARD__ENUMIMPL_SIEVE_PROFILING
                        Profiler<TimedDataCollector> p(dc_reduce);
#endif
                        typename std::list<std::pair<linalg::math_rowvector<typename IntTypeContext::Integer>, typename RealIntTypeContext::Type> >::iterator w;
                        bool found;
                        typename IntTypeContext::Integer ii;
                        do
                        {
                            // Look for vectors in L which can reduce v
                            found = false;
                            for (w = L.begin(); w != L.end(); ++w)
                            {
                                if (w->second > n)
                                    // We won't find anything shorter, so stop
                                    break;
                                dp(nn, w->first); // nn = dot(v, w->first)
                                computeRoundedQuotient(ii, nn, w->second, d_ic, d_ric); // round(ii, nn / w->second);
                                if (!isZero(ii))
                                {
                                    found = true;
                                    // Reduce v by w->first
                                    vec -= ii * w->first;
                                    // Update dot product linear form
                                    dp.reset_dotprods(vec, dot_products);
                                    dp(n, vec);
                                }
                            }
                        } while (found);
                        if (!isZero(vec))
                        {
                            // Insert element at the right place
                            L.insert(w, std::make_pair(vec, n));
                            // Look for vectors in L which can be reduced by v
                            while (w != L.end())
                            {
                                dp(nn, w->first); // nn = dot(v, w->first)
                                computeRoundedQuotient(ii, nn, n, d_ic, d_ric); // round(ii, nn / n);
                                if (!isZero(ii))
                                {
                                    // Reduce w->first by v
                                    w->first -= ii * vec;
                                    // Add to stack
                                    S.push_back(w->first);
                                    w = L.erase(w);
                                }
                                else
                                    ++w;
                            }
                        }
                    }
                    
                    // Collision?
                    if (isZero(vec))
                        ++collisions; // yes
                    
                    // Keep track of maximal list size
                    max_size = std::max(max_size, L.size());
                    
                    // Found new shortest?
                    if (!L.empty())
                        if (isZero(shortest_n) || (L.front().second < shortest_n))
                        {
                            shortest_n = L.front().second;
                            // Callback available?
                            if (!d_cf.empty() && lattice.canGetIntegerVector())
                            {
                                std::pair<linalg::math_matrix<typename IntTypeContext::Integer> *, unsigned> ll = lattice.getMatrixIndex(begin);
                                if (ll.first)
                                    d_cf(*ll.first, ll.second, L.front().first);
                            }
                        }
                }
                while (10 * collisions < max_size + 2000); // Heuristic: if there are too many collisions, terminate
            }
            catch (LatticeReduction::stop_enumeration &)
            {
                d_verbose(LatticeReduction::VL_Chatter) << "Stopping enumeration";
            }
            d_verbose(LatticeReduction::VL_Information) << "Total of " << iterations << " iterations and " << collisions << " collisions, with a maximal list size of " << max_size << ".";
        
            // Find shortest list entry
            if (!L.empty())
            {
                result.resize(dim);
                result = L.front().first;
                // Found something!
                Verbose::VerboseStream vv = d_verbose(LatticeReduction::VL_Information);
                vv << "Found vector";
                if (lattice.canGetIntegerVector())
                {
                    linalg::math_rowvector<typename IntTypeContext::Integer> res;
                    lattice.getLinearCombination(res, begin, result);
                    vv << " " << res << " of squared norm " << normSq(res) << " and";
                }
                vv << " of projected squared norm " << L.front().second << "!";
            }
        }
        
        RealTypeContext & rc()
        {
            return d_rc;
        }
    };
    
    template<class SourceRealTypeContext,
             class RealTypeContext,
             class IntTypeContext>
    class GaussSieveImpl
    /*
      Implements the GaussSieve algorithm described in D. Micciancio, P. Voulgaris: "Faster exponential
      time algorithms for the shortest vector problem", 2009.
    */
    {
    public:
        typedef boost::function<void(const linalg::math_matrix<typename IntTypeContext::Integer> & basis,
                                     int p,
                                     const linalg::math_rowvector<typename IntTypeContext::Integer> & vec)> CallbackFunction;
        
    private:
        CallbackFunction d_cf;
        Verbose & d_verbose;
        
        template<class SRTC, class RTC, class ITC> // No need for the general case.
        class InitFPDPs;
        
        template<class RTC, class ITC>
        class InitFPDPs<arithmetic::RationalContext, RTC, ITC>
        {
        public:
            void operator() (Lattice<arithmetic::RationalContext, ITC> & lattice,
                             unsigned begin, unsigned end, linalg::math_matrix<typename RTC::Real> & dps, RTC & rc) const
            {
                typename SourceRealTypeContext::Real t(lattice.rc());
                for (unsigned i = 0; i < dps.rows(); ++i)
                    for (unsigned j = 0; j <= i; ++j)
                    {
                        lattice.computeDotProductProjected(t, begin, begin + i, begin + j);
                        arithmetic::convert(dps(i, j), t, rc);
                        if (j < i)
                            dps(j, i) = dps(i, j);
                    }
            }
        };
        
        template<class RTC, class ITC>
        class InitFPDPs<RTC, RTC, ITC>
        {
        public:
            void operator() (Lattice<RTC, ITC> & lattice,
                             unsigned begin, unsigned end, linalg::math_matrix<typename RTC::Real> & dps, RTC & rc) const
            {
                for (unsigned i = 0; i < dps.rows(); ++i)
                    for (unsigned j = 0; j <= i; ++j)
                    {
                        lattice.computeDotProductProjected(dps(i, j), begin, begin + i, begin + j);
                        if (j < i)
                            dps(j, i) = dps(i, j);
                    }
            }
        };
        
    public:
        GaussSieveImpl(Verbose & v, GaussianFactorComputer &, unsigned)
            : d_cf(NULL), d_verbose(v)
        {
        }
        
        void enumerate_GS(Lattice<SourceRealTypeContext, IntTypeContext> & lattice, unsigned begin, unsigned end,
                          linalg::math_rowvector<typename IntTypeContext::Integer> & result,
                          typename SourceRealTypeContext::Real & bound)
        /*
          Implements the GaussSieve algorithm described in D. Micciancio, P. Voulgaris: "Faster
          exponential time algorithms for the shortest vector problem", 2009.
        */
        {
            bool do_int = lattice.canGetIntegerVector() && (begin == 0);
            if (do_int)
            {
                GaussSieveImpl2<SourceRealTypeContext, RealTypeContext, IntTypeContext, IntTypeContext>
                    gs(lattice.ic(), lattice.ic(), d_verbose, end - begin + 1, d_cf);
                for (unsigned i = 0; i < gs.dot_products.rows(); ++i)
                    for (unsigned j = 0; j <= i; ++j)
                    {
                        lattice.getDotProduct(gs.dot_products(i, j), begin + i, begin + j);
                        if (j < i)
                            gs.dot_products(j, i) = gs.dot_products(i, j);
                    }
                gs.enumerate_GS(lattice, begin, end, result, bound);
            }
            else
            {
                GaussSieveImpl2<SourceRealTypeContext, RealTypeContext, IntTypeContext, RealTypeContext>
                    gs(lattice.ic(), d_verbose, end - begin + 1, d_cf);
                InitFPDPs<SourceRealTypeContext, RealTypeContext, IntTypeContext>()(lattice, begin, end, gs.dot_products, gs.rc());
                gs.enumerate_GS(lattice, begin, end, result, bound);
            }
        }
        
        void setCallback(CallbackFunction cf)
        {
            d_cf = cf;
        }
    };
    
    template<class RealTypeContext, class IntTypeContext>
    class Sieve
    {
    public:
        enum Algorithm { A_ListSieve, A_ListSieveBirthday, A_GaussSieve };
        
    private:
        typedef typename helper::SelectFirstType<RealTypeContext::has_uniform_rng && RealTypeContext::has_squareroot,
                                                 ListSieveImpl<RealTypeContext, RealTypeContext, IntTypeContext>,
                                                 ListSieveImpl<RealTypeContext, arithmetic::RealContext, IntTypeContext> >::result LSImplementation;
        
        typedef typename helper::SelectFirstType<RealTypeContext::has_uniform_rng && RealTypeContext::has_squareroot,
                                                 GaussSieveImpl<RealTypeContext, RealTypeContext, IntTypeContext>,
                                                 GaussSieveImpl<RealTypeContext, arithmetic::RealContext, IntTypeContext> >::result GSImplementation;
        
        typedef boost::function<void(const linalg::math_matrix<typename IntTypeContext::Integer> & basis,
                                     int p,
                                     const linalg::math_rowvector<typename IntTypeContext::Integer> & vec)> CallbackFunction;
        
        LSImplementation d_ls_impl;
        GSImplementation d_gs_impl;
        
    public:
        Sieve(Verbose & v, GaussianFactorComputer & gf, unsigned enumdimension)
            : d_ls_impl(v, gf, enumdimension), d_gs_impl(v, gf, enumdimension)
        {
        }
        
        bool enumerate(Algorithm alg, Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned begin, unsigned end,
                       linalg::math_rowvector<typename IntTypeContext::Integer> & result, typename RealTypeContext::Real & bound)
        // Finds a shortest vector in the lattice generated by the orthogonal projections of the vectors
        // A.row(begin) to A.row(end) into the orthogonal complement of the vectors A.row(0) to
        // A.row(begin-1).
        {
            result.resize(0);
            switch (alg)
            {
            case A_ListSieve:
                d_ls_impl.enumerate_LS(lattice, begin, end, result, bound);
                break;
            case A_ListSieveBirthday:
                d_ls_impl.enumerate_LSB(lattice, begin, end, result, bound);
                break;
            default:
            case A_GaussSieve:
                d_gs_impl.enumerate_GS(lattice, begin, end, result, bound);
                break;
            }
            return result.size() != 0;
        }
        
        void setCallback(CallbackFunction cf)
        {
            d_ls_impl.setCallback(cf);
            d_gs_impl.setCallback(cf);
        }
    };
}

#endif

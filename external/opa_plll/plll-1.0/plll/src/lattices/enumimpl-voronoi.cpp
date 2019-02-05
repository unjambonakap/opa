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

#ifndef PLLL_INCLUDE_GUARD__ENUMIMPL_VORONOI_CPP
#define PLLL_INCLUDE_GUARD__ENUMIMPL_VORONOI_CPP

#include "enumimpl-sieve-impl.cpp"
#include <map>
#include <set>

namespace plll
{
    template<class RealTypeContext, class IntTypeContext, class RealIntTypeContext>
    class VoronoiCell
    {
    public:
        typedef boost::function<void(const linalg::math_matrix<typename IntTypeContext::Integer> & basis,
                                     int p,
                                     const linalg::math_rowvector<typename IntTypeContext::Integer> & vec)> CallbackFunction;
        
    private:
        class DPLinForm // Dot product linear form with basis
        {
        private:
            linalg::math_rowvector<typename RealIntTypeContext::Type> d_coeffs;
        
        public:
            void reset(RealIntTypeContext & c, const linalg::math_rowvector<typename IntTypeContext::Integer> & u, const linalg::math_matrix<typename RealIntTypeContext::Type> & dot_products)
            {
                d_coeffs.resize(u.size(), linalg::Initialize(c));
                for (unsigned i = 0; i < u.size(); ++i)
                {
                    setZero(d_coeffs[i]);
                    for (unsigned j = 0; j < u.size(); ++j)
                        d_coeffs[i] += arithmetic::convert(u[j], c) * dot_products(i, j);
                }
            }
        
            void operator() (RealIntTypeContext & c, typename RealIntTypeContext::Type & result, const linalg::math_rowvector<typename IntTypeContext::Integer> & v)
            // Computes dot product of <u> and <v> w.r.t <basis>.
            {
                setZero(result);
                for (unsigned i = 0; i < d_coeffs.size(); ++i)
                    result += arithmetic::convert(v[i], c) * d_coeffs[i];
            }
        };
    
        /*
          Note that we only store *half* the Voronoi cell, i.e. we leave away -v if v is included in it.
        */
    
        template<class IntTypeContext_, class RealIntTypeContext_, bool hasFractions>
        class CVPP;
    
        template<class IntTypeContext_, class RealIntTypeContext_>
        class CVPP<IntTypeContext_, RealIntTypeContext_, true>
        {
        private:
            DPLinForm dp;
            bool min_sign_val, sign_val;
            bool first;
            typename RealIntTypeContext_::Type min_val, val, vala;
            typename IntTypeContext_::Integer alpha;
        
        public:
            CVPP(IntTypeContext_ & ic, RealIntTypeContext_ & c)
                : min_val(c), val(c), vala(c)
            {
            }
        
            void operator()(IntTypeContext_ & ic, RealIntTypeContext_ & c, unsigned dim,
                            linalg::math_rowvector<typename IntTypeContext::Integer> & v,
                            const std::list<std::pair<linalg::math_rowvector<typename IntTypeContext::Integer>, typename RealIntTypeContext_::Type> > & V,
                            const linalg::math_matrix<typename RealIntTypeContext_::Type> & dot_products)
            {
                while (true)
                {
                    first = true;
                    dp.reset(c, v, dot_products);
                    typename std::list<std::pair<linalg::math_rowvector<typename IntTypeContext::Integer>, typename RealIntTypeContext_::Type> >::const_iterator min_i = V.end();
                    for (typename std::list<std::pair<linalg::math_rowvector<typename IntTypeContext::Integer>, typename RealIntTypeContext_::Type> >::const_iterator i = V.begin(); i != V.end(); ++i)
                    {
                        dp(c, val, i->first);
                        if (sign(val) < 0)
                        {
                            neg(val, val);
                            sign_val = true;
                        }
                        else
                            sign_val = false;
                        val /= i->second;
                        if (first || (val < min_val))
                        {
                            min_val = val;
                            min_sign_val = sign_val;
                            min_i = i;
                            first = false;
                        }
                    }
                    setOne(val);
                    if (min_val <= val)
                        break;
                    arithmetic::convert_ceil(alpha, min_val, ic);
                    // Note that alpha is always positive
                    if (isOne(alpha))
                    {
                        if (min_sign_val)
                            v += min_i->first;
                        else
                            v -= min_i->first;
                    }
                    else
                    {
                        if (min_sign_val)
                            v += alpha * min_i->first;
                        else
                            v -= alpha * min_i->first;
                    }
                }
            }
        };
    
        template<class IntTypeContext_, class RealIntTypeContext_>
        class CVPP<IntTypeContext_, RealIntTypeContext_, false>
        {
        private:
            DPLinForm dp;
            bool min_sign_val, sign_val;
            typename RealIntTypeContext_::Type min_val, min_val_d, val, vala;
            typename IntTypeContext_::Integer alpha, tmp;
        
        public:
            CVPP(IntTypeContext_ & ic, RealIntTypeContext_ & c)
                : min_val(c), min_val_d(c), val(c), vala(c)
            {
            }
        
            void operator()(IntTypeContext_ & ic, RealIntTypeContext_ & c, unsigned dim,
                            linalg::math_rowvector<typename IntTypeContext::Integer> & v,
                            const std::list<std::pair<linalg::math_rowvector<typename IntTypeContext::Integer>, typename RealIntTypeContext_::Type> > & V,
                            const linalg::math_matrix<typename RealIntTypeContext_::Type> & dot_products)
            {
                while (true)
                {
                    setZero(min_val_d);
                    setOne(min_val);
                    dp.reset(c, v, dot_products);
                    typename std::list<std::pair<linalg::math_rowvector<typename IntTypeContext::Integer>, typename RealIntTypeContext_::Type> >::const_iterator min_i = V.end();
                    for (typename std::list<std::pair<linalg::math_rowvector<typename IntTypeContext::Integer>, typename RealIntTypeContext_::Type> >::const_iterator i = V.begin(); i != V.end(); ++i)
                    {
                        dp(c, val, i->first);
                        if (sign(val) < 0)
                        {
                            neg(val, val);
                            sign_val = true;
                        }
                        else
                            sign_val = false;
                        if (val * min_val_d < min_val * i->second)
                        {
                            min_val = val;
                            min_val_d = i->second;
                            min_sign_val = sign_val;
                            min_i = i;
                        }
                    }
                    if (min_val <= min_val_d)
                        break;
                    euclideanDivision(alpha, tmp, min_val, min_val_d);
                    if (!isZero(tmp))
                        ++alpha;
                    // Note that alpha is always positive
                    if (isOne(alpha))
                    {
                        if (min_sign_val)
                            v += min_i->first;
                        else
                            v -= min_i->first;
                    }
                    else
                    {
                        if (min_sign_val)
                            v += alpha * min_i->first;
                        else
                            v -= alpha * min_i->first;
                    }
                }
            }
        };
    
        static unsigned getIndex(unsigned dim, const linalg::math_rowvector<typename IntTypeContext::Integer> & v)
        {
            unsigned res = 0;
            for (unsigned i = 0; i < dim; ++i)
            {
                res <<= 1;
                if (bit(v[i], 0))
                    res |= 1;
            }
            return res;
        }
    
        class SortIteratorSecond
        {
        public:
            template<class It>
            bool operator() (const It & x, const It & y) const
            {
                return x->second < y->second;
            }
        };
    
        class VecIntPairComparatorNeg
        {
        private:
            const std::vector<std::pair<linalg::math_rowvector<typename IntTypeContext::Integer>, typename RealIntTypeContext::Type> > & d_A;
        
        public:
            VecIntPairComparatorNeg(const std::vector<std::pair<linalg::math_rowvector<typename IntTypeContext::Integer>, typename RealIntTypeContext::Type> > & A)
                : d_A(A)
            {
            }
        
            bool operator() (unsigned x, unsigned y) const
            {
                int c = compare(d_A[x].second, d_A[y].second);
                return (c == 0) ? x < y : c > 0;
            }
        };
    
        void VEnum(IntTypeContext & ic, RealIntTypeContext & c, unsigned dim,
                   std::vector<std::pair<linalg::math_rowvector<typename IntTypeContext::Integer>, typename RealIntTypeContext::Type> > & A, std::vector<bool> & Aused,
                   const std::list<std::pair<linalg::math_rowvector<typename IntTypeContext::Integer>, typename RealIntTypeContext::Type> > & V, const std::vector<unsigned> & Vidx,
                   const std::pair<linalg::math_rowvector<typename IntTypeContext::Integer>, typename RealIntTypeContext::Type> & t)
        {
            for (unsigned i = 0; i < Aused.size(); ++i)
                Aused[i] = false;
        
            // Solve CVPP
            std::pair<linalg::math_rowvector<typename IntTypeContext::Integer>, typename RealIntTypeContext::Type> y0 = t;
            d_cvpp(ic, c, dim, y0.first, V, dot_products);
            DPLinForm dp;
            dp.reset(c, y0.first, dot_products);
            dp(c, y0.second, y0.first);
        
            // Store into A and initialize priority queue
            unsigned idx = getIndex(dim, y0.first);
            A[idx] = y0;
            Aused[idx] = true;
            std::set<unsigned, VecIntPairComparatorNeg>
                Q = std::set<unsigned, VecIntPairComparatorNeg>(VecIntPairComparatorNeg(A));
            Q.insert(idx);
        
            // Do while loop
            typename RealIntTypeContext::Type yn(c), n(c);
            while (!Q.empty())
            {
                // Get shortest element from queue
                typename std::set<unsigned, VecIntPairComparatorNeg>::iterator i = Q.end();
                --i;
                unsigned ui = *i;
                // std::pair<linalg::math_rowvector<Integer>, Integer> u = A[ui];
                // std::cout << "u = " << A[ui].first << " (" << A[ui].second << ")\n";
                Q.erase(i);
            
                // Iterate through all elements of u+V and u-V
                dp.reset(c, A[ui].first, dot_products);
                unsigned yii = 0;
                for (typename std::list<std::pair<linalg::math_rowvector<typename IntTypeContext::Integer>, typename RealIntTypeContext::Type> >::const_iterator yi = V.begin(); yi != V.end(); ++yi, ++yii)
                {
                    // Compute norm of y = u +- *yi
                    yn = A[ui].second + yi->second;
                    dp(c, n, yi->first);
                    n <<= 1;
                    bool wasneg = sign(n) < 0;
                    if (!wasneg)
                        // We start with smaller of the two
                        neg(n, n);
                    // Compare norms
                    if (A[ui].second <= yn + n)
                    {
                        yn += n;
                        // Trick 17: compute index for y.first = A[ui].first + yi->first using indices of
                        // A[ui].first (which is ui) and yi->first (which is Vidx[yii])
                        unsigned idx = ui ^ Vidx[yii];
                        // Compare to table
                        if (!Aused[idx])
                        {
                            // Update table and priority queue
                            Aused[idx] = true;
                            if (wasneg)
                                A[idx].first = A[ui].first + yi->first;
                            else
                                A[idx].first = A[ui].first - yi->first;
                            A[idx].second = yn;
                            Q.insert(idx);
                        }
                        else if (yn < A[idx].second)
                        {
                            // Update table and priority queue
                            i = Q.find(idx);
                            if (i != Q.end())
                                Q.erase(i);
                            if (wasneg)
                                A[idx].first = A[ui].first + yi->first;
                            else
                                A[idx].first = A[ui].first - yi->first;
                            A[idx].second = yn;
                            Q.insert(idx);
                        }
                    }
                    else
                    {
                        neg(n, n);
                        yn += n;
                        if (A[ui].second <= yn)
                        {
                            // Trick 17: compute index for y.first = A[ui].first + yi->first using indices of
                            // A[ui].first (which is ui) and yi->first (which is Vidx[yii])
                            unsigned idx = ui ^ Vidx[yii];
                            // Compare to table
                            if (!Aused[idx])
                            {
                                // Update table and priority queue
                                Aused[idx] = true;
                                if (wasneg)
                                    A[idx].first = A[ui].first - yi->first;
                                else
                                    A[idx].first = A[ui].first + yi->first;
                                A[idx].second = yn;
                                Q.insert(idx);
                            }
                            else if (yn < A[idx].second)
                            {
                                // Update table and priority queue
                                i = Q.find(idx);
                                if (i != Q.end())
                                    Q.erase(i);
                                if (wasneg)
                                    A[idx].first = A[ui].first - yi->first;
                                else
                                    A[idx].first = A[ui].first + yi->first;
                                A[idx].second = yn;
                                Q.insert(idx);
                            }
                        }
                    }
                }
            }
        }
    
        void VRelevant(IntTypeContext & ic, RealTypeContext & rc, RealIntTypeContext & c, unsigned dim,
                       std::list<std::pair<linalg::math_rowvector<typename IntTypeContext::Integer>, typename RealIntTypeContext::Type> > & dest,
                       const std::list<std::pair<linalg::math_rowvector<typename IntTypeContext::Integer>, typename RealIntTypeContext::Type> > & src,
                       const std::pair<linalg::math_rowvector<typename IntTypeContext::Integer>, typename RealIntTypeContext::Type> & b,
                       const typename RealTypeContext::Real & gs_norm_sq)
        {
            std::vector<std::pair<linalg::math_rowvector<typename IntTypeContext::Integer>, typename RealIntTypeContext::Type> > T, A; // "T[i].first.size() == 0" means "infinity"
            std::vector<bool> Tused, Aused;
            T.resize(1u << dim);
            A.resize(1u << dim);
            Tused.resize(1u << dim);
            Aused.resize(1u << dim);
            for (unsigned i = 0; i < T.size(); ++i)
            {
                T[i].first.resize(dim);
                T[i].second.setContext(c);
            }
            for (unsigned i = 0; i < T.size(); ++i)
            {
                A[i].first.resize(dim);
                A[i].second.setContext(c);
            }
            for (unsigned i = 0; i < T.size(); ++i)
                Tused[i] = false;
            unsigned infinity_count = 1u << dim;
            bool recompute_max = false;
            typename RealIntTypeContext::Type max(c), mm(c);
            setZero(max);
            std::vector<unsigned> srcidx;
            srcidx.resize(src.size());
            unsigned ii = 0;
            for (typename std::list<std::pair<linalg::math_rowvector<typename IntTypeContext::Integer>, typename RealIntTypeContext::Type> >::const_iterator i = src.begin(); i != src.end(); ++i, ++ii)
                srcidx[ii] = getIndex(dim, i->first);
            for (unsigned m = 0; ; ++m)
            {
                if (infinity_count == 0)
                {
                    arithmetic::RealContext rc;
                    if (arithmetic::convert(m, rc) * arithmetic::convert(gs_norm_sq, rc) >= arithmetic::convert(max, rc))
                        // work with arithmetic::Real instead of RealTypeContext::Real
                        break;
                }
                std::pair<linalg::math_rowvector<typename IntTypeContext::Integer>, typename RealIntTypeContext::Type> t = b;
                t.first *= arithmetic::convert(m, ic);
                arithmetic::convert(mm, m, c);
                mm = square(mm);
                t.second *= mm;
                VEnum(ic, c, dim, A, Aused, src, srcidx, t);
                for (unsigned i = 0; i < A.size(); ++i)
                    if (Aused[i])
                    {
                        unsigned idx = getIndex(dim, A[i].first);
                        if (!Tused[idx])
                        {
                            Tused[idx] = true;
                            if (infinity_count == (1u << dim))
                                max = A[i].second;
                            else
                                if (max < A[i].second)
                                    max = A[i].second;
                            --infinity_count;
                            T[idx] = A[i];
                            Tused[idx] = true;
                        }
                        else if (T[idx].second > A[i].second)
                        {
                            if (T[idx].second == max)
                                recompute_max = true;
                            T[idx] = A[i];
                        }
                    }
                if (recompute_max)
                {
                    recompute_max = false;
                    setZero(max);
                    for (unsigned i = 0; i < T.size(); ++i)
                        if (Tused[i])
                            if (max < T[i].second)
                                max = T[i].second;
                }
            }
            dest.clear();
            for (unsigned i = 0; i < T.size(); ++i)
                if (Tused[i])
                    dest.push_back(T[i]);
        }
    
        void VFilter(IntTypeContext & ic, RealIntTypeContext & c,
                     std::list<std::pair<linalg::math_rowvector<typename IntTypeContext::Integer>, typename RealIntTypeContext::Type> > & dest,
                     const std::list<std::pair<linalg::math_rowvector<typename IntTypeContext::Integer>, typename RealIntTypeContext::Type> > & src)
        {
            dest.clear();
            // Create sorted version of src by norm
            std::vector<typename std::list<std::pair<linalg::math_rowvector<typename IntTypeContext::Integer>, typename RealIntTypeContext::Type> >::const_iterator> sorted;
            sorted.reserve(src.size());
            for (typename std::list<std::pair<linalg::math_rowvector<typename IntTypeContext::Integer>, typename RealIntTypeContext::Type> >::const_iterator i = src.begin(); i != src.end(); ++i)
                if (!isZero(i->second))
                    sorted.push_back(i);
            std::sort(sorted.begin(), sorted.end(), SortIteratorSecond());
            // Callback for shortest element
            if (!sorted.empty())
            {
                // Callback available?
                if (!d_cf.empty() && d_ll.first)
                    d_cf(*d_ll.first, d_ll.second, sorted.front()->first);
            }
            // Go through sorted list
            typename RealIntTypeContext::Type n1(c), n2(c);
            DPLinForm dp;
            for (unsigned i = 0; i < sorted.size(); ++i)
            {
                dp.reset(c, sorted[i]->first, dot_products);
                bool fail1 = false, fail2 = false;
                for (typename std::list<std::pair<linalg::math_rowvector<typename IntTypeContext::Integer>, typename RealIntTypeContext::Type> >::iterator j = dest.begin(); j != dest.end(); ++j)
                {
                    // Note that ||2*sorted[i]->first +- j->first||^2 = 4 * sorted[i]->second + j->second +-
                    // 4 * <sorted[i]->first, j->first>. Therefore, comparing this norm to ||j->first||^2
                    // only involves comparing sorted[i]->second +- <sorted[i]->first, j->first> to 0.
                    n1 = sorted[i]->second;
                    dp(c, n2, j->first);
                    // Check whether ||2*sorted[i]->first - j->first||^2 > ||j->first||^2 is violated
                    if (n1 <= n2)
                    {
                        fail1 = true;
                        if (fail2)
                            break;
                    }
                    // Check whether ||2*sorted[i]->first + j->first||^2 > ||j->first||^2 is violated
                    if (n1 <= -n2)
                    {
                        fail2 = true;
                        if (fail1)
                            break;
                    }
                }
                if (!fail1 || !fail2)
                    dest.push_back(*sorted[i]);
            }
        }
    
        Verbose & d_v;
        CallbackFunction d_cf;
        std::list<std::pair<linalg::math_rowvector<typename IntTypeContext::Integer>, typename RealIntTypeContext::Type> > d_V, d_U;
        IntTypeContext & d_ic;
        RealIntTypeContext & d_c;
        RealTypeContext & d_rc;
        CVPP<IntTypeContext, RealIntTypeContext, RealIntTypeContext::is_realtype> d_cvpp;
    
    public:
        linalg::math_matrix<typename RealIntTypeContext::Type> dot_products;
        linalg::math_rowvector<typename RealTypeContext::Real> squared_GS_norms;
        std::pair<linalg::math_matrix<typename IntTypeContext::Integer> *, unsigned> d_ll;
        
        VoronoiCell(Verbose & v, const CallbackFunction & cf, unsigned dimension, RealIntTypeContext & context, RealTypeContext & rc, IntTypeContext & ic)
            : d_v(v), d_cf(cf), d_ic(ic), d_c(context), d_rc(rc), d_cvpp(d_ic, d_c), d_ll(NULL, 0)
        {
            dot_products.resize(dimension, dimension, linalg::Initialize(context));
            squared_GS_norms.resize(dimension, linalg::Initialize(rc));
        }
        
        void compute()
        {
            linalg::math_rowvector<typename IntTypeContext::Integer> e;
            e.resize(squared_GS_norms.size(), linalg::Initialize(d_ic));
            setOne(e[0]);
            d_V.push_back(std::make_pair(e, dot_products(0, 0)));
            d_v(LatticeReduction::VL_Information) << "Starting Voronoi Cell Computation";
        
            d_v(LatticeReduction::VL_Information) << 1 << ": " << d_V.size();
            for (unsigned k = 1; k < squared_GS_norms.size(); ++k)
            {
                for (unsigned i = 0; i < e.size(); ++i)
                    setZero(e[i]);
                setOne(e[k]);
                VRelevant(d_ic, d_rc, d_c, k + 1, d_U, d_V, std::make_pair(e, dot_products(k, k)), squared_GS_norms[k]);
                if (k + 1 == squared_GS_norms.size())
                {
                    d_v(LatticeReduction::VL_Information) << k + 1 << ": " << d_U.size() << " (unfiltered)";
                    break;
                }
                VFilter(d_ic, d_c, d_V, d_U);
                d_v(LatticeReduction::VL_Information) << k + 1 << ": " << d_V.size();
            }
            // Release part of the memory
            d_V.clear();
        }
        
        bool retrieveShortest(linalg::math_rowvector<typename IntTypeContext::Integer> & result)
        {
            typename std::list<std::pair<linalg::math_rowvector<typename IntTypeContext::Integer>, typename RealIntTypeContext::Type> >::iterator dest_i = d_U.end();
            for (typename std::list<std::pair<linalg::math_rowvector<typename IntTypeContext::Integer>, typename RealIntTypeContext::Type> >::iterator i = d_U.begin(); i != d_U.end(); ++i)
                if (((dest_i == d_U.end()) || (i->second < dest_i->second)) && !isZero(i->second))
                    dest_i = i;
            if (dest_i != d_U.end())
            {
                d_v(LatticeReduction::VL_Information) << "Found vector " << dest_i->first << " (squared norm " << dest_i->second << ")";
                result = dest_i->first;
                // Callback available?
                if (!d_cf.empty() && d_ll.first)
                {
                    try
                    {
                        d_cf(*d_ll.first, d_ll.second, result);
                    }
                    catch (LatticeReduction::stop_enumeration &)
                    {
                        d_v(LatticeReduction::VL_Chatter) << "Stopping enumeration";
                    }
                }
                return true;
            }
            else
                return false;
        }
    };
    
    template<class RealTypeContext, class IntTypeContext>
    class VoronoiCellWrapper
    {
    public:
        typedef boost::function<void(const linalg::math_matrix<typename IntTypeContext::Integer> & basis,
                                     int p,
                                     const linalg::math_rowvector<typename IntTypeContext::Integer> & vec)> CallbackFunction;
        
    private:
        bool d_int;
        STD_AUTO_PTR<VoronoiCell<RealTypeContext, IntTypeContext, RealTypeContext> > d_vc_real;
        STD_AUTO_PTR<VoronoiCell<RealTypeContext, IntTypeContext, IntTypeContext> > d_vc_int;
    
        // Disable copying and copy-constructing
        VoronoiCellWrapper(const VoronoiCellWrapper &);
        VoronoiCellWrapper & operator = (const VoronoiCellWrapper &);
    
    public:
        inline VoronoiCellWrapper(Verbose & v, const CallbackFunction & cf,
                                  Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned begin, unsigned end)
            : d_vc_real(), d_vc_int()
        {
            d_int = lattice.canGetIntegerVector() && (begin == 0);
            if (d_int)
            {
                d_vc_int.reset(new VoronoiCell<RealTypeContext, IntTypeContext, IntTypeContext>(v, cf, end - begin + 1, lattice.ic(), lattice.rc(), lattice.ic()));
                d_vc_int->d_ll = lattice.getMatrixIndex(begin);
                for (unsigned i = 0; i < d_vc_int->dot_products.rows(); ++i)
                {
                    for (unsigned j = 0; j <= i; ++j)
                    {
                        lattice.getDotProduct(d_vc_int->dot_products(i, j), begin + i, begin + j);
                        if (j < i)
                            d_vc_int->dot_products(j, i) = d_vc_int->dot_products(i, j);
                    }
                    d_vc_int->squared_GS_norms[i] = lattice.getNormSq(i);
                }
                try
                {
                    d_vc_int->compute();
                }
                catch (LatticeReduction::stop_enumeration &)
                {
                    v(LatticeReduction::VL_Chatter) << "Stopping enumeration";
                }
            }
            else
            {
                d_vc_real.reset(new VoronoiCell<RealTypeContext, IntTypeContext, RealTypeContext>(v, cf, end - begin + 1, lattice.rc(), lattice.rc(), lattice.ic()));
                d_vc_int->d_ll = lattice.getMatrixIndex(begin);
                for (unsigned i = 0; i < d_vc_real->dot_products.rows(); ++i)
                {
                    for (unsigned j = 0; j <= i; ++i)
                    {
                        lattice.computeDotProductProjected(d_vc_real->dot_products(i, j), begin, begin + i, begin + j);
                        if (j < i)
                            d_vc_real->dot_products(j, i) = d_vc_real->dot_products(i, j);
                    }
                    d_vc_real->squared_GS_norms[i] = lattice.getNormSq(i);
                }
                try
                {
                    d_vc_real->compute();
                }
                catch (LatticeReduction::stop_enumeration &)
                {
                    v(LatticeReduction::VL_Chatter) << "Stopping enumeration";
                }
            }
        }
        
        inline bool retrieveShortest(linalg::math_rowvector<typename IntTypeContext::Integer> & result, typename RealTypeContext::Real & /*bound*/)
        {
            return d_int ? d_vc_int->retrieveShortest(result) : d_vc_real->retrieveShortest(result);
        }
    };

    template<class RealTypeContext, class IntTypeContext>
    class VoronoiCellComputer
    {
    private:
        Verbose & d_verbose;
        typename VoronoiCellWrapper<RealTypeContext, IntTypeContext>::CallbackFunction d_cf;
        
    public:
        VoronoiCellComputer(Verbose & v, GaussianFactorComputer & /*gf*/, unsigned /*enumdimension*/)
            : d_verbose(v), d_cf(NULL)
        {
        }
        
        bool enumerate(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned begin, unsigned end,
                       linalg::math_rowvector<typename IntTypeContext::Integer> & result,
                       typename RealTypeContext::Real & bound)
        {
            VoronoiCellWrapper<RealTypeContext, IntTypeContext> vc(d_verbose, d_cf, lattice, begin, end);
            return vc.retrieveShortest(result, bound);
        }
        
        void setCallback(typename VoronoiCellWrapper<RealTypeContext, IntTypeContext>::CallbackFunction cf)
        {
            d_cf = cf;
        }
    };
}

#endif

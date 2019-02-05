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

#ifndef PLLL_INCLUDE_GUARD__ENUMIMPL_SCHNORRENUM_HPP
#define PLLL_INCLUDE_GUARD__ENUMIMPL_SCHNORRENUM_HPP

#include <queue>
#include "timer.hpp"

namespace plll
{
    template<class RealTypeContext, class IntTypeContext>
    class SchnorrEnumeratorImpl
    {
    public:
        typedef boost::function<void(const linalg::math_matrix<typename IntTypeContext::Integer> & basis,
                                     int p,
                                     const linalg::math_rowvector<typename IntTypeContext::Integer> & vec)> CallbackFunction;
        
    private:
        Updater<RealTypeContext, IntTypeContext> d_update;
        CallbackFunction d_cf;
        
        static inline void prepareBound(const Lattice<RealTypeContext, IntTypeContext> & lattice,
                                        unsigned begin, unsigned end, typename RealTypeContext::Real & bound)
        {
            const unsigned dim = end - begin + 1;
            
            // Prepare bound
            setOne(bound);
            for (unsigned i = begin; i <= end; ++i)
                bound *= lattice.getNormSq(i);
            power(bound, bound, arithmetic::convert(1, lattice.rc()) / arithmetic::convert(dim, lattice.rc()));
            bound *= arithmetic::convert(dim, lattice.rc());
            bound >>= 2;
        }
        
        struct Stage
        {
            long stage;
            unsigned t;
            linalg::math_rowvector<typename IntTypeContext::Integer> x;
            typename RealTypeContext::Real y;
            typename RealTypeContext::Real c;
            long delta;
            int delta2;
            
            Stage(long _stage, unsigned _t, const linalg::math_rowvector<typename IntTypeContext::Integer> _x,
                  const typename RealTypeContext::Real & _y, const typename RealTypeContext::Real & _c, long _delta, int _delta2)
                : stage(_stage), t(_t), x(_x), y(_y), c(_c), delta(_delta), delta2(_delta2)
            {
            }
            
            bool operator < (const Stage & s) const
            {
                return stage < s.stage;
            }
        };
        
        inline long getStage(const Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned t, const typename RealTypeContext::Real & beta)
        {
            typename RealTypeContext::Real tmp(lattice.rc());
            tmp = beta;
            tmp /= arithmetic::convert(t, lattice.rc());
            tmp = log2(tmp);
            tmp = -tmp;
            return arithmetic::convert_ceil<long>(tmp);
        }
        
        std::priority_queue<Stage> d_queue;
        
        void delayStage(const Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned t, const typename RealTypeContext::Real & beta,
                        const linalg::math_rowvector<typename IntTypeContext::Integer> & x, const typename RealTypeContext::Real & y,
                        const typename RealTypeContext::Real & c, long delta, int delta2)
        {
            long s = getStage(lattice, t, beta);
            if (d_verbose.yieldsOutput(LatticeReduction::VL_Chatter))
                d_verbose(LatticeReduction::VL_Chatter) << "Adding stage with s = " << s << ", t = " << t;
            d_queue.push(Stage(s, t, x, y, c, delta, delta2));
        }
        
        Verbose & d_verbose;
        GaussianFactorComputer & d_gaussianfactors;
        
    public:
        SchnorrEnumeratorImpl(Verbose & v, GaussianFactorComputer & gf, unsigned enumdimension)
            : d_update(), d_cf(NULL), d_verbose(v), d_gaussianfactors(gf)
        {
        }
        
        static inline unsigned minimalDimension()
        // Returns the minimal dimension for enumeration
        {
            return 0;
        }
        
        bool enumerate(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned begin, unsigned end,
                       linalg::math_rowvector<typename IntTypeContext::Integer> & result, typename RealTypeContext::Real & bound)
        // Finds a shortest vector in the lattice generated by the orthogonal projections of the vectors
        // A.row(begin) to A.row(end) into the orthogonal complement of the vectors A.row(0) to
        // A.row(begin-1). Uses Schnorr's fast (though usually not accurate) SVP solver.
        {
            enum { CRITICALSTAGES = 20 }; // if depth is less than this, just enumerate full subtree and
            // don't store things...
            
            const unsigned dim = end - begin + 1;
            d_gaussianfactors.computeHPF(dim);
            RealTypeContext & rc = lattice.rc();
            IntTypeContext & ic = lattice.ic();
            
            // Initialize updater
            prepareBound(lattice, begin, end, bound);
            d_update.initialize(d_verbose, lattice, begin, end, bound);
            
            // Prepare enumeration
            linalg::math_rowvector<typename IntTypeContext::Integer> x(dim);
            linalg::math_rowvector<typename RealTypeContext::Real> x_real(dim);
            for (unsigned i = 0; i < dim; ++i)
                setZero(x_real[i]);
            linalg::base_rowvector<long> delta(dim);
            linalg::base_rowvector<int> delta2(dim);
            setOne(x[0]);
            setOne(x_real[0]);
            delta[0] = 1;
            delta2[0] = 1;
            for (unsigned i = 1; i < dim; ++i)
                delta2[i] = -1;
            linalg::math_rowvector<typename RealTypeContext::Real> c(dim + 1), y(dim);
            for (unsigned i = 0; i <= dim; ++i)
                setZero(c[i]);
            for (unsigned i = 0; i < dim; ++i)
                setZero(y[i]);
            linalg::base_rowvector<int> r(dim + 1);
            linalg::math_matrix<typename RealTypeContext::Real> sigma(dim, dim);
            for (unsigned i = 0; i < dim; ++i)
            {
                r[i] = dim - 1;
                for (unsigned j = 0; j < dim; ++j)
                {
                    sigma(i, j).setContext(rc);
                    setZero(sigma(i, j));
                }
            }
            typename RealTypeContext::Real rhosq(rc), beta(rc);
            
            // Do enumeration
            long s = std::numeric_limits<long>::min()/2;
            unsigned t = 0, tmax = 0, tend = dim;
            bool noSolutionYet = true;
            typename RealTypeContext::Real tmp(rc);
            bool go_up;
            try
            {
                CPUTimer cputimer;
                cputimer.start();
                while (true)
                {
                    // Process current stage
                    while (true)
                    {
                        tmp = x_real[t] - y[t];
                        square(c[t], tmp);
                        c[t] *= lattice.getNormSq(begin + t);
                        c[t] += c[t + 1];
                        if (c[t] <= bound)
                        {
                            // Compute constants needed for Gaussian heuristic
                            rhosq = bound - c[t];
                            power(beta, rhosq, (long)t);
                            for (unsigned i = 0; i < t; ++i)
                                beta /= sqrt(lattice.getNormSq(begin + i));
                            sqrt(beta, beta);
                            beta *= arithmetic::convert(d_gaussianfactors.HPF(t), rc);
                            if (t == 0)
                            {
                                if (d_update(lattice, begin, end, result, bound, x, c[0], noSolutionYet))
                                {
                                    // Found a new shortest vector!
                                    if (!d_cf.empty())
                                    {
                                        std::pair<linalg::math_matrix<typename IntTypeContext::Integer> *, unsigned> l = lattice.getMatrixIndex(begin);
                                        if (l.first)
                                            d_cf(*l.first, l.second, result);
                                    }
                                    if (lattice.canGetIntegerVector())
                                    {
                                        Verbose::VerboseStream vv = d_verbose(LatticeReduction::VL_Information);
                                        vv << "Found new shortest vector: ";
                                        linalg::math_rowvector<typename IntTypeContext::Integer> vec;
                                        lattice.getLinearCombination(vec, begin, result);
                                        vv << vec;
                                        vv << ", bound = " << bound << ", norm = " << sqrt(c[0])
                                           << " [" << (long double)cputimer.elapsed() / (long double)cputimer.d_one_timeunit << " " << cputimer.d_timeunit_name << "]";
                                    }
                                    go_up = true;
                                }
                                else
                                    go_up = false;
                            }
                            else
                            {
                                // Decide whether to go down or put piece of tree into backlog
                                if ((t <= CRITICALSTAGES) || (beta >= (arithmetic::convert(t + 1, rc) >> s)))
                                {
                                    // Go stage down
                                    t = t - 1;
                                    if (t > 0)
                                    {
                                        if (r[t - 1] < r[t])
                                            r[t - 1] = r[t];
                                    }
                                    for (unsigned j = r[t]; j > t; --j) // this is the summation order considered in Pujol-Stehle
                                    {
                                        tmp = x_real[j] * lattice.getCoeff(begin + j, begin + t);
                                        sigma(t, j - 1) = sigma(t, j) + tmp;
                                    }
                                    y[t] = -sigma(t, t);
                                    // Find out x coordinate
                                    bool ru;
                                    arithmetic::convert_round(x[t], y[t], ru, ic);
                                    arithmetic::convert(x_real[t], x[t], rc);
                                    delta[t] = 0;
                                    delta2[t] = ru ? 1 : -1;
                                    continue;
                                }
                                else
                                {
                                    // Store stage
                                    delayStage(lattice, t, beta, x, y[t], c[t], delta[t], delta2[t]);
                                    go_up = false;
                                }
                            }
                        }
                        else
                            go_up = true;
                        // Possibly go stage up
                        if (go_up)
                        {
                            // Go stage up
                            if (++t >= tend)
                                // In this case, we are done!
                                break;
                            r[t - 1] = t;
                            if (t > tmax)
                                tmax = t;
                        }
                        // Modify x[t]
                        if (t == tmax)
                        {
                            // Just go straight up
                            ++x[t];
                            arithmetic::convert(x_real[t], x[t], rc);
                        }
                        else
                        {
                            // Go zig-zag course to next x[t]
                            delta2[t] = -delta2[t];
                            delta[t] = -delta[t] + delta2[t];
                            x[t] += arithmetic::convert(delta[t], lattice.ic());
                            arithmetic::convert(x_real[t], x[t], rc);
                        }
                    }
                    
                    // Continue with another stage
                    if (d_queue.empty())
                        break;
                    if (d_verbose.yieldsOutput(LatticeReduction::VL_Chatter))
                        d_verbose(LatticeReduction::VL_Chatter) << "Popping stage [s=" << d_queue.top().stage << ", t=" << d_queue.top().t << "]";
                    s = d_queue.top().stage;
                    tend = t = d_queue.top().t;
                    x = d_queue.top().x;
                    y[t] = d_queue.top().y;
                    c[t] = d_queue.top().c;
                    delta[t] = d_queue.top().delta;
                    delta2[t] = d_queue.top().delta2;
                    d_queue.pop();
                    // Prepare x_real
                    for (unsigned i = t; i <= tmax; ++i)
                        arithmetic::convert(x_real[i], x[i], rc);
                    // Go stage down
                    t = t - 1;
                    r[t] = tmax;
                    if (t > 0)
                    {
                        if (r[t - 1] < r[t])
                            r[t - 1] = r[t];
                    }
                    for (unsigned j = r[t]; j > t; --j) // this is the summation order considered in Pujol-Stehle
                    {
                        tmp = x_real[j] * lattice.getCoeff(begin + j, begin + t);
                        sigma(t, j - 1) = sigma(t, j) + tmp;
                    }
                    y[t] = -sigma(t, t);
                    // Find out x coordinate
                    bool ru;
                    arithmetic::convert_round(x[t], y[t], ru, ic);
                    arithmetic::convert(x_real[t], x[t], rc);
                    delta[t] = 0;
                    delta2[t] = ru ? 1 : -1;
                }
            }
            catch (LatticeReduction::stop_enumeration &)
            {
                d_verbose(LatticeReduction::VL_Chatter) << "Stopping enumeration";
            }
            
            // Nothing found?
            return !noSolutionYet;
        }
        
        void setCallback(CallbackFunction cf)
        {
            d_cf = cf;
        }
    };
}

#include "enumimpl-simpleenum.cpp"

namespace plll
{
    template<class RealTypeContext, class IntTypeContext>
    class SchnorrEnumerator
    {
    private:
        typedef typename helper::SelectFirstType<RealTypeContext::has_full_power,
                                                 SchnorrEnumeratorImpl<RealTypeContext, IntTypeContext>,
                                                 SimpleEnumerator<RealTypeContext, IntTypeContext> >::result EnumImpl;
        
        EnumImpl d_impl;
        
    public:
        SchnorrEnumerator(Verbose & v, GaussianFactorComputer & gf, unsigned enumdimension)
            : d_impl(v, gf, enumdimension)
        {
        }
        
        bool enumerate(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned begin, unsigned end,
                       linalg::math_rowvector<typename IntTypeContext::Integer> & result, typename RealTypeContext::Real & bound)
        // Finds a shortest vector in the lattice generated by the orthogonal projections of the vectors
        // A.row(begin) to A.row(end) into the orthogonal complement of the vectors A.row(0) to
        // A.row(begin-1).
        {
            return d_impl.enumerate(lattice, begin, end, result, bound);
        }
        
        void setCallback(typename EnumImpl::CallbackFunction cf)
        {
            d_impl.setCallback(cf);
        }
    };
}

#endif

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

#ifndef PLLL_INCLUDE_GUARD__ENUM_IMPL_CPP
#define PLLL_INCLUDE_GUARD__ENUM_IMPL_CPP

#include "myalloc.hpp"
#include <vector>
#include <algorithm>
#include <stack>
#include <limits>
#include <list>

#include "lll2-internal.hpp"
#include "lattice.cpp"
#include "enumimpl.hpp"
#include "enumimpl-boundsel.cpp"
#include "enumimpl-preproc.cpp"
#include "enumimpl-updater.cpp"

namespace plll
{
    template<class RealTypeContext, class IntTypeContext>
    class PruningHelper
    /*
      Computes bounds for various pruning settings
    */
    {
    private:
        unsigned d_begin;
        linalg::math_rowvector<typename RealTypeContext::Real> d_coeffs;
        StandardEnumSettings::PruningMethod d_method;
        double d_parameter;
    
        template<class RTC, class ITC, bool has_full_power>
        class ComputeSchnorrHoerner2;
    
        template<class RTC, class ITC>
        class ComputeSchnorrHoerner2<RTC, ITC, true>
        // coeff = 1/pi * [ 2^{-2 d_parameter} vol (L(b_1, ..., b_{d_dimension-k-1})^2) * (Gamma((d_dimension-k-1)/2 + 1))^2 ]^{1/(d_dimension-k-1)}
        {
        private:
            unsigned d_begin;
            unsigned d_dimension;
            double d_parameter;
            typename RTC::Real d_det, d_oopi;
        
        public:
            ComputeSchnorrHoerner2(unsigned begin, unsigned dimension, double parameter, const Lattice<RTC, ITC> & lattice)
                : d_begin(begin), d_dimension(dimension), d_parameter(parameter), d_det(lattice.rc()), d_oopi(lattice.rc())
            {
                setOne(d_det);
                typename RTC::Real t(lattice.rc());
                lattice.rc().getPi(t);
                setOne(d_oopi);
                d_oopi /= t;
            }
        
            void operator() (unsigned k, typename RTC::Real & coeff, const Lattice<RTC, ITC> & lattice)
            // Garuanteed to be called with k = 0, 1, ..., dimension-2, dimension-1 (in this order).
            {
                if (k == 0)
                {
                    setZero(coeff);
                    return;
                }
                // Extend determinant
                d_det *= lattice.getNormSq(d_begin + k - 1);
                // Compute coefficient
                setOne(coeff);
                power(coeff, arithmetic::convert(2, lattice.rc()), arithmetic::convert(-2 * d_parameter, lattice.rc()));
                coeff *= d_det;
                coeff *= square(arithmetic::convert(tgammal(0.5 * k + 1), lattice.rc()));
                power(coeff, coeff, arithmetic::convert(1.0 / k, lattice.rc()));
                coeff *= d_oopi;
            }
        };
    
        template<class RTC, class ITC>
        class ComputeSchnorrHoerner2<RTC, ITC, false>
        {
        private:
            unsigned d_begin;
            unsigned d_dimension;
            double d_parameter;
        
        public:
            ComputeSchnorrHoerner2(unsigned begin, unsigned dimension, double parameter, const Lattice<RTC, ITC> & lattice)
                : d_begin(begin), d_dimension(dimension), d_parameter(parameter)
            {
            }
        
            void operator() (unsigned k, typename RTC::Real & coeff, const Lattice<RTC, ITC> & lattice)
            // Garuanteed to be called with k = 0, 1, ..., dimension-2, dimension-1 (in this order).
            {
                // coeff = 1/pi * [ 2^{-d_parameter} vol L(b_1, ..., b_{d_dimension-k-1}) * Gamma((d_dimension-k-1)/2 + 1) ]^{2/(d_dimension-k-1)}
                setZero(coeff); // !!! ??? ...
            }
        };
        
        inline void getEpsOrZero(const RealTypeContext & rc, typename RealTypeContext::Real & epsilon, helper::BoolToType<true>)
        {
            rc.getEpsilon(epsilon);
        }
        
        inline void getEpsOrZero(const RealTypeContext & rc, typename RealTypeContext::Real & epsilon, helper::BoolToType<false>)
        {
            setZero(epsilon);
        }
        
    public:
        PruningHelper()
            : d_begin(0), d_method(StandardEnumSettings::PM_None), d_parameter(1.0)
        {
        }
    
        void setup(unsigned begin, unsigned dimension, const Lattice<RealTypeContext, IntTypeContext> & lattice,
                   StandardEnumSettings::PruningMethod pruning, double pruning_parameter)
        // Call this on beginning of enumerate. It will do some precomputations to enable fast update of the pruning bounds.
        {
            d_begin = begin;
            d_method = pruning;
            d_parameter = pruning_parameter;
            if (d_method == StandardEnumSettings::PM_None)
                return;
            d_coeffs.resize(dimension);
            for (unsigned i = 0; i < dimension; ++i)
                d_coeffs[i].setContext(lattice.rc());
            unsigned mid;
            switch (d_method)
            {
            default:
            case StandardEnumSettings::PM_Linear:
                for (unsigned i = 0; i < dimension; ++i)
                    arithmetic::convert(d_coeffs[i], (double)(dimension - i) / (double)dimension, lattice.rc());
                break;
            case StandardEnumSettings::PM_Step:
                mid = (dimension - 1) / 2;
                assert((0 < d_parameter) && (d_parameter <= 1));
                arithmetic::convert(d_coeffs[mid], d_parameter, lattice.rc());
                for (unsigned i = mid + 1; i < dimension; ++i)
                    d_coeffs[i] = d_coeffs[mid];
                for (unsigned i = 0; i < mid; ++i)
                    setOne(d_coeffs[i]);
                break;
            case StandardEnumSettings::PM_Piecewise:
                mid = (dimension - 1) / 2;
                assert((0 < d_parameter) && (d_parameter < 0.25));
                for (unsigned i = 0; i < mid; ++i)
                    arithmetic::convert(d_coeffs[i], 1.0 - 2.0 * (double)i * (1.0 - d_parameter) / (double)dimension, lattice.rc());
                for (unsigned i = mid; i < dimension; ++i)
                    arithmetic::convert(d_coeffs[i], 2.0 * (1.0 - (double)i / (double)dimension) * d_parameter, lattice.rc());
                break;
            case StandardEnumSettings::PM_SchnorrHoerner1:
                for (unsigned i = 0; i < dimension; ++i)
                    arithmetic::convert(d_coeffs[i], std::min(1.0, 1.05 * (double)(dimension - i) / (double)dimension), lattice.rc());
                break;
            case StandardEnumSettings::PM_SchnorrHoerner2:
            {
                ComputeSchnorrHoerner2<RealTypeContext, IntTypeContext, RealTypeContext::has_full_power>
                    cc(d_begin, d_coeffs.size(), d_parameter, lattice);
                for (unsigned i = 0; i < dimension; ++i)
                    cc(i, d_coeffs[i], lattice);
            }
            break;
            case StandardEnumSettings::PM_GNR110_Poly8:
            {
                typename RealTypeContext::Real x(lattice.rc()), res(lattice.rc()), one(lattice.rc());
                setOne(one);
                // Prepare polynomial
                linalg::math_rowvector<typename RealTypeContext::Real> poly;
                poly.resize(9, linalg::Initialize(lattice.rc()));
                /* These values are taken from the GPU implementation available at
                   http://homes.esat.kuleuven.be/~jhermans/gpuenum/index.html */
                arithmetic::convert(poly[0],  0.000914465, lattice.rc());
                arithmetic::convert(poly[1],  0.0400812, lattice.rc());
                arithmetic::convert(poly[2], -0.00424356, lattice.rc());
                arithmetic::convert(poly[3],  0.00022931, lattice.rc());
                arithmetic::convert(poly[4], -6.91288e-06, lattice.rc());
                arithmetic::convert(poly[5],  1.21218e-07, lattice.rc());
                arithmetic::convert(poly[6], -1.20165e-09, lattice.rc());
                arithmetic::convert(poly[7],  6.20066e-12, lattice.rc());
                arithmetic::convert(poly[8], -1.29185e-14, lattice.rc());
            
                // Evaluate polynomial using (Horner's method)
                for (unsigned i = 0; i < dimension; ++i)
                {
                    // Determine evaluation point
                    arithmetic::convert(x, (dimension - i) * 110, lattice.rc());
                    x /= arithmetic::convert(dimension, lattice.rc());
                    // Evaluate
                    res = poly[8];
                    for (int j = 7; j >= 0; --j)
                    {
                        res *= x;
                        res += poly[j];
                    }
                    // Clip result to range [eps, 1]
                    if (res > one)
                        res = one;
                    typename RealTypeContext::Real eps(lattice.rc());
                    getEpsOrZero(lattice.rc(), eps, helper::BoolToType<RealTypeContext::has_constants>());
                    if (res < eps)
                        res = eps;
                    // Store result
                    d_coeffs[i] = res;
                }
            }
            break;
            }
        }
    
        void computeBounds(linalg::math_rowvector<typename RealTypeContext::Real> & bounds, unsigned dimension,
                           const typename RealTypeContext::Real & bound) const
        // Given a global bound, will compute the pruning bounds for all stages.
        {
            switch (d_method)
            {
            case StandardEnumSettings::PM_None:
                for (unsigned i = 0; i < dimension; ++i)
                    bounds[i] = bound;
                break;
            case StandardEnumSettings::PM_Step:
            {
                unsigned mid = (dimension - 1) / 2;
                bounds[0] = bound * d_coeffs[0];
                for (unsigned i = 1; i <= mid; ++i)
                    bounds[i] = bounds[0];
                for (unsigned i = mid + 1; i < dimension; ++i)
                    bounds[i] = bound;
            }
            break;
            case StandardEnumSettings::PM_Linear:
            case StandardEnumSettings::PM_Piecewise:
            case StandardEnumSettings::PM_SchnorrHoerner1:
            case StandardEnumSettings::PM_GNR110_Poly8:
                for (unsigned i = 0; i < dimension; ++i)
                    bounds[i] = bound * d_coeffs[i];
                break;
            case StandardEnumSettings::PM_SchnorrHoerner2:
                for (unsigned i = 0; i < dimension; ++i)
                    bounds[i] = bound - d_coeffs[i];
                break;
            }
        }
    };
}

#include "enumimpl-parallelserial.cpp"
#include "enumimpl-simpleenum.cpp"
#include "enumimpl-schnorrenum.cpp"
#include "enumimpl-voronoi.cpp"
#include "enumimpl-sieve.cpp"

namespace plll
{
    template<class RealTypeContext, class IntTypeContext>
    void Enumerator<RealTypeContext, IntTypeContext>::fullInit()
    {
        d_serial_enum = new SerialEnumerator<RealTypeContext, IntTypeContext>(d_verbose, d_gaussian_factors, d_enumdim);
        if (d_maxthreads != 1)
            d_parallel_enum = new ParallelEnumerator<RealTypeContext, IntTypeContext>(d_verbose, d_gaussian_factors, d_enumdim, d_maxthreads);
        d_schnorr_enum = new SchnorrEnumerator<RealTypeContext, IntTypeContext>(d_verbose, d_gaussian_factors, d_enumdim);
        d_sieve = new Sieve<RealTypeContext, IntTypeContext>(d_verbose, d_gaussian_factors, d_enumdim);
        d_voronoi = new VoronoiCellComputer<RealTypeContext, IntTypeContext>(d_verbose, d_gaussian_factors, d_enumdim);
        
        if (d_ecf || d_ecf2)
            setCallback(d_ecf, d_ecf2);
    }
    
    template<class RealTypeContext, class IntTypeContext>
    template<class Enum, class CallbackF>
    Enumerator<RealTypeContext, IntTypeContext>::Enumerator(Workspace<RealTypeContext, IntTypeContext, Enum, CallbackF> & workspace,
                                                            unsigned enumdimension, unsigned max_threads,
                                                            LatticeReduction::EnumCallbackFunction ecf, LatticeReduction::EnumCallbackFunction_LI ecf2)
        : d_enumdim(enumdimension), d_maxthreads(max_threads), d_ecf(NULL), d_ecf2(NULL),
          d_parallel_enum(NULL), d_serial_enum(NULL), d_sieve(NULL), d_simple_enum(NULL),
        d_schnorr_enum(NULL), d_voronoi(NULL),
        d_verbose(workspace.verbose()), d_gaussian_factors(workspace.gaussianFactors())
    {
        // Check concurrency
        if (d_maxthreads == 0)
            d_maxthreads = boost::thread::hardware_concurrency();
        if (d_maxthreads == 0)
            d_maxthreads = 1;
        
        // Create enumerators (if necessary)
        d_simple_enum = new SimpleEnumerator<RealTypeContext, IntTypeContext>(d_verbose, d_gaussian_factors, d_enumdim);
        if (enumdimension >= c_enum_threshold)
            fullInit();
        
        if (ecf || ecf2)
            setCallback(ecf, ecf2);
    }
    
    template<class RealTypeContext, class IntTypeContext>
    bool Enumerator<RealTypeContext, IntTypeContext>::enumerate(Lattice<RealTypeContext, IntTypeContext> & lattice,
                                                                unsigned begin, unsigned end, linalg::math_rowvector<typename IntTypeContext::Integer> & result,
                                                                typename RealTypeContext::Real & bound, LatticeReduction::SVPMode algorithm,
                                                                const StandardEnumSettings & settings, bool DontFallback)
    // Finds a shortest vector in the lattice generated by the orthogonal projections of the vectors
    // A.row(begin) to A.row(end) into the orthogonal complement of the vectors A.row(0) to
    // A.row(begin-1). Uses the Kannan-Schnorr-Euchner enumeration method.
    {
        ++lattice.getStatistics().enumcalls;
        // Check dimension threshold
        unsigned dim = end - begin + 1;
        bool retval = false;
        if ((dim < c_enum_threshold) && !DontFallback)
            // Fallback: use simple enumerator
            retval = d_simple_enum->enumerate(lattice, begin, end, result, bound);
        else
        {
            // Check if enumeration was initialized
            if (d_serial_enum == NULL)
                fullInit();
            // Run
            switch (algorithm)
            {
            case LatticeReduction::SVP_ParallelKannanSchnorrEuchner:
                if (d_parallel_enum)
                {
                    retval = d_parallel_enum->enumerate(lattice, begin, end, result, bound, settings);
                    break;
                }
                // else: d_serial_enum
            case LatticeReduction::SVP_KannanSchnorrEuchner:
                retval = d_serial_enum->enumerate(lattice, begin, end, result, bound, settings);
                break;
            case LatticeReduction::SVP_SchnorrFast:
                retval = d_schnorr_enum->enumerate(lattice, begin, end, result, bound);
                break;
            case LatticeReduction::SVP_VoronoiCellSVP:
                retval = d_voronoi->enumerate(lattice, begin, end, result, bound);
                break;
            case LatticeReduction::SVP_ListSieve:
                retval = d_sieve->enumerate(Sieve<RealTypeContext, IntTypeContext>::A_ListSieve, lattice, begin, end, result, bound);
                break;
            case LatticeReduction::SVP_ListSieveBirthday:
                retval = d_sieve->enumerate(Sieve<RealTypeContext, IntTypeContext>::A_ListSieveBirthday, lattice, begin, end, result, bound);
                break;
            case LatticeReduction::SVP_GaussSieve:
                retval = d_sieve->enumerate(Sieve<RealTypeContext, IntTypeContext>::A_GaussSieve, lattice, begin, end, result, bound);
                break;
            default:
                assert(!"Invalid enumeration method!");
            }
        }
        if (!retval)
        {
            if (d_verbose.yieldsOutput(LatticeReduction::VL_Chatter))
                d_verbose(LatticeReduction::VL_Chatter) << "Enumeration failed!";
            ++lattice.getStatistics().enumfails;
        }
        else
        {
            if (d_verbose.yieldsOutput(LatticeReduction::VL_Chatter))
                d_verbose(LatticeReduction::VL_Chatter) << "Enumeration yielded " << result;
        }
        return retval;
    }
    
    template<typename Enumerator>
    void enumSetCallback_impl(Enumerator & e, LatticeReduction::EnumCallbackFunction cf,
                              LatticeReduction::EnumCallbackFunction_LI cf2,
                              arithmetic::IntegerContext *)
    {
        e.setCallback(cf);
    }
    
    template<typename Enumerator>
    void enumSetCallback_impl(Enumerator & e, LatticeReduction::EnumCallbackFunction cf,
                              LatticeReduction::EnumCallbackFunction_LI cf2,
                              arithmetic::NIntContext<long int> *)
    {
        e.setCallback(cf2);
    }
    
    template<typename IntTypeContext, typename Enumerator>
    void enumSetCallback(Enumerator & e, LatticeReduction::EnumCallbackFunction cf,
                         LatticeReduction::EnumCallbackFunction_LI cf2)
    {
        enumSetCallback_impl(e, cf, cf2, static_cast<IntTypeContext *>(NULL));
    }
    
    template<class RealTypeContext, class IntTypeContext>
    void Enumerator<RealTypeContext, IntTypeContext>::setCallback(LatticeReduction::EnumCallbackFunction cf,
                                                                  LatticeReduction::EnumCallbackFunction_LI cf2)
    {
        if ((cf != NULL) || (cf2 != NULL))
        {
            if ((cf == NULL) && helper::CompareTypes<IntTypeContext, arithmetic::IntegerContext>::equal)
                cf = boost::bind(EnumCallbackAdaptor, _1, _2, _3, cf2);
            if ((cf2 == NULL) && helper::CompareTypes<IntTypeContext, arithmetic::NIntContext<long int> >::equal)
                cf2 = boost::bind(EnumCallbackAdaptor_LI, _1, _2, _3, cf);
        }
        d_ecf = cf;
        d_ecf2 = cf2;
        if (d_parallel_enum)
            enumSetCallback<IntTypeContext>(*d_parallel_enum, cf, cf2);
        if (d_serial_enum)
            enumSetCallback<IntTypeContext>(*d_serial_enum, cf, cf2);
        if (d_sieve)
            enumSetCallback<IntTypeContext>(*d_sieve, cf, cf2);
        if (d_simple_enum)
            enumSetCallback<IntTypeContext>(*d_simple_enum, cf, cf2);
        if (d_schnorr_enum)
            enumSetCallback<IntTypeContext>(*d_schnorr_enum, cf, cf2);
    }
}

#endif

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

#ifndef PLLL_INCLUDE_GUARD__ENUMIMPL_UPDATER_CPP
#define PLLL_INCLUDE_GUARD__ENUMIMPL_UPDATER_CPP

namespace plll
{
    template<class RealTypeContext, class IntTypeContext, bool exact, bool has_constants>
    class UpdaterImpl;

    template<class RealTypeContext, class IntTypeContext, bool has_constants>
    class UpdaterImpl<RealTypeContext, IntTypeContext, true, has_constants>
    {
    public:
        void initialize(Verbose & v, Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned begin, unsigned end,
                        typename RealTypeContext::Real & bound)
        {
        }
    
        bool operator() (const Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned begin, unsigned end,
                         linalg::math_rowvector<typename IntTypeContext::Integer> & result, typename RealTypeContext::Real & bound,
                         const linalg::math_rowvector<typename IntTypeContext::Integer> & vec, const typename RealTypeContext::Real & len,
                         bool & noSolutionYet) const
        {
            noSolutionYet = false;
            result = vec;
            bound = len;
            return true;
        }
    };
    
    template<class RealTypeContext, class IntTypeContext>
    class UpdaterImpl<RealTypeContext, IntTypeContext, false, true>
    {
    private:
        typename RealTypeContext::Real d_A, d_B;
    
    public:
        void initialize(Verbose & v, Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned begin, unsigned end,
                        typename RealTypeContext::Real & bound)
        {
            RealTypeContext & rc = lattice.rc();

            while (true)
            {
                typename RealTypeContext::Real eps(rc);
                rc.getEpsilon(eps);
            
                // Assuming (eta,delta)-LLL-reduced with: eta >= 1/5, eta^2 + 1 / (1 + eta)^2 < delta < 1
                // We assume eta = 0.5 and delta = 0.999. This is "usually" correct.
            
                // Assuming that exact Gram-Schmidt coefficients as rational numbers were converted, we get
                // kappa <= 3.
            
                // R = (1 + kappa * eps) * max_i ||b_i^*||^2
                // R = (1 + 3 * eps)
                typename RealTypeContext::Real R(lattice.getNormSq(begin));
                for (unsigned i = begin + 1; i <= end; ++i)
                    if (R < lattice.getNormSq(i))
                        R = lattice.getNormSq(i);
                R *= arithmetic::convert(1, rc) + arithmetic::convert(3, rc) * eps;
            
                // alpha = 1 / sqrt(delta - eta^2)
                typename RealTypeContext::Real alpha(rc);
                alpha = arithmetic::convert(1, rc) / sqrt(arithmetic::convert(0.999 - 0.25, rc));
            
                // rho = (1 + eta) * alpha
                typename RealTypeContext::Real rho(rc);
                rho = arithmetic::convert(1 + 0.5, rc) * alpha;
            
                // C_2 = (kappa + 2) / (alpha - 1) + (kappa + 4) / (rho - 1)
                typename RealTypeContext::Real C2(rc);
                C2 = arithmetic::convert(3 + 2, rc) / (alpha - arithmetic::convert(1, rc)) + arithmetic::convert(3 + 4, rc) / (rho - arithmetic::convert(1, rc));
            
                // C_3 = 2 * alpha * (2 + kappa + 2 * C_2) / (1 + eta - alpha)
                typename RealTypeContext::Real C3(rc);
                C3 = arithmetic::convert(2, rc) * alpha * (arithmetic::convert(2 + 3, rc) + arithmetic::convert(2, rc) * C2) / (arithmetic::convert(1 + 0.5, rc) - alpha);
            
                // lambda = \lambda_1(\Lambda)
                // We use: lambda = sqrt(bound)
            
                // Assume C_2 * rho^dim * eps <= 0.01
                unsigned long dim = end - begin + 1;
                if (C2 * power(rho, dim) * eps > arithmetic::convert(0.01, rc))
                {
                    v(LatticeReduction::VL_Warning) << "Warning: Floating-Point KFP Assertion 1 is invalid!";
                    if (RealTypeContext::is_variable_precision)
                    {
                        typename RealTypeContext::Real prec_r(rc);
                        prec_r = log2(C2) + log2(rho) * arithmetic::convert(dim, rc) + arithmetic::convert(7, rc); // 7 > log_2(100) = -log_2(0.01)
                        long prec = arithmetic::convert_ceil<long>(prec_r) + 8; // add some backup
                        prec = (prec + 63) & ~63; // round to multiple of 64
                        v(LatticeReduction::VL_Warning) << "Increasing precision from " << rc.getRealPrecision() << " to " << prec << "...";
                        rc.setRealPrecision(prec);
                        lattice.changeOfPrecision();
                        continue; // restart loop
                    }
                }
            
                // Now bound should be at least (1 + 2 * dim * eps) * lambda^2 + C_3 * rho^dim * eps * R
                d_A.setContext(rc);
                d_A = (arithmetic::convert(1, rc) + arithmetic::convert(2 * dim, rc) * eps);
                d_B.setContext(rc);
                d_B = C3 * power(rho, dim) * eps * R;
                bool verbose = (((end - begin) > 52) && v.yieldsOutput(LatticeReduction::VL_Information)) || v.yieldsOutput(LatticeReduction::VL_Chatter);
                if ((d_B << 7) > bound)
                {
                    typename RealTypeContext::Real old = d_B;
                    d_B = bound >> 7;
                    v(LatticeReduction::VL_Information) << "Restricting B=" << old << " to " << d_B << "!";
                }
                if (verbose)
                {
                    v(LatticeReduction::VL_Information) << "===== FP-KFP ===== [" << begin << "," << end << "]";
                    v(LatticeReduction::VL_Information) << "A = " << d_A;
                    v(LatticeReduction::VL_Information) << "B = " << d_B;
                    v(LatticeReduction::VL_Information) << "bound_old = " << bound;
                }
                bound = d_A * bound + d_B;
                if (verbose)
                    v(LatticeReduction::VL_Information) << "bound_new = " << bound;
                break;
            }
        }
    
        bool operator() (const Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned begin, unsigned end,
                         linalg::math_rowvector<typename IntTypeContext::Integer> & result, typename RealTypeContext::Real & bound,
                         const linalg::math_rowvector<typename IntTypeContext::Integer> & vec, const typename RealTypeContext::Real & len,
                         bool & noSolutionYet) const
        {
            typename RealTypeContext::Real bound2(lattice.rc());
            lattice.computeProjectionLength(bound2, begin, begin, vec);
//        std::cout << vec << " " << bound2 << " " << len << "\n";
            if (!noSolutionYet)
            {
                typename RealTypeContext::Real bound1(lattice.rc());
                lattice.computeProjectionLength(bound1, begin, begin, result);
//            std::cout << "  already had: " << result << " " << bound1 << "\n";
                if (bound2 > bound1)
                    return false;
            }
            noSolutionYet = false;
            result = vec;
//        std::cerr << "  NEW SHORTEST!\n";
        
            /*
              Because Pujol and Stehle did not provide a good choice of eps'' in their paper, we will
              not update the bound here. This screws up the running time, but at least currently it
              makes our live simpler... !!! ???
          
              NEW: do a "stupid" (possibly incorrect!!!) approach: use same constants as before!
            */
            bound2 = d_A * bound2 + d_B;
            if (bound2 < bound)
                bound = bound2;
            return true;
        }
    };
    
    // Currently no implementation for non-exact real types which don't have constants (like
    // epsilon, the machine constant).
    
    template<class RealTypeContext, class IntTypeContext>
    class Updater
    {
    private:
        UpdaterImpl<RealTypeContext, IntTypeContext, RealTypeContext::is_exact, RealTypeContext::has_constants> d_impl;
        
    public:
        // Might adjust the bound. Return true in case precision was changed
        inline void initialize(Verbose & v, Lattice<RealTypeContext, IntTypeContext> & lattice,
                               unsigned begin, unsigned end, typename RealTypeContext::Real & bound)
        {
            d_impl.initialize(v, lattice, begin, end, bound);
        }
    
        inline bool operator() (const Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned begin, unsigned end,
                                linalg::math_rowvector<typename IntTypeContext::Integer> & result, typename RealTypeContext::Real & bound,
                                const linalg::math_rowvector<typename IntTypeContext::Integer> & vec, const typename RealTypeContext::Real & len,
                                bool & noSolutionYet) const
        {
            return d_impl(lattice, begin, end, result, bound, vec, len, noSolutionYet);
        }
    };
}

#endif

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

#ifndef PLLL_INCLUDE_GUARD__GRAMSCHMIDT_CLASSIC_INTEGER_CPP
#define PLLL_INCLUDE_GUARD__GRAMSCHMIDT_CLASSIC_INTEGER_CPP

namespace plll
{
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Gram-Schmidt computers: exact (integer arithmetic, rounding only in last step) classical Gram-Schmidt
    
    template<class RealTypeContext, class IntTypeContext>
    class GSClassicInteger
    // classical Gram-Schmidt orthogonalization using integer arithmetic, rounding is only done in the last step
    {
    private:
        RealTypeContext & d_rc;
        IntTypeContext & d_ic;
        linalg::math_matrix<typename IntTypeContext::Integer> & d_A;
    
        unsigned d_level, d_vectors_used;
        linalg::math_matrix<typename RealTypeContext::Real> d_coeffs;
        linalg::math_rowvector<typename RealTypeContext::Real> d_sqnorms;
        linalg::base_rowvector<bool> d_dependent;
    
        template<class RTC, class ITC, bool exact>
        class Storage;
    
        template<class RTC, class ITC>
        class Storage<RTC, ITC, true>
        {
        private:
            GSClassicInteger<RTC, ITC> & d_base;
        
        public:
            typedef typename RTC::Real Real;
        
            RTC & d_rc;
        
            static void copyFrom(const Storage &)
            {
            }
        
            Storage(GSClassicInteger<RTC, ITC> & base)
                : d_base(base), d_rc(d_base.d_rc)
            {
            }
        
            Storage(GSClassicInteger<RTC, ITC> & base, const Storage &)
                : d_base(base), d_rc(d_base.d_rc)
            {
            }
        
            static inline void reassign(unsigned up_to)
            {
            }
        
            static inline void resize(unsigned)
            {
            }
        
            inline void swapCoeffs(unsigned i, unsigned j, unsigned k, unsigned l)
            {
                arithmetic::swap(d_base.d_coeffs(i, j), d_base.d_coeffs(k, l));
            }
        
            inline void swapCoeffsRow(unsigned i, unsigned j)
            {
                linalg::swap(d_base.d_coeffs.row(i), d_base.d_coeffs.row(j));
            }
        
            inline void swapCoeffsCol(unsigned i, unsigned j)
            {
                linalg::swap(d_base.d_coeffs.col(i), d_base.d_coeffs.col(j));
            }
        
            inline void swapSqNorms(unsigned i, unsigned j)
            {
                arithmetic::swap(d_base.d_sqnorms[i], d_base.d_sqnorms[j]);
            }
        
            inline void setCoeffZero(unsigned i, unsigned j)
            {
                setZero(d_base.d_coeffs(i, j));
            }
        
            inline void setCoeff(unsigned i, unsigned j, const arithmetic::Integer & num, const arithmetic::Integer & denom)
            {
                convert(d_base.d_coeffs(i, j), num, denom, d_rc);
            }
        
            inline void incCoeff(unsigned i, unsigned j)
            {
                ++d_base.d_coeffs(i, j);
            }
        
            inline void decCoeff(unsigned i, unsigned j)
            {
                --d_base.d_coeffs(i, j);
            }
        
            inline void setSqNorm(unsigned i, const arithmetic::Integer & num, const arithmetic::Integer & denom)
            {
                convert(d_base.d_sqnorms[i], num, denom, d_rc);
            }
        
            inline void setCoeff(unsigned i, unsigned j, const Real & v)
            {
                d_base.d_coeffs(i, j) = v;
            }
        
            inline void setSqNorm(unsigned i, const Real & v)
            {
                d_base.d_sqnorms[i] = v;
            }
        
            inline void getRoundedCoeff(arithmetic::Integer & res, unsigned i, unsigned j) const
            {
                convert_round(res, d_base.d_coeffs(i, j), arithmetic::IntegerContext());
            }
        
            inline const Real & getCoeff(unsigned i, unsigned j) const
            {
                return d_base.d_coeffs(i, j);
            }
        
            inline const Real & getSqNorm(unsigned i) const
            {
                return d_base.d_sqnorms[i];
            }
        };
    
        template<class RTC, class ITC>
        class Storage<RTC, ITC, false>
        {
        private:
            GSClassicInteger<RTC, ITC> & d_base;
        
            linalg::math_matrix<arithmetic::Rational> d_coeffs;
            linalg::math_rowvector<arithmetic::Rational> d_sqnorms;
        
        public:
            typedef arithmetic::Rational Real;
        
            arithmetic::RationalContext d_rc;
        
            void copyFrom(const Storage & s)
            {
                d_coeffs = s.d_coeffs;
                d_sqnorms = s.d_sqnorms;
            }
        
            Storage(GSClassicInteger<RTC, ITC> & base)
                : d_base(base)
            {
                d_sqnorms.resize(d_base.d_A.rows());
                d_coeffs.resize(d_base.d_A.rows(), d_base.d_A.rows());
            }
        
            Storage(GSClassicInteger<RTC, ITC> & base, const Storage & storage)
                : d_base(base), d_coeffs(storage.d_coeffs), d_sqnorms(storage.d_sqnorms)
            {
            }
        
            inline void reassign(unsigned up_to)
            {
                for (unsigned i = 0; i < up_to; ++i)
                {
                    for (unsigned j = 0; j < i; ++j)
                        convert(d_base.d_coeffs(i, j), d_coeffs(i, j), d_base.d_rc);
                    convert(d_base.d_sqnorms[i], d_sqnorms[i], d_base.d_rc);
                }
            }
        
            inline void resize(unsigned z)
            {
                d_coeffs.resize(z, z);
                d_sqnorms.resize(z);
            }
        
            inline void swapCoeffs(unsigned i, unsigned j, unsigned k, unsigned l)
            {
                arithmetic::swap(d_coeffs(i, j), d_coeffs(k, l));
                arithmetic::swap(d_base.d_coeffs(i, j), d_base.d_coeffs(k, l));
            }
        
            inline void swapCoeffsRow(unsigned i, unsigned j)
            {
                linalg::swap(d_coeffs.row(i), d_coeffs.row(j));
                linalg::swap(d_base.d_coeffs.row(i), d_base.d_coeffs.row(j));
            }
        
            inline void swapCoeffsCol(unsigned i, unsigned j)
            {
                linalg::swap(d_coeffs.col(i), d_coeffs.col(j));
                linalg::swap(d_base.d_coeffs.col(i), d_base.d_coeffs.col(j));
            }
        
            inline void swapSqNorms(unsigned i, unsigned j)
            {
                arithmetic::swap(d_sqnorms[i], d_sqnorms[j]);
                arithmetic::swap(d_base.d_sqnorms[i], d_base.d_sqnorms[j]);
            }
        
            inline void setCoeffZero(unsigned i, unsigned j)
            {
                setZero(d_coeffs(i, j));
                setZero(d_base.d_coeffs(i, j));
            }
        
            inline void setCoeff(unsigned i, unsigned j, const arithmetic::Integer & num, const arithmetic::Integer & denom)
            {
                d_coeffs(i, j) = arithmetic::Rational(num, denom);
                convert(d_base.d_coeffs(i, j), d_coeffs(i, j), d_base.d_rc);
            }
        
            inline void incCoeff(unsigned i, unsigned j)
            {
                ++d_coeffs(i, j);
                convert(d_base.d_coeffs(i, j), d_coeffs(i, j), d_base.d_rc);
            }
        
            inline void decCoeff(unsigned i, unsigned j)
            {
                --d_coeffs(i, j);
                convert(d_base.d_coeffs(i, j), d_coeffs(i, j), d_base.d_rc);
            }
        
            inline void setSqNorm(unsigned i, const arithmetic::Integer & num, const arithmetic::Integer & denom)
            {
                d_sqnorms[i] = arithmetic::Rational(num, denom);
                convert(d_base.d_sqnorms[i], d_sqnorms[i], d_base.d_rc);
            }
        
            inline void setCoeff(unsigned i, unsigned j, const Real & v)
            {
                d_coeffs(i, j) = v;
                convert(d_base.d_coeffs(i, j), d_coeffs(i, j), d_base.d_rc);
            }
        
            inline void setSqNorm(unsigned i, const Real & v)
            {
                d_sqnorms[i] = v;
                convert(d_base.d_sqnorms[i], d_sqnorms[i], d_base.d_rc);
            }
        
            inline void getRoundedCoeff(arithmetic::Integer & res, unsigned i, unsigned j) const
            {
                convert_round(res, d_coeffs(i, j), arithmetic::IntegerContext());
            }
        
            inline const Real & getCoeff(unsigned i, unsigned j) const
            {
                return d_coeffs(i, j);
            }
        
            inline const Real & getSqNorm(unsigned i) const
            {
                return d_sqnorms[i];
            }
        };
    
        typedef Storage<RealTypeContext, IntTypeContext, RealTypeContext::is_exact> StorageType;
        StorageType d_storage;
    
        void doGS(unsigned start, unsigned stop)
        {
            typename StorageType::Real t(d_storage.d_rc), t2(d_storage.d_rc), norm(d_storage.d_rc);
            
            for (unsigned i = start; i <= stop; ++i)
            {
                // Start with norm (unprojected)
                convert(norm, normSq(d_A.row(i)), d_storage.d_rc);
                // Compute coefficients and projected norm
                for (unsigned j = 0; j < i; ++j)
                {
                    if (d_dependent[j])
                        d_storage.setCoeffZero(i, j);
                    else
                    {
                        // Compute coefficient
                        convert(t, dot(d_A.row(i), d_A.row(j)), d_storage.d_rc);
                        for (unsigned k = 0; k < j; ++k)
                        {
                            t2 = d_storage.getCoeff(j, k) * d_storage.getCoeff(i, k);
                            t -= t2 * d_storage.getSqNorm(k);
                        }
                        t /= d_storage.getSqNorm(j);
                        d_storage.setCoeff(i, j, t);
                        
                        // Continue with computation of squared norm
                        square(t, t);
                        t *= d_storage.getSqNorm(j);
                        norm -= t;
                    }
                }
                
                // Store norm
                d_storage.setSqNorm(i, norm);
                // Is dependent on the previous vectors?
                d_dependent[i] = isZero(norm);
            }
        }
        
        std::deque<std::pair<typename RealTypeContext::Real *, typename RealTypeContext::Real *> > d_precupdater_stack;
        
        mutable typename IntTypeContext::Integer d_tmp;
        
        Verbose * d_verbose;
        
        void do_add_row()
        {
            // Add row
            d_coeffs.resize(d_coeffs.rows() + 1, d_coeffs.cols() + 1, linalg::Initialize(d_rc));
            d_sqnorms.resize(d_sqnorms.size() + 1, linalg::Initialize(d_rc));
            d_dependent.resize(d_dependent.size() + 1);
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
            DuplicateStorage(const GSClassicInteger & gs, linalg::math_matrix<typename IntTypeContext::Integer> & m, RealTypeContext & rc, IntTypeContext & ic)
                : d_gs(gs, rc, ic, m)
            {
            }
            
            GSClassicInteger d_gs;
        };
        
        void copyFrom(const GSClassicInteger & gs)
        {
            d_level = gs.d_level;
            d_vectors_used = gs.d_vectors_used;
            d_coeffs = gs.d_coeffs;
            d_sqnorms = gs.d_sqnorms;
            d_dependent = gs.d_dependent;
            d_storage.copyFrom(gs.d_storage);
            d_precupdater_stack = gs.d_precupdater_stack;
        }
        
        void copyFrom(const DuplicateStorage & storage)
        {
            copyFrom(storage.d_gs);
        }
        
        GSClassicInteger(const GSClassicInteger & source, RealTypeContext & rc, IntTypeContext & ic, linalg::math_matrix<typename IntTypeContext::Integer> & A)
            : d_rc(rc), d_ic(ic), d_A(A), d_level(source.d_level), d_vectors_used(source.d_vectors_used),
              d_coeffs(source.d_coeffs), d_sqnorms(source.d_sqnorms), d_dependent(source.d_dependent),
              d_storage(*this, source.d_storage), d_precupdater_stack(source.d_precupdater_stack),
              d_verbose(source.d_verbose)
        {
        }
        
        GSClassicInteger(Verbose & verbose, RealTypeContext & rc, IntTypeContext & ic, linalg::math_matrix<typename IntTypeContext::Integer> & A)
            : d_rc(rc), d_ic(ic), d_A(A), d_level(0), d_vectors_used(d_A.rows()), d_storage(*this), d_verbose(&verbose)
        {
            d_coeffs.resize(A.rows(), A.rows(), linalg::Initialize(d_rc));
            d_sqnorms.resize(A.rows(), linalg::Initialize(d_rc));
            d_dependent.resize(A.rows());
            for (unsigned i = 0; i < A.rows(); ++i)
                d_dependent[i] = false;
        }
        
        linalg::math_matrix<typename RealTypeContext::Real> * getCoeffsStorage()
        {
            return &d_coeffs;
        }
        
        linalg::math_rowvector<typename RealTypeContext::Real> * getSqNormsStorage()
        {
            return &d_sqnorms;
        }
        
        void registerZDFAlpha(typename RealTypeContext::Real & zerodotfive, typename RealTypeContext::Real & ralpha)
        // Informs the Gram-Schmidt computer about the current values of zerodotfive and ralpha used. In
        // case the precision is adjusted, the RealTypeContext's of these values will be reset as well.
        {
            d_precupdater_stack.push_back(std::make_pair(&zerodotfive, &ralpha));
        }
        
        void unregisterZDFAlpha()
        // Removes the previously added values of zerodotfive and ralpha from the precision updating stack.
        {
            d_precupdater_stack.pop_back();
        }
        
        static void adjustZerodotfive(typename RealTypeContext::Real &, const typename RealTypeContext::Real &)
        // Adjust the value of "zerodotfive", which should be 0.5 for exact arithmetic, and slightly
        // larger for floating point approximations
        {
            // No adjustment needed.
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
            typename StorageType::Real t(d_storage.d_rc);
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
            if (d_dependent[j] || (d_dependent[i] && isZero(d_coeffs(i, j))))
            {
                // In this case, swapping is easy, since vector i (=j+1) does not uses j for reduction.
                d_storage.swapSqNorms(i, j);
                for (unsigned k = 0; k < j; ++k)
                    d_storage.swapCoeffs(j, k, i, k);
                for (unsigned k = i + 1; k < d_level; ++k)
                    d_storage.swapCoeffs(k, i, k, j);
                std::swap(d_dependent[i], d_dependent[j]);
            }
            else if (d_dependent[i])
            {
                assert(!isZero(d_coeffs(i, j))); // In this case, the j-th vector has to be independent
                                                 // and d_coeffs(i, j) not equal to zero
                /*
                  This case is a bit more tricky: the i-th vector depends on the j-th vector. After
                  swapping, this will still be the case.
                  
                  As b_i^* = b_i - \sum_{k=0}^j \mu_{ik} b_k^* = 0 with \mu_{ik} \neq 0, we have
                             b_i = \mu_{ij} b_j^* + \sum_{k=0}^{j-1} \mu_{ik} b_k^*
                        and  b_j =          b_j^* + \sum_{k=0}^{j-1} \mu_{jk} b_k^*.
                  
                  Therefore, if b_i' = b_j and b_j' = b_i,
                             b_j'   = \mu_{ij} b_j^* + \sum_{k=0}^{j-1} \mu_{ik} b_k^*,
                             b_j'^* = \mu_{ij} b_j^*,
                        and  b_i'   = 1/\mu_{ij} b_i'^* + \sum_{k=0}^{j-1} \mu_{jk} b_k^*,
                             b_i'^* = 0.
                  
                  For k > i,
                             b_k = b_k^* + \sum_{l=0}^{k-1} \mu_{kl} b_l^*
                                 = b_k^* + \sum_{l=0}^{j-1} \mu_{kl} b_l^* + \sum_{l=i+1}^{k-1} \mu_{kl} b_l^*
                                 + \mu_{kj} b_j^* + 0 * b_i^*
                                 = b_k^* + \sum_{l=0}^{j-1} \mu_{kl} b_l^* + \sum_{l=i+1}^{k-1} \mu_{kl} b_l^*
                                 + \mu_{kj}*\mu_{ij}/\mu_{ij} b_j^* + 0
                                 = b_k^* + \sum_{l=0}^{j-1} \mu_{kl} b_l^* + \sum_{l=i+1}^{k-1} \mu_{kl} b_l^* 
                                 + \mu_{kj}/\mu_{ij} b_j'^* + 0 * b_i'^*.
                */
                typename StorageType::Real t(d_storage.d_rc), tt(d_storage.d_rc);
                square(t, d_storage.getCoeff(i, j));
                d_storage.setSqNorm(j, d_storage.getSqNorm(j) * t);
                setOne(t);
                t /= d_storage.getCoeff(i, j);
                d_storage.setCoeff(i, j, t);
                for (unsigned k = 0; k < j; ++k)
                    d_storage.swapCoeffs(j, k, i, k);
                for (unsigned k = i + 1; k < d_level; ++k)
                {
                    tt = d_storage.getCoeff(k, j) * t;
                    d_storage.setCoeff(k, j, tt);
                }
            }
            else
            {
                // Modify the GS coefficients and the squared norms involved for i and j
                typename StorageType::Real t(d_storage.d_rc), tt(d_storage.d_rc);
                typename StorageType::Real muij(d_storage.d_rc), lambdai(d_storage.d_rc), lambdaj(d_storage.d_rc);
                muij = d_storage.getCoeff(i, j);
                lambdai = d_storage.getSqNorm(i);
                lambdaj = d_storage.getSqNorm(j);
                square(t, muij);
                t *= lambdaj;
                t += lambdai;
                d_storage.setSqNorm(j, t);
                t = muij * lambdaj;
                t /= d_storage.getSqNorm(j);
                d_storage.setCoeff(i, j, t);
                typename StorageType::Real A(d_storage.d_rc);
                A = lambdai / d_storage.getSqNorm(j);
                t = A * lambdaj;
                d_storage.setSqNorm(i, t);
                
                // "Swap" the GS vectors
                for (unsigned k = 0; k < j; ++k)
                    d_storage.swapCoeffs(j, k, i, k);
                
                // Update the other coefficients
                if (i + 1 < d_level)
                {
                    for (unsigned k = i + 1; k < d_level; ++k)
                    {
                        t = d_storage.getCoeff(k, j);
                        tt = A * d_storage.getCoeff(k, i);
                        tt += t * d_storage.getCoeff(i, j);
                        d_storage.setCoeff(k, j, tt);
                        tt = muij * d_storage.getCoeff(k, i);
                        t -= tt;
                        d_storage.setCoeff(k, i, t);
                    }
                }
            }
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
            if (isPMOne(m))
            {
                if (sign(m) > 0)
                {
                    for (unsigned k = 0; k < j; ++k)
                        if (!d_dependent[k])
                            d_storage.setCoeff(i, k, d_storage.getCoeff(i, k) + d_storage.getCoeff(j, k));
                    d_storage.incCoeff(i, j);
                }
                else
                {
                    for (unsigned k = 0; k < j; ++k)
                        if (!d_dependent[k])
                            d_storage.setCoeff(i, k, d_storage.getCoeff(i, k) - d_storage.getCoeff(j, k));
                    d_storage.decCoeff(i, j);
                }
            }
            else
            {
                for (unsigned k = 0; k < j; ++k)
                    if (!d_dependent[k])
                        d_storage.setCoeff(i, k, d_storage.getCoeff(i, k) + convert(m, d_storage.d_rc) * d_storage.getCoeff(j, k));
                d_storage.setCoeff(i, j, d_storage.getCoeff(i, j) + convert(m, d_storage.d_rc));
            }
        }
        
        void flip(unsigned i)
        {
            if (d_level <= i)
                return;
            for (unsigned k = 0; k < i; ++k)
                d_storage.setCoeff(i, k, -d_storage.getCoeff(i, k));
            for (unsigned k = i + 1; k < d_level; ++k)
                d_storage.setCoeff(k, i, -d_storage.getCoeff(k, i));
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
        
            typename StorageType::Real A(d_storage.d_rc), B(d_storage.d_rc), C(d_storage.d_rc), D(d_storage.d_rc);
            convert(A, B00, d_storage.d_rc); // A = B00
            convert(B, B01, d_storage.d_rc); // B = B01
            convert(C, B10, d_storage.d_rc); // C = B10
            convert(D, B11, d_storage.d_rc); // D = B11
            typename StorageType::Real E(d_storage.d_rc);
            E = B;
            E *= d_storage.getCoeff(ell, k);
            E += A;
            // E = (A + B \mu_{\ell k})
            typename StorageType::Real F(d_storage.d_rc), G(d_storage.d_rc);
            F = E * E;
            F *= d_storage.getSqNorm(k);
            G = B * B;
            G *= d_storage.getSqNorm(ell);
            F += G;
            // \lambda_k' = F
            G = D * d_storage.getCoeff(ell, k);
            G += C;
            // G = C + D \mu_{\ell k}
            typename StorageType::Real H(d_storage.d_rc);
            H = G * E;
            H *= d_storage.getSqNorm(k);
            H += B * D * d_storage.getSqNorm(ell);
            H /= F;
            // \mu_{\ell k}' = H
            typename StorageType::Real I(d_storage.d_rc), J(d_storage.d_rc);
            I = H * E;
            I = G - I;
            // I = G - H E
            J = H * B;
            J = D - J;
            // J = D - H B
            typename StorageType::Real K(d_storage.d_rc), L(d_storage.d_rc);
            K = I * I;
            K *= d_storage.getSqNorm(k);
            L = J * J;
            L *= d_storage.getSqNorm(ell);
            K += L;
            // \lambda_\ell' = K
            typename StorageType::Real M(d_storage.d_rc);
            L = E * d_storage.getSqNorm(k);
            L /= F;
            M = B * d_storage.getSqNorm(ell);
            M /= F;
            typename StorageType::Real N(d_storage.d_rc), O(d_storage.d_rc);
            N = I * d_storage.getSqNorm(k);
            N /= K;
            O = J * d_storage.getSqNorm(ell);
            O /= K;

            typename StorageType::Real T(d_storage.d_rc);
            d_storage.setSqNorm(k, F);
            d_storage.setCoeff(ell, k, H);
            d_storage.setSqNorm(ell, K);
            for (unsigned j = 0; j < k; ++j)
            {
                T = d_storage.getCoeff(k, j);
                d_storage.setCoeff(k, j, T * A);
                d_storage.setCoeff(k, j, d_storage.getCoeff(k, j) + B * d_storage.getCoeff(ell, j));
                d_storage.setCoeff(ell, j, d_storage.getCoeff(ell, j) * D);
                d_storage.setCoeff(ell, j, d_storage.getCoeff(ell, j) + C * T);
            }
            for (unsigned i = ell + 1; i < d_level; ++i)
            {
                T = d_storage.getCoeff(i, k);
                d_storage.setCoeff(i, k, T * L);
                d_storage.setCoeff(i, k, d_storage.getCoeff(i, k) + M * d_storage.getCoeff(i, ell));
                d_storage.setCoeff(i, ell, d_storage.getCoeff(i, ell) * O);
                d_storage.setCoeff(i, ell, d_storage.getCoeff(i, ell) + N * T);
            }
        }
        
        void print()
        {
            std::cout << "Computed: " << d_level << "\n";
            for (unsigned i = 0; i < d_level; ++i)
            {
                std::cout << i << ":";
                for (unsigned j = 0; j < i; ++j)
                {
                    std::cout << " " << d_coeffs(i, j);
                }
                std::cout << " -> " << d_sqnorms[i] << " [" << d_storage.getSqNorm(i) << "]\n";
            }
            std::cout << d_A << "\n";
        }
        
        void update(unsigned level)
        {
            assert(level <= d_vectors_used);
            if (level <= d_level)
                return;
            
            // Apply exact Gram-Schmidt
            doGS(d_level, level - 1);
            
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
            arithmetic::Integer lambda;
            typename IntTypeContext::Integer lambdap;
            bool changed = false;
            // Do size reduction with end down to begin
            for (signed i = end; i >= (signed)begin; --i)
            {
                if (arithmetic::abs(d_coeffs(dest, i)) > zerodotfive)
                {
                    d_storage.getRoundedCoeff(lambda, dest, i);
                    neg(lambda, lambda);
                    convert(lambdap, lambda, d_ic);
                    lattice.add(dest, lambdap, i);
                    changed = true;
                }
            }
            return changed;
        }
        
        void computeDotProductProjected(typename RealTypeContext::Real & r, unsigned k, unsigned i, unsigned j) const
        {
            // FIX ME: USE INTEGER ARITHMETIC !!! ??? ...
            GSGeneric::computeDotProductProjected(d_rc, d_ic, d_coeffs, d_sqnorms, d_dependent, d_level, r, k, i, j);
        }
    
        template<class RTC, class ITC, bool exact> class computeProjections;
    
        template<class RTC, class ITC> class computeProjections<RTC, ITC, true>
        {
        public:
            static void computeProjectionLength(const GSClassicInteger & gs, const Storage<RTC, ITC, true> & storage, RTC & rc, ITC & ic,
                                                typename RTC::Real & result, unsigned k, unsigned b,
                                                const linalg::math_rowvector<typename ITC::Integer> & vec)
            // Computes the squared length of the projection of \sum_{i=0}^{m-1} vec[i] B_{b+i} onto the
            // orthogonal complement of B_0, ..., B_{k-1}. (Here, B_0, ..., B_{n-1} is the basis of the
            // lattice.)
            {
//    std::cout << k << " " << b << " " << &vec << " " << std::flush;
//    std::cout << vec.size() << " " << vec << "\n";
                typename RTC::Real t(rc);
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
                        t += convert(vec[j - b], rc) * gs.d_coeffs(j, i);
                    square(t, t);
                    result += t * gs.d_sqnorms[i];
                }
                for (unsigned i = beg; i <= end; ++i)
                {
                    convert(t, vec[i - b], rc);
                    for (unsigned j = i + 1; j <= end; ++j)
                        t += convert(vec[j - b], rc) * gs.d_coeffs(j, i);
                    square(t, t);
                    result += t * gs.d_sqnorms[i];
                }
            }
        
            static void computeProjectionLengthBV(const GSClassicInteger & gs, const Storage<RTC, ITC, true> & storage, RTC & rc, ITC & ic,
                                                  typename RTC::Real & result, unsigned k, unsigned b)
            // Computes the squared length of the projection of B_b onto the orthogonal complement of B_0, ...,
            // B_{k-1}. (Here, B_0, ..., B_{n-1} is the basis of the lattice.)
            {
                typename RTC::Real t(rc);
                if (b < k)
                {
                    // Projection is 0
                    setZero(result);
                    return;
                }
                result = gs.d_sqnorms[b];
                for (unsigned i = k; i < b; ++i)
                {
                    square(t, gs.d_coeffs(b, i));
                    result += t * gs.d_sqnorms[i];
                }
            }
        
            static void computeProjectionLength(const GSClassicInteger & gs, const Storage<RTC, ITC, true> & storage, RTC & rc, ITC & ic,
                                                typename RTC::Real & result, unsigned k, const linalg::math_rowvector<typename ITC::Integer> & vec)
            // Computes the squared length of the projection of vec onto the orthogonal complement of B_0,
            // ..., B_{k-1}. (Here, B_0, ..., B_{n-1} is the basis of the lattice.)
            {
                typename RTC::Real t(rc);
                // Compute squared norm (unprojected)
                convert(result, normSq(vec), rc);
                if (k == 0)
                    // In case we just need the norm, return
                    return;
                // Orthogonalize vector
                linalg::math_rowvector<typename RTC::Real> r(k);
                for (unsigned j = 0; j < k; ++j)
                    if (!gs.d_dependent[j])
                    {
                        r[j].setContext(rc); // r[j] = coeffs(new_vector, j) * d_sqnorms[j]
                        // Compute dot product <B_j, vec>
                        convert(r[j], dot(gs.d_A.row(j), vec), rc);
                        for (unsigned k = 0; k < j; ++k)
                        {
                            t = r[k] * gs.d_coeffs(j, k);
                            r[j] -= t;
                        }
                        // Compute length of projection iteratively
                        t = r[j] / gs.d_sqnorms[j];
                        t *= r[j];
                        result = result - t;
                    }
            }
        };
    
        template<class RTC, class ITC> class computeProjections<RTC, ITC, false>
        {
        public:
            static void computeProjectionLength(const GSClassicInteger & gs, const Storage<RTC, ITC, false> & storage, RTC & rc, ITC & ic,
                                                typename RTC::Real & result, unsigned k, unsigned b,
                                                const linalg::math_rowvector<typename ITC::Integer> & vec)
            // Computes the squared length of the projection of \sum_{i=0}^{m-1} vec[i] B_{b+i} onto the
            // orthogonal complement of B_0, ..., B_{k-1}. (Here, B_0, ..., B_{n-1} is the basis of the
            // lattice.)
            {
//    std::cout << k << " " << b << " " << &vec << " " << std::flush;
//    std::cout << vec.size() << " " << vec << "\n";
                typename StorageType::Real t(storage.d_rc), res(storage.d_rc);
                unsigned beg = k < b ? b : k;
                unsigned end = b + vec.size() - 1;
                if ((end < k) || (vec.size() == 0))
                {
                    // Projection is 0
                    setZero(result);
                    return;
                }
                setZero(res);
                for (unsigned i = k; i < beg; ++i)
                {
                    setZero(t);
                    unsigned begg = i < beg ? beg : i + 1;
                    for (unsigned j = begg; j <= end; ++j)
                        t += convert(vec[j - b], storage.d_rc) * storage.getCoeff(j, i);
                    square(t, t);
                    res += t * storage.getSqNorm(i);
                }
                for (unsigned i = beg; i <= end; ++i)
                {
                    convert(t, vec[i - b], storage.d_rc);
                    for (unsigned j = i + 1; j <= end; ++j)
                        t += convert(vec[j - b], storage.d_rc) * storage.getCoeff(j, i);
                    square(t, t);
                    res += t * storage.getSqNorm(i);
                }
                convert(result, res, rc);
            }
        
            static void computeProjectionLengthBV(const GSClassicInteger & gs, const Storage<RTC, ITC, false> & storage, RTC & rc, ITC & ic,
                                                  typename RTC::Real & result, unsigned k, unsigned b)
            // Computes the squared length of the projection of B_b onto the orthogonal complement of B_0, ...,
            // B_{k-1}. (Here, B_0, ..., B_{n-1} is the basis of the lattice.)
            {
                if (b < k)
                {
                    // Projection is 0
                    setZero(result);
                    return;
                }
                typename StorageType::Real t(storage.d_rc), res(storage.d_rc);
                res = storage.getSqNorm(b);
                for (unsigned i = k; i < b; ++i)
                {
                    square(t, storage.getCoeff(b, i));
                    res += t * storage.getSqNorm(i);
                }
                convert(result, res, rc);
            }
        
            static void computeProjectionLength(const GSClassicInteger & gs, const Storage<RTC, ITC, false> & storage, RTC & rc, ITC & ic,
                                                typename RTC::Real & result, unsigned k, const linalg::math_rowvector<typename ITC::Integer> & vec)
            // Computes the squared length of the projection of vec onto the orthogonal complement of B_0,
            // ..., B_{k-1}. (Here, B_0, ..., B_{n-1} is the basis of the lattice.)
            {
                // Compute squared norm (unprojected)
                if (k == 0)
                {
                    convert(result, normSq(vec), rc);
                    // In case we just need the norm, return
                    return;
                }
                typename StorageType::Real t(storage.d_rc), res(storage.d_rc);
                convert(res, normSq(vec), storage.d_rc);
                // Orthogonalize vector
                linalg::math_rowvector<typename StorageType::Real> r(k);
                for (unsigned j = 0; j < k; ++j)
                    if (!gs.d_dependent[j])
                    {
                        r[j].setContext(storage.d_rc); // r[j] = coeffs(new_vector, j) * d_sqnorms[j]
                        // Compute dot product <B_j, vec>
                        convert(r[j], dot(gs.d_A.row(j), vec), storage.d_rc);
                        for (unsigned k = 0; k < j; ++k)
                        {
                            t = r[k] * storage.getCoeff(j, k);
                            r[j] -= t;
                        }
                        // Compute length of projection iteratively
                        t = r[j] / storage.getSqNorm(j);
                        t *= r[j];
                        res = res - t;
                    }
                convert(result, res, rc);
            }
        };
    
        void computeProjectionLength(typename RealTypeContext::Real & result, unsigned k, unsigned b,
                                     const linalg::math_rowvector<typename IntTypeContext::Integer> & vec) const
        // Computes the squared length of the projection of \sum_{i=0}^{m-1} vec[i] B_{b+i} onto the
        // orthogonal complement of B_0, ..., B_{k-1}. (Here, B_0, ..., B_{n-1} is the basis of the
        // lattice.)
        {
            computeProjections<RealTypeContext, IntTypeContext, RealTypeContext::is_exact>::
                computeProjectionLength(*this, d_storage, d_rc, d_ic, result, k, b, vec);
        }
    
        void computeProjectionLengthBV(typename RealTypeContext::Real & result, unsigned k, unsigned b) const
        // Computes the squared length of the projection of B_b onto the orthogonal complement of B_0, ...,
        // B_{k-1}. (Here, B_0, ..., B_{n-1} is the basis of the lattice.)
        {
            computeProjections<RealTypeContext, IntTypeContext, RealTypeContext::is_exact>::
                computeProjectionLengthBV(*this, d_storage, d_rc, d_ic, result, k, b);
        }
    
        void computeProjectionLength(typename RealTypeContext::Real & result, unsigned k,
                                     const linalg::math_rowvector<typename IntTypeContext::Integer> & vec) const
        // Computes the squared length of the projection of vec onto the orthogonal complement of B_0,
        // ..., B_{k-1}. (Here, B_0, ..., B_{n-1} is the basis of the lattice.)
        {
            computeProjections<RealTypeContext, IntTypeContext, RealTypeContext::is_exact>::
                computeProjectionLength(*this, d_storage, d_rc, d_ic, result, k, vec);
        }
    
        void insertVectorLC(unsigned ofs, const linalg::math_rowvector<typename IntTypeContext::Integer> & result)
        // Inserts a new vector at index ofs, which will equal \sum_{i=0}^{result.size()-1} result[i] * b[ofs+i]
        {
            bool add_row = d_vectors_used == d_coeffs.rows();
            if (d_level > ofs)
                d_level = ofs;
            if (add_row)
                do_add_row();
            ++d_vectors_used;
        }
    
        void insertVector(unsigned ofs, const linalg::math_rowvector<typename IntTypeContext::Integer> & result)
        // Inserts a new vector <result> at index ofs
        {
            bool add_row = d_vectors_used == d_coeffs.rows();
            if (d_level > ofs)
                d_level = ofs;
            if (add_row)
                do_add_row();
            ++d_vectors_used;
        }
    
        void removeZeroVector(unsigned ofs)
        // Removes the zero vector at position ofs
        {
            if (ofs < d_level)
            {
                // Remove row
                for (unsigned i = ofs + 1; i < d_level; ++i)
                    d_storage.swapCoeffsRow(i, i - 1);
                // Remove column
                for (unsigned i = ofs + 1; i < d_level; ++i)
                    d_storage.swapCoeffsCol(i, i - 1);
                // Remove norms
                for (unsigned i = ofs + 1; i < d_level; ++i)
                    d_storage.swapSqNorms(i, i - 1);
                // Remove dependence info
                for (unsigned i = ofs + 1; i < d_level; ++i)
                    d_dependent[i - 1] = d_dependent[i];
                // Decrease level
                --d_level;
            }
            --d_vectors_used;
        }

        static void compactify()
        {
            /*
              if (d_vectors_used < d_coeffs.rows())
              {
              d_coeffs.resize(d_vectors_used, d_vectors_used, linalg::Initialize(d_rc));
              d_sqnorms.resize(d_vectors_used, linalg::Initialize(d_rc));
              d_dependent.resize(d_vectors_used, linalg::Initialize(d_rc));
              d_storage.resize(d_vectors_used, linalg::Initialize(d_rc));
              }
            */
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
            d_storage.reassign(d_level);
            for (typename std::deque<std::pair<typename RealTypeContext::Real *, typename RealTypeContext::Real *> >::iterator i = d_precupdater_stack.begin(); i != d_precupdater_stack.end(); ++i)
            {
                i->first->setContext(d_rc);
                i->second->setContext(d_rc);
            }
        }
    };
}

#endif

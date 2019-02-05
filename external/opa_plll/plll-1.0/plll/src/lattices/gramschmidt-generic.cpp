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

#ifndef PLLL_INCLUDE_GUARD__GRAMSCHMIDT_GENERIC_CPP
#define PLLL_INCLUDE_GUARD__GRAMSCHMIDT_GENERIC_CPP

namespace plll
{
    namespace GSGeneric
    {
        /*
          Contains generic Gram-Schmid decomposition updating methods.
        */
    
        template<class RealTypeContext, class IntTypeContext>
        void swap(RealTypeContext & d_rc, IntTypeContext & d_ic,
                  linalg::math_matrix<typename RealTypeContext::Real> & d_coeffs,
                  linalg::math_rowvector<typename RealTypeContext::Real> & d_sqnorms,
                  linalg::base_rowvector<bool> & d_dependent,
                  unsigned d_number_of_computed_values,
                  unsigned i, unsigned j)
        // Assumes i - j = 1.
        {
            assert(j + 1 == i);
        
            if (d_dependent[j] || (d_dependent[i] && isZero(d_coeffs(i, j))))
            {
                // In this case, swapping is easy, since vector i (=j+1) does not uses j for reduction.
                arithmetic::swap(d_sqnorms[i], d_sqnorms[j]);
                for (unsigned k = 0; k < j; ++k)
                    arithmetic::swap(d_coeffs(j, k), d_coeffs(i, k));
                for (unsigned k = i + 1; k < d_number_of_computed_values; ++k)
                    arithmetic::swap(d_coeffs(k, i), d_coeffs(k, j));
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
                typename RealTypeContext::Real t(d_rc);
                square(t, d_coeffs(i, j));
                d_sqnorms[j] *= t;
                setOne(t);
                d_coeffs(i, j) = t / d_coeffs(i, j);
                for (unsigned k = 0; k < j; ++k)
                    arithmetic::swap(d_coeffs(j, k), d_coeffs(i, k));
                for (unsigned k = i + 1; k < d_number_of_computed_values; ++k)
                    d_coeffs(k, j) *= d_coeffs(i, j);
            }
            else
            {
                typename RealTypeContext::Real t(d_rc);
                // Modify the GS coefficients and the squared norms involved for i and j
                typename RealTypeContext::Real muij(d_rc), lambdai(d_rc), lambdaj(d_rc), A(d_rc);
                muij = d_coeffs(i, j);
                lambdai = d_sqnorms[i];
                lambdaj = d_sqnorms[j];
                d_sqnorms[j] = muij * muij;
                d_sqnorms[j] *= lambdaj;
                d_sqnorms[j] += lambdai;
                d_coeffs(i, j) = muij * lambdaj;
                d_coeffs(i, j) /= d_sqnorms[j];
                A = lambdai / d_sqnorms[j];
                d_sqnorms[i] = A * lambdaj;
            
                // "Swap" the GS vectors
                for (unsigned k = 0; k < j; ++k)
                    arithmetic::swap(d_coeffs(j, k), d_coeffs(i, k));
            
                // Update the other coefficients
                for (unsigned k = i + 1; k < d_number_of_computed_values; ++k)
                {
                    t = d_coeffs(k, j);
                    d_coeffs(k, j) *= d_coeffs(i, j);
                    d_coeffs(k, j) += A * d_coeffs(k, i);
                    d_coeffs(k, i) *= muij;
                    d_coeffs(k, i) = -d_coeffs(k, i);
                    d_coeffs(k, i) += t;
                }
            }
        }
    
        template<class RealTypeContext, class IntTypeContext>
        void add(RealTypeContext & d_rc, IntTypeContext & d_ic,
                 linalg::math_matrix<typename RealTypeContext::Real> & d_coeffs,
                 linalg::math_rowvector<typename RealTypeContext::Real> & d_sqnorms,
                 linalg::base_rowvector<bool> & d_dependent,
                 unsigned d_number_of_computed_values,
                 unsigned i, const typename IntTypeContext::Integer & m, unsigned j)
        // Assumes i > j.
        {
            assert(i > j);
            if (isPMOne(m))
            {
                if (sign(m) > 0)
                {
                    for (unsigned k = 0; k < j; ++k)
                        if (!d_dependent[k])
                            d_coeffs(i, k) += d_coeffs(j, k);
                    ++d_coeffs(i, j);
                }
                else
                {
                    for (unsigned k = 0; k < j; ++k)
                        if (!d_dependent[k])
                            d_coeffs(i, k) -= d_coeffs(j, k);
                    --d_coeffs(i, j);
                }
            }
            else
            {
                for (unsigned k = 0; k < j; ++k)
                    if (!d_dependent[k])
                        d_coeffs(i, k) += convert(m, d_rc) * d_coeffs(j, k);
                d_coeffs(i, j) += convert(m, d_rc);
            }
        }
    
        template<class RealTypeContext, class IntTypeContext>
        void flip(RealTypeContext & d_rc, IntTypeContext & d_ic,
                  linalg::math_matrix<typename RealTypeContext::Real> & d_coeffs,
                  linalg::math_rowvector<typename RealTypeContext::Real> & d_sqnorms,
                  linalg::base_rowvector<bool> & d_dependent,
                  unsigned d_number_of_computed_values,
                  unsigned i)
        {
            for (unsigned k = 0; k < i; ++k)
                neg(d_coeffs(i, k), d_coeffs(i, k));
            for (unsigned k = i + 1; k < d_number_of_computed_values; ++k)
                neg(d_coeffs(k, i), d_coeffs(k, i));
        }
        
        template<class RealTypeContext, class IntTypeContext>
        void trans(RealTypeContext & d_rc, IntTypeContext & d_ic,
                   linalg::math_matrix<typename RealTypeContext::Real> & d_coeffs,
                   linalg::math_rowvector<typename RealTypeContext::Real> & d_sqnorms,
                   linalg::base_rowvector<bool> & d_dependent,
                   unsigned d_number_of_computed_values,
                   unsigned k, unsigned ell,
                   const typename IntTypeContext::Integer & B00, const typename IntTypeContext::Integer & B01,
                   const typename IntTypeContext::Integer & B10, const typename IntTypeContext::Integer & B11)
        // Assumes k - ell = 1, and that both vectors are not dependent.
        {
            assert(k + 1 == ell);
            assert(!d_dependent[k] && !d_dependent[ell]);
            
            typename RealTypeContext::Real A(d_rc), B(d_rc), C(d_rc), D(d_rc);
            convert(A, B00, d_rc); // A = B00
            convert(B, B01, d_rc); // B = B01
            convert(C, B10, d_rc); // C = B10
            convert(D, B11, d_rc); // D = B11
            typename RealTypeContext::Real E(d_rc);
            E = B;
            E *= d_coeffs(ell, k);
            E += A;
            // E = (A + B \mu_{\ell k})
            typename RealTypeContext::Real F(d_rc), G(d_rc);
            F = E * E;
            F *= d_sqnorms[k];
            G = B * B;
            G *= d_sqnorms[ell];
            F += G;
            // \lambda_k' = F
            G = D * d_coeffs(ell, k);
            G += C;
            // G = C + D \mu_{\ell k}
            typename RealTypeContext::Real H(d_rc);
            H = G * E;
            H *= d_sqnorms[k];
            H += B * D * d_sqnorms[ell];
            H /= F;
            // \mu_{\ell k}' = H
            typename RealTypeContext::Real I(d_rc), J(d_rc);
            I = H * E;
            I = G - I;
            // I = G - H E
            J = H * B;
            J = D - J;
            // J = D - H B
            typename RealTypeContext::Real K(d_rc), L(d_rc);
            K = I * I;
            K *= d_sqnorms[k];
            L = J * J;
            L *= d_sqnorms[ell];
            K += L;
            // \lambda_\ell' = K
            typename RealTypeContext::Real M(d_rc);
            L = E * d_sqnorms[k];
            L /= F;
            M = B * d_sqnorms[ell];
            M /= F;
            typename RealTypeContext::Real N(d_rc), O(d_rc);
            N = I * d_sqnorms[k];
            N /= K;
            O = J * d_sqnorms[ell];
            O /= K;
            
            typename RealTypeContext::Real T(d_rc);
            d_sqnorms[k] = F;
            d_coeffs(ell, k) = H;
            d_sqnorms[ell] = K;
            for (unsigned j = 0; j < k; ++j)
            {
                T = d_coeffs(k, j);
                d_coeffs(k, j) *= A;
                d_coeffs(k, j) += B * d_coeffs(ell, j);
                d_coeffs(ell, j) *= D;
                d_coeffs(ell, j) += C * T;
            }
            for (unsigned i = ell + 1; i < d_number_of_computed_values; ++i)
            {
                T = d_coeffs(i, k);
                d_coeffs(i, k) *= L;
                d_coeffs(i, k) += M * d_coeffs(i, ell);
                d_coeffs(i, ell) *= O;
                d_coeffs(i, ell) += N * T;
            }
        }
        
        template<class RealTypeContext, class IntTypeContext, class Lattice>
        bool sizereduce(RealTypeContext & d_rc, IntTypeContext & d_ic,
                        linalg::math_matrix<typename RealTypeContext::Real> & d_coeffs,
                        linalg::math_rowvector<typename RealTypeContext::Real> & d_sqnorms,
                        linalg::base_rowvector<bool> & d_dependent,
                        unsigned d_number_of_computed_values,
                        Lattice & lattice, typename RealTypeContext::Real & zerodotfive,
                        unsigned begin, unsigned end, unsigned dest)
        // zerodotfive: something in [0.5, 1)
        //
        // Applies size reduction: A.row(begin) to A.row(end) will be used to size-reduce A.row(dest).
        // Assumes begin <= end < dest.
        {
            typename IntTypeContext::Integer lambda(d_ic);
            bool changed = false;
            // Find maximal exponent
            long exp = d_coeffs(dest, end).getApproxExponent();
            for (signed i = end - 1; i >= (signed)begin; --i)
            {
                long e = d_coeffs(dest, i).getApproxExponent();
                if (exp < e)
                    exp = e;
            }
            
            bool keep_repeating;
            long badloops = 0;
            // Loop until done (or until no decrease for too long time)
            do
            {
                keep_repeating = false;
                lattice.update(dest + 1);
                // Do size reduction with end down to begin
                long newexp = d_coeffs(dest, end).getApproxExponent();
                for (signed i = end; i >= (signed)begin; --i)
                {
                    long newe = d_coeffs(dest, i).getApproxExponent();
                    if (newexp < newe)
                        newexp = newe;
//        std::cout << dest << " " << i << " " << arithmetic::convert<arithmetic::Real>(d_coeffs(dest, i)) << "\n";
                    if (arithmetic::abs(d_coeffs(dest, i)) > zerodotfive)
                    {
                        convert_round(lambda, d_coeffs(dest, i), d_ic);
                        neg(lambda, lambda);
                        lattice.add(dest, lambda, i);
                        keep_repeating = true;
//            std::cout << "    " << lambda << " " << arithmetic::convert<arithmetic::Real>(d_coeffs(dest, i)) << "\n";
                    }
                }
                if (keep_repeating)
                    changed = true;
                // Watch out for infinite loops
                if (newexp > exp - 5)
                    ++badloops;
                else
                    badloops = 0; // reset counter
                exp = newexp;
                // Terminate
                if (badloops == 64)
                    throw reduction_error("Too long size reduction loop detected. Precision probably not sufficient.");
            }
            while (keep_repeating);
            return changed;
        }
        
        template<class RealTypeContext, class IntTypeContext>
        void computeDotProductProjected(RealTypeContext & d_rc, IntTypeContext & d_ic,
                                        const linalg::math_matrix<typename RealTypeContext::Real> & d_coeffs,
                                        const linalg::math_rowvector<typename RealTypeContext::Real> & d_sqnorms,
                                        const linalg::base_rowvector<bool> & d_dependent,
                                        unsigned d_number_of_computed_values,
                                        typename RealTypeContext::Real & result,
                                        unsigned k, unsigned i, unsigned j)
        {
            if ((i < k) || (j < k))
            {
                // Projection is 0
                setZero(result);
                return;
            }
        
            // Ensure that j <= i
            if (i < j)
                std::swap(i, j);
        
            // Use that   \pi_k(b_j) = b_j^* + \sum_{\ell=k}^j \mu_{j,\ell} b_\ell^*
            // and        \pi_k(b_i) = b_i^* + \sum_{\ell=k}^i \mu_{i,\ell} b_\ell^*.
            // Therefore, the dot product of \pi_k(b_i) and \pi_k(b_j) equals
            //            \mu_{ij} \langle b_j^*, b_j^* \rangle + \sum_{\ell=k}^{j-1} \mu_{i,\ell} \mu_{j,\ell} \langle b_\ell^*, b_\ell^* \rangle
            // in case of j < i, and
            //            \langle b_i^*, b_i^* \rangle + \sum_{\ell=k}^{i-1} \mu_{i,\ell}^2 \langle b_\ell^*, b_\ell^* \rangle
            // in case of j = i.
            typename RealTypeContext::Real t(d_rc);
            if (i == j)
            {
                result = d_sqnorms[j];
                for (unsigned ell = k; ell < j; ++ell)
                {
                    square(t, d_coeffs(i, ell));
                    result += t * d_sqnorms[ell];
                }
            }
            else
            {
                result = d_coeffs(i, j) * d_sqnorms[j];
                for (unsigned ell = k; ell < j; ++ell)
                {
                    t = d_coeffs(i, ell) * d_coeffs(j, ell);
                    result += t * d_sqnorms[ell];
                }
            }
        }
    
        template<class RealTypeContext, class IntTypeContext>
        void computeProjectionLength(RealTypeContext & d_rc, IntTypeContext & d_ic,
                                     const linalg::math_matrix<typename RealTypeContext::Real> & d_coeffs,
                                     const linalg::math_rowvector<typename RealTypeContext::Real> & d_sqnorms,
                                     const linalg::base_rowvector<bool> & d_dependent,
                                     unsigned d_number_of_computed_values,
                                     typename RealTypeContext::Real & result,
                                     unsigned k, unsigned b, const linalg::math_rowvector<typename IntTypeContext::Integer> & vec)
        // Computes the squared length of the projection of \sum_{i=0}^{m-1} vec[i] B_{b+i} onto the
        // orthogonal complement of B_0, ..., B_{k-1}. (Here, B_0, ..., B_{n-1} is the basis of the
        // lattice.)
        {
//    std::cout << k << " " << b << " " << &vec << " " << std::flush;
//    std::cout << vec.size() << " " << vec << "\n";
            setZero(result);
            unsigned beg = k < b ? b : k;
            unsigned end = b + vec.size() - 1;
            if ((k > end) || (vec.size() == 0))
                // Projection is 0
                return;
            typename RealTypeContext::Real t(d_rc);
            for (unsigned i = k; i < beg; ++i)
            {
                setZero(t);
                unsigned begg = i < beg ? beg : i + 1;
                for (unsigned j = begg; j <= end; ++j)
                    t += convert(vec[j - b], d_rc) * d_coeffs(j, i);
                square(t, t);
                result += t * d_sqnorms[i];
            }
            for (unsigned i = beg; i <= end; ++i)
            {
                convert(t, vec[i - b], d_rc);
                for (unsigned j = i + 1; j <= end; ++j)
                    t += convert(vec[j - b], d_rc) * d_coeffs(j, i);
                square(t, t);
                result += t * d_sqnorms[i];
            }
        }
    
        template<class RealTypeContext, class IntTypeContext>
        void computeProjectionLengthBV(RealTypeContext & d_rc, IntTypeContext & d_ic,
                                       const linalg::math_matrix<typename RealTypeContext::Real> & d_coeffs,
                                       const linalg::math_rowvector<typename RealTypeContext::Real> & d_sqnorms,
                                       const linalg::base_rowvector<bool> & d_dependent,
                                       unsigned d_number_of_computed_values,
                                       typename RealTypeContext::Real & result, unsigned k, unsigned b)
        // Computes the squared length of the projection of B_b onto the orthogonal complement of B_0, ...,
        // B_{k-1}. (Here, B_0, ..., B_{n-1} is the basis of the lattice.)
        {
            if (b < k)
            {
                // Projection is 0
                setZero(result);
                return;
            }
            result = d_sqnorms[b];
            if (b > k)
            {
                typename RealTypeContext::Real t(d_rc);
                for (unsigned i = k; i < b; ++i)
                {
                    square(t, d_coeffs(b, i));
                    result += t * d_sqnorms[i];
                }
            }
        }
    
        template<class RealTypeContext, class IntTypeContext>
        void computeProjectionLength(RealTypeContext & d_rc, IntTypeContext & d_ic,
                                     const linalg::math_matrix<typename RealTypeContext::Real> & d_coeffs,
                                     const linalg::math_rowvector<typename RealTypeContext::Real> & d_sqnorms,
                                     const linalg::base_rowvector<bool> & d_dependent,
                                     const linalg::math_matrix<typename IntTypeContext::Integer> & d_A,
                                     unsigned d_number_of_computed_values,
                                     typename RealTypeContext::Real & result, unsigned k,
                                     const linalg::math_rowvector<typename IntTypeContext::Integer> & vec)
        // Computes the squared length of the projection of vec onto the orthogonal complement of B_0,
        // ..., B_{k-1}. (Here, B_0, ..., B_{n-1} is the basis of the lattice.)
        {
            typename RealTypeContext::Real t(d_rc);
            // Compute squared norm (unprojected)
            convert(result, normSq(vec), d_rc);
            if (k == 0)
                // In case we just need the norm, return
                return;
            // Orthogonalize vector
            linalg::math_rowvector<typename RealTypeContext::Real> r(k);
            for (unsigned j = 0; j < k; ++j)
                if (!d_dependent[j])
                {
                    r[j].setContext(d_rc); // r[j] = coeffs(new_vector, j) * d_sqnorms[j]
                    // Compute dot product <B_j, vec>
                    convert(r[j], dot(d_A.row(j), vec), d_rc);
                    for (unsigned k = 0; k < j; ++k)
                    {
                        t = r[k] * d_coeffs(j, k);
                        r[j] -= t;
                    }
                    // Compute length of projection iteratively
                    t = r[j] / d_sqnorms[j];
                    t *= r[j];
                    result = result - t;
                }
        }
    
        template<class RealTypeContext, class IntTypeContext>
        bool insertVector(RealTypeContext & d_rc, IntTypeContext & d_ic,
                          linalg::math_matrix<typename RealTypeContext::Real> & d_coeffs,
                          linalg::math_rowvector<typename RealTypeContext::Real> & d_sqnorms,
                          linalg::base_rowvector<bool> & d_dependent,
                          unsigned & d_vectors_used,
                          unsigned ofs,
                          bool shift)
        // Makes space for a new vector at index ofs. Returns true in case a row/column was added.
        {
            bool add_row = d_vectors_used == d_coeffs.rows();
            if (add_row)
            {
                d_coeffs.resize(d_coeffs.rows() + 1, d_coeffs.cols() + 1, linalg::Initialize(d_rc));
                d_sqnorms.resize(d_sqnorms.size() + 1, linalg::Initialize(d_rc));
                d_dependent.resize(d_dependent.size() + 1);
            }
            ++d_vectors_used;
            if (shift)
            {
                // Set new entries to zero
                for (unsigned i = 0; i < d_vectors_used; ++i)
                {
                    setZero(d_coeffs(i, d_vectors_used - 1));
                    setZero(d_coeffs(d_vectors_used - 1, i));
                }
                setZero(d_sqnorms[d_vectors_used - 1]);
                // Shift row
                for (unsigned i = d_vectors_used - 1; i > ofs; --i)
                    linalg::swap(d_coeffs.row(i), d_coeffs.row(i - 1));
                // Shift column
                for (unsigned i = d_vectors_used - 1; i > ofs; --i)
                    linalg::swap(d_coeffs.col(i), d_coeffs.col(i - 1));
                // Shift norms
                for (unsigned i = d_vectors_used - 1; i > ofs; --i)
                    arithmetic::swap(d_sqnorms[i], d_sqnorms[i - 1]);
                // Shift dependence info
                for (unsigned i = d_vectors_used - 1; i > ofs; --i)
                    d_dependent[i] = d_dependent[i - 1];
                d_dependent[ofs] = true; // the new vector is obviously dependent
            }
            return add_row;
        }
    
        template<class RealTypeContext, class IntTypeContext>
        void removeVector(RealTypeContext & d_rc, IntTypeContext & d_ic,
                          linalg::math_matrix<typename RealTypeContext::Real> & d_coeffs,
                          linalg::math_rowvector<typename RealTypeContext::Real> & d_sqnorms,
                          linalg::base_rowvector<bool> & d_dependent,
                          unsigned & d_vectors_used,
                          unsigned d_level,
                          unsigned ofs)
        // Removes a vector at position ofs.
        {
            // Remove row
            for (unsigned i = ofs + 1; i < d_level; ++i)
                linalg::swap(d_coeffs.row(i), d_coeffs.row(i - 1));
            // Remove column
            for (unsigned i = ofs + 1; i < d_level; ++i)
                linalg::swap(d_coeffs.col(i), d_coeffs.col(i - 1));
            // Remove norms
            for (unsigned i = ofs + 1; i < d_level; ++i)
                arithmetic::swap(d_sqnorms[i], d_sqnorms[i - 1]);
            // Remove dependence info
            for (unsigned i = ofs + 1; i < d_level; ++i)
                d_dependent[i - 1] = d_dependent[i];
            // Decrease number of rows/cols
            --d_vectors_used;
        }
    
        template<class RealTypeContext, class IntTypeContext>
        void compactify(RealTypeContext & d_rc, IntTypeContext & d_ic,
                        linalg::math_matrix<typename RealTypeContext::Real> & d_coeffs,
                        linalg::math_rowvector<typename RealTypeContext::Real> & d_sqnorms,
                        linalg::base_rowvector<bool> & d_dependent,
                        unsigned & d_vectors_used)
        // Makes sure no unused cols/rows are there
        {
            if (d_vectors_used < d_coeffs.rows())
            {
                d_coeffs->resize(d_vectors_used, d_vectors_used, linalg::Initialize(d_rc));
                d_sqnorms->resize(d_vectors_used, linalg::Initialize(d_rc));
                d_dependent.resize(d_vectors_used);
            }
        }
    }
}

#endif

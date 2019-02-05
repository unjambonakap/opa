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

#include <plll/config.hpp>
#include <plll/linalg.hpp>
#include "primes.hpp"
#include "factor.hpp"
#include <cmath>
#include <vector>

namespace plll
{
    namespace linalg
    {
        static unsigned long ModInv(unsigned long x, unsigned long p)
        // Find x^{-1} modulo p.
        {
            signed long xx = x, yy = p, a = 0, b = 1, aa = 1, bb = 0, t;
            while (yy != 0)
            {
                signed long q = xx / yy, r = xx % yy;
                xx = yy;
                yy = r;
                t = aa - q * a;
                aa = a;
                a = t;
                t = bb - q * b;
                bb = b;
                b = t;
            }
            // Now xx = gcd(x, p) = aa*x + bb*p
            if (aa < 0)
                aa += p;
            return aa;
        }
        
        template<class T>
        std::ostream& operator << (std::ostream& s, const std::vector<T> & v)
        {
            s << "[";
            bool first = true;
            for (typename std::vector<T>::const_iterator i = v.begin(); i != v.end(); ++i)
            {
                if (first)
                    first = false;
                else
                    s << ", ";
                s << *i;
            }
            s << "]";
            return s;
        }
        
        void computeRankProfile(std::vector<std::pair<unsigned, unsigned> > & rp, const math_matrix<arithmetic::Integer> & A, unsigned long p)
        // Computes the rank profile of A modulo p and stores it into rp. Returns false in case p is not a
        // prime. Assumes that p squared plus some more fits into unsigned long.
        {
            const unsigned m = A.rows(), n = A.cols();
            rp.clear();
            rp.reserve(std::min(m, n));
            
            arithmetic::Integer pp(p), qq, rr;
            math_matrix<unsigned long> Acpy;
            Acpy.resize(m, n);
            for (unsigned i = 0; i < m; ++i)
                for (unsigned j = 0; j < n; ++j)
                {
                    euclideanDivision(qq, rr, A(i, j), pp);
                    Acpy(i, j) = arithmetic::convert<long>(rr);
                }
            std::vector<unsigned> row_indices(m);
            for (unsigned i = 0; i < m; ++i)
                row_indices[i] = i;
            // Transform Acyp to row echelon form using row operations
            unsigned r = 0; // current rank
            for (unsigned i = 0; i < n; ++i)
            {
                // Iterate over columns
                unsigned nz = r;
                for (; nz < m; ++nz)
                    if (Acpy(nz, i) != 0)
                        break;
                if (nz == m)
                {
                    // No non-zero entry!
                    continue;
                }
                // Swap rows nz and r
                if (nz != r)
                {
                    swap(Acpy.row(r), Acpy.row(nz));
                    std::swap(row_indices[r], row_indices[nz]);
                }
                // Add entry (row_indices[r], i) to rank profile
                rp.push_back(std::make_pair(row_indices[r], i));
                // Make entry (r, i) one
                unsigned long a = ModInv(Acpy(r, i), p);
                Acpy(r, i) = 1;
                for (unsigned j = i + 1; j < n; ++j)
                    Acpy(r, j) = (Acpy(r, j) * a) % p;
                // Use the r-th row to clear out the rows below
                for (unsigned j = r + 1; j < m; ++j)
                {
                    if (Acpy(j, i) != 0)
                    {
                        unsigned long mul = p - Acpy(j, i);
                        Acpy(j, i) = 0;
                        for (unsigned k = i + 1; k < n; ++k)
                            Acpy(j, k) = (Acpy(j, k) + ((Acpy(r, k) * mul) % p)) % p;
                    }
                }
                // Increase rank
                ++r;
                if (r == m)
                    break;
            }
        }
        
        inline void apply2x2row(math_matrix<arithmetic::Integer> & B, unsigned i, unsigned j,
                                const arithmetic::Integer & A00, const arithmetic::Integer & A01,
                                const arithmetic::Integer & A10, const arithmetic::Integer & A11, unsigned c = 0)
        // [ B.row(i) ]   [ A00 A01 ]   [ B.row(i) ]
        // [ B.row(j) ] = [ A10 A11 ] * [ B.row(j) ]   (assumes that columns < c are zero in rows i and j)
        {
            arithmetic::Integer t;
            for (; c < B.cols(); ++c)
            {
                t = A00 * B(i, c);
                t += A01 * B(j, c);
                B(j, c) *= A11;
                B(j, c) += A10 * B(i, c);
                B(i, c) = t;
            }
        }
        
        inline void apply2x2row(math_matrix<arithmetic::Integer> & B, unsigned i, unsigned j,
                                const arithmetic::Integer & A00, const arithmetic::Integer & A01,
                                const arithmetic::Integer & A10, const arithmetic::Integer & A11, const arithmetic::Integer & mod, unsigned c = 0)
        // [ B.row(i) ]   [ A00 A01 ]   [ B.row(i) ]   (computations done modulo <mod>)
        // [ B.row(j) ] = [ A10 A11 ] * [ B.row(j) ]   (assumes that columns < c are zero in rows i and j)
        {
            arithmetic::Integer t;
            for (; c < B.cols(); ++c)
            {
                t = A00 * B(i, c);
                t %= mod;
                t += A01 * B(j, c);
                t %= mod;
                B(j, c) *= A11;
                B(j, c) %= mod;
                B(j, c) += A10 * B(i, c);
                B(j, c) %= mod;
                B(i, c) = t;
            }
        }
        
        unsigned modularHNF(math_matrix<arithmetic::Integer> & B, const arithmetic::Integer & h)
        // Applies row operations to B (by multiplying unimodular matrices from the left) to put B into
        // Hermite Normal Form. Works modulo h, where h should be a positive strictly larger multiple of the
        // determinant of B. Might not work if matrix is not square or not invertible or if h equals the
        // absolute value of the determinant/product of invariant factors.
        {
            assert(sign(h) > 0);
            unsigned r = 0;
            arithmetic::Integer d, a, b, aa, bb, quot, res;
            B %= h;
            for (unsigned c = 0; c < B.cols(); ++c)
            {
                // Combine rows r to B.rows()-1
                unsigned nz = r;
                for (; nz < B.rows(); ++nz)
                    if (!isZero(B(nz, c)))
                        break;
                if (nz == B.rows())
                    continue;
                
                // Work with row nz
                for (unsigned rr = nz + 1; rr < B.rows(); ++rr)
                    if (!isZero(B(rr, c)))
                    {
                        XGCD(d, a, b, B(nz, c), B(rr, c));
                        // Now d = a * B(nz, c) + b * B(rr, c).
                        if (c + 1 < B.cols())
                        {
                            aa = B(rr, c) / d;
                            bb = B(nz, c) / d;
                            bb = -bb;
                            apply2x2row(B, nz, rr, a, b, aa, bb, h, c + 1);
                        }
                        B(nz, c) = d;
                        setZero(B(rr, c));
                    }
                
                // Make B(nz, c) as small as possible
                XGCD(d, a, b, B(nz, c), h);
                if (!isOne(a))
                {
                    B.row(nz) *= a;
                    B.row(nz) %= h;
                }
                
                // Make B(nz, c) "positive"
                for (unsigned j = c; j < B.cols(); ++j)
                    if (sign(B(nz, j)) < 0)
                        B(nz, j) += h;
                if ((B(nz, c) << 1) > h)
                {
                    for (unsigned j = c; j < B.cols(); ++j)
                        if (!isZero(B(nz, j)))
                            B(nz, j) = h - B(nz, j);
                }
                
                // Swap row nz to r
                if (nz > r)
                    swap(B.row(r), B.row(nz));
                
                // Clean up above
                for (unsigned rr = 0; rr < r; ++rr)
                    if (!isZero(B(rr, c)))
                    {
                        euclideanDivision(quot, res, B(rr, c), B(r, c));
                        if (!isZero(quot))
                        {
                            if (isPMOne(quot))
                            {
                                if (sign(quot) > 0)
                                    B.row(rr) -= B.row(r);
                                else
                                    B.row(rr) += B.row(r);
                            }
                            else
                                B.row(rr) -= quot * B.row(r);
                            
                            B.row(rr) %= h;
                        }
                    }
                
                // Increase rank
                ++r;
            }
            return r;
        }
        
        bool addRow(math_matrix<arithmetic::Integer> & AA,
                    const std::vector<std::pair<unsigned, unsigned> > & rank_profile,
                    unsigned dst, unsigned src,
                    const math_matrix<arithmetic::Integer> & B)
        // Insert B.row(src) to AA.row(dst).
        {
            math_rowvector<arithmetic::Integer> b;
            b = B.row(src);
            arithmetic::Integer t, d, x, y, xx, yy;
            for (unsigned i = 0; i < rank_profile.size(); ++i)
            {
                if (isZero(b[rank_profile[i].second]))
                    continue;
                if (isZero(b[rank_profile[i].second] % AA(i, rank_profile[i].second)))
                {
                    t = b[rank_profile[i].second] / AA(i, rank_profile[i].second);
                    b -= t * AA.row(i);
                }
                else
                {
                    XGCD(d, x, y, AA(i, rank_profile[i].second), b[rank_profile[i].second]);
                    // Now d = x * B(nz, c) + y * B(rr, c).
                    xx = b[rank_profile[i].second] / d;
                    yy = AA(i, rank_profile[i].second) / d;
                    yy = -yy;
                    for (unsigned j = 0; j < b.size(); ++j)
                    {
                        t = b[j] * y;
                        t += AA(i, j) * x;
                        b[j] *= yy;
                        b[j] += AA(i, j) * xx;
                        AA(i, j) = t;
                    }
                    // Reduce rows above
                    for (unsigned k = 0; k < i; ++k)
                    {
                        if (!isZero(AA(k, rank_profile[i].second)))
                        {
                            euclideanDivision(x, y, AA(k, rank_profile[i].second), d);
                            if (!isZero(x))
                                AA.row(k) -= x * AA.row(i);
                        }
                    }
                }
            }
            // Check if b is zero
            bool zero = true;
            for (unsigned i = 0; i < b.size(); ++i)
                if (!isZero(b[i]))
                {
                    zero = false;
                    break;
                }
            if (zero)
                // If zero, there's nothing left to do
                return true;
            // This should not happen! Apparently, the rank profile was incorrect...
            return false;
        }
        
        unsigned hnf(math_matrix<arithmetic::Integer> & B)
        // Applies row operations to B (by multiplying unimodular matrices from the left) to put B into
        // Hermite Normal Form. Uses the 1996 algorithm by Storjohann and Labahn.
        {
            // Test for zero matrix
            bool zero = true;
            for (unsigned i = 0; (i < B.rows()) && zero; ++i)
                for (unsigned j = 0; j < B.cols(); ++j)
                    if (!isZero(B(i, j)))
                    {
                        zero = false;
                        break;
                    }
            if (zero)
                // In this case, the rank is 0
                return 0;
            
            // Now do the real work
            std::vector<std::pair<unsigned, unsigned> > rank_profile;
            arithmetic::RandomNumberGenerator rng;
            rng.randomizeTime(); // we don't need good quality random numbers
            unsigned long bound = std::sqrt(std::numeric_limits<unsigned long>::max()) / 2;
            unsigned minrank = 0;
            while (true)
            {
                unsigned long p = rng.random(bound);
                p = nextPrime(p);
                computeRankProfile(rank_profile, B, p);
                unsigned rank = rank_profile.size();
                if (rank == 0)
                    // p divides all entries; since we tested whether B is zero, we have to choose
                    // another prime
                    continue;
                if ((rank == B.rows()) && (B.rows() == B.cols()))
                    break;
                // If rank is smaller than before, try another prime p
                if (rank < minrank)
                    continue;
                minrank = rank;
                
//            std::cout << "Trying rank profile with rank " << rank_profile.size() << " obtained from prime " << p << "\n";
                
                // Get submatrix
                math_matrix<arithmetic::Integer> A;
                A.resize(rank, rank);
                for (unsigned i = 0; i < rank; ++i)
                    for (unsigned j = 0; j < rank; ++j)
                        A(i, j) = B(rank_profile[i].first, rank_profile[j].second);
                
                // Compute its determinant
                arithmetic::Integer h = det(A);
                abs(h, h);
                
                // Compute HNF of A (modulo 2*h)
                h <<= 1;
                modularHNF(A, h);
                
                // Extend to full matrix
                math_matrix<arithmetic::Integer> AA;
                AA.resize(B.rows(), B.cols());
                for (unsigned i = 0; i < rank; ++i)
                    for (unsigned j = 0; j < rank; ++j)
                        AA(i, rank_profile[j].second) = A(i, j);
                
                // Add columns
                if (rank < B.cols())
                {
                    // Here, we are in the nice situation that the missing column entries in the first
                    // <rank> rows are uniquely determined from the corresponding entries of the pivot rows
                    // of the original matrix. Therefore, we can recover them using a p-adic method. As a
                    // prime, we can re-use p.
                    
                    // Get submatrix
                    math_matrix<arithmetic::Integer> BB;
                    BB.resize(rank, rank);
                    for (unsigned i = 0; i < rank; ++i)
                        for (unsigned j = 0; j < rank; ++j)
                            BB(i, j) = B(rank_profile[i].first, rank_profile[j].second);
                    
                    // Solve
                    math_matrix<arithmetic::Integer> res = solveUniqInt(BB.transpose(), A.transpose(), arithmetic::convert<arithmetic::Integer>(p));
                    if (res.rows() == 0)
                        continue;
                    
                    // Compute missing entries
                    for (unsigned k = 0; k < rank_profile[0].second; ++k)
                    {
                        for (unsigned i = 0; i < rank; ++i)
                            for (unsigned ell = 0; ell < rank; ++ell)
                                AA(i, k) += B(rank_profile[ell].first, k) * res(ell, i);
                    }
                    for (unsigned i = 1; i < rank; ++i)
                        for (unsigned k = rank_profile[i - 1].second + 1; k < rank_profile[i].second; ++k)
                        {
                            for (unsigned i = 0; i < rank; ++i)
                                for (unsigned ell = 0; ell < rank; ++ell)
                                    AA(i, k) += B(rank_profile[ell].first, k) * res(ell, i);
                        }
                    for (unsigned k = rank_profile[rank - 1].second + 1; k < B.cols(); ++k)
                    {
                        for (unsigned i = 0; i < rank; ++i)
                            for (unsigned ell = 0; ell < rank; ++ell)
                                AA(i, k) += B(rank_profile[ell].first, k) * res(ell, i);
                    }
                }
                
                // Add rows
                if (rank < B.rows())
                {
                    unsigned r = rank;
                    bool try_again = false;
                    for (unsigned k = 0; k < rank_profile[0].first; ++k)
                    {
                        if (!addRow(AA, rank_profile, r, k, B))
                        {
                            try_again = true;
                            break;
                        }
                        ++r;
                    }
                    if (try_again)
                        continue;
                    for (unsigned i = 1; (i < rank) && !try_again; ++i)
                        for (unsigned k = rank_profile[i - 1].first + 1; k < rank_profile[i].first; ++k)
                        {
                            if (!addRow(AA, rank_profile, r, k, B))
                            {
                                try_again = true;
                                break;
                            }
                            ++r;
                        }
                    if (try_again)
                        continue;
                    for (unsigned k = rank_profile[rank - 1].first + 1; k < B.rows(); ++k)
                    {
                        if (!addRow(AA, rank_profile, r, k, B))
                        {
                            try_again = true;
                            break;
                        }
                        ++r;
                    }
                    if (try_again)
                        continue;
                }
                
                swap(AA, B);
                return rank;
            }
            // If this line is reached, then the matrix is square and has full rank.
            
            // Compute determinant
            arithmetic::Integer h = det(B);
            abs(h, h);
            
            // Compute HNF of B (modulo 2*h)
            h <<= 1;
            return modularHNF(B, h);
        }
        
        unsigned hnf(math_matrix<arithmetic::Integer> & B, math_matrix<arithmetic::Integer> & T)
        // Applies row operations to B (by multiplying unimodular matrices from the left) to put B into
        // Hermite Normal Form.
        {
            math_matrix<arithmetic::Integer> A;
            A.resize(B.rows(), B.cols() + B.rows());
            for (unsigned i = 0; i < B.rows(); ++i)
            {
                for (unsigned j = 0; j < B.cols(); ++j)
                    A(i, j) = B(i, j);
                setOne(A(i, B.cols() + i));
            }
            // Compute classical HNF
            hnf(A);
            // Restore B from left part of A, and retrieve rank
            unsigned rank = 0;
            for (unsigned i = 0; i < B.rows(); ++i)
            {
                bool zero = true;
                for (unsigned j = 0; j < B.cols(); ++j)
                {
                    B(i, j) = A(i, j);
                    if (!isZero(B(i, j)))
                        zero = false;
                }
                if (!zero)
                    ++rank;
            }
            // Create T from right part of A
            T.resize(B.rows(), B.rows());
            for (unsigned i = 0; i < B.rows(); ++i)
                for (unsigned j = 0; j < B.rows(); ++j)
                    T(i, j) = A(i, B.cols() + j);
            // Return rank
            return rank;
        }
    
        arithmetic::Integer det(const math_matrix<arithmetic::Integer> & A, const arithmetic::Integer & p)
        // Assumes that A is a square matrix and p a prime. Computes the determinant of A modulo p.
        {
            assert(A.rows() == A.cols());
            math_matrix<arithmetic::Integer> Acpy;
            Acpy = A % p;
            arithmetic::Integer res;
            setOne(res);
            // Transform Acyp to identity using row operations
            for (unsigned i = 0; i < A.rows(); ++i)
            {
                unsigned nz = i;
                for (; nz < A.rows(); ++nz)
                    if (!isZero(Acpy(nz, i)))
                        break;
                if (nz == A.rows())
                {
                    // No non-zero entry!
                    return arithmetic::Integer(); // determinant is 0
                }
                // Swap rows nz and i
                if (nz != i)
                {
                    swap(Acpy.row(i), Acpy.row(nz));
                    res = -res;
                }
                // Make entry (i, i) one
                arithmetic::Integer r, a, b;
                XGCD(r, a, b, Acpy(i, i), p); // r == a * Acpy(i, i) + b * p;
                if (!isOne(r))
                    return arithmetic::Integer(); // can only happen if r is zero (or if p is not a prime!)
                res *= Acpy(i, i);
                res %= p;
                Acpy.row(i) *= a;
                Acpy.row(i) %= p;
                // Use the i-th row to clear out all other rows
                for (unsigned j = 0; j < A.rows(); ++j)
                    if (j != i)
                    {
                        arithmetic::Integer m = Acpy(j, i);
                        Acpy.row(j) -= Acpy.row(i) * m;
                        Acpy.row(j) %= p;
                    }
            }
            // Move result to right range
            if (sign(res) < 0)
                res += p;
            return res;
        }
    
        arithmetic::Integer det(const math_matrix<arithmetic::Integer> & A)
        // Computes the determinant of A. Uses a CRT approach
        {
            assert(A.rows() == A.cols());
        
            arithmetic::Integer B;
            // Use two methods to compute a bound on the determinant. Choose the minimum of the two
            // bounds.
            {
                // Firstly, compute bound on modulus using Hadamard's inequality on the rows
                arithmetic::Integer B1;
                setOne(B1);
                for (unsigned i = 0; i < A.rows(); ++i)
                    B1 *= normSq(A.row(i));
                B1 = sqrtCeil(B1);
                B = B1;
            }
            {
                // Secondly, compute bound on modulus using Hadamard's inequality on the columns
                arithmetic::Integer B1;
                setOne(B1);
                for (unsigned i = 0; i < A.cols(); ++i)
                    B1 *= normSq(A.col(i));
                B1 = sqrtCeil(B1);
                if (B > B1)
                    B = B1;
            }
            // times 2 plus 1
            B <<= 1;
            ++B;
        
            // Next, collect primes
            FactorizedNumber modulus;
            arithmetic::Integer p;
//    std::cout << "Using primes";
            while (modulus.number() < B)
            {
                p = nextPrime(p);
                modulus.mulCoprimeLarger(FactorizedNumber(p, 1));
//        std::cout << " " << p;
            }
//    std::cout << "\n";
//    std::cout << "modulus: " << modulus.number() << "\n";
//    std::cout << "bound:   " << B << "\n";
        
            // Prepare Chinese Remaindering
            ChineseRemainder CR(modulus);
            for (ChineseRemainder::Iterator i = CR.begin(); i != CR.end(); ++i)
            {
                arithmetic::Integer v = det(A, i.prime());
//        std::cout << "det = " << v << " (mod " << i.prime() << ")\n";
                i.setValue(v);
            }
        
            return CR.process(true);
        }
    
        bool modularInvert(const math_matrix<arithmetic::Integer> & A,
                           math_matrix<arithmetic::Integer> & Ainv,
                           const arithmetic::Integer & p)
        // Assumes that A is a square matrix and p a prime. Tries to invert A modulo p; if the inversion is
        // successful, the result will be stored in Ainv and true is returned. Otherwise, false is returned.
        {
            assert(A.rows() == A.cols());
            // Set Ainv = [I A], where I is the identity matrix (and make everything modulo p)
            Ainv.resize(0, 0);
            Ainv.resize(A.rows(), A.cols() * 2);
            for (unsigned i = 0; i < A.rows(); ++i)
            {
                setOne(Ainv(i, i));
                for (unsigned j = 0; j < A.cols(); ++j)
                    Ainv(i, j + A.cols()) = A(i, j) % p;
            }
            // Transform A part of Ainv to identity using row operations
            for (unsigned i = 0; i < A.rows(); ++i)
            {
                unsigned nz = i;
                for (; nz < A.rows(); ++nz)
                    if (!isZero(Ainv(nz, A.cols() + i)))
                        break;
                if (nz == A.rows())
                {
                    // No non-zero entry!
                    return false;
                }
                // Swap rows nz and i
                if (nz != i)
                    swap(Ainv.row(i), Ainv.row(nz));
                // Make entry (i, A.cols() + i) one
                arithmetic::Integer r, a, b;
                XGCD(r, a, b, Ainv(i, A.cols() + i), p); // r == a * Ainv(i, A.cols() + i) + b * p;
                if (!isOne(r))
                {
                    Ainv.resize(0, 0);
                    return false;
                }
                Ainv.row(i) *= a;
                Ainv.row(i) %= p;
                // Use the i-th row to clear out all other rows
                for (unsigned j = 0; j < A.rows(); ++j)
                    if (j != i)
                    {
                        arithmetic::Integer m = Ainv(j, A.cols() + i);
                        Ainv.row(j) -= Ainv.row(i) * m;
                        Ainv.row(j) %= p;
                    }
            }
            // Get rid of identity part
            Ainv.resize(A.rows(), A.cols());
            // Make sure everything is in range [0, p)
            for (unsigned i = 0; i < A.rows(); ++i)
                for (unsigned j = 0; j < A.cols(); ++j)
                    if (sign(Ainv(i, j)) < 0)
                        Ainv(i, j) += p;
            return true;
        }
    
        math_matrix<arithmetic::Integer> solveUniqInt(const math_matrix<arithmetic::Integer> & A,
                                                      const math_matrix<arithmetic::Integer> & w,
                                                      const arithmetic::Integer & startPrime)
        // Uses p-adic iterative solver; assumes that A is square and that the equation A * res = v has
        // a unique integral solution res.
        {
            assert(A.rows() == A.cols());
            assert(A.rows() == w.rows());
            /*
              Assume that  A * v = w  holds modulo p^e. Write v' = v + p^e v''; then
              A v' = w                    (modulo p^{e+1})
              <=> p^e A v'' = w - A v         (modulo p^{e+1})
              <=> A v'' = (w - A v) / p^e     (modulo p)
              yields
              v'' = A^{-1} ([w - A v] / p^e)  (modulo p).
            */
            math_matrix<arithmetic::Integer> res;
            // Find a p such that A is invertible modulo p.
            arithmetic::Integer p = startPrime;
            if (p < arithmetic::convert<arithmetic::Integer>(2))
                p = arithmetic::convert<arithmetic::Integer>(2);
            math_matrix<arithmetic::Integer> Ainvp; // inverse of A modulo p
            while (true)
            {
//            std::cout << "Inverting modulo " << p << "...\n";
                if (modularInvert(A, Ainvp, p))
                    break;
                p = nextPrime(p);
            }
//        std::cout << "A is invertible modulo " << p << "\n";
            // Estimate maximal exponent needed (using Cramer's rule to bound solution)
            arithmetic::Integer maxentry;
            setZero(maxentry);
            for (unsigned i = 0; i < A.rows(); ++i)
                for (unsigned j = 0; j < A.cols(); ++j)
                    if (compareAbsValues(maxentry, A(i, j)) < 0)
                        maxentry = abs(A(i, j));
            unsigned maxe = 0;
            while (!isZero(maxentry))
            {
                ++maxe;
                maxentry /= p;
            }
            maxe *= A.cols() - 1;
            setZero(maxentry);
            for (unsigned i = 0; i < w.rows(); ++i)
                for (unsigned j = 0; j < w.cols(); ++j)
                    if (compareAbsValues(maxentry, w(i, j)) < 0)
                        maxentry = abs(w(i, j));
            unsigned e = 0;
            while (!isZero(maxentry))
            {
                ++e;
                maxentry /= p;
            }
            maxe += e;
            // Estimate size of A.cols()!
            setOne(maxentry);
            for (unsigned i = 2; i <= A.cols(); ++i)
                maxentry *= arithmetic::convert<arithmetic::Integer>(i);
            e = 0;
            while (!isZero(maxentry))
            {
                ++e;
                maxentry /= p;
            }
            maxe += e;
//        std::cout << "Maximal exponent is " << maxe << "\n";
            // Iteratively find a solution
            e = 0;
            arithmetic::Integer pe; // p^e
            setOne(pe);
            res.resize(A.cols(), w.cols()); // res is the all-zero vector
            while (true)
            {
                math_matrix<arithmetic::Integer> err;
                err = w - A * res;
                bool zero = true;
                for (unsigned i = 0; (i < err.rows()) && zero; ++i)
                    for (unsigned j = 0; j < err.cols(); ++j)
                        if (!isZero(err(i, j)))
                        {
                            zero = false;
                            break;
                        }
                if (zero)
                    break;
                if (e > maxe)
                {
//                std::cout << "No solution found!\n";
                    res.resize(0, 0);
                    break;
                }
                arithmetic::Integer pee = pe * p;
//            std::cout << "Computing solution modulo " << p << "^" << e+1 << "\n";
                // Divide by p^e
                err /= pe;
                // Multiply with Ainvp modulo p, and then scale with p^e
                err = Ainvp * err;
                for (unsigned i = 0; i < err.rows(); ++i)
                    for (unsigned j = 0; j < err.cols(); ++j)
                    {
                        err(i, j) %= p;
                        err(i, j) *= pe;
                    }
                // Add to solution
                arithmetic::Integer peeh = pee >> 1; // pee/2
                res += err;
                for (unsigned i = 0; i < err.rows(); ++i)
                    for (unsigned j = 0; j < err.cols(); ++j)
                    {
                        res(i, j) %= pee;
                        if (sign(res(i, j)) < 0)
                            res(i, j) += pee;
                        if (res(i, j) > peeh)
                            res(i, j) -= pee;
                    }
                // Increase exponent
                pe = pee;
                ++e;
                // Now A * res = w modulo p^e.
            }
//        std::cout << "Found solution\n";
            return res;
        }
    
        void invert(math_matrix<arithmetic::Integer> & res,
                    const math_matrix<arithmetic::Integer> & A,
                    const arithmetic::Integer & startPrime)
        // Uses p-adic iterative solver; assumes that A is square and that the equation A * res = v has
        // a unique integral solution res.
        {
            assert(A.rows() == A.cols());
            /*
              Assume that  A * v = I  holds modulo p^e. Write v' = v + p^e v''; then
              A v' = I                    (modulo p^{e+1})
              <=> p^e A v'' = I - A v         (modulo p^{e+1})
              <=> A v'' = (I - A v) / p^e     (modulo p)
              yields
              v'' = A^{-1} ([I - A v] / p^e)  (modulo p).
            */
            // Find a p such that A is invertible modulo p.
            arithmetic::Integer p = startPrime;
            if (p < arithmetic::convert<arithmetic::Integer>(2))
                p = arithmetic::convert<arithmetic::Integer>(2);
            math_matrix<arithmetic::Integer> Ainvp; // inverse of A modulo p
            while (true)
            {
//            std::cout << "Inverting modulo " << p << "...\n";
                if (modularInvert(A, Ainvp, p))
                    break;
                p = nextPrime(p);
            }
//        std::cout << "A is invertible modulo " << p << "\n";
            // Estimate maximal exponent needed (using Cramer's rule to bound solution)
            arithmetic::Integer maxentry;
            setZero(maxentry);
            for (unsigned i = 0; i < A.rows(); ++i)
                for (unsigned j = 0; j < A.cols(); ++j)
                    if (compareAbsValues(maxentry, A(i, j)) < 0)
                        maxentry = abs(A(i, j));
            unsigned maxe = 0;
            while (!isZero(maxentry))
            {
                ++maxe;
                maxentry /= p;
            }
            maxe *= A.cols() - 1;
            ++maxe; // for entry 1
            // Estimate size of A.cols()!
            setOne(maxentry);
            for (unsigned i = 2; i <= A.cols(); ++i)
                maxentry *= arithmetic::convert<arithmetic::Integer>(i);
            unsigned e = 0;
            while (!isZero(maxentry))
            {
                ++e;
                maxentry /= p;
            }
            maxe += e;
//        std::cout << "Maximal exponent is " << maxe << "\n";
            // Iteratively find a solution
            e = 1;
            if (maxe < 2)
                maxe = 2;
            arithmetic::Integer pe; // p^e
            setOne(pe);
            res = Ainvp;
            while (true)
            {
                math_matrix<arithmetic::Integer> err;
                err = A * res;
                err = -err;
                for (unsigned i = 0; i < err.rows(); ++i)
                    ++err(i, i);
                bool zero = true;
                for (unsigned i = 0; (i < err.rows()) && zero; ++i)
                    for (unsigned j = 0; j < err.cols(); ++j)
                        if (!isZero(err(i, j)))
                        {
                            zero = false;
                            break;
                        }
                if (zero)
                    break;
                if (e > maxe)
                {
//                std::cout << "No solution found!\n";
                    res.resize(0, 0);
                    break;
                }
                arithmetic::Integer pee = pe * p;
//            std::cout << "Computing solution modulo " << p << "^" << e+1 << "\n";
                // Divide by p^e
                err /= pe;
                // Multiply with Ainvp modulo p, and then scale with p^e
                err = Ainvp * err;
                for (unsigned i = 0; i < err.rows(); ++i)
                    for (unsigned j = 0; j < err.cols(); ++j)
                    {
                        err(i, j) %= p;
                        err(i, j) *= pe;
                    }
                // Add to solution
                arithmetic::Integer peeh = pee >> 1; // pee/2
                res += err;
                for (unsigned i = 0; i < err.rows(); ++i)
                    for (unsigned j = 0; j < err.cols(); ++j)
                    {
                        res(i, j) %= pee;
                        if (sign(res(i, j)) < 0)
                            res(i, j) += pee;
                        if (res(i, j) > peeh)
                            res(i, j) -= pee;
                    }
                // Increase exponent
                pe = pee;
                ++e;
                // Now A * res = I modulo p^e.
            }
//        std::cout << "Found solution\n";
        }
    
        bool solve(math_colvector<arithmetic::Integer> & sol,
                   arithmetic::Integer & g,
                   const math_matrix<arithmetic::Integer> & A,
                   const math_colvector<arithmetic::Integer> & w)
        // Computes a minimal t > 0 such that A x = t w is solvable, and returns a solution together
        // with t.  In case no solution exists, false is returned.
        {
            assert(A.rows() == w.size());
        
//    std::cout << "A = " << A << "\n";
//    std::cout << "w = " << w << "\n";
        
            math_matrix<arithmetic::Integer> AA(A.cols() + 1, A.rows()), T;
            for (unsigned i = 0; i < A.rows(); ++i)
            {
                for (unsigned j = 0; j < A.cols(); ++j)
                    AA(j, i) = A(i, j);
                AA(A.cols(), i) = -w[i];
            }
//    std::cout << "(A, -w)^T = " << AA << "\n";
            unsigned r = hnf(AA, T);
            if (r == AA.rows())
                return false;
//    std::cout << "hnf((A, -b)^T) = " << AA << "\n";
        
            // Now the last A.cols()+1-r = T.rows()-r rows of T are a Z-basis of the "solution
            // space". We have to find a vector in the solution space with minimal last coefficient. The
            // minimal entry equals the GCD of all last entries of the basis vectors.

            // Start with last row
            g = T(T.rows() - 1, T.cols() - 1);
            sol.resize(A.cols());
            for (unsigned i = 0; i < A.cols(); ++i)
                sol[i] = T(T.rows() - 1, i);
        
            if (isPMOne(g))
            {
                // We are done instantly!
                if (isNegative(g))
                {
                    // Make sure g is positive
                    makeAbs(g);
                    for (unsigned j = 0; j < A.cols(); ++j)
                        sol[j] = -sol[j];
                }
                return true;
            }
        
            // Work on the last T.rows()-r rows (minus the last one)
            for (unsigned i = r; i < T.rows() - 1; ++i)
            {
                if (!arithmetic::isZero(T(i, T.cols() - 1)))
                {
                    arithmetic::Integer u, v;
                    arithmetic::XGCD(g, u, v, g, T(i, T.cols() - 1)); // g_new = u*g_old + v*T(i, T.cols() - 1)
                    for (unsigned j = 0; j < A.cols(); ++j)
                    {
                        sol[j] *= u;
                        sol[j] += v * T(i, j);
                    }
                }
                if (isPMOne(g))
                    break;
            }
        
            // Make sure g is positive
            if (isNegative(g))
            {
                makeAbs(g);
                for (unsigned j = 0; j < A.cols(); ++j)
                    sol[j] = -sol[j];
            }
        
            // Shortcut: g == 1
            if (isOne(g))
                return true;
        
            // Finally, compute GCD of g and solution entries
            arithmetic::Integer G = g;
            for (unsigned i = 0; i < sol.size(); ++i)
                if (!isZero(sol[i]))
                {
                    arithmetic::GCD(G, G, sol[i]);
                    if (isOne(G))
                        break;
                }
            if (!isOne(G))
            {
                // Divide out GCD!
                g /= G;
                sol /= G;
            }
            return true;
        }
    
        std::pair<math_colvector<arithmetic::Integer>, arithmetic::Integer> solve(const math_matrix<arithmetic::Integer> & A,
                                                                                  const math_colvector<arithmetic::Integer> & w)
        // Computes a minimal t > 0 such that A x = t w is solvable, and returns a solution together with t.
        // In case no solution exists, a random vector together with 0 is returned.
        {
            std::pair<math_colvector<arithmetic::Integer>, arithmetic::Integer> res;
            solve(res.first, res.second, A, w);
            return res;
        }
    
        math_colvector<arithmetic::Integer> solveInt(const math_matrix<arithmetic::Integer> & A,
                                                     const math_colvector<arithmetic::Integer> & w)
        // Assumes that A*v = w has a solution v in the integers. Computes v.
        {
            math_colvector<arithmetic::Integer> res;
            arithmetic::Integer g;
            if (solve(res, g, A, w))
                if (!isOne(g))
                    res.resize(0);
            return res;
        }
    
        math_matrix<arithmetic::Integer> kernel(const math_matrix<arithmetic::Integer> & A)
        // The columns of the returned matrix contain a Z-basis for the kernel of A, i.e. every vector v
        // with A v = 0 has a unique representation as a Z-linear combination of the columns of the returned
        // matrix.
        {
            math_matrix<arithmetic::Integer> AA, T;
            AA.resize(A.cols(), A.rows());
            transpose(AA, A);
            int r = hnf(AA, T);
            // Extract the last A.cols() - r rows from T as the columns of the result
            AA.resize(A.cols(), A.cols() - r);
            for (unsigned c = 0; c < A.cols() - r; ++c)
                for (unsigned i = 0; i < A.cols(); ++i)
                    AA(i, c) = T(r + c, i);
            return AA;
        }
    }
}

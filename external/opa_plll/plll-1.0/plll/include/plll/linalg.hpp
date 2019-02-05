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

#ifndef PLLL_INCLUDE_GUARD__LINALG_HPP
#define PLLL_INCLUDE_GUARD__LINALG_HPP

#include <plll/arithmetic.hpp>
#include <plll/matrix.hpp>
#include <utility>

/**
   \file
   \brief Basic linear algebra algorithms.
   
   This header provides basic linear algebra algorithms over `plll` integers
   (`plll::arithmetic::Integer`), such as solving systems of equations, computing determinants and
   Hermite Normal Form.
*/
namespace plll
{
    namespace linalg
    {
        /**@{
           \name Linear algebra over rationals and integers.
        */
        
        /**
           \brief Computes the Hermite Normal Form (HNF) of the given matrix `A` by applying row
                  operations.
           
           The resulting HNF is an upper triangular matrix, with all zero rows moved to the bottom.
           Note that the matrix in `A` is modified. The algorithm returns the rank of `A`.
           
           Implements the 1996 algorithm by Storjohann and Labahn \cite storjohann-labahn-hnf.
           
           \param A The matrix to compute the HNF of. This matrix is modified.
           \return The rank of `A`.
         */
        unsigned hnf(math_matrix<arithmetic::Integer> & A);
        /**
           \brief Computes the Hermite Normal Form (HNF) of the given matrix `A` by applying row
                  operations, together with a transformation matrix `T`.

           The resulting HNF is an upper triangular matrix, with all zero rows moved to the bottom.
           Note that the matrices in `A` and `T` are modified. The algorithm returns the rank of
           `A`, and the resulting `T` satisfies \f$T \cdot A_{original} = A_{HNF}\f$.
           
           Implements the 1996 algorithm by Storjohann and Labahn \cite storjohann-labahn-hnf.
           
           \param A The matrix to compute the HNF of. This matrix is modified.
           \param T The matrix where the transformation matrix is stored in.
           \return The rank of `A`.
         */
        unsigned hnf(math_matrix<arithmetic::Integer> & A, math_matrix<arithmetic::Integer> & T);

        /**
           \brief Computes the determinant of an integer matrix `A`.

           It does this by first finding an upper bound on it, then computing the determinant modulo
           many small primes, and finally recovering the result using an effective version of the
           Chinese Remainder Theorem.
           
           \param A The matrix whose determinant should be computed. Must be square.
           \return The determinant of the matrix.
         */
        arithmetic::Integer det(const math_matrix<arithmetic::Integer> & A);
        
        /**
           \brief Assumes that `A` is a square matrix and that the equation `A * res == w` has a
                  unique solutions in the integers. Returns the
           
           Note that the assumption is the case if and only if `det(A) != 0`. solution `res`.
           
           In case no solution exists, the algorithm might not terminate. The algorithm uses p-adic
           methods by first finding a prime `p` not dividing the determinant of `A`. If special
           knowledge of the determinant is available, a prime not dividing it can be specified as
           the `startPrime` parameter.
           
           Note that `w` and thus also `res` need not to be vectors, but can also be matrices with
           any number of columns.
           
           \param A The left-hand side of the linear equations.
           \param w The right-hand side of the linear equations.
           \param startPrime The prime to start at. If given, should not be much less than a prime
                             not dividing the determinant of `A`.
           \return A solution `res` to `A * res == w`.
         */
        math_matrix<arithmetic::Integer> solveUniqInt(const math_matrix<arithmetic::Integer> & A,
                                                      const math_matrix<arithmetic::Integer> & w,
                                                      const arithmetic::Integer & startPrime = arithmetic::Integer());
        
        /**
           \brief Assumes that `A` is invertible. Inverts the matrix `A` and stores the result in
                  `inverse`.

           In case no solution exists, the algorithm might not terminate. The algorithm uses p-adic
           methods by first finding a prime `p` not dividing the determinant of `A`. If special
           knowledge of the determinant is available, a prime not dividing it can be specified as
           the `startPrime` parameter.
           
           \param inverse A matrix where the inverse is stored in.
           \param A The matrix to invert.
           \param startPrime The prime to start at. If given, should not be much less than a prime
                             not dividing the determinant of `A`.
         */
        void invert(math_matrix<arithmetic::Integer> & inverse,
                    const math_matrix<arithmetic::Integer> & A,
                    const arithmetic::Integer & startPrime = arithmetic::Integer());
        
        /**
           \brief Solves the linear system `A * x == b` over the rationals.

           Determines the minimal positive integer `t` such that `A * x == t * b` is solvable with
           an integral vector `x`, and returns `t` together with an arbitrary such solution `x`; the
           quotient `x/t` is then a "minimal" solution.
           
           In case no solution exists, it returns a random vector together with the integer `0`.
           
           \param A The left-hand side of the linear equations.
           \param b The right-hand side of the linear equations.
           \return A pair `(x, t)` such that `A * x == t * b` with `t > 0` minimal such that such an
                   `x` exists.
         */
        std::pair<math_colvector<arithmetic::Integer>, arithmetic::Integer> solve(const math_matrix<arithmetic::Integer> & A,
                                                                                  const math_colvector<arithmetic::Integer> & b);
        
        /**
           \brief Assumes that the linear system `A * v == w` has a solution in the integers. Finds
                  this solution.
           
           This algorithm essentially uses `plll::linalg::solveUniqInt()` and might also not
           terminate in case no solution exists.
           
           \param A The left-hand side of the linear equations.
           \param w The right-hand side of the linear equations.
           \return A vector `res` such that `A * res == b`.
         */
        math_colvector<arithmetic::Integer> solveInt(const math_matrix<arithmetic::Integer> & A,
                                                     const math_colvector<arithmetic::Integer> & w);

        /**
           \brief Computes a \f$\mathbb{Z}\f$-basis of the right-kernel of `A`.
           
           The right-kernel of `A` is the \f$\mathbb{Z}\f$-module \f$\{ v \mid A v = 0 \}\f$. A
           basis is a set of vectors from the kernel such that every vector in the kernel can be
           uniquely written as a linear combination of these vectors.
           
           \param A The matrix to compute the kernel of.
           \return A matrix whose columns form a basis of the kernel of `A`.
         */
        math_matrix<arithmetic::Integer> kernel(const math_matrix<arithmetic::Integer> & A);
        
        ///@}
    }
}

#endif

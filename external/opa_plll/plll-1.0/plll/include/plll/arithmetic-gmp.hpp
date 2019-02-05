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

#ifndef PLLL_INCLUDE_GUARD__ARITHMETIC_GMP_HPP
#define PLLL_INCLUDE_GUARD__ARITHMETIC_GMP_HPP

#include <plll/config.hpp>

/**
   \file
   \brief Header for arbitrary precision integer and floating point arithmetic, provided by GMP and MPFR.
   
   This header provides arbitrary precision integer and floating point arithmetic for `plll`. It provides
   `plll::arithmetic::IntegerContext` and `plll::arithmetic::RealContext`, two arithmetic contexts as described in \ref
   arithctxts.

   Concrete integers are represented by `plll::arithmetic::Integer` objects and concrete floating
   point numbers are represented by `plll::arithmetic::Real` objects. All expressions are
   represented by `plll::arithmetic::expressions::Expression<>` templates.
*/
namespace plll
{
    namespace arithmetic
    {
#ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
        namespace internal
        {
            extern bool Real_precision_check_enabled;
        }
        
        inline void Real_precision_check_disable() { arithmetic::internal::Real_precision_check_enabled = false; }
        inline void Real_precision_check_enable() { arithmetic::internal::Real_precision_check_enabled = true; }
#else
        inline void Real_precision_check_disable() { }
        inline void Real_precision_check_enable() { }
#endif
    }
}

#include <stdint.h>
#define MPFR_USE_INTMAX_T 1

#include <gmp.h>
#include <mpfr.h>
#include <cmath>
#include <math.h>
#include <limits>
#include <iosfwd>
#include <cassert>
#include <utility>
#include <algorithm>
#include <plll/helper.hpp>

namespace plll
{
    namespace arithmetic
    {
#if (MPFR_VERSION <= 0x0204FF)
        // Versions 2.4.x or before: define MPFR_RNDN and others
        #define MPFR_RNDN GMP_RNDN
        #define MPFR_RNDZ GMP_RNDZ
        #define MPFR_RNDU GMP_RNDU
        #define MPFR_RNDD GMP_RNDD
        typedef mp_exp_t mpfr_exp_t;
#endif
#if (__GMP_MP_RELEASE < 50000)
        // Before version 5.0.0: define mp_bitcnt_t
        typedef unsigned long mp_bitcnt_t;
#endif
        
        class Integer;
        class IntegerContext;
        class Real;
        class RealContext;

        /**@{
           \name Swap functions.
        */
        /**
           \brief Swaps two `plll::arithmetic::Integer` objects.
         */
        void swap(plll::arithmetic::Integer &, plll::arithmetic::Integer &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        
        /**
           \brief Swaps two `plll::arithmetic::Real` objects.
         */
        void swap(plll::arithmetic::Real &, plll::arithmetic::Real &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        ///@}
        
        /**
           \brief Initializes the local thread allocator for the current thread, and sets the GMP
                  and MPFR allocators to the local thread allocator.
           
           This function must be called at the beginning of each thread (and, most importantly, of
           the process itself!) before any `plll::arithmetic::Integer` or `pll::arithmetic::Real`
           variable was initialized.
           
           \warning In case the program uses GMP and/or MPFR directly, this function must also be
                    called before any GMP respectively MPFR method is used directly.
                    
                    Violation of this will quite surely lead to crashes and/or segmentation faults.
         */
        void initArithmeticThreadAllocators();
        
        /**@{
           \name Predicates.
        */
        /**
           \brief Tests the given `plll::arithmetic::Integer` object for being zero.
           
           \return Returns `true` if and only if the argument is zero.
           \sa \ref arithctxts_preds
        */
        inline bool isZero(const Integer &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        
        /**
           \brief Tests the given `plll::arithmetic::Real` object for being zero.
           
           \return Returns `true` if and only if the argument is zero.
           \sa \ref arithctxts_preds
        */
        inline bool isZero(const Real &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        
        /**
           \brief Tests the given `plll::arithmetic::Integer` object for being one.
           
           \return Returns `true` if and only if the argument is one.
           \sa \ref arithctxts_preds
        */
        inline bool isOne(const Integer &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        
        /**
           \brief Tests the given `plll::arithmetic::Integer` object for being one or minus one.
           
           \return Returns `true` if and only if the argument is \f$\pm 1\f$.
           \sa \ref arithctxts_preds
        */
        inline bool isPMOne(const Integer &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        
        /**
           \brief Tests the given `plll::arithmetic::Integer` object for being two or minus two.
           
           \return Returns `true` if and only if the argument is \f$\pm 2\f$.
           \sa \ref arithctxts_preds
        */
        inline bool isPMTwo(const Integer &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        
        /**
           \brief Tests the given `plll::arithmetic::Real` object for being one.
           
           \return Returns `true` if and only if the argument is one.
           \sa \ref arithctxts_preds
        */
        inline bool isOne(const Real &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        
        /**
           \brief Tests the given `plll::arithmetic::Integer` object for being strictly positive.
           
           \return Returns `true` if and only if the argument is strictly positive.
           \sa \ref arithctxts_preds
        */
        inline bool isPositive(const Integer &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        
        /**
           \brief Tests the given `plll::arithmetic::Real` object for being strictly positive.
           
           \return Returns `true` if and only if the argument is strictly positive.
           \sa \ref arithctxts_preds
        */
        inline bool isPositive(const Real &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        
        /**
           \brief Tests the given `plll::arithmetic::Integer` object for being positive or zero.
           
           \return Returns `true` if and only if the argument is positive or zero.
           \sa \ref arithctxts_preds
        */
        inline bool isNonNegative(const Integer &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        
        /**
           \brief Tests the given `plll::arithmetic::Real` object for being positive or zero.
           
           \return Returns `true` if and only if the argument is positive or zero.
           \sa \ref arithctxts_preds
        */
        inline bool isNonNegative(const Real &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        
        /**
           \brief Tests the given `plll::arithmetic::Integer` object for being strictly negative.
           
           \return Returns `true` if and only if the argument is strictly negative.
           \sa \ref arithctxts_preds
        */
        inline bool isNegative(const Integer &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        
        /**
           \brief Tests the given `plll::arithmetic::Real` object for being strictly negative.
           
           \return Returns `true` if and only if the argument is strictly negative.
           \sa \ref arithctxts_preds
        */
        inline bool isNegative(const Real &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        
        /**
           \brief Tests the given `plll::arithmetic::Integer` object for being negative or zero.
           
           \return Returns `true` if and only if the argument is negative or zero.
           \sa \ref arithctxts_preds
        */
        inline bool isNonPositive(const Integer &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        
        /**
           \brief Tests the given `plll::arithmetic::Real` object for being negative or zero.
           
           \return Returns `true` if and only if the argument is negative or zero.
           \sa \ref arithctxts_preds
        */
        inline bool isNonPositive(const Real &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        ///@}
        
        /**@{
           \name Euclidean ring functions.
        */
        /**
           \brief Computes an Euclidean Division of `a` by `b`.
           
           \param q The quotient of `a` divided by `b`.
           \param r The remainder of `a` divided by `b`, i.e. `a` modulo `b`.
           \param a The first operand.
           \param b The second operand. Must be non-zero.
           
           Computes `q` and `r` such that `a == q * b + r` and that \f$0 \le |r| \le |b|\f$ and \f$b
           \cdot r \ge 0\f$.
         */
        inline void euclideanDivision(Integer & q, Integer & r, const Integer & a, const Integer & b);
        
        /**
           \brief Computes an Euclidean Division of `a` by `b`.
           
           \param q The quotient of `a` divided by `b`.
           \param r The remainder of `a` divided by `b`, i.e. `a` modulo `b`.
           \param a The first operand.
           \param b The second operand. Must be non-zero.
           
           Computes `q` and `r` such that `a == q * b + r` and that \f$0 \le r \le |b|\f$.
         */
        inline void euclideanDivisionPos(Integer & q, Integer & r, const Integer & a, const Integer & b);
        
        /**
           \brief Computes the non-negative Greatest Common Divisior `r` of `x` and `y`.
           
           \param r The result is stored in here.
           \param x The first operand.
           \param y The first operand.
           \sa \ref arithctxts_ints
         */
        inline void GCD(Integer & r, const Integer & x, const Integer & y);
        
        /**
           \brief Computes the non-negative extended Greatest Common Divisior `r` of `x` and `y`.
           
           \param r The result is stored in here.
           \param a The Bezout coefficient of `x`.
           \param b The Bezout coefficient of `y`.
           \param x The first operand.
           \param y The first operand.
           
           Afterwards, the variables `r`, `a` and `b` satisfy `r == a * x + b * y`.
           
           \sa \ref arithctxts_ints
         */
        inline void XGCD(Integer & r, Integer & a, Integer & b, const Integer & x, const Integer & y);
        
        /**
           \brief Computes the non-negative Least Common Multiple `r` of `x` and `y`.
           
           \param r The result is stored in here.
           \param x The first operand.
           \param y The first operand.
           \sa \ref arithctxts_ints
         */
        inline void LCM(Integer & r, const Integer & x, const Integer & y);
        ///@}
        
        /**@{
           \name Stream input/output.
        */
        /**
           \brief Outputs the integer on the given output stream.
         */
        std::ostream & operator << (std::ostream &, const Integer &);
        /**
           \brief Reads the integer from the given input stream.
         */
        std::istream & operator >> (std::istream &, Integer &);
        /**
           \brief Outputs the floating point number on the given output stream.
         */
        std::ostream & operator << (std::ostream &, const Real &);
        /**
           \brief Reads the floating point number from the given input stream.
         */
        std::istream & operator >> (std::istream &, Real &);
        ///@}
        
        /**@{
           \name Setting to specific constants.
        */
        /**
           \brief Sets the given floating point number to `Not a Number`.
         */
        inline void setNaN(Real &);
        /**
           \brief Sets the given floating point number to \f$\pm \infty\f$.
           
           \param r The floating point variable to set to infinity.
           \param sign If `true`, the number is set to \f$+\infty\f$, and otherwise to
                       \f$-\infty\f$. The default value is `true`.
         */
        inline void setInfinity(Real & r, bool sign = true);
        /**
           \brief Sets the given integer to zero.
         */
        inline void setZero(Integer &);
        /**
           \brief Sets the given floating point number to \f$\pm 0\f$.
           
           \param r The floating point variable to set to zero.
           \param sign If `true`, the number is set to \f$+0\f$, and otherwise to \f$-0\f$. The
                       default value is `true`.
         */
        inline void setZero(Real & r, bool sign = true);
        /**
           \brief Sets the given integer to one.
         */
        inline void setOne(Integer &);
        /**
           \brief Sets the given floating point number to one.
         */
        inline void setOne(Real &);
        ///@}
        
        /**@{
           \name Comparisons.
        */
        /**
           \brief Compares the two integers.

           \param a The first operand.
           \param b The second operand.
           \return Returns a negative number if \f$a < b\f$, zero if \f$a = b\f$ and a positive
                   number if \f$a > b\f$.
         */
        inline int compare(const Integer & a, const Integer & b);
        /**
           \brief Compares the two floating point numbers.

           \param a The first operand.
           \param b The second operand.
           \return Returns a negative number if \f$a < b\f$, zero if \f$a = b\f$ and a positive
                   number if \f$a > b\f$.
         */
        inline int compare(const Real & a, const Real & b);
        /**
           \brief Compares the two integers in absolute value.

           \param a The first operand.
           \param b The second operand.
           \return Returns a negative number if \f$|a| < |b|\f$, zero if \f$|a| = |b|\f$ and a
                   positive number if \f$|a| > |b|\f$.
         */
        inline int compareAbsValues(const Integer & a, const Integer & b);
        /**
           \brief Compares the two floating point numbers in absolute value.

           \param a The first operand.
           \param b The second operand.
           \return Returns a negative number if \f$|a| < |b|\f$, zero if \f$|a| = |b|\f$ and a
                   positive number if \f$|a| > |b|\f$.
         */
        inline int compareAbsValues(const Real & a, const Real & b);
        ///@}
        
        /**@{
           \name Sign querying/modification.
        */
        /**
           \brief Returns the sign of the given integer.
           
           \return Returns a negative value if the value is negative, 0 if it is zero, and a
                   positive value if it is positive.
         */
        inline int sign(const Integer &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Returns the sign of the given floating point number.
           
           \return Returns a negative value if the value is negative, 0 if it is zero, and a
                   positive value if it is positive.
           \sa \ref arithctxts_preds
         */
        inline int sign(const Real &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        ///@}
        /**@{
           \name Bit manipulation.
        */
        /**
           \brief Returns the `n` bit of \f$|x|\f$ in the usual binary representation.
           \param x The integer whose bit to query.
           \param n The index of the bit to query.
         */
        inline int bit(const Integer & x, long n) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Sets the `n`-th bit of \f$|x|\f$ to `value`.
           \param x The integer whose bits to modify.
           \param n The index of the bit to set or clear.
           \param value The new value of the `n`-th bit. The default value is `true`.
           \sa \ref arithctxts_ints
         */
        inline void setbit(Integer & x, long n, bool value = true);
        ///@}

        /**@{
           \name Functional versions of operators.
        */
        /**
           \brief Increments `a` by one and stores the result in `r`.
           
           \param r The result.
           \param a The operand.
           \sa \ref arithctxts_ints
         */
        inline void increment(Integer & r, const Integer & a);
        /**
           \brief Decrements `a` by one and stores the result in `r`.
           
           \param r The result.
           \param a The operand.
         */
        inline void decrement(Integer & r, const Integer & a);
        /**
           \brief Adds `a` and `b` and stores the result in `r`.
           
           \param r The result.
           \param a The first operand.
           \param b The second operand.
         */
        inline void add(Integer & r, const Integer & a, const Integer & b);
        /**
           \brief Subtracts `b` from `a` and stores the result in `r`.
           
           \param r The result.
           \param a The first operand.
           \param b The second operand.
         */
        inline void sub(Integer & r, const Integer & a, const Integer & b);
        /**
           \brief Multiplies `a` with `b` and stores the result in `r`.
           
           \param r The result.
           \param a The first operand.
           \param b The second operand.
         */
        inline void mul(Integer & r, const Integer & a, const Integer & b);
        /**
           \brief Negates `a` and stores the result in `r`.
           
           \param r The result.
           \param a The operand.
         */
        inline void neg(Integer & r, const Integer & a);
        /**
           \brief Divides `a` by `b` and stores the result in `r`.
           
           \param r The result.
           \param a The first operand.
           \param b The second operand.
         */
        inline void div(Integer & r, const Integer & a, const Integer & b);
        /**
           \brief Takes the remainder of the division of `a` by `b` and stores it in `r`.
           
           \param r The result.
           \param a The first operand.
           \param b The second operand.
         */
        inline void mod(Integer & r, const Integer & a, const Integer & b);
        /**
           \brief Stores quotient and remainder of the division of `a` by `b` in `q` respectively
                  `r`.
           
           \param q The quotient.
           \param r The remainder.
           \param a The first operand.
           \param b The second operand.
         */
        inline void divmod(Integer & q, Integer & r, const Integer & a, const Integer & b);
        /**
           \brief Takes the absolute value of `a` and stores the result in `r`.
           
           \param r The result.
           \param a The operand.
         */
        inline void abs(Integer & r, const Integer & a);
        /**
           \brief Multiplies `a` and `b` and adds the result to `r`.
           
           \param r The accumulator.
           \param a The first operand.
           \param b The second operand.
         */
        inline void addmul(Integer & r, const Integer & a, const Integer & b);
        /**
           \brief Multiplies `a` and `b` and subtracts the result from `r`.
           
           \param r The accumulator.
           \param a The first operand.
           \param b The second operand.
         */
        inline void submul(Integer & r, const Integer & a, const Integer & b);
        /**
           \brief Computes the bitwise and of `a` and `b` and stores the result in `r`.
           
           \param r The result.
           \param a The first operand.
           \param b The second operand.
         */
        inline void band(Integer & r, const Integer & a, const Integer & b);
        /**
           \brief Computes the bitwise or of `a` and `b` and stores the result in `r`.
           
           \param r The result.
           \param a The first operand.
           \param b The second operand.
         */
        inline void bor(Integer & r, const Integer & a, const Integer & b);
        /**
           \brief Computes the bitwise exclusive or of `a` and `b` and stores the result in `r`.
           
           \param r The result.
           \param a The first operand.
           \param b The second operand.
         */
        inline void bxor(Integer & r, const Integer & a, const Integer & b);
        /**
           \brief Takes the bitwise complement of `a` and stores the result in `r`.
           
           \param r The result.
           \param a The operand.
         */
        inline void bneg(Integer & r, const Integer & a);
        /**
           \brief Shifts `a` by `b` bits to the left and stores the result in `r`.
           
           \param r The result.
           \param a The first operand.
           \param b The second operand.
         */
        inline void shl(Integer & r, const Integer & a, long b);
        /**
           \brief Shifts `a` by `b` bits to the right and stores the result in `r`.
           
           \param r The result.
           \param a The first operand.
           \param b The second operand.
         */
        inline void shl(Integer & r, const Integer & a, const Integer & b);
        /**
           \brief Shifts `a` by `b` bits to the left and stores the result in `r`.
           
           \param r The result.
           \param a The first operand.
           \param b The second operand.
         */
        inline void shr(Integer & r, const Integer & a, long b);
        /**
           \brief Shifts `a` by `b` bits to the right and stores the result in `r`.
           
           \param r The result.
           \param a The first operand.
           \param b The second operand.
         */
        inline void shr(Integer & r, const Integer & a, const Integer & b);
        /**
           \brief Computes the square of `a` and stores the result in `r`.
           
           \param r The result.
           \param a The operand.
         */
        inline void square(Integer & r, const Integer & a);
    
        /**
           \brief Adds `a` and `b` and stores the result in `r`.
           
           \param r The result.
           \param a The first operand.
           \param b The second operand.
         */
        inline void add(Real & r, const Real & a, const Real & b);
        /**
           \brief Subtracts `b` from `a` and stores the result in `r`.
           
           \param r The result.
           \param a The first operand.
           \param b The second operand.
         */
        inline void sub(Real & r, const Real & a, const Real & b);
        /**
           \brief Multiplies `a` with `b` and stores the result in `r`.
           
           \param r The result.
           \param a The first operand.
           \param b The second operand.
         */
        inline void mul(Real & r, const Real & a, const Real & b);
        /**
           \brief Divides `a` by `b` and stores the result in `r`.
           
           \param r The result.
           \param a The first operand.
           \param b The second operand.
         */
        inline void div(Real & r, const Real & a, const Real & b);
        /**
           \brief Takes the remainder of the division of `a` by `b` and stores it in `r`.
           
           \param r The result.
           \param a The first operand.
           \param b The second operand.
         */
        inline void mod(Real & r, const Real & a, const Real & b);
        /**
           \brief Stores quotient and remainder of the division of `a` by `b` in `q` respectively
                  `r`.
           
           \param q The quotient.
           \param r The remainder.
           \param a The first operand.
           \param b The second operand.
         */
        inline void divmod(Real & q, Real & r, const Real & a, const Real & b);
        /**
           \brief Multiplies `a` by \f$2^b\f$ and stores the result in `r`.
           
           \param r The result.
           \param a The first operand.
           \param b The second operand.
         */
        inline void shl(Real & r, const Real & a, const Real & b);
        /**
           \brief Multiplies `a` by \f$2^{-b}\f$ and stores the result in `r`.
           
           \param r The result.
           \param a The first operand.
           \param b The second operand.
           
           \sa \ref arithctxts_sqrtfp
         */
        inline void shr(Real & r, const Real & a, const Real & b);
        /**
           \brief Multiplies `a` by \f$2^b\f$ and stores the result in `r`.
           
           \param r The result.
           \param a The first operand.
           \param b The second operand.
           
           \sa \ref arithctxts_sqrtfp
         */
        inline void shl(Real & r, const Real & a, long b);
        /**
           \brief Multiplies `a` by \f$2^{-b}\f$ and stores the result in `r`.
           
           \param r The result.
           \param a The first operand.
           \param b The second operand.
         */
        inline void shr(Real & r, const Real & a, long b);
        /**
           \brief Increments `a` by one and stores the result in `r`.
           
           \param r The result.
           \param a The operand.
         */
        inline void increment(Real & r, const Real & a);
        /**
           \brief Decrements `a` by one and stores the result in `r`.
           
           \param r The result.
           \param a The operand.
         */
        inline void decrement(Real & r, const Real & a);
        /**
           \brief Negates `a` and stores the result in `r`.
           
           \param r The result.
           \param a The operand.
         */
        inline void neg(Real & r, const Real & a);
        /**
           \brief Takes the absolute value of `a` and stores the result in `r`.
           
           \param r The result.
           \param a The operand.
         */
        inline void abs(Real & r, const Real & a);
        /**
           \brief Multiplies `a` and `b` and adds the result to `r`.
           
           \param r The accumulator.
           \param a The first operand.
           \param b The second operand.
         */
        inline void addmul(Real & r, const Real & a, const Real & b);
        /**
           \brief Multiplies `a` and `b` and subtracts the result from `r`.
           
           \param r The accumulator.
           \param a The first operand.
           \param b The second operand.
         */
        inline void submul(Real & r, const Real & a, const Real & b);
        /**
           \brief Computes the square of `a` and stores the result in `r`.
           
           \param r The result.
           \param a The operand.
         */
        inline void square(Real & r, const Real & a);
        ///@}
        /**@{
           \name Sign querying/modification.
        */
        /**
           \brief Makes the operand non-negative.
           
           \param a The operand.
         */
        inline void makeAbs(Integer & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Makes the operand non-negative.
           
           \param a The operand.
         */
        inline void makeAbs(Real & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        ///@}
        
        /**@{
           \name Exponentiation.
        */
        /**
           \brief Raises `a` to the power `b` and stores the result in `r`.
           
           \param r The result.
           \param a The base.
           \param b The exponent.
         */
        inline void power(Integer & r, const Integer & a, long b);
        /**
           \brief Raises `a` to the power `b` and stores the result in `r`.
           
           \param r The result.
           \param a The base.
           \param b The exponent.
         */
        inline void power(Integer & r, const Integer & a, const Integer & b);
        ///@}
        
        /**@{
           \name Integer approximation.
        */
        /**
           \brief Computes \f$\lceil\sqrt{a}\rceil\f$ and stores the result in `r`.
           
           \param r The result.
           \param a The operand.
         */
        inline void sqrtCeil(Integer & r, const Integer & a);
        /**
           \brief Computes \f$\lfloor\sqrt{a}\rfloor\f$ and stores the result in `r`.
           
           \param r The result.
           \param a The operand.
         */
        inline void sqrtFloor(Integer & r, const Integer & a);
        
        /**
           \brief Quickly approximates \f$\log_2 |x|\f$ and returns the approximation.
           
           \param x A non-zero integer.
           \return \f$\approx \log_2 |x|\f$.
         */
        long approxLog2(const Integer & x) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Computes and returns \f$\lceil \log_2 |x| \rceil\f$.
           
           \param x A non-zero integer.
           \return \f$n \in \mathbb{N}\f$ such that `\f$2^{n-1} < |x| \le 2^n\f$.
           \sa \ref arithctxts_ints
         */
        long ceilOfLog2(const Integer & x) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Computes and returns \f$\lfloor \log_2 |x| \rfloor\f$.
           
           \param x A non-zero integer.
           \return \f$n \in \mathbb{N}\f$ such that `\f$2^n \le |x| < 2^{n+1}\f$.
           \sa \ref arithctxts_ints
         */
        long floorOfLog2(const Integer & x) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Computes and returns `n` such that \f$2^{n-1} \le |x| < 2^n\f$.
           
           \param x A non-zero integer.
           \return \f$n \in \mathbb{N}\f$ such that `\f$2^{n-1} \le |x| < 2^n\f$.
           \sa \ref arithctxts_ints
         */
        inline long bitLength(const Integer & x) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        
        /**
           \brief Computes \f$\lfloor \tfrac{a}{b} \rfloor\f$ and stores the result in `r`.
           
           \param r The result.
           \param a The divident.
           \param b The divisor.
           \sa \ref arithctxts_ints
         */
        void floorDiv(Integer & r, const Integer & a, const Integer & b);
        /**
           \brief Computes \f$\lceil \tfrac{a}{b} \rceil\f$ and stores the result in `r`.
           
           \param r The result.
           \param a The divident.
           \param b The divisor.
           \sa \ref arithctxts_ints
         */
        void ceilDiv(Integer & r, const Integer & a, const Integer & b);
        /**
           \brief Computes \f$\lfloor \tfrac{a}{b} \rceil\f$ (rounding to the next integer) and
                  stores the result in `r`.
           
           \param r The result.
           \param a The divident.
           \param b The divisor.
           \sa \ref arithctxts_ints
         */
        void roundDiv(Integer & r, const Integer & a, const Integer & b);
        ///@}
        
        /**@{
           \name Trigonometric functions.
        */
        /**
           \brief Computes the sine of `a` and stores the result in `res`.
           
           \param res The result.
           \param a The value the function is evaluated at.
           \sa \ref arithctxts_trig
         */
        inline void sin(Real & res, const Real & a);
        /**
           \brief Computes the cosine of `a` and stores the result in `res`.
           
           \param res The result.
           \param a The value the function is evaluated at.
           \sa \ref arithctxts_trig
         */
        inline void cos(Real & res, const Real & a);
        /**
           \brief Computes the tangent of `a` and stores the result in `res`.
           
           \param res The result.
           \param a The value the function is evaluated at.
           \sa \ref arithctxts_trig
         */
        inline void tan(Real & res, const Real & a);
        /**
           \brief Computes the arcsine of `a` and stores the result in `res`.
           
           \param res The result.
           \param a The value the function is evaluated at.
           \sa \ref arithctxts_trig
         */
        inline void asin(Real & res, const Real & a);
        /**
           \brief Computes the arccosine of `a` and stores the result in `res`.
           
           \param res The result.
           \param a The value the function is evaluated at.
           \sa \ref arithctxts_trig
         */
        inline void acos(Real & res, const Real & a);
        /**
           \brief Computes the arctangent of `a` and stores the result in `res`.
           
           \param res The result.
           \param a The value the function is evaluated at.
           \sa \ref arithctxts_trig
         */
        inline void atan(Real & res, const Real & a);
        /**
           \brief Computes the arctangent of \f$\tfrac{y}{x}\f$` and stores the result in `res`. The
                  signs of `x` and `y` are used to determine the quadrant and yield a result in
                  \f$[-\pi, \pi]\f$.
           
           \param res The result.
           \param y The numerator of the fraction the function is evaluated at.
           \param x The denominator of the fraction the function is evaluated at.
           \sa \ref arithctxts_trig
         */
        inline void atan2(Real & res, const Real & y, const Real & x);
        ///@}
        
        /**@{
           \name Exponential and logarithmic functions.
        */
        /**
           \brief Computes the exponential function at `a` and stores the result in `res`.
           
           \param res The result.
           \param a The value the function is evaluated at.
           \sa \ref arithctxts_trig
         */
        inline void exp(Real & res, const Real & a);
        /**
           \brief Computes the natural logarithm of `a` and stores the result in `res`.
           
           \param res The result.
           \param a The value the function is evaluated at.
         */
        inline void log(Real & res, const Real & a);
        /**
           \brief Computes the logarithm of `a` to base 2 and stores the result in `res`.
           
           \param res The result.
           \param a The value the function is evaluated at.
         */
        inline void log2(Real & res, const Real & a);
        /**
           \brief Computes the logarithm of `a` to base 10 and stores the result in `res`.
           
           \param res The result.
           \param a The value the function is evaluated at.
         */
        inline void log10(Real & res, const Real & a);
        /**
           \brief Computes the square root of `a` and stores the result in `res`.
           
           \param res The result.
           \param a The value the function is evaluated at.
         */
        inline void sqrt(Real & res, const Real & a);
        ///@}
        
        /**@{
           \name Special functions.
        */
        /**
           \brief Computes the Gamma function at `a` and stores the result in `res`.
           
           \param res The result.
           \param a The value the function is evaluated at.
           \sa \ref arithctxts_sqrtfp
         */
        inline void gamma(Real & res, const Real & a);
        /**
           \brief Computes the logarithm of the absolute value of the Gamma function at `a` and
                  stores the result in `res`.
           
           \param res The result.
           \param a The value the function is evaluated at.
           \sa \ref arithctxts_specfns
         */
        inline void lgamma(Real & res, const Real & a);
        /**
           \brief Computes the logarithm of the absolute value of the Gamma function at `a` and
                  stores the result in `res`. The sign of the Gamma function at `a` is stored in
                  `sign`.

           \param res The result.
           \param sign Will be set to 1 if \f$\Gamma(a) \ge 0\f$ and -1 otherwise.
           \sa \ref arithctxts_specfns
           \param a The value the function is evaluated at.
         */
        inline void lgamma(Real & res, int & sign, const Real & a);
        /**
           \brief Raises `a` to the power `b` and stores the result in `res`.
           
           \param res The result.
           \param a The base.
           \param b The exponent.
           \sa \ref arithctxts_specfns
         */
        ///@}
        
        /**@{
           \name Exponentiation.
        */
        inline void power(Real & res, const Real & a, long b);
        /**
           \brief Raises `a` to the power `b` and stores the result in `res`.
           
           \param res The result.
           \param a The base.
           \param b The exponent.
         */
        inline void power(Real & res, const Real & a, const Integer & b);
        /**
           \brief Raises `a` to the power `b` and stores the result in `res`.
           
           \param res The result.
           \param a The base.
           \param b The exponent.
           \sa \ref arithctxts_sqrtfp
         */
        inline void power(Real & res, const Real & a, const Real & b);
        ///@}
    
        /**
           \brief Represents an arithmetic context for arbitrary precision integer.
           
           \sa \ref arithctxts_contexts.
         */
        class IntegerContext
        {
            friend class arithmetic::Integer;
        
        public:
            /**
               \brief The integer type.
             */
            typedef arithmetic::Integer Integer;
            /**
               \brief The integer type.
             */
            typedef arithmetic::Integer Type;
            
            /**
               \brief The properties of arbitrary precision integers.
             */
            enum { is_cputype = false, is_realtype = false, is_inttype = true, is_exact = true,
                   is_modulo = false, has_infinity = false, has_uniform_rng = true };
            
            /**
               \brief A uniform random number generator frontend.
               
               \sa \ref arithctxts_urng
             */
            class UniformRNG;
        };
        
        /**
           \brief Represents an arbitrary precision integer.
         */
        class Integer
        {
            friend class Real;
            friend class RandomNumberGenerator;
        
#ifndef PLLL_INTERNAL_NO_TEMPLATE_FRIENDS
            template<class X, class Y>
            friend class implementation::conversion_impl;
            
            friend class implementation::nativeconversion_impl<Integer>;
            friend class implementation::from_string_conversion<IntegerContext>;
            friend class implementation::to_string_conversion<Integer>;
            
        private:
#else
        public:
#endif
            mpz_t d_value;
        
            static void mpz_init_set_ld(mpz_t &, long double);
            static void mpz_set_ld(mpz_t &, long double);
            static long double mpz_get_ld(const mpz_t &);
            static void mpz_set_ll(mpz_t &, long long);
            static long long mpz_get_ll(const mpz_t &);

        public:
            /**
               \brief The context type.
             */
            typedef IntegerContext Context;
            
            /**
               \brief Creates a new integer. Default value is zero.
             */
            inline Integer()
            {
                mpz_init(d_value);
            }
            
            /**
               \brief Creates a new integer. Default value is zero.
               
               \param c An integer context.
             */
            inline explicit Integer(const IntegerContext & c)
            {
                mpz_init(d_value);
            }
            
            /**
               \brief Creates a copy of the given integer.
               
               \param i The integer to be copied.
             */
            inline Integer(const Integer & i)
            {
                mpz_init_set(d_value, i.d_value);
            }
            
            /**
               \brief Creates a copy of the given integer.
               
               \param i The integer to be copied.
               \param ic An integer context.
             */
            inline Integer(const Integer & i, const IntegerContext & ic)
            {
                mpz_init_set(d_value, i.d_value);
            }
            
            /**
               \brief Creates a arbitrary precision integer from the given native integer.
               
               \param i The native integer.
             */
            inline explicit Integer(signed int i)
            {
                mpz_init_set_si(d_value, i);
            }
            
            /**
               \brief Creates a arbitrary precision integer from the given native integer.
               
               \param i The native integer.
             */
            inline explicit Integer(unsigned int i)
            {
                mpz_init_set_ui(d_value, i);
            }
            
            /**
               \brief Creates a arbitrary precision integer from the given native integer.
               
               \param i The native integer.
             */
            inline explicit Integer(signed long i)
            {
                mpz_init_set_si(d_value, i);
            }
            
            /**
               \brief Creates a arbitrary precision integer from the given native integer.
               
               \param i The native integer.
             */
            inline explicit Integer(unsigned long i)
            {
                mpz_init_set_ui(d_value, i);
            }
            
            /**
               \brief Creates a arbitrary precision integer from the given native integer.
               
               \param i The native integer.
             */
            inline explicit Integer(long long i)
            {
                mpz_init(d_value);
                mpz_set_ll(d_value, i);
            }
            
            /**
               \brief Creates a arbitrary precision integer from the given native floating point
                      number.
               
               \param d The native floating point number.
             */
            inline explicit Integer(double d)
            {
                mpz_init_set_d(d_value, d);
            }
            
            /**
               \brief Creates a arbitrary precision integer from the given native floating point
                      number.
               
               \param d The native floating point number.
             */
            inline explicit Integer(long double d)
            {
                mpz_init_set_ld(d_value, d);
            }
            
            /**
               \brief Creates an integer from the given integer expression.
               
               \param E An integer expression.
             */
            template<class A, template<typename, typename> class O>
            inline Integer(const expressions::Expression<IntegerContext, A, O> & E)
            {
                mpz_init(d_value);
                E.assignTo(*this);
            }
            
            /**
               \brief Creates an integer from the given integer expression.
               
               \param E An integer expression.
               \param ic An integer context.
             */
            template<class A, template<typename, typename> class O>
            inline Integer(const expressions::Expression<IntegerContext, A, O> & E, const IntegerContext & ic)
            {
                mpz_init(d_value);
                E.assignTo(*this);
            }
            
            /**
               \brief Releases the memory used by the arbitrary precision integer.
             */
            inline ~Integer() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
#if __cplusplus >= 201103L
                if (d_value[0]._mp_d != NULL)
#endif
                    mpz_clear(d_value);
            }
            
#if __cplusplus >= 201103L
            /**
               \brief Creates a new integer and moves the given integer into this one.
               
               \param r A integer to be moved into this one.
             */
            inline Integer(Integer && i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_value{i.d_value[0]}
            {
                i.d_value[0]._mp_d = NULL;
            }
            
            /**
               \brief Moves the given integer into this one.
               
               \param r A integer to be moved into this one.
             */
            inline Integer & operator = (Integer && i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                mpz_clear(d_value);
                d_value[0] = i.d_value[0];
                i.d_value[0]._mp_d = NULL;
                return *this;
            }
#endif
            
            /**
               \brief Creates an integer from the given arbitrary precision floating point number.
               
               \param r An arbitrary precision floating point number.
             */
            explicit Integer(const Real & r);
            
            /**
               \brief Sets the integer context `c`.
               
               \param c The integer context.
             */
            static inline void setContext(const IntegerContext & c)
            {
            }
            
            /**
               \brief Assigns the integer `i` to the current integer.
             */
            inline Integer & operator = (const Integer & i)
            {
                mpz_set(d_value, i.d_value);
                return *this;
            }
            
            // Assign generic expression
            template<class A, template<typename, typename> class O>
            inline Integer & operator = (const expressions::Expression<IntegerContext, A, O> & E)
            {
                E.assignTo(*this);
                return *this;
            }
            
            inline friend void increment(Integer & r, const Integer & a) // r = a + 1
            {
                mpz_add_ui(r.d_value, a.d_value, 1);
            }
            
            inline friend void decrement(Integer & r, const Integer & a) // r = a - 1
            {
                mpz_sub_ui(r.d_value, a.d_value, 1);
            }
            
            inline friend void add(Integer & r, const Integer & a, const Integer & b)
            {
                mpz_add(r.d_value, a.d_value, b.d_value);
            }
            
            inline friend void add_ui(Integer & r, const Integer & a, unsigned long b) // r = a + b
            // Internal for expression templates
            {
                mpz_add_ui(r.d_value, a.d_value, b);
            }
            
            inline friend void add_si(Integer & r, const Integer & a, signed long b) // r = a + b
            // Internal for expression templates
            {
                if (b < 0)
                    mpz_sub_ui(r.d_value, a.d_value, -b);
                else
                    mpz_add_ui(r.d_value, a.d_value, b);
            }
            
            inline friend void sub(Integer & r, const Integer & a, const Integer & b)
            {
                mpz_sub(r.d_value, a.d_value, b.d_value);
            }
            
            inline friend void sub_ui(Integer & r, const Integer & a, unsigned long b) // r = a - b
            // Internal for expression templates
            {
                mpz_sub_ui(r.d_value, a.d_value, b);
            }
            
            inline friend void sub_si(Integer & r, const Integer & a, signed long b) // r = a - b
            // Internal for expression templates
            {
                if (b < 0)
                    mpz_add_ui(r.d_value, a.d_value, -b);
                else
                    mpz_sub_ui(r.d_value, a.d_value, b);
            }
            
            inline friend void ui_sub(Integer & r, unsigned long a, const Integer & b) // r = a - b
            // Internal for expression templates
            {
                mpz_ui_sub(r.d_value, a, b.d_value);
            }
            
            inline friend void si_sub(Integer & r, signed long a, const Integer & b) // r = a - b
            // Internal for expression templates
            {
                if (a < 0)
                {
                    mpz_add_ui(r.d_value, b.d_value, -a);
                    mpz_neg(r.d_value, r.d_value);
                }
                else
                    mpz_ui_sub(r.d_value, a, b.d_value);
            }
            
            inline friend void mul(Integer & r, const Integer & a, const Integer & b)
            {
                mpz_mul(r.d_value, a.d_value, b.d_value);
            }
            
            inline friend void mul_ui(Integer & r, const Integer & a, unsigned long b) // r = a * b
            // Internal for expression templates
            {
                mpz_mul_ui(r.d_value, a.d_value, b);
            }
            
            inline friend void mul_si(Integer & r, const Integer & a, signed long b) // r = a * b
            // Internal for expression templates
            {
                mpz_mul_si(r.d_value, a.d_value, b);
            }
            
            inline friend void neg(Integer & r, const Integer & a)
            {
                mpz_neg(r.d_value, a.d_value);
            }
            
            inline friend void div(Integer & r, const Integer & a, const Integer & b)
            {
                mpz_tdiv_q(r.d_value, a.d_value, b.d_value);
            }
            
            inline friend void div_ui(Integer & r, const Integer & a, unsigned long b) // r = a / b
            // Internal for expression templates
            {
                mpz_tdiv_q_ui(r.d_value, a.d_value, b);
            }
            
            inline friend void div_si(Integer & r, const Integer & a, signed long b) // r = a / b
            // Internal for expression templates
            {
                if (b < 0)
                {
                    mpz_tdiv_q_ui(r.d_value, a.d_value, -b);
                    mpz_neg(r.d_value, r.d_value);
                }
                else
                    mpz_tdiv_q_ui(r.d_value, a.d_value, b);
            }
            
            inline friend void mod(Integer & r, const Integer & a, const Integer & b)
            {
                mpz_tdiv_r(r.d_value, a.d_value, b.d_value);
            }
            
            inline friend void mod_ui(Integer & r, const Integer & a, unsigned long b) // r = a % b
            // Internal for expression templates
            {
                mpz_tdiv_r_ui(r.d_value, a.d_value, b);
            }
            
            inline friend void mod_si(Integer & r, const Integer & a, signed long b) // r = a % b
            // Internal for expression templates
            {
                mpz_tdiv_r_ui(r.d_value, a.d_value, (b < 0) ? -b : b);
            }
            
            inline friend void divmod(Integer & q, Integer & r, const Integer & a, const Integer & b)
            {
                mpz_tdiv_qr(q.d_value, r.d_value, a.d_value, b.d_value);
            }
            
            inline friend void divmod_ui(Integer & q, Integer & r, const Integer & a, unsigned long b)
            // Internal for expression templates
            {
                mpz_tdiv_qr_ui(q.d_value, r.d_value, a.d_value, b);
            }
            
            inline friend void divmod_si(Integer & q, Integer & r, const Integer & a, signed long b)
            // Internal for expression templates
            {
                mpz_tdiv_qr_ui(q.d_value, r.d_value, a.d_value, (b < 0) ? -b : b);
                if (b < 0)
                    mpz_neg(q.d_value, q.d_value);
            }
            
            inline friend void abs(Integer & r, const Integer & a)
            {
                mpz_abs(r.d_value, a.d_value);
            }
            
            inline friend void addmul(Integer & r, const Integer & a, const Integer & b) // r += a * b
            {
                mpz_addmul(r.d_value, a.d_value, b.d_value);
            }
            
            inline friend void addmul_ui(Integer & r, const Integer & a, unsigned long b) // r += a * b
            // Internal for expression templates
            {
                mpz_addmul_ui(r.d_value, a.d_value, b);
            }
            
            inline friend void addmul_si(Integer & r, const Integer & a, signed long b) // r += a * b
            // Internal for expression templates
            {
                if (b < 0)
                    mpz_submul_ui(r.d_value, a.d_value, -b);
                else
                    mpz_addmul_ui(r.d_value, a.d_value, b);
            }
            
            inline friend void submul(Integer & r, const Integer & a, const Integer & b) // r -= a * b
            {
                mpz_submul(r.d_value, a.d_value, b.d_value);
            }
            
            inline friend void submul_ui(Integer & r, const Integer & a, unsigned long b) // r -= a * b
            // Internal for expression templates
            {
                mpz_submul_ui(r.d_value, a.d_value, b);
            }
            
            inline friend void submul_si(Integer & r, const Integer & a, signed long b) // r -= a * b
            // Internal for expression templates
            {
                if (b < 0)
                    mpz_addmul_ui(r.d_value, a.d_value, -b);
                else
                    mpz_submul_ui(r.d_value, a.d_value, b);
            }
            
            inline friend void power(Integer & r, const Integer & a, long b)
            {
                if (b < 0)
                    mpz_set_ui(r.d_value, 0);
                else
                    mpz_pow_ui(r.d_value, a.d_value, b);
            }
            
            inline friend void power(Integer & r, const Integer & a, const Integer & b)
            {
                if (sign(b) < 0)
                    mpz_set_ui(r.d_value, 0);
                else
                    mpz_pow_ui(r.d_value, a.d_value, mpz_get_ui(b.d_value));
            }
            
            inline friend void floorDiv(Integer & r, const Integer & a, const Integer & b) // Computes r = floor(a/b), where b \neq 0.
            {
                mpz_fdiv_q(r.d_value, a.d_value, b.d_value);
            }
            
            inline friend void ceilDiv(Integer & r, const Integer & a, const Integer & b) // Computes r = ceil(a/b), where b \neq 0.
            {
                mpz_cdiv_q(r.d_value, a.d_value, b.d_value);
            }
            
            inline friend void roundDiv(Integer & r, const Integer & a, const Integer & b) // Computes r = round(a/b), where b \neq 0.
            {
                mpz_t rem;
                mpz_init(rem);
                mpz_fdiv_qr(r.d_value, rem, a.d_value, b.d_value);
                mpz_mul_2exp(rem, rem, 1);
                if (mpz_cmp(rem, b.d_value) > 0)
                    mpz_add_ui(r.d_value, r.d_value, 1);
                mpz_clear(rem);
            }
            
            inline friend void sqrtCeil(Integer & r, const Integer & a)
            {
                mpz_t rr;
                mpz_init(rr);
                mpz_sqrtrem(r.d_value, rr, a.d_value);
                if (mpz_sgn(rr) != 0)
                    increment(r, r);
                mpz_clear(rr);
            }
            
            inline friend void sqrtFloor(Integer & r, const Integer & a)
            {
                mpz_sqrt(r.d_value, a.d_value);
            }
            
            inline friend void band(Integer & r, const Integer & a, const Integer & b)
            {
                mpz_and(r.d_value, a.d_value, b.d_value);
            }
            
            inline friend void bor(Integer & r, const Integer & a, const Integer & b)
            {
                mpz_ior(r.d_value, a.d_value, b.d_value);
            }
            
            inline friend void bxor(Integer & r, const Integer & a, const Integer & b)
            {
                mpz_xor(r.d_value, a.d_value, b.d_value);
            }
            
            inline friend void bneg(Integer & r, const Integer & a)
            {
                mpz_com(r.d_value, a.d_value);
            }
            
            inline friend void shl(Integer & r, const Integer & a, long b)
            {
                if (b < 0)
                    mpz_tdiv_q_2exp(r.d_value, a.d_value, -b);
                else
                    mpz_mul_2exp(r.d_value, a.d_value, b);
            }
            
            inline friend void shl(Integer & r, const Integer & a, const Integer & b)
            {
                mpz_mul_2exp(r.d_value, a.d_value, mpz_get_ui(b.d_value));
            }
            
            inline friend void shr(Integer & r, const Integer & a, long b)
            {
                if (b < 0)
                    mpz_mul_2exp(r.d_value, a.d_value, -b);
                else
                    mpz_tdiv_q_2exp(r.d_value, a.d_value, b);
            }
            
            inline friend void shr(Integer & r, const Integer & a, const Integer & b)
            {
                mpz_tdiv_q_2exp(r.d_value, a.d_value, mpz_get_ui(b.d_value));
            }
            
            inline friend bool isZero(const Integer & i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return mpz_sgn(i.d_value) == 0;
            }
            
            inline friend bool isOne(const Integer & i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return mpz_cmp_ui(i.d_value, 1) == 0;
            }
            
            inline friend bool isPMOne(const Integer & i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return mpz_cmpabs_ui(i.d_value, 1) == 0;
            }
            
            inline friend bool isPMTwo(const Integer & i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return mpz_cmpabs_ui(i.d_value, 2) == 0;
            }
            
            inline friend bool isPositive(const Integer & i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return mpz_sgn(i.d_value) > 0;
            }
            
            inline friend bool isNonNegative(const Integer & i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return mpz_sgn(i.d_value) >= 0;
            }
            
            inline friend bool isNegative(const Integer & i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return mpz_sgn(i.d_value) < 0;
            }
            
            inline friend bool isNonPositive(const Integer & i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return mpz_sgn(i.d_value) <= 0;
            }
            
            inline friend void makeAbs(Integer & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                abs(a, a);
            }
            
            inline friend void euclideanDivision(Integer & q, Integer & r, const Integer & a, const Integer & b)
            // Given b \neq 0, computes q, r such that a = q b + r, 0 <= |r| < |b| and b*r >= 0.
            {
                mpz_tdiv_qr(q.d_value, r.d_value, a.d_value, b.d_value);
            }
            
            inline friend void euclideanDivisionPos(Integer & q, Integer & r, const Integer & a, const Integer & b)
            // Given b \neq 0, computes q, r such that a = q b + r, 0 <= r < |b|
            {
                mpz_fdiv_qr(q.d_value, r.d_value, a.d_value, b.d_value);
            }
            
            inline friend void GCD(Integer & r, const Integer & x, const Integer & y)
            // Computes the gcd of x and y.
            {
                mpz_gcd(r.d_value, x.d_value, y.d_value);
            }
            
            inline friend void GCD_ui(Integer & r, const Integer & x, unsigned long y)
            // Internal for expression templates
            {
                mpz_gcd_ui(r.d_value, x.d_value, y);
            }
            
            inline friend void GCD_si(Integer & r, const Integer & x, signed long y)
            // Internal for expression templates
            {
                mpz_gcd_ui(r.d_value, x.d_value, y < 0 ? -y : y);
            }
            
            inline friend void XGCD(Integer & r, Integer & a, Integer & b, const Integer & x, const Integer & y)
            // Computes the gcd of x and y. Also computes a, b such that   gcd = a*x + b*y.
            {
                mpz_gcdext(r.d_value, a.d_value, b.d_value, x.d_value, y.d_value);
            }
            
            inline friend void LCM(Integer & r, const Integer & x, const Integer & y)
            // Computes the lcm of x and y.
            {
                mpz_lcm(r.d_value, x.d_value, y.d_value);
            }
            
            inline friend void LCM_ui(Integer & r, const Integer & x, unsigned long y)
            // Internal for expression templates
            {
                mpz_lcm_ui(r.d_value, x.d_value, y);
            }
            
            inline friend void LCM_si(Integer & r, const Integer & x, signed long y)
            // Internal for expression templates
            {
                mpz_lcm_ui(r.d_value, x.d_value, y < 0 ? -y : y);
            }
            
            friend std::ostream & operator << (std::ostream & s, const Integer & i);
            // Output to stream
            
            friend std::istream & operator >> (std::istream & s, Integer & i);
            // Input from stream
            
            inline friend void setZero(Integer & i)
            // Set to zero
            {
                mpz_set_ui(i.d_value, 0);
            }
            
            inline friend void setOne(Integer & i)
            // Set to one
            {
                mpz_set_ui(i.d_value, 1);
            }
            
            inline friend int compare(const Integer & a, const Integer & b)
            // Tests whether the first integer is < (-1), = (0) or > (1) than the second.
            {
                return mpz_cmp(a.d_value, b.d_value);
            }
            
            inline friend int compare_d(const Integer & a, double b)
            // Internal for expression templates
            {
                return mpz_cmp_d(a.d_value, b);
            }
            
            inline friend int compare_ui(const Integer & a, unsigned long b)
            // Internal for expression templates
            {
                return mpz_cmp_ui(a.d_value, b);
            }
            
            inline friend int compare_si(const Integer & a, signed long b)
            // Internal for expression templates
            {
                return mpz_cmp_si(a.d_value, b);
            }
            
            inline friend int compareAbsValues(const Integer & a, const Integer & b)
            // Tests whether the absolute value of the first integer is < (-1), = (0) or > (1) than the
            // absolute value of the second integer.
            {
                return mpz_cmpabs(a.d_value, b.d_value);
            }
            
            inline friend int compareAbsValues_d(const Integer & a, double b)
            // Internal for expression templates
            {
                return mpz_cmpabs_d(a.d_value, b);
            }
            
            inline friend int compareAbsValues_si(const Integer & a, signed long b)
            // Internal for expression templates
            {
                return mpz_cmpabs_ui(a.d_value, b < 0 ? -b : b);
            }
            
            inline friend int compareAbsValues_ui(const Integer & a, unsigned long b)
            // Internal for expression templates
            {
                return mpz_cmpabs_ui(a.d_value, b);
            }
            
            inline friend int sign(const Integer & i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            // Returns the sign (i.e. -1, 0, 1).
            {
                return mpz_sgn(i.d_value);
            }
            
            inline friend int bit(const Integer & x, long n) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            // Returns the n-th bit of the given integer.
            {
                return mpz_tstbit(x.d_value, n);
            }
            
            inline friend void setbit(Integer & x, long n, bool value)
            // Sets the n-th bit of the given integer.
            {
                if (value)
                    mpz_setbit(x.d_value, n);
                else
                    mpz_clrbit(x.d_value, n);
            }
            
            inline friend void square(Integer & r, const Integer & a)
            {
                mpz_mul(r.d_value, a.d_value, a.d_value);
            }
            
            inline friend long approxLog2(const Integer & x) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE // Returns   ~log_2(|x|).
            {
                return mpz_sizeinbase(x.d_value, 2);
            }
            
            inline friend long ceilOfLog2(const Integer & x) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE // Returns   ceil(log_2(|x|)).
            {
                mp_bitcnt_t r = mpz_sizeinbase(x.d_value, 2);
                mp_bitcnt_t lsb = mpz_scan1(x.d_value, 0);
                return (lsb < r - 1) ? r : r - 1;
            }
            
            inline friend long floorOfLog2(const Integer & x) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE // Returns   floor(log_2(|x|)).
            {
                return (long)mpz_sizeinbase(x.d_value, 2) - 1;
            }
            
            inline friend long bitLength(const Integer & x) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE // Returns n such that 2^{n-1} <= |x| < 2^n
            {
                return (long)mpz_sizeinbase(x.d_value, 2);
            }
            
            inline friend void arithmetic::swap(plll::arithmetic::Integer &, plll::arithmetic::Integer &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            
            inline friend void power(Real & res, const Real & a, const Integer & b);
            
            // These ones should only be used in very, very special and rare cases!
            inline const mpz_t & getInternal() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE { return d_value; }
            inline mpz_t & getInternal() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE { return d_value; }
            
            inline friend bool operator == (const Integer & a, const Integer & b);
            inline friend bool operator != (const Integer & a, const Integer & b);
            inline friend bool operator <= (const Integer & a, const Integer & b);
            inline friend bool operator >= (const Integer & a, const Integer & b);
            inline friend bool operator < (const Integer & a, const Integer & b);
            inline friend bool operator > (const Integer & a, const Integer & b);
        };
        
        // Internal for expression templates
        inline void add_ui(Integer & r, const Integer & a, unsigned long b); // r = a + b
        inline void add_si(Integer & r, const Integer & a, signed long b); // r = a + b
        inline void sub_ui(Integer & r, const Integer & a, unsigned long b); // r = a - b
        inline void sub_si(Integer & r, const Integer & a, signed long b); // r = a - b
        inline void ui_sub(Integer & r, unsigned long a, const Integer & b); // r = a - b
        inline void si_sub(Integer & r, signed long a, const Integer & b); // r = a - b
        inline void mul_ui(Integer & r, const Integer & a, unsigned long b); // r = a * b
        inline void mul_si(Integer & r, const Integer & a, signed long b); // r = a * b
        inline void div_ui(Integer & r, const Integer & a, unsigned long b); // r = a / b
        inline void div_si(Integer & r, const Integer & a, signed long b); // r = a / b
        inline void mod_ui(Integer & r, const Integer & a, unsigned long b); // r = a % b
        inline void mod_si(Integer & r, const Integer & a, signed long b); // r = a % b
        inline void addmul_ui(Integer & r, const Integer & a, unsigned long b); // r += a * b
        inline void addmul_si(Integer & r, const Integer & a, signed long b); // r += a * b
        inline void submul_ui(Integer & r, const Integer & a, unsigned long b); // r -= a * b
        inline void submul_si(Integer & r, const Integer & a, signed long b); // r -= a * b
        inline void GCD_ui(Integer & r, const Integer & x, unsigned long y);
        inline void GCD_si(Integer & r, const Integer & x, signed long y);
        inline void LCM_ui(Integer & r, const Integer & x, unsigned long y);
        inline void LCM_si(Integer & r, const Integer & x, signed long y);
        inline int compare_d(const Integer & a, double b);
        inline int compare_ui(const Integer & a, unsigned long b);
        inline int compare_si(const Integer & a, signed long b);
        inline int compareAbsValues_d(const Integer & a, double b);
        inline int compareAbsValues_si(const Integer & a, signed long b);
        inline int compareAbsValues_ui(const Integer & a, unsigned long b);
    }
}

#include "arithmetic-gmp-iops.hpp"

namespace plll
{
    namespace arithmetic
    {
        /**
           \brief Represents an arithmetic context for arbitrary precision floating point values.
           
           \sa \ref arithctxts_contexts.
         */
        class RealContext
        {
            friend class arithmetic::Real;
            friend RealContext & getThreadRealContext();
        
            // stores information like precision (for Real)
        
            long d_prec;
            bool d_global;
        
        private:
            inline explicit RealContext(bool global)
                : d_prec(mpfr_get_default_prec()), d_global(true)
            {
            }
        
        public:
            /**
               \brief The floating point type.
             */
            typedef arithmetic::Real Real;
            /**
               \brief The floating point type.
             */
            typedef arithmetic::Real Type;
        
            /**
               \brief The properties of arbitrary precision floating point numbers.
             */
            enum { is_cputype = false, is_realtype = true, is_inttype = false, is_exact = false, is_variable_precision = true,
                   has_squareroot = true, has_full_power = true, has_special_fns = true, has_huge_exponent = true,
                   has_infinity = true, has_uniform_rng = true, has_constants = true, has_trigonometric = true };
            
            /**
               \brief Creates a new real context with the current global default precision.
             */
            inline RealContext() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_prec(mpfr_get_default_prec()), d_global(false)
            {
            }
            
            /**
               \brief Creates a new real context with the given precision.
             */
            inline explicit RealContext(long prec) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_prec(prec < MPFR_PREC_MIN ? MPFR_PREC_MIN : (prec > MPFR_PREC_MAX ? MPFR_PREC_MAX : prec)), d_global(false)
            {
            }
            
            /**
               \brief Sets the precision of this context.
               
               If the precision is too small or too large, it will be clipped to the permissive
               interval bounds.
               
               \param prec The new precision.
             */
            inline void setRealPrecision(unsigned long prec) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                if (prec < MPFR_PREC_MIN) prec = MPFR_PREC_MIN;
                if (prec > MPFR_PREC_MAX) prec = MPFR_PREC_MAX;
                if (d_global)
                    mpfr_set_default_prec(prec);
                else
                    d_prec = prec;
            }
            
            /**
               \brief Returns the precision of this context.
             */
            inline unsigned long getRealPrecision() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return d_global ? mpfr_get_default_prec() : d_prec;
            }
            
            /**
               \brief Returns the minimal possible precision for this context.
             */
            static inline unsigned long getMinRealPrecision() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE // returns minimal value for precision
            {
                return MPFR_PREC_MIN;
            }
        
            /**
               \brief Returns the maximal possible precision for this context.
             */
            static inline unsigned long getMaxRealPrecision() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE // returns maximal value for precision
            {
                return MPFR_PREC_MAX;
            }
            
            /**
               \brief Returns the machine constant for the current context.
               
               \sa \ref arithctxts_consts
             */
            inline Real getEpsilon() const;
            /**
               \brief Sets the argument to the machine constant for the current context.
               
               \sa \ref arithctxts_consts
             */
            inline void getEpsilon(Real &) const;
        
            /**
               \brief Returns an approximation of \f$\pi\f$ for the current context.
               
               \sa \ref arithctxts_consts
             */
            inline Real getPi() const;
            /**
               \brief Sets the argument to an approximation of \f$\pi\f$ for the current context.
               
               \sa \ref arithctxts_consts
             */
            inline void getPi(Real &) const;
        
            /**
               \brief Returns an approximation of the Euler number \f$\exp(1)\f$ for the current
                      context.
               
               \sa \ref arithctxts_consts
             */
            inline Real getEuler() const;
            /**
               \brief Sets the argument to an approximation of the Euler number \f$\exp(1)\f$ for
                      the current context.
               
               \sa \ref arithctxts_consts
             */
            inline void getEuler(Real &) const;
        
            /**
               \brief Returns an approximation of the natural logarithm \f$\log 2\f$ for the current
                      context.
               
               \sa \ref arithctxts_consts
             */
            inline Real getLog2() const;
            /**
               \brief Sets the argument to an approximation of the natural logarithm \f$\log 2\f$
                      for the current context.
               
               \sa \ref arithctxts_consts
             */
            inline void getLog2(Real &) const;
        
            /**
               \brief A uniform random number generator frontend.
               
               \sa \ref arithctxts_urng
             */
            class UniformRNG;
        };
        
        /**
           \brief Retrieves a context for the current thread. The context is thread local and cannot
                  be accessed by this method from other threads.
         */
        RealContext & getThreadRealContext(); // returns context for current thread
        
        /**
           \brief Represents an arbitrary precision floating point value.
         */
        class Real
        {
            friend class RandomNumberGenerator;
            friend class RealContext;
            
#ifndef PLLL_INTERNAL_NO_TEMPLATE_FRIENDS
            template<class X, class Y>
            friend class implementation::conversion_impl;
            
            friend class implementation::nativeconversion_impl<Real>;
            friend class implementation::from_string_conversion<RealContext>;
            friend class implementation::to_string_conversion<Real>;
            
        private:
#else
        public:
#endif
            mpfr_t d_value;
            
            inline Real(bool, unsigned long p)
            {
                mpfr_init2(d_value, p);
            }
            
            static void mpfr_set_ll(mpfr_t &, long long, mpfr_rnd_t);
            static long long mpfr_get_ll(const mpfr_t &, mpfr_rnd_t);
            static long long mpfr_get_ll(const mpfr_t &, mpfr_rnd_t, bool & roundUp);
            
        public:
            // Internal helper for the expression template classes (we don't want to friend them
            // all explicitly!).
            struct PrecisionInit
            {
                unsigned long precision;
                
                inline PrecisionInit(unsigned long prec) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    : precision(prec)
                {
                }
            };
            
            inline Real(PrecisionInit p)
            {
                mpfr_init2(d_value, p.precision);
            }
            
        public:
            /**
               \brief The context type.
             */
            typedef RealContext Context;
            
            /**
               \brief Creates a new floating point number.
               
               \warning Default value is `NaN` and not zero!
             */
            inline Real()
            {
                mpfr_init(d_value);
            }
            
            /**
               \brief Releases the memory used by the arbitrary precision floating point number.
             */
            inline ~Real() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
#if __cplusplus >= 201103L
                if (d_value[0]._mpfr_d != NULL)
#endif
                    mpfr_clear(d_value);
            }
            
#if __cplusplus >= 201103L
            /**
               \brief Creates a new real and moves the given real into this one.
               
               \param r A real to be moved into this one.
             */
            inline Real(Real && r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_value{r.d_value[0]}
            {
                r.d_value[0]._mpfr_d = NULL;
            }
            
            /**
               \brief Moves the given real into this one.
               
               \param r A real to be moved into this one.
             */
            inline Real & operator = (Real && r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                mpfr_clear(d_value);
                d_value[0] = r.d_value[0];
                r.d_value[0]._mpfr_d = NULL;
                return *this;
            }
#endif
            
            /**
               \brief Creates a real with the precision given by the context.
               
               \param rc A context.
             */
            inline explicit Real(const RealContext & rc)
            {
                mpfr_init2(d_value, rc.getRealPrecision());
            }
            
            /**
               \brief Creates a copy of the given floating point number.
               
               \param r The floating point number to copy.
               \param clone If `true`, creates a copy with the same precision. If `false`, creates a
                            copy of the default precision.
             */
            inline Real(const Real & r, bool clone = true)
            {
                if (clone)
                {
                    mpfr_init2(d_value, mpfr_get_prec(r.d_value));
                    mpfr_set(d_value, r.d_value, MPFR_RNDN);
                }
                else
                    mpfr_init_set(d_value, r.d_value, MPFR_RNDN);
            }
            
            /**
               \brief Creates a copy of the given floating point number with the given context's
                      precision.

               \param r The floating point number to copy.
               \param rc The context whose precision to use for the result.
             */
            inline explicit Real(const Real & r, const RealContext & rc)
            {
                mpfr_init2(d_value, rc.getRealPrecision());
                mpfr_set(d_value, r.d_value, MPFR_RNDN);
            }
            
            /**
               \brief Retrieves the precision of this floating point number.
               
               \return The precision in mantissa bits of the current floating point number.
             */
            inline unsigned long precision() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return mpfr_get_prec(d_value);
            }
            
            /**
               \brief Adjusts the precision of this floating point number to the given context's
                      precision.
               
               \param rc A context whose precision to use.
             */
            inline void setContext(const RealContext & rc)
            {
                if (rc.getRealPrecision() != precision())
                {
                    mpfr_t tmp;
                    mpfr_init2(tmp, rc.getRealPrecision());
                    mpfr_set(tmp, d_value, MPFR_RNDN);
                    mpfr_swap(tmp, d_value);
                    mpfr_clear(tmp);
                }
            }
            
            /**
               \brief Creates a new floating point number from the given native integer.
               
               \param i The native integer whose value to take.
             */
            inline explicit Real(long i)
            {
                mpfr_init_set_si(d_value, i, MPFR_RNDN);
            }
        
            /**
               \brief Creates a new floating point number from the given native integer with the
                      given context's precision.
               
               \param i The native integer whose value to take.
               \param rc The context whose precision to take.
             */
            inline explicit Real(long i, const RealContext & rc)
            {
                mpfr_init2(d_value, rc.getRealPrecision());
                mpfr_set_si(d_value, i, MPFR_RNDN);
            }
        
            /**
               \brief Creates a new floating point number from the given native integer.
               
               \param i The native integer whose value to take.
             */
            inline explicit Real(unsigned long i)
            {
                mpfr_init_set_ui(d_value, i, MPFR_RNDN);
            }
        
            /**
               \brief Creates a new floating point number from the given native integer with the
                      given context's precision.
               
               \param i The native integer whose value to take.
               \param rc The context whose precision to take.
             */
            inline explicit Real(unsigned long i, const RealContext & rc)
            {
                mpfr_init2(d_value, rc.getRealPrecision());
                mpfr_set_ui(d_value, i, MPFR_RNDN);
            }
        
            /**
               \brief Creates a new floating point number from the given native integer.
               
               \param i The native integer whose value to take.
             */
            inline explicit Real(long long i)
            {
                mpfr_init(d_value);
                mpfr_set_ll(d_value, i, MPFR_RNDN);
            }
        
            /**
               \brief Creates a new floating point number from the given native integer with the
                      given context's precision.
               
               \param i The native integer whose value to take.
               \param rc The context whose precision to take.
             */
            inline explicit Real(long long i, const RealContext & rc)
            {
                mpfr_init2(d_value, rc.getRealPrecision());
                mpfr_set_ll(d_value, i, MPFR_RNDN);
            }
        
            /**
               \brief Creates a new floating point number from the given native floating point
                      number.
               
               \param d The native floating point number whose value to take.
             */
            inline explicit Real(double d)
            {
                mpfr_init_set_d(d_value, d, MPFR_RNDN);
            }
        
            /**
               \brief Creates a new floating point number from the given native floating point
                      number with the given context's precision.
               
               \param d The native floating point number whose value to take.
               \param rc The context whose precision to take.
             */
            inline explicit Real(double d, const RealContext & rc)
            {
                mpfr_init2(d_value, rc.getRealPrecision());
                mpfr_set_d(d_value, d, MPFR_RNDN);
            }
        
            /**
               \brief Creates a new floating point number from the given native floating point
                      number.
               
               \param d The native floating point number whose value to take.
             */
            inline explicit Real(long double d)
            {
                mpfr_init_set_ld(d_value, d, MPFR_RNDN);
            }
        
            /**
               \brief Creates a new floating point number from the given native floating point
                      number with the given context's precision.
               
               \param d The native floating point number whose value to take.
               \param rc The context whose precision to take.
             */
            inline explicit Real(long double d, const RealContext & rc)
            {
                mpfr_init2(d_value, rc.getRealPrecision());
                mpfr_set_ld(d_value, d, MPFR_RNDN);
            }
        
            /**
               \brief Creates a new floating point number from the given arbitrary precision
                      integer.
               
               \param i The arbitrary precision integer whose value to take.
             */
            inline explicit Real(const Integer & i)
            {
                mpfr_init_set_z(d_value, i.d_value, MPFR_RNDN);
            }
        
            /**
               \brief Creates a new floating point number from the given arbitrary precision integer
                      with the given context's precision.
               
               \param i The arbitrary precision integer whose value to take.
               \param rc The context whose precision to take.
             */
            inline explicit Real(const Integer & i, const RealContext & rc)
            {
                mpfr_init2(d_value, rc.getRealPrecision());
                mpfr_set_z(d_value, i.d_value, MPFR_RNDN);
            }
            
            /**
               \brief Creates a floating point number from the given floating point expression.
               
               \param E The expression used to create this floating point number.
             */
            template<class A, template<typename, typename> class O>
            inline Real(const expressions::Expression<RealContext, A, O> & E)
            {
                mpfr_init2(d_value, E.precision());
                E.assignTo(*this);
            }
            
            /**
               \brief Creates a floating point number from the given floating point expression and
                      Real context.
               
               \param E The expression used to create this floating point number.
               \param rc The context whose precision to take.
             */
            template<class A, template<typename, typename> class O>
            inline Real(const expressions::Expression<RealContext, A, O> & E, const RealContext & rc)
            {
                mpfr_init2(d_value, rc.getRealPrecision());
                E.assignTo(*this);
            }
            
            /**
               \brief Assigns the given floating point number `r` to this floating point number.
               
               \param r The floating point number to assign to this floating point number.
               \return A reference to this floating point number.
             */
            inline Real & operator = (const Real & r)
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(d_value) == mpfr_get_prec(r.d_value));
              #endif
                mpfr_set(d_value, r.d_value, MPFR_RNDN);
                return *this;
            }
            
            /**
               \brief Assigns the floating point expression `E` to this floating point number.
               
               \param E The floating point expression.
               \return A reference to this floating point number.
             */
            template<class A, template<typename, typename> class O>
            inline Real & operator = (const expressions::Expression<RealContext, A, O> & E)
            {
                E.assignTo(*this);
                return *this;
            }
            
            inline friend bool operator == (const Real & a, const Real & b);
            inline friend bool operator != (const Real & a, const Real & b);
            inline friend bool operator <= (const Real & a, const Real & b);
            inline friend bool operator >= (const Real & a, const Real & b);
            inline friend bool operator < (const Real & a, const Real & b);
            inline friend bool operator > (const Real & a, const Real & b);
            
            inline friend bool isZero(const Real & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            // Tests whether the given value is = 0.
            {
                return mpfr_zero_p(r.d_value);
            }
            
            inline friend bool isOne(const Real & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            // Tests whether the given value is = 1.
            {
                return mpfr_cmp_ui(r.d_value, 1) == 0;
            }
            
            inline friend bool isPositive(const Real & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            // Tests whether the given value is > 0.
            {
                return mpfr_sgn(r.d_value) > 0;
            }
            
            inline friend bool isNonNegative(const Real & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            // Tests whether the given value is >= 0.
            {
                return mpfr_sgn(r.d_value) >= 0;
            }
            
            inline friend bool isNegative(const Real & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            // Tests whether the given value is < 0.
            {
                return mpfr_sgn(r.d_value) < 0;
            }
            
            inline friend bool isNonPositive(const Real & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            // Tests whether the given value is <= 0.
            {
                return mpfr_sgn(r.d_value) <= 0;
            }
            
            inline friend void add(Real & r, const Real & a, const Real & b)
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                {
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(a.d_value));
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(b.d_value));
                }
              #endif
                mpfr_add(r.d_value, a.d_value, b.d_value, MPFR_RNDN);
            }
            
            inline friend void add_ui(Real & r, const Real & a, unsigned long b)
            // Internal for expression templates
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(a.d_value));
              #endif
                mpfr_add_ui(r.d_value, a.d_value, b, MPFR_RNDN);
            }
            
            inline friend void add_si(Real & r, const Real & a, signed long b)
            // Internal for expression templates
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(a.d_value));
              #endif
                mpfr_add_si(r.d_value, a.d_value, b, MPFR_RNDN);
            }
            
            inline friend void add_d(Real & r, const Real & a, double b)
            // Internal for expression templates
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(a.d_value));
              #endif
                mpfr_add_d(r.d_value, a.d_value, b, MPFR_RNDN);
            }
            
            inline friend void add_z(Real & r, const Real & a, const Integer & b)
            // Internal for expression templates
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(a.d_value));
              #endif
                mpfr_add_z(r.d_value, a.d_value, b.getInternal(), MPFR_RNDN);
            }
            
            inline friend void sub(Real & r, const Real & a, const Real & b)
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                {
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(a.d_value));
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(b.d_value));
                }
              #endif
                mpfr_sub(r.d_value, a.d_value, b.d_value, MPFR_RNDN);
            }
            
            inline friend void sub_ui(Real & r, const Real & a, unsigned long b)
            // Internal for expression templates
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(a.d_value));
              #endif
                mpfr_sub_ui(r.d_value, a.d_value, b, MPFR_RNDN);
            }
            
            inline friend void sub_si(Real & r, const Real & a, signed long b)
            // Internal for expression templates
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(a.d_value));
              #endif
                mpfr_sub_si(r.d_value, a.d_value, b, MPFR_RNDN);
            }
            
            inline friend void sub_d(Real & r, const Real & a, double b)
            // Internal for expression templates
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(a.d_value));
              #endif
                mpfr_sub_d(r.d_value, a.d_value, b, MPFR_RNDN);
            }
            
            inline friend void sub_z(Real & r, const Real & a, const Integer & b)
            // Internal for expression templates
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(a.d_value));
              #endif
                mpfr_sub_z(r.d_value, a.d_value, b.getInternal(), MPFR_RNDN);
            }
            
            inline friend void ui_sub(Real & r, unsigned long a, const Real & b)
            // Internal for expression templates
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(b.d_value));
              #endif
                mpfr_ui_sub(r.d_value, a, b.d_value, MPFR_RNDN);
            }
            
            inline friend void si_sub(Real & r, signed long a, const Real & b)
            // Internal for expression templates
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(b.d_value));
              #endif
                mpfr_si_sub(r.d_value, a, b.d_value, MPFR_RNDN);
            }
            
            inline friend void d_sub(Real & r, double a, const Real & b)
            // Internal for expression templates
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(b.d_value));
              #endif
                mpfr_d_sub(r.d_value, a, b.d_value, MPFR_RNDN);
            }
            
            inline friend void z_sub(Real & r, const Integer & a, const Real & b)
            // Internal for expression templates
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(b.d_value));
              #endif
                mpfr_z_sub(r.d_value, a.getInternal(), b.d_value, MPFR_RNDN);
            }
            
            inline friend void addmul(Real & r, const Real & a, const Real & b) // r += a * b
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                {
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(a.d_value));
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(b.d_value));
                }
              #endif
                mpfr_fma(r.d_value, a.d_value, b.d_value, r.d_value, MPFR_RNDN);
            }
            
            inline friend void addmul4(Real & r, const Real & a, const Real & b, const Real & c) // r = a * b + c
            // Internal for expression templates
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                {
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(a.d_value));
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(b.d_value));
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(c.d_value));
                }
              #endif
                mpfr_fma(r.d_value, a.d_value, b.d_value, c.d_value, MPFR_RNDN);
            }
            
            inline friend void submul(Real & r, const Real & a, const Real & b) // r -= a * b
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                {
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(a.d_value));
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(b.d_value));
                }
              #endif
                mpfr_fms(r.d_value, a.d_value, b.d_value, r.d_value, MPFR_RNDN);
                mpfr_neg(r.d_value, r.d_value, MPFR_RNDN);
            }
            
            inline friend void submul4(Real & r, const Real & a, const Real & b, const Real & c) // r = a * b - c
            // Internal for expression templates
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                {
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(a.d_value));
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(b.d_value));
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(c.d_value));
                }
              #endif
                mpfr_fms(r.d_value, a.d_value, b.d_value, c.d_value, MPFR_RNDN);
            }
            
            inline friend void mul(Real & r, const Real & a, const Real & b)
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                {
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(a.d_value));
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(b.d_value));
                }
              #endif
                mpfr_mul(r.d_value, a.d_value, b.d_value, MPFR_RNDN);
            }
            
            inline friend void mul_ui(Real & r, const Real & a, unsigned long b)
            // Internal for expression templates
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(a.d_value));
              #endif
                mpfr_mul_ui(r.d_value, a.d_value, b, MPFR_RNDN);
            }
            
            inline friend void mul_si(Real & r, const Real & a, signed long b)
            // Internal for expression templates
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(a.d_value));
              #endif
                mpfr_mul_si(r.d_value, a.d_value, b, MPFR_RNDN);
            }
            
            inline friend void mul_d(Real & r, const Real & a, double b)
            // Internal for expression templates
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(a.d_value));
              #endif
                mpfr_mul_d(r.d_value, a.d_value, b, MPFR_RNDN);
            }
            
            inline friend void mul_z(Real & r, const Real & a, const Integer & b)
            // Internal for expression templates
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(a.d_value));
              #endif
                mpfr_mul_z(r.d_value, a.d_value, b.getInternal(), MPFR_RNDN);
            }
            
            inline friend void div(Real & r, const Real & a, const Real & b)
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                {
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(a.d_value));
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(b.d_value));
                }
              #endif
                mpfr_div(r.d_value, a.d_value, b.d_value, MPFR_RNDN);
            }
            
            inline friend void div_ui(Real & r, const Real & a, unsigned long b)
            // Internal for expression templates
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(a.d_value));
              #endif
                mpfr_div_ui(r.d_value, a.d_value, b, MPFR_RNDN);
            }
            
            inline friend void div_si(Real & r, const Real & a, signed long b)
            // Internal for expression templates
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(a.d_value));
              #endif
                mpfr_div_si(r.d_value, a.d_value, b, MPFR_RNDN);
            }
            
            inline friend void div_d(Real & r, const Real & a, double b)
            // Internal for expression templates
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(a.d_value));
              #endif
                mpfr_div_d(r.d_value, a.d_value, b, MPFR_RNDN);
            }
            
            inline friend void div_z(Real & r, const Real & a, const Integer & b)
            // Internal for expression templates
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(a.d_value));
              #endif
                mpfr_div_z(r.d_value, a.d_value, b.getInternal(), MPFR_RNDN);
            }
            
            inline friend void ui_div(Real & r, unsigned long a, const Real & b)
            // Internal for expression templates
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(b.d_value));
              #endif
                mpfr_ui_div(r.d_value, a, b.d_value, MPFR_RNDN);
            }
            
            inline friend void si_div(Real & r, signed long a, const Real & b)
            // Internal for expression templates
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(b.d_value));
              #endif
                mpfr_si_div(r.d_value, a, b.d_value, MPFR_RNDN);
            }
            
            inline friend void d_div(Real & r, double a, const Real & b)
            // Internal for expression templates
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(b.d_value));
              #endif
                mpfr_d_div(r.d_value, a, b.d_value, MPFR_RNDN);
            }
            
            inline friend void mod(Real & r, const Real & a, const Real & b)
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                {
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(a.d_value));
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(b.d_value));
                }
              #endif
                mpfr_fmod(r.d_value, a.d_value, b.d_value, MPFR_RNDN);
            }
            
            inline friend void divmod(Real & q, Real & r, const Real & a, const Real & b)
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                {
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(a.d_value));
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(b.d_value));
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(q.d_value));
                }
              #endif
                if ((&r == &a) || (&r == &b))
                {
                    if ((&q == &a) || (&q == &b) || (mpfr_get_prec(q.d_value) != mpfr_get_prec(r.d_value)))
                    {
                        mpfr_t tmp;
                        mpfr_init2(tmp, mpfr_get_prec(r.d_value));
                        mpfr_fmod(tmp, a.d_value, b.d_value, MPFR_RNDN);
                        mpfr_fms(q.d_value, tmp, b.d_value, a.d_value, MPFR_RNDN);
                        mpfr_neg(q.d_value, q.d_value, MPFR_RNDN);
                        mpfr_swap(r.d_value, tmp);
                        mpfr_clear(tmp);
                    }
                    else
                    {
                        mpfr_fmod(q.d_value, a.d_value, b.d_value, MPFR_RNDN);
                        mpfr_fms(r.d_value, q.d_value, b.d_value, a.d_value, MPFR_RNDN);
                        mpfr_neg(r.d_value, r.d_value, MPFR_RNDN);
                        mpfr_swap(q.d_value, r.d_value);
                    }
                }
                else
                {
                    mpfr_fmod(r.d_value, a.d_value, b.d_value, MPFR_RNDN);
                    mpfr_fms(q.d_value, r.d_value, b.d_value, a.d_value, MPFR_RNDN);
                    mpfr_neg(q.d_value, q.d_value, MPFR_RNDN);
                }
            }
                
            inline friend void shl(Real & r, const Real & a, const Real & b)
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                {
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(a.d_value));
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(b.d_value));
                }
              #endif
                mpfr_t tmp;
                mpfr_init2(tmp, mpfr_get_prec(r.d_value));
                mpfr_ui_pow(tmp, 2, b.d_value, MPFR_RNDN);
                mpfr_mul(r.d_value, a.d_value, tmp, MPFR_RNDN);
                mpfr_clear(tmp);
            }
            
            inline friend void shr(Real & r, const Real & a, const Real & b)
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                {
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(a.d_value));
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(b.d_value));
                }
              #endif
                mpfr_t tmp;
                mpfr_init2(tmp, mpfr_get_prec(r.d_value));
                mpfr_ui_pow(tmp, 2, b.d_value, MPFR_RNDN);
                mpfr_div(r.d_value, a.d_value, tmp, MPFR_RNDN);
                mpfr_clear(tmp);
            }
            
            inline friend void shl(Real & r, const Real & a, long b)
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(a.d_value));
              #endif
                mpfr_mul_2si(r.d_value, a.d_value, b, MPFR_RNDN);
            }
            
            inline friend void shr(Real & r, const Real & a, long b)
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(a.d_value));
              #endif
                mpfr_div_2si(r.d_value, a.d_value, b, MPFR_RNDN);
            }
            
            inline friend void increment(Real & r, const Real & a)
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(a.d_value));
              #endif
                mpfr_add_ui(r.d_value, a.d_value, 1, MPFR_RNDN);
            }
            
            inline friend void decrement(Real & r, const Real & a)
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(a.d_value));
              #endif
                mpfr_sub_ui(r.d_value, a.d_value, 1, MPFR_RNDN);
            }
            
            inline friend void neg(Real & r, const Real & a)
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(a.d_value));
              #endif
                mpfr_neg(r.d_value, a.d_value, MPFR_RNDN);
            }
            
            inline friend void abs(Real & r, const Real & a)
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(a.d_value));
              #endif
                mpfr_abs(r.d_value, a.d_value, MPFR_RNDN);
            }
            
            inline friend void makeAbs(Real & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                abs(a, a);
            }
            
            friend std::ostream & operator << (std::ostream & s, const Real & r);
            // Output to stream
            
            friend std::istream & operator >> (std::istream & s, Real & r);
            // Input from stream
            
            inline friend void setNaN(Real & r)
            {
                mpfr_set_nan(r.d_value);
            }
            
            inline friend void setInfinity(Real & r, bool sign)
            {
                mpfr_set_inf(r.d_value, sign ? 1 : -1);
            }
            
            inline friend void setZero(Real & r, bool sign)
            // Set to zero
            {
              #if (MPFR_VERSION <= 0x0204FF)
                mpfr_set_ui(r.d_value, 0, MPFR_RNDN);
                if (!sign)
                    mpfr_neg(r.d_value, r.d_value, MPFR_RNDN);
              #else
                mpfr_set_zero(r.d_value, sign ? 1 : -1);
              #endif
            }
            
            inline friend void setOne(Real & r)
            // Set to one
            {
                mpfr_set_ui(r.d_value, 1, MPFR_RNDN);
            }
            
            inline friend int compare(const Real & a, const Real & b)
            // Tests whether the first real is < (-1), = (0) or > (1) than the second.
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(a.d_value) == mpfr_get_prec(b.d_value));
              #endif
                return mpfr_cmp(a.d_value, b.d_value);
            }
            
            inline friend int compare_ui(const Real & a, unsigned long b)
            // Internal for expression templates
            {
                return mpfr_cmp_ui(a.d_value, b);
            }
            
            inline friend int compare_si(const Real & a, signed long b)
            // Internal for expression templates
            {
                return mpfr_cmp_si(a.d_value, b);
            }
            
            inline friend int compare_d(const Real & a, double b)
            // Internal for expression templates
            {
                return mpfr_cmp_d(a.d_value, b);
            }
            
            inline friend int compare_ld(const Real & a, long double b)
            // Internal for expression templates
            {
                return mpfr_cmp_ld(a.d_value, b);
            }
            
            inline friend int compare_z(const Real & a, const Integer & b)
            // Internal for expression templates
            {
                return mpfr_cmp_z(a.d_value, b.getInternal());
            }
            
            inline friend int compare_ui2(const Real & a, unsigned long b, long exp) // compares with b*2^exp
            // Internal for expression templates
            {
                return mpfr_cmp_ui_2exp(a.d_value, b, exp);
            }
            
            inline friend int compare_si2(const Real & a, signed long b, long exp) // compares with b*2^exp
            // Internal for expression templates
            {
                return mpfr_cmp_si_2exp(a.d_value, b, exp);
            }
            
            inline friend int compareAbsValues(const Real & a, const Real & b)
            // Tests whether the absolute value of the first Real is < (-1), = (0) or > (1) than the
            // absolute value of the second Real.
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(a.d_value) == mpfr_get_prec(b.d_value));
              #endif
                return mpfr_cmpabs(a.d_value, b.d_value);
            }
            
            inline friend int sign(const Real & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            // Returns the sign (i.e. -1, 0, 1).
            {
                return mpfr_sgn(r.d_value);
            }
            
            inline friend void square(Real & r, const Real & a)
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(r.d_value) == mpfr_get_prec(a.d_value));
              #endif
                mpfr_sqr(r.d_value, a.d_value, MPFR_RNDN);
            }
            
            inline friend void sin(Real & res, const Real & a)
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(a.d_value) == mpfr_get_prec(res.d_value));
              #endif
                mpfr_sin(res.d_value, a.d_value, MPFR_RNDN);
            }
            
            inline friend void cos(Real & res, const Real & a)
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(a.d_value) == mpfr_get_prec(res.d_value));
              #endif
                mpfr_cos(res.d_value, a.d_value, MPFR_RNDN);
            }
            
            inline friend void tan(Real & res, const Real & a)
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(a.d_value) == mpfr_get_prec(res.d_value));
              #endif
                mpfr_tan(res.d_value, a.d_value, MPFR_RNDN);
            }
            
            inline friend void asin(Real & res, const Real & a)
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(a.d_value) == mpfr_get_prec(res.d_value));
              #endif
                mpfr_asin(res.d_value, a.d_value, MPFR_RNDN);
            }
            
            inline friend void acos(Real & res, const Real & a)
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(a.d_value) == mpfr_get_prec(res.d_value));
              #endif
                mpfr_acos(res.d_value, a.d_value, MPFR_RNDN);
            }
            
            inline friend void atan(Real & res, const Real & a)
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(a.d_value) == mpfr_get_prec(res.d_value));
              #endif
                mpfr_atan(res.d_value, a.d_value, MPFR_RNDN);
            }
            
            inline friend void atan2(Real & res, const Real & y, const Real & x)
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                {
                    assert(mpfr_get_prec(y.d_value) == mpfr_get_prec(res.d_value));
                    assert(mpfr_get_prec(x.d_value) == mpfr_get_prec(res.d_value));
                }
              #endif
                mpfr_atan2(res.d_value, y.d_value, x.d_value, MPFR_RNDN);
            }
            
            inline friend void exp(Real & res, const Real & a)
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(a.d_value) == mpfr_get_prec(res.d_value));
              #endif
                mpfr_exp(res.d_value, a.d_value, MPFR_RNDN);
            }
            
            inline friend void log(Real & res, const Real & a)
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(a.d_value) == mpfr_get_prec(res.d_value));
              #endif
                mpfr_log(res.d_value, a.d_value, MPFR_RNDN);
            }
            
            inline friend void log2(Real & res, const Real & a)
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(a.d_value) == mpfr_get_prec(res.d_value));
              #endif
                mpfr_log2(res.d_value, a.d_value, MPFR_RNDN);
            }
            
            inline friend void log10(Real & res, const Real & a)
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(a.d_value) == mpfr_get_prec(res.d_value));
              #endif
                mpfr_log10(res.d_value, a.d_value, MPFR_RNDN);
            }
            
            inline friend void sqrt(Real & res, const Real & a)
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(a.d_value) == mpfr_get_prec(res.d_value));
              #endif
                mpfr_sqrt(res.d_value, a.d_value, MPFR_RNDN);
            }
            
            inline friend void gamma(Real & res, const Real & a)
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(a.d_value) == mpfr_get_prec(res.d_value));
              #endif
                mpfr_gamma(res.d_value, a.d_value, MPFR_RNDN);
            }
            
            inline friend void lgamma(Real & res, const Real & a)
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(a.d_value) == mpfr_get_prec(res.d_value));
              #endif
                int sign;
                mpfr_lgamma(res.d_value, &sign, a.d_value, MPFR_RNDN);
            }
            
            inline friend void lgamma(Real & res, int & sign, const Real & a)
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(a.d_value) == mpfr_get_prec(res.d_value));
              #endif
                mpfr_lgamma(res.d_value, &sign, a.d_value, MPFR_RNDN);
            }
            
            inline friend void power(Real & res, const Real & a, signed long b)
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(res.d_value) == mpfr_get_prec(a.d_value));
              #endif
                mpfr_pow_si(res.d_value, a.d_value, b, MPFR_RNDN);
            }
            
            inline friend void power(Real & res, const Real & a, unsigned long b)
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(res.d_value) == mpfr_get_prec(a.d_value));
              #endif
                mpfr_pow_ui(res.d_value, a.d_value, b, MPFR_RNDN);
            }
            
            inline friend void power(Real & res, const Real & a, const Integer & b)
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                    assert(mpfr_get_prec(res.d_value) == mpfr_get_prec(a.d_value));
              #endif
                mpfr_pow_z(res.d_value, a.d_value, b.d_value, MPFR_RNDN);
            }
            
            inline friend void power(Real & res, const Real & a, const Real & b)
            {
              #ifdef PLLL_INTERNAL__BIGINT__HELP_FINDING_PRECISION_BUGS
                if (arithmetic::internal::Real_precision_check_enabled)
                {
                    assert(mpfr_get_prec(res.d_value) == mpfr_get_prec(a.d_value));
                    assert(mpfr_get_prec(res.d_value) == mpfr_get_prec(b.d_value));
                }
              #endif
                mpfr_pow(res.d_value, a.d_value, b.d_value, MPFR_RNDN);
            }
            
            inline friend void arithmetic::swap(plll::arithmetic::Real &, plll::arithmetic::Real &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            
            /**
               \brief Returns approximate `e` such that |x| is approximately \f$2^e\f$.
               
               The value of `e` should be OK up to +-2.
               
               \return \f$\approx \log_2 |x|\f$.
             */
            inline long getApproxExponent() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return (mpfr_zero_p(d_value) != 0) ? std::numeric_limits<long>::min() : mpfr_get_exp(d_value);
            }
            
            // These ones should only be used in very, very special and rare cases!
            const mpfr_t & getInternal() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE { return d_value; }
            mpfr_t & getInternal() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE { return d_value; }
        };
        
        inline void swap(Integer & a, Integer & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            mpz_swap(a.d_value, b.d_value);
        }
        
        // Internal for expression templates
        inline void add_ui(Real & r, const Real & a, unsigned long b);
        inline void add_si(Real & r, const Real & a, signed long b);
        inline void add_d(Real & r, const Real & a, double b);
        inline void add_z(Real & r, const Real & a, const Integer & b);
        inline void sub_ui(Real & r, const Real & a, unsigned long b);
        inline void sub_si(Real & r, const Real & a, signed long b);
        inline void sub_d(Real & r, const Real & a, double b);
        inline void sub_z(Real & r, const Real & a, const Integer & b);
        inline void ui_sub(Real & r, unsigned long a, const Real & b);
        inline void si_sub(Real & r, signed long a, const Real & b);
        inline void d_sub(Real & r, double a, const Real & b);
        inline void z_sub(Real & r, const Integer & a, const Real & b);
        inline void addmul4(Real & r, const Real & a, const Real & b, const Real & c); // r = a * b + c
        inline void submul4(Real & r, const Real & a, const Real & b, const Real & c); // r = a * b - c
        inline void mul_ui(Real & r, const Real & a, unsigned long b);
        inline void mul_si(Real & r, const Real & a, signed long b);
        inline void mul_d(Real & r, const Real & a, double b);
        inline void mul_z(Real & r, const Real & a, const Integer & b);
        inline void div_ui(Real & r, const Real & a, unsigned long b);
        inline void div_si(Real & r, const Real & a, signed long b);
        inline void div_d(Real & r, const Real & a, double b);
        inline void div_z(Real & r, const Real & a, const Integer & b);
        inline void ui_div(Real & r, unsigned long a, const Real & b);
        inline void si_div(Real & r, signed long a, const Real & b);
        inline void d_div(Real & r, double a, const Real & b);
        inline int compare_ui(const Real & a, unsigned long b);
        inline int compare_si(const Real & a, signed long b);
        inline int compare_d(const Real & a, double b);
        inline int compare_ld(const Real & a, long double b);
        inline int compare_z(const Real & a, const Integer & b);
        inline int compare_ui2(const Real & a, unsigned long b, long exp); // compares with b*2^exp
        inline int compare_si2(const Real & a, signed long b, long exp); // compares with b*2^exp
    }
}

#include "arithmetic-gmp-rops.hpp"

namespace plll
{
    namespace arithmetic
    {
        /**
           \brief Represents a random number generator.
           
           As the underlying generator, the <a
           href="https://en.wikipedia.org/wiki/Mersenne_twister">Mersenne Twister</a> is used.

           \warning Note that while this random number generator is perfectly fine for most
                    purposes, it is _not_ for cryptographical purposes.
         */
        class RandomNumberGenerator
        {
            // stores information on a random number generator
        
        private:
            void * d_state; // we hide the implementation
            Integer d_seed;
        
        public:
            /**
               \brief Creates a new default initialized random number generator.
             */
            RandomNumberGenerator();
            
            /**
               \brief Creates a clone of the given random number generator.
             */
            RandomNumberGenerator(const RandomNumberGenerator &);
            
            /**
               \brief Releases the state for this random number generator.
             */
            ~RandomNumberGenerator();
            
            /**
               \brief Copies the state of the given random number generator to the current one.
               
               \param rng The random number generator whose state should be copied.
               \return A reference to the current random number generator.
             */
            RandomNumberGenerator & operator = (const RandomNumberGenerator & rng);
            
            /**
               \brief Retrieves the last set seed of the generator.
               
               \return The previously set seed.
             */
            Integer getSeed() const;
            /**
               \brief Sets a new seed for the generator.
               
               The state is generated from this seed, so a short or small seed does not lead to a
               bad state.
               
               \param seed The integer to use as a seed.
             */
            void setSeed(const Integer & seed);
            
            /**
               \brief Uses the system's time to seed this random number generator.
               
               \warning This is implementation dependent and might behave completely different even
                        on the same system if two programs are run at the precisely same time.
             */
            void randomizeTime();

            /**
               \brief Uses the operating system's random number generator to initialize this random
                      number generator.
               
               \param goodness Specify how much of the state should be taken from /dev/random
                               instead of /dev/urandom. The default fraction is very small, to
                               minimize the amount of blocking when trying to read too much from
                               /dev/random. A value of 0.0 reads everything from /dev/urandom, and
                               thus might be completely predictable to a watcher, while a value of
                               1.0 reads everything from /dev/random, resulting in a state
                               completely independent from previously read random values.
                               
                               For most purposes, the default value 0.01 is totally fine, and even
                               0.0 is fine in many cases. Note that only a value of 0.0 garuantees
                               that absolutely no blocking occurs.
             */
            void randomizeDevRandom(double goodness = 0.01);
            /**
               \brief Creates a seed using an internal random number generator.
               
               The internal random number generator is initialized using /dev/random and
               /dev/urandom on the first call to this function.
               
               \return A seed for use with `setSeed()`.
               
               \warning Might block on the first call. (But only then, and this can be avoided by
                        calling `RandomNumberGenerator::randomizeSeed()` earlier.)
               
               \sa `RandomNumberGenerator::randomizeSeed()`
             */
            static Integer createSeed();
            /**
               \brief Initializes this random number generator with a random seed.
               
               Equivalent to `setSeed(RandomNumberGenerator::createSeed())`.
             */
            void randomizeSeed();
            /**
               \brief Calls `randomizeDevRandom(goodness)` for the internal random number generator
                      which is used for `createSeed()`.
               
               Allows to initialize the internal random number generator at this point of time. If
               `RandomNumberGenerator::createSeed()` is then called later during the program, no
               blocking will occur at that point.
               
               \warning This call might block until enough entropy is available to the operating
                        system.
               
               \sa `RandomNumberGenerator::createSeed()`
             */
            static void initializeSeeder(double goodness = 0.01);
            
            /**
               \brief Creates a random arbitrary precision integer in the range \f$[0, bound)\f$.
               
               \param res Where to store the result.
               \param bound A bound on the maximal integer returned.
             */
            void random(Integer & res, const Integer & bound);
            
            /**
               \brief Creates and returns a random arbitrary precision integer in the range \f$[0, bound)\f$.
               
               \param bound A bound on the maximal integer returned.
               \return The result.
             */
            Integer random(const Integer & bound);
            
            /**
               \brief Creates and returns a native integer in the range \f$[0, bound)\f$.
               
               \param bound A bound on the maximal integer returned.
               \return The result.
             */
            unsigned long random(unsigned long bound);
            
            /**
               \brief Writes a given number of random bytes to the memory location pointed to.
               
               \param ptr A pointer to the memory location.
               \param count The number of bytes to write there.
             */
            void randomBits(void * ptr, unsigned long count);
            
            /**
               \brief Creates random bits.
               
               \param res Will be filled with a random integer in range \f$[0, 2^{bits})\f$.
               \param bits The number of bits.
             */
            void randomBits(Integer & res, unsigned long bits);
            
            /**
               \brief Creates random bits.
               
               \param bits The number of bits.
               \return A random integer in range \f$[0, 2^{bits})\f$.
             */
            Integer randomBits(unsigned long bits);
            
            /**
               \brief Creates a random integer of a fixed bit length.
               
               \param res Will be filled with a random integer in range \f$[2^{bits-1}, 2^{bits})\f$.
               \param bits The number of bits.
             */
            void randomLen(Integer & res, unsigned long bits);
            
            /**
               \brief Creates a random integer of a fixed bit length.
               
               \param bits The number of bits.
               \return A random integer in range \f$[2^{bits-1}, 2^{bits})\f$.
             */
            Integer randomLen(unsigned long bits);
            
            /**
               \brief Creates a uniformly distributed floating point number in the interval \f$[0,
                      1)\f$.
               
               \param r Where to store the result.
             */
            void randomUniform(Real & r);
            
            /**
               \brief Creates a uniformly distributed floating point number in the interval \f$[0,
                      1)\f$.
               
               \param r Where to store the result.
               \param rc The context whose precision to use for the result.
             */
            void randomUniform(Real & r,const RealContext & rc);
        
            /**
               \brief Creates and returns a uniformly distributed floating point number in the
                      interval \f$[0, 1)\f$.
               
               \result The random number.
             */
            inline Real randomUniform()
            {
                Real r;
                randomUniform(r);
                return r;
            }
        
            /**
               \brief Creates and returns a uniformly distributed floating point number in the
                      interval \f$[0, 1)\f$.
               
               \param rc The context whose precision to use for the result.
               \return The random number.
             */
            inline Real randomUniform(const RealContext & rc)
            {
                Real r(rc);
                randomUniform(r, rc);
                return r;
            }
        };
        
        class RealContext::UniformRNG
        {
        private:
            RandomNumberGenerator & d_rng;
        
        public:
            /**
               \brief Creates a new uniform random number generator based on the given random number
                      generator.
             */
            UniformRNG(RandomNumberGenerator & rng) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_rng(rng)
            {
            }
            
            /**
               \brief Creates and returns a uniformly distributed floating point number in the
                      interval \f$[0, 1)\f$.
               
               \param rc The context whose precision to use for the result.
               \return The random number.
             */
            inline Real randomUniform(const RealContext & rc)
            {
                Real r(rc);
                randomUniform(r);
                return r;
            }
            
            /**
               \brief Creates a uniformly distributed floating point number in the interval \f$[0,
                      1)\f$.
               
               \param r Where to store the result.
             */
            inline void randomUniform(Real & r)
            {
                d_rng.randomUniform(r);
            }
        };
    
        inline Real RealContext::getEpsilon() const
        {
            Real r(*this);
            mpfr_set_ui_2exp(r.d_value, 1, 1 - getRealPrecision(), MPFR_RNDN);
            return r;
        }
    
        inline void RealContext::getEpsilon(Real & x) const
        {
            if ((unsigned long)mpfr_get_prec(x.d_value) != getRealPrecision())
            {
                mpfr_clear(x.d_value);
                mpfr_init2(x.d_value, getRealPrecision());
            }
            mpfr_set_ui_2exp(x.d_value, 1, 1 - getRealPrecision(), MPFR_RNDN);
        }
     
        inline Real RealContext::getPi() const
        {
            Real x(*this);
            getPi(x);
            return x;
        }
    
        inline void RealContext::getPi(Real & x) const
        {
            mpfr_const_pi(x.d_value, MPFR_RNDN);
        }
    
        inline Real RealContext::getEuler() const
        {
            Real x(*this);
            getEuler(x);
            return x;
        }
    
        inline void RealContext::getEuler(Real & x) const
        {
            mpfr_const_euler(x.d_value, MPFR_RNDN);
        }
    
        inline Real RealContext::getLog2() const
        {
            Real x(*this);
            getLog2(x);
            return x;
        }
    
        inline void RealContext::getLog2(Real & x) const
        {
            mpfr_const_log2(x.d_value, MPFR_RNDN);
        }
    
        class IntegerContext::UniformRNG
        {
        private:
            RandomNumberGenerator & d_rng;
            
        public:
            /**
               \brief Creates a new uniform random number generator based on the given random number
                      generator.
             */
            inline UniformRNG(RandomNumberGenerator & rng) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_rng(rng)
            {
            }
            
            /**
               \brief Creates a random arbitrary precision integer in the range \f$[0, bound)\f$.
               
               \param res Where to store the result.
               \param bound A bound on the maximal integer returned.
             */
            inline void random(Integer & res, const Integer & bound)
            {
                d_rng.random(res, bound);
            }
            
            /**
               \brief Creates and returns a random arbitrary precision integer in the range \f$[0, bound)\f$.
               
               \param bound A bound on the maximal integer returned.
               \param ic An integer context.
               \return The result.
             */
            inline Integer random(const Integer & bound, const IntegerContext & ic)
            {
                Integer res;
                random(res, bound);
                return res;
            }
            
            /**
               \brief Creates random bits.
               
               \param res Will be filled with a random integer in range \f$[0, 2^{bits})\f$.
               \param bits The number of bits.
             */
            inline void randomBits(Integer & res, unsigned long bits)
            {
                d_rng.randomBits(res, bits);
            }
            
            /**
               \brief Creates random bits.
               
               \param bits The number of bits.
               \param ic An integer context.
               \return A random integer in range \f$[0, 2^{bits})\f$.
             */
            inline Integer randomBits(unsigned long bits, const IntegerContext & ic)
            {
                Integer res;
                randomBits(res, bits);
                return res;
            }
            
            /**
               \brief Creates a random integer of a fixed bit length.
               
               \param res Will be filled with a random integer in range \f$[2^{bits-1}, 2^{bits})\f$.
               \param bits The number of bits.
             */
            inline void randomLen(Integer & res, unsigned long bits) // returns a random integer x with 2^(bits-1) <= x < 2^bits
            {
                d_rng.randomLen(res, bits);
            }
            
            /**
               \brief Creates a random integer of a fixed bit length.
               
               \param bits The number of bits.
               \param ic An integer context.
               \return A random integer in range \f$[2^{bits-1}, 2^{bits})\f$.
             */
            inline Integer randomLen(unsigned long bits, const IntegerContext & ic)
            {
                Integer res;
                randomLen(res, bits);
                return res;
            }
        };
        
        inline void swap(Real & a, Real & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            mpfr_swap(a.d_value, b.d_value);
        }
    }
}

#include "arithmetic-gmp-conv.hpp"

#endif

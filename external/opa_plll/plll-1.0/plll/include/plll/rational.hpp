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

#ifndef PLLL_INCLUDE_GUARD__RATIONAL_HPP
#define PLLL_INCLUDE_GUARD__RATIONAL_HPP

#include <plll/arithmetic.hpp>
#include "helper.hpp"
#include <sstream>
#include <cassert>
#include <cmath>

/**
   \file
   \brief Header for rational number arithmetic.
   
   This header provides rational number arithmetic for `plll`. It provides
   `plll::arithmetic::RationalContext`, an arithmetic context as described in \ref
   arithctxts. Concrete rational numbers are represented by `plll::arithmetic::Rational` objects,
   which essentially consist of two `plll::arithmetic::Integer` objects, and expressions are
   represented by `plll::arithmetic::expressions::expression<>` templates.
*/
namespace plll
{
    namespace arithmetic
    {
        class Rational;
        class RationalContext;
        
        /**@{
           \name Predicates.
        */
        /**
           \brief Tests the given rational number for being zero.
           
           \return Returns `true` if and only if the argument is zero.
           \sa \ref arithctxts_preds
        */
        inline bool isZero(const Rational & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        
        /**
           \brief Tests the given `plll::arithmetic::Real` object for being zero.
           
           \return Returns `true` if and only if the argument is zero.
           \sa \ref arithctxts_preds
        */
        inline bool isOne(const Rational & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        
        /**
           \brief Tests the given rational number for being strictly positive.
           
           \return Returns `true` if and only if the argument is strictly positive.
           \sa \ref arithctxts_preds
        */
        inline bool isPositive(const Rational & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        
        /**
           \brief Tests the given rational number for being positive or zero.
           
           \return Returns `true` if and only if the argument is positive or zero.
           \sa \ref arithctxts_preds
        */
        inline bool isNonNegative(const Rational & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        
        /**
           \brief Tests the given rational number for being strictly negative.
           
           \return Returns `true` if and only if the argument is strictly negative.
           \sa \ref arithctxts_preds
        */
        inline bool isNegative(const Rational & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        
        /**
           \brief Tests the given rational number for being negative or zero.
           
           \return Returns `true` if and only if the argument is negative or zero.
           \sa \ref arithctxts_preds
        */
        inline bool isNonPositive(const Rational & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        ///@}
        
        /**@{
           \name Functional versions of operators.
        */
        /**
           \brief Adds `a` and `b` and stores the result in `r`.
           
           \param r The result.
           \param a The first operand.
           \param b The second operand.
        */
        inline void add(Rational & r, const Rational & a, const Rational & b);
        
        /**
           \brief Subtracts `b` from `a` and stores the result in `r`.
           
           \param r The result.
           \param a The first operand.
           \param b The second operand.
        */
        inline void sub(Rational & r, const Rational & a, const Rational & b);
        
        /**
           \brief Multiplies `a` and `b` and adds the result to `r`.
           
           \param r The accumulator.
           \param a The first operand.
           \param b The second operand.
        */
        inline void addmul(Rational & r, const Rational & a, const Rational & b);
        
        /**
           \brief Multiplies `a` and `b` and subtracts the result from `r`.
           
           \param r The accumulator.
           \param a The first operand.
           \param b The second operand.
        */
        inline void submul(Rational & r, const Rational & a, const Rational & b);
        
        /**
           \brief Multiplies `a` with `b` and stores the result in `r`.
           
           \param r The result.
           \param a The first operand.
           \param b The second operand.
        */
        inline void mul(Rational & r, const Rational & a, const Rational & b);
        
        /**
           \brief Divides `a` by `b` and stores the result in `r`.
           
           \param r The result.
           \param a The first operand.
           \param b The second operand.
        */
        inline void div(Rational & r, const Rational & a, const Rational & b);
        
        /**
           \brief Takes the remainder of the division of `a` by `b` and stores it in `r`.
           
           \param r The result.
           \param a The first operand.
           \param b The second operand.
        */
        inline void mod(Rational & r, const Rational & a, const Rational & b);
        
        /**
           \brief Multiplies `a` by \f$2^b\f$ and stores the result in `r`.
           
           \param r The result.
           \param a The first operand.
           \param b The second operand.
        */
        inline void shl(Rational & r, const Rational & a, long b);
        
        /**
           \brief Multiplies `a` by \f$2^b\f$ and stores the result in `r`.
           
           \param r The result.
           \param a The first operand.
           \param b The second operand.
        */
        inline void shr(Rational & r, const Rational & a, long b);
        
        /**
           \brief Increments `a` by one and stores the result in `r`.
           
           \param r The result.
           \param a The operand.
        */
        inline void increment(Rational & r, const Rational & a);
        
        /**
           \brief Negates `a` and stores the result in `r`.
           
           \param r The result.
           \param a The operand.
        */
        inline void decrement(Rational & r, const Rational & a);
        
        /**
           \brief Negates `a` and stores the result in `r`.
           
           \param r The result.
           \param a The operand.
        */
        inline void neg(Rational & r, const Rational & a);
        
        /**
           \brief Takes the absolute value of `a` and stores the result in `r`.
           
           \param r The result.
           \param a The operand.
        */
        inline void abs(Rational & r, const Rational & a);
        
        /**
           \brief Computes the square of `a` and stores the result in `r`.
           
           \param r The result.
           \param a The operand.
        */
        inline void square(Rational & r, const Rational & a);
        ///@}
        
        /**@{
           \name Sign querying/modification.
        */
        /**
           \brief Makes the operand non-negative.
           
           \param a The operand.
        */
        inline void makeAbs(Rational & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        ///@}
        
        /**@{
           \name Stream input/output.
        */
        /**
           \brief Outputs the rational number on the given output stream.
        */
        std::ostream & operator << (std::ostream & s, const Rational & r);
        
        /**
           \brief Reads the rational number from the given input stream.
        */
        std::istream & operator >> (std::istream & s, Rational & r);
        ///@}
        
        /**@{
           \name Setting to specific constants.
        */
        /**
           \brief Sets the given rational number to \f$\pm 0\f$.
           
           Note that the sign is ignored, as there is no distinction between the rational +0 and
           -0.
           
           \param r The rational number variable to set to zero.
           \param sign This is ignored. The default value is `true`.
        */
        inline void setZero(Rational & r, bool sign = true);
        
        /**
           \brief Sets the given rational number to one.
           
           \param r The rational number variable to set to one.
        */
        inline void setOne(Rational & r);
        ///@}
        
        /**@{
           \name Comparisons.
        */
        /**
           \brief Compares the two rational numbers.
           
           \param a The first operand.
           \param b The second operand.
           \return Returns a negative number if \f$a < b\f$, zero if \f$a = b\f$ and a positive
                   number if \f$a > b\f$.
        */
        inline int compare(const Rational & a, const Rational & b);
        
        /**
           \brief Compares the two rational numbers in absolute value.
           
           \param a The first operand.
           \param b The second operand.
           \return Returns a negative number if \f$|a| < |b|\f$, zero if \f$|a| = |b|\f$ and a
                   positive number if \f$|a| > |b|\f$.
        */
        inline int compareAbsValues(const Rational & a, const Rational & b);
        ///@}
        
        /**@{
           \name Sign querying/modification.
        */
        /**
           \brief Returns the sign of the given rational number.
           
           \return Returns a negative value if the value is negative, 0 if it is zero, and a
                   positive value if it is positive.
           \sa \ref arithctxts_preds
         */
        inline int sign(const Rational & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        ///@}
        
        /**@{
           \name Exponentiation.
        */
        /**
           \brief Raises `a` to the power `e` and stores the result in `r`.
           
           \param r The result.
           \param a The base.
           \param e The exponent.
        */
        inline void power(Rational & r, const Rational & a, signed long e);
        
        /**
           \brief Raises `a` to the power `e` and stores the result in `r`.
           
           \param r The result.
           \param a The base.
           \param e The exponent.
        */
        inline void power(Rational & r, const Rational & a, unsigned long e);
        
        /**
           \brief Raises `a` to the power `e` and stores the result in `r`.
           
           \param r The result.
           \param a The base.
           \param e The exponent.
        */
        void power(Rational & r, const Rational & a, const Integer & e);
        ///@}
        
        /**@{
           \name Swap functions.
        */
        /**
           \brief Swaps two `plll::arithmetic::Rational` objects.
        */
        inline void swap(Rational &, Rational &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        ///@}
        
        /**
           \brief Represents an arithmetic context for rational numbers.
           
           \sa \ref arithctxts_contexts.
         */
        class RationalContext
        // Emulates RealContext
        {
            friend class Rational;
            
            enum { d_precision = 21474836479 }; // should be std::numeric_limits<unsigned long>::max(), but that is not working in enums.
            
        public:
            /**
               \brief The rational type.
             */
            typedef Rational Real;
            /**
               \brief The rational type.
             */
            typedef Rational Type;
            
            /**
               \brief The properties of rational numbers.
             */
            enum { is_cputype = false, is_realtype = true, is_inttype = false, is_exact = true, is_variable_precision = false,
                   has_squareroot = false, has_full_power = false, has_special_fns = false, has_huge_exponent = true,
                   has_infinity = false, has_uniform_rng = false, has_constants = false, has_trigonometric = false };
            
            /**
               \brief Creates a new rational context.
             */
            inline RationalContext() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
            }
            
            /**
               \brief Sets the precision of this context. (Will do nothing.)
               
               Since rationals always have infinite precision, this simply does nothing.
               
               \param prec The precision to be set.
             */
            static inline void setRealPrecision(unsigned long prec) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
            }
            
            /**
               \brief Returns the precision of this context.
               
               Will return a very large number.
             */
            static inline unsigned long getRealPrecision() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return d_precision;
            }
            
            /**
               \brief Returns the minimal possible precision for this context.
               
               Will return a very large number.
             */
            static inline unsigned long getMinRealPrecision() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return d_precision;
            }
            
            /**
               \brief Returns the maximal possible precision for this context.
               
               Will return a very large number.
             */
            static inline unsigned long getMaxRealPrecision() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE // returns maximal value for precision
            {
                return d_precision;
            }
        };
        
        /**
           \brief Represents a rational number as a quotient of two arbitrary precision integers.
         */
        class Rational
        {
#ifndef PLLL_INTERNAL_NO_TEMPLATE_FRIENDS
            template<class X, class Y>
            friend class implementation::conversion_impl;
            
            friend class implementation::nativeconversion_impl<Rational>;
            
        private:
#else
        public:
#endif
            
            Integer d_num, d_denom; // the denominator is always positive (and in particular, non-zero!)
            
            inline void normalize()
            {
                Integer d = GCD(d_num, d_denom);
                if (!isOne(d))
                {
                    d_num /= d;
                    d_denom /= d;
                }
            }
        
            void initFromDouble(double d)
            {
                const int MantissaBitsInDouble = std::numeric_limits<double>::digits + 2;
                int e;
                d = std::frexp(d, &e);
                // Now 0.5 <= d < 1
                convert(d_num, std::ldexp(d, MantissaBitsInDouble), IntegerContext());
                setOne(d_denom);
                // Add exponent
                shl(*this, *this, e - MantissaBitsInDouble);
            }
        
            void initFromDouble(long double d)
            {
                const int MantissaBitsInDouble = std::numeric_limits<long double>::digits + 2;
                int e;
                d = std::frexp(d, &e);
                // Now 0.5 <= d < 1
                convert(d_num, std::ldexp(d, MantissaBitsInDouble), IntegerContext());
                setOne(d_denom);
                // Add exponent
                shl(*this, *this, e - MantissaBitsInDouble);
            }
        
        public:
            typedef RationalContext Context;
            
            /**
               \brief Retrieves the numerator of this rational number.

               \return An `arithmetic::Integer` object representing the numerator.
             */
            inline const Integer & numerator() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return d_num;
            }
            
            /**
               \brief Retrieves the denominator of this rational number.

               \return An `arithmetic::Integer` object representing the denominator.
             */
            inline const Integer & denominator() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return d_denom;
            }
            
            /**
               \brief Creates a new rational number representing 0.
             */
            inline Rational()
            {
                setZero(d_num);
                setOne(d_denom);
            }
            
#if __cplusplus >= 201103L
            /**
               \brief Creates a new rational number and moves the given rational number into this
                      one.
               
               \param r A rational number to be moved into this one.
             */
            inline Rational(Rational && r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_num(std::move(r.d_num)), d_denom(std::move(r.d_denom))
            {
            }
            
            /**
               \brief Moves the given rational number into this one.
               
               \param r A rational number to be moved into this one.
             */
            Rational & operator = (Rational && r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                d_num = std::move(r.d_num);
                d_denom = std::move(r.d_denom);
                return *this;
            }
#endif
            
            /**
               \brief Creates a new rational number representing 0.
               
               \param rc A rational context.
             */
            inline Rational(const RationalContext & rc)
            {
                setZero(d_num);
                setOne(d_denom);
            }
            
            /**
               \brief Creates a copy of the given rational number.
               
               \param r The rational number to copy.
             */
            inline Rational(const Rational & r)
                : d_num(r.d_num), d_denom(r.d_denom)
            {
            }
            
            /**
               \brief Creates a copy of the given rational number.
               
               \param r The rational number to copy.
               \param rc A rational context.
             */
            inline Rational(const Rational & r, const RationalContext & rc)
                : d_num(r.d_num), d_denom(r.d_denom)
            {
            }
        
            /**
               \brief Evaluates the given rational expression into a new rational number.
               
               \param E A rational expression.
             */
            template<class A, template<typename, typename> class O>
            inline Rational(const expressions::Expression<RationalContext, A, O> & E)
            {
                E.assignTo(*this);
            }
            
            /**
               \brief Retrieves the precision of the current rational number.
               
               Will return a large number equal to the one returned by
               `RationalContext::getRealPrecision()`.
             */
            static inline unsigned long precision() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return RationalContext::d_precision;
            }
            
            /**
               \brief Applies the given rational context to this number.
               
               Does nothing.
               
               \param rc A rational context.
             */
            static inline void setContext(const RationalContext & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
            }
            
            /**
               \brief Creates a rational number from the given native floating point value.
               
               \param d The native floating point value.
             */
            inline explicit Rational(double d)
            {
                initFromDouble(d);
            }
        
            /**
               \brief Creates a rational number from the given native floating point value.
               
               \param d The native floating point value.
               \param rc A rational context.
             */
            inline explicit Rational(double d, const RationalContext & rc)
            {
                initFromDouble(d);
            }
            
            /**
               \brief Creates a rational number from the given native floating point value.
               
               \param d The native floating point value.
             */
            inline explicit Rational(long double d)
            {
                initFromDouble(d);
            }
            
            /**
               \brief Creates a rational number from the given native floating point value.
               
               \param d The native floating point value.
               \param rc A rational context.
             */
            inline explicit Rational(long double d, const RationalContext & rc)
            {
                initFromDouble(d);
            }
            
            /**
               \brief Creates a rational number from the given native integer.
               
               \param i The native integer.
             */
            inline explicit Rational(long i)
                : d_num(i)
            {
                setOne(d_denom);
            }
            
            /**
               \brief Creates a rational number from the given native integer.
               
               \param i The native integer.
               \param rc A rational context.
             */
            inline explicit Rational(long i, const RationalContext & rc)
                : d_num(i)
            {
                setOne(d_denom);
            }
            
            /**
               \brief Creates a rational number from the given native integer.
               
               \param i The native integer.
             */
            inline explicit Rational(unsigned long i)
                : d_num(i)
            {
                setOne(d_denom);
            }
            
            /**
               \brief Creates a rational number from the given native integer.
               
               \param i The native integer.
               \param rc A rational context.
             */
            inline explicit Rational(unsigned long i, const RationalContext & rc)
                : d_num(i)
            {
                setOne(d_denom);
            }
            
            /**
               \brief Creates a rational number from the given native integer.
               
               \param i The native integer.
             */
            inline explicit Rational(long long i)
                : d_num(i)
            {
                setOne(d_denom);
            }
            
            /**
               \brief Creates a rational number from the given native integer.
               
               \param i The native integer.
               \param rc A rational context.
             */
            inline explicit Rational(long long i, const RationalContext & rc)
                : d_num(i)
            {
                setOne(d_denom);
            }
            
            /**
               \brief Creates a rational number from the given arbitrary precision integer pair of
                      numerator and denominator.
               
               \param n The numerator.
               \param d The denominator.
             */
            inline explicit Rational(const Integer & n, const Integer & d)
                : d_num(n), d_denom(d)
            {
                if (sign(d) < 0)
                {
                    neg(d_denom, d_denom);
                    neg(d_num, d_num);
                }
                normalize();
            }
            
            /**
               \brief Creates a rational number from the given arbitrary precision integer
                      expression.
               
               \param i The arbitrary precision integer.
             */
            inline explicit Rational(const Integer & i)
                : d_num(i)
            {
                setOne(d_denom);
            }
            
            template<class Data, template<typename, typename> class Op>
            inline explicit Rational(const expressions::Expression<IntegerContext, Data, Op> & i)
                : d_num(i)
            {
                setOne(d_denom);
            }
            
            /**
               \brief Creates a rational number from the given arbitrary precision integer.
               
               \param i The arbitrary precision integer.
               \param rc A rational context.
             */
            inline explicit Rational(const Integer & i, const RationalContext & rc)
                : d_num(i)
            {
                setOne(d_denom);
            }
            
            template<class Data, template<typename, typename> class Op>
            inline explicit Rational(const expressions::Expression<IntegerContext, Data, Op> & i, const RationalContext & rc)
                : d_num(i)
            {
                setOne(d_denom);
            }
            
            /**
               \brief Assigns the given rational number `r` to this rational number.
               
               \param r The rational number to assign to this rational number.
               \return A reference to this rational number.
             */
            inline Rational & operator = (const Rational & r)
            {
                d_num = r.d_num;
                d_denom = r.d_denom;
                return *this;
            }
            
            /**
               \brief Assigns the rational expression `E` to this rational number.
               
               \param E The rational expression.
               \return A reference to this rational number.
             */
            template<class A, template<typename, typename> class O>
            inline Rational & operator = (const expressions::Expression<RationalContext, A, O> & E)
            {
                E.assignTo(*this);
                return *this;
            }
            
            /**
               \brief Assigns the integer expression `E` to this rational number.
               
               \param E The integer expression.
               \return A reference to this rational number.
             */
            template<class A, template<typename, typename> class O>
            inline Rational & operator = (const expressions::Expression<IntegerContext, A, O> & E)
            {
                E.assignTo(d_num);
                setOne(d_denom);
                return *this;
            }
            
            inline friend bool isZero(const Rational & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return isZero(r.d_num);
            }
            
            inline friend bool isOne(const Rational & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return isOne(r.d_num) && isOne(r.d_denom);
            }
            
            inline friend bool isPositive(const Rational & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return isPositive(r.d_num);
            }
            
            inline friend bool isNonNegative(const Rational & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return isNonNegative(r.d_num);
            }
            
            inline friend bool isNegative(const Rational & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return isNegative(r.d_num);
            }
            
            inline friend bool isNonPositive(const Rational & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return isNonPositive(r.d_denom);
            }
            
        private:
            inline void add_r_neq_1st(const Rational & a, const Rational & b)
            // Adds a and b and stores result in *this. Assumes that this != &a.
            {
                Integer g1;
                GCD(g1, a.d_denom, b.d_denom);
                if (isOne(g1))
                {
                    d_num = b.d_num * a.d_denom;
                    d_num += a.d_num * b.d_denom;
                    d_denom = a.d_denom * b.d_denom;
                }
                else
                {
                    Integer t = a.d_denom / g1;
                    d_num = b.d_num * t;
                    t = b.d_denom / g1;
                    d_num += a.d_num * t;
                    d_denom = a.d_denom * t;
                }
                normalize();
            }
            
            inline void sub_r_neq_1st(const Rational & a, const Rational & b, bool compute_negative)
            // Subtracts b from a and stores result in *this. Assumes that this != &a.
            {
                Integer g1;
                GCD(g1, a.d_denom, b.d_denom);
                if (isOne(g1))
                {
                    d_num = b.d_num * a.d_denom;
                    d_num -= a.d_num * b.d_denom;
                    d_denom = a.d_denom * b.d_denom;
                }
                else
                {
                    Integer t = a.d_denom / g1;
                    d_num = b.d_num * t;
                    t = b.d_denom / g1;
                    d_num -= a.d_num * t;
                    d_denom = a.d_denom * t;
                }
                if (!compute_negative)
                    d_num = -d_num;
                normalize();
            }
            
        public:
            inline friend void add(Rational & r, const Rational & a, const Rational & b)
            {
                if (&r == &a)
                {
                    if (&r == &b)
                    {
                        // r == a == b, i.e. multiply r by 2
                        if (bit(r.d_denom, 0))
                            r.d_denom >>= 1;
                        else
                            r.d_num <<= 1;
                    }
                    else
                        r.add_r_neq_1st(b, a);
                }
                else
                    r.add_r_neq_1st(a, b);
            }
            
            template<typename IntIntExpr>
            inline friend void add_z(Rational & r, const Rational & a, const IntIntExpr & b)
            {
                if (&r == &a)
                    r.d_num += b * r.d_denom;
                else
                {
                    r.d_num = a.d_num + b * a.d_denom;
                    r.d_denom = a.d_denom;
                }
            }
            
            inline friend void sub(Rational & r, const Rational & a, const Rational & b)
            {
                if (&r == &a)
                {
                    if (&r == &b)
                    {
                        // r == a == b, i.e. set r to zero
                        setZero(r.d_num);
                        setOne(r.d_denom);
                    }
                    else
                        r.sub_r_neq_1st(b, a, true);
                }
                else
                    r.sub_r_neq_1st(a, b, false);
            }
            
            template<typename IntIntExpr>
            inline friend void sub_z(Rational & r, const Rational & a, const IntIntExpr & b)
            {
                if (&r == &a)
                    r.d_num -= b * r.d_denom;
                else
                {
                    r.d_num = a.d_num;
                    r.d_num -= b * a.d_denom;
                    r.d_denom = a.d_denom;
                }
            }
            
            template<typename IntIntExpr>
            inline friend void z_sub(Rational & r, const IntIntExpr & a, const Rational & b)
            {
                if (&r == &b)
                {
                    r.d_num -= a * r.d_denom;
                    neg(r.d_num, r.d_num);
                }
                else
                {
                    r.d_num = a * b.d_denom;
                    r.d_num -= b.d_num;
                    r.d_denom = b.d_denom;
                }
            }
            
            inline friend void addmul(Rational & r, const Rational & a, const Rational & b);
            inline friend void submul(Rational & r, const Rational & a, const Rational & b);
            
            template<typename IntIntExpr>
            inline friend void addmul_z(Rational & r, const Rational & a, const IntIntExpr & b)
            {
                if (&r == &a)
                {
                    r.d_denom *= r.d_num;
                    increment(r.d_num, b);
                    r.d_num = b + convert(1u, IntegerContext());
                    r.d_num *= r.d_denom;
                }
                else
                {
                    Integer t = r.d_denom * a.d_num;
                    t *= b;
                    r.d_num *= a.d_denom;
                    r.d_num += t * b;
                    r.d_denom *= a.d_denom;
                }
            }
            
            template<typename IntIntExpr>
            inline friend void submul_z(Rational & r, const Rational & a, const IntIntExpr & b)
            {
                if (&r == &a)
                {
                    r.d_denom *= r.d_num;
                    decrement(r.d_num, b);
                    neg(r.d_num, r.d_num);
                    r.d_num = b + convert(1u, IntegerContext());
                    r.d_num *= r.d_denom;
                }
                else
                {
                    Integer t = r.d_denom * a.d_num;
                    t *= b;
                    r.d_num *= a.d_denom;
                    r.d_num -= t * b;
                    r.d_denom *= a.d_denom;
                }
            }
            
            inline friend void mul(Rational & r, const Rational & a, const Rational & b)
            {
                if (isZero(a.d_num) || isZero(b.d_num))
                {
                    setZero(r.d_num);
                    setOne(r.d_denom);
                }
                else
                {
                    Integer g1, g2, t;
                    GCD(g1, a.d_num, b.d_denom);
                    GCD(g2, b.d_num, a.d_denom);
                    if (isOne(g1))
                    {
                        if (isOne(g2))
                        {
                            r.d_num = a.d_num * b.d_num;
                            r.d_denom = a.d_denom * b.d_denom;
                        }
                        else
                        {
                            t = b.d_num / g2;
                            r.d_num = a.d_num * t;
                            t = a.d_denom / g2;
                            r.d_denom = b.d_denom * t;
                        }
                    }
                    else
                    {
                        if (isOne(g2))
                        {
                            t = a.d_num / g1;
                            r.d_num = b.d_num * t;
                            t = b.d_denom / g1;
                            r.d_denom = a.d_denom * t;
                        }
                        else
                        {
                            t = b.d_num / g2;
                            r.d_num = a.d_num / g1;
                            r.d_num *= t;
                            t = a.d_denom / g2;
                            r.d_denom = b.d_denom / g1;
                            r.d_denom *= t;
                        }
                    }
//                r.normalize();  -- the above ensures that numerator and denominator are coprime
                }
            }
            
            template<typename IntIntExpr>
            inline friend void mul_z(Rational & r, const Rational & a, const IntIntExpr & b)
            {
                if (isZero(a.d_num))
                {
                    setZero(r.d_num);
                    setOne(r.d_denom);
                }
                else
                {
                    Integer g;
                    GCD(g, a.d_denom, b);
                    if (isOne(g))
                        r.d_num = a.d_num * b;
                    else
                    {
                        g = b / g;
                        r.d_num = a.d_num * g;
                    }
                    r.d_denom = a.d_denom;
                }
            }
            
            inline friend void div(Rational & r, const Rational & a, const Rational & b)
            {
                if (&r == &b)
                {
                    if (&r == &a)
                    {
                        // Result is 1 (except if r == a == b == 0, but we ignore this case)
                        setOne(r);
                        return;
                    }
                    // r == a, but r != b
                    swap(r.d_num, r.d_denom);
                    if (sign(r.d_denom) < 0)
                        r.d_denom = -r.d_denom;
                    r.d_num *= a.d_num;
                    r.d_denom *= a.d_denom;
                }
                else
                {
                    if (&r == &a)
                    {
                        // r is not b, but a
                        if (sign(b.d_num) < 0)
                        {
                            r.d_num *= b.d_denom;
                            r.d_num = -r.d_num;
                            r.d_denom *= b.d_num;
                            r.d_denom = -r.d_denom;
                        }
                        else
                        {
                            r.d_num *= b.d_denom;
                            r.d_denom *= b.d_num;
                        }
                    }
                    else
                    {
                        // r is not b and not a
                        if (sign(b.d_num) < 0)
                        {
                            r.d_num = a.d_num * b.d_denom;
                            r.d_num = -r.d_num;
                            r.d_denom = a.d_denom * b.d_num;
                            r.d_denom = -r.d_denom;
                        }
                        else
                        {
                            r.d_num = a.d_num * b.d_denom;
                            r.d_denom = a.d_denom * b.d_num;
                        }
                    }
                }
                r.normalize();
            }
            
            template<typename IntIntExpr>
            inline friend void div_z(Rational & r, const Rational & a, const IntIntExpr & b)
            {
                if (isZero(a.d_num))
                {
                    setZero(r.d_num);
                    setOne(r.d_denom);
                }
                else
                {
                    Integer g;
                    GCD(g, a.d_num, b);
                    if (isOne(g))
                        r.d_denom = a.d_denom * b;
                    else
                    {
                        g = b / g;
                        r.d_denom = a.d_denom * g;
                    }
                    r.d_num = a.d_num;
                }
            }
            
            template<typename IntIntExpr>
            inline friend void z_div(Rational & r, const IntIntExpr & a, const Rational & b)
            {
                div_z(r, b, a);
                swap(r.d_num, r.d_denom);
                if (sign(r.d_denom) < 0)
                {
                    neg(r.d_num, r.d_num);
                    makeAbs(r.d_denom);
                }
            }
            
            inline friend void mod(Rational & r, const Rational & a, const Rational & b)
            {
                arithmetic::mod(r.d_num, a.d_num * b.d_denom, b.d_num * a.d_denom);
                r.d_denom = a.d_denom * b.d_denom;
                r.normalize();
            }
            
            template<typename IntIntExpr>
            inline friend void mod_z(Rational & r, const Rational & a, const IntIntExpr & b)
            {
                arithmetic::mod(r.d_num, a.d_num, b * a.d_denom);
                r.d_denom = a.d_denom;
                r.normalize();
            }
            
            template<typename IntIntExpr>
            inline friend void z_mod(Rational & r, const IntIntExpr & a, const Rational & b)
            {
                arithmetic::mod(r.d_num, a * b.d_denom, b.d_num);
                r.d_denom = b.d_denom;
                r.normalize();
            }
            
            inline friend void shl(Rational & r, const Rational & a, long b)
            {
                if (b > 0)
                {
                    r.d_num = a.d_num << b;
                    r.d_denom = a.d_denom;
                }
                else
                {
                    r.d_num = a.d_num;
                    r.d_denom = a.d_denom << (-b);
                }
                r.normalize();
            }
            
            template<typename IntIntExpr>
            inline friend void shl_z(Rational & r, const IntIntExpr & a, long b)
            {
                if (b > 0)
                {
                    r.d_num = a << b;
                    setOne(r.d_denom);
                }
                else
                {
                    r.d_num = a;
                    setOne(r.d_denom);
                    r.d_denom <<= -b;
                    r.normalize();
                }
            }
            
            inline friend void shr(Rational & r, const Rational & a, long b)
            {
                if (b > 0)
                {
                    r.d_num = a.d_num;
                    r.d_denom = a.d_denom << b;
                }
                else
                {
                    r.d_num = a.d_num << (-b);
                    r.d_denom = a.d_denom;
                }
                r.normalize();
            }
            
            template<typename IntIntExpr>
            inline friend void shr_z(Rational & r, const IntIntExpr & a, long b)
            {
                if (b < 0)
                {
                    r.d_num = a.d_num << -b;
                    setOne(r.d_denom);
                }
                else
                {
                    r.d_num = a.d_num;
                    setOne(r.d_denom);
                    r.d_denom <<= b;
                    r.normalize();
                }
            }
            
            inline friend void increment(Rational & r, const Rational & a)
            {
                r.d_num = a.d_num + a.d_denom;
                r.d_denom = a.d_denom;
                r.normalize();
            }
            
            inline friend void decrement(Rational & r, const Rational & a)
            {
                r.d_num = a.d_num - a.d_denom;
                r.d_denom = a.d_denom;
                r.normalize();
            }
            
            inline friend void neg(Rational & r, const Rational & a)
            {
                neg(r.d_num, a.d_num);
                r.d_denom = a.d_denom;
            }
            
            template<typename IntIntExpr>
            inline friend void neg_z(Rational & r, const IntIntExpr & a)
            {
                neg(r.d_num, a);
                setOne(r.d_denom);
            }
            
            inline friend void abs(Rational & r, const Rational & a)
            {
                abs(r.d_num, a.d_num);
                r.d_denom = a.d_denom;
            }
            
            template<typename IntIntExpr>
            inline friend void abs_z(Rational & r, const IntIntExpr & a)
            {
                abs(r.d_num, a);
                setOne(r.d_denom);
            }
            
            inline friend void makeAbs(Rational & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                makeAbs(a.d_num);
            }
            
            friend std::ostream & operator << (std::ostream & s, const Rational & r);
            friend std::istream & operator >> (std::istream & s, Rational & r);
            
            inline friend void setZero(Rational & r, bool sign)
            // Set to zero
            {
                setZero(r.d_num);
                setOne(r.d_denom);
            }
            
            inline friend void setOne(Rational & r)
            // Set to one
            {
                setOne(r.d_num);
                setOne(r.d_denom);
            }
            
            inline friend int compare(const Rational & a, const Rational & b)
            {
                return compare(a.d_num * b.d_denom, b.d_num * a.d_denom);
            }
            
            template<typename IntIntExpr>
            inline friend int compare_z(const Rational & a, const IntIntExpr & b)
            {
                return compare(a.d_num, b * a.d_denom);
            }
            
            inline friend int compareAbsValues(const Rational & a, const Rational & b)
            {
                return compareAbsValues(a.d_num * b.d_denom, b.d_num * a.d_denom);
            }
            
            template<typename IntIntExpr>
            inline friend int compareAbsValues_z(const Rational & a, const IntIntExpr & b)
            {
                return compareAbsValues(a.d_num, b * a.d_denom);
            }
            
            inline friend int sign(const Rational & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            // Returns the sign (i.e. -1, 0, 1).
            {
                return sign(r.d_num);
            }
            
            inline friend void square(Rational & r, const Rational & a)
            {
                square(r.d_num, a.d_num);
                square(r.d_denom, a.d_denom);
            }
            
            template<typename IntIntExpr>
            inline friend void square_z(Rational & r, const IntIntExpr & a)
            {
                square(r.d_num, a);
                setOne(r.d_denom);
            }
            
            inline friend void power(Rational & r, const Rational & a, signed long e)
            {
                power(r, a, convert(e, IntegerContext()));
            }
            
            template<typename IntIntExpr>
            inline friend void power_z(Rational & r, const IntIntExpr & a, signed long e)
            {
                if (e >= 0)
                {
                    power(r.d_num, a, e);
                    setOne(r.d_denom);
                }
                else
                {
                    power(r.d_denom, a, -e);
                    setOne(r.d_num);
                }
            }
            
            inline friend void power(Rational & r, const Rational & a, unsigned long e)
            {
                power(r, a, convert(e, IntegerContext()));
            }
            
            template<typename IntIntExpr>
            inline friend void power_z(Rational & r, const IntIntExpr & a, unsigned long e)
            {
                power(r.d_num, a, e);
                setOne(r.d_denom);
            }
            
            friend void power(Rational & r, const Rational & a, const Integer & e);
            
            template<typename IntIntExpr>
            inline friend void power_z(Rational & r, const IntIntExpr & a, const Integer & e)
            {
                if (isNonNegative(e))
                {
                    power(r.d_num, a, e);
                    setOne(r.d_denom);
                }
                else
                {
                    power(r.d_denom, a, -e);
                    setOne(r.d_num);
                }
            }
            
            inline friend void arithmetic::swap(Rational &, Rational &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            
            /**
               \brief Returns approximate `e` such that |x| is approximately \f$2^e\f$.
               
               The value of `e` should be OK up to +-2.
               
               \return \f$\approx \log_2 |x|\f$.
             */
            inline long getApproxExponent() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return isZero(d_num) ? std::numeric_limits<long>::min() : (ceilOfLog2(d_num) - ceilOfLog2(d_denom));
            }
        };
        
        template<typename IntIntExpr>
        inline void add_z(Rational & r, const Rational & a, const IntIntExpr & b);
        
        template<typename IntIntExpr>
        inline void sub_z(Rational & r, const Rational & a, const IntIntExpr & b);
        
        template<typename IntIntExpr>
        inline void z_sub(Rational & r, const IntIntExpr & a, const Rational & b);
        
        template<typename IntIntExpr>
        inline void addmul_z(Rational & r, const Rational & a, const IntIntExpr & b);
        
        template<typename IntIntExpr>
        inline void submul_z(Rational & r, const Rational & a, const IntIntExpr & b);
        
        template<typename IntIntExpr>
        inline void mul_z(Rational & r, const Rational & a, const IntIntExpr & b);
        
        template<typename IntIntExpr>
        inline void div_z(Rational & r, const Rational & a, const IntIntExpr & b);
        
        template<typename IntIntExpr>
        inline void z_div(Rational & r, const IntIntExpr & a, const Rational & b);
        
        template<typename IntIntExpr>
        inline void mod_z(Rational & r, const Rational & a, const IntIntExpr & b);
        
        template<typename IntIntExpr>
        inline void z_mod(Rational & r, const IntIntExpr & a, const Rational & b);
        
        template<typename IntIntExpr>
        inline void shl_z(Rational & r, const IntIntExpr & a, long b);
        
        template<typename IntIntExpr>
        inline void shr_z(Rational & r, const IntIntExpr & a, long b);
        
        template<typename IntIntExpr>
        inline void neg_z(Rational & r, const IntIntExpr & a);
        
        template<typename IntIntExpr>
        inline void abs_z(Rational & r, const IntIntExpr & a);
        
        template<typename IntIntExpr>
        inline int compare_z(const Rational & a, const IntIntExpr & b);
        
        template<typename IntIntExpr>
        inline int compareAbsValues_z(const Rational & a, const IntIntExpr & b);
        
        template<typename IntIntExpr>
        inline void square_z(Rational & r, const IntIntExpr & a);
        
        template<typename IntIntExpr>
        inline void power_z(Rational & r, const IntIntExpr & a, signed long e);
        
        template<typename IntIntExpr>
        inline void power_z(Rational & r, const IntIntExpr & a, unsigned long e);
        
        template<typename IntIntExpr>
        inline void power_z(Rational & r, const IntIntExpr & a, const Integer & e);
        
        inline void swap(Rational & a, Rational & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            swap(a.d_num, b.d_num);
            swap(a.d_denom, b.d_denom);
        }
    }
}

#include "rational-ops.hpp"
#include "rational-conv.hpp"

#endif

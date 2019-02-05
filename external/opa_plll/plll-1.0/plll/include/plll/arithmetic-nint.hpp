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

#ifndef PLLL_INCLUDE_GUARD__NINT_WRAPPER
#define PLLL_INCLUDE_GUARD__NINT_WRAPPER

/**
   \file
   \brief Header for native integer arithmetic in the `plll` arithmetic context framework.
   
   This header provides a wrapper for native integer types such as `int`, `long int` and `long
   long`. It provides the arithmetic context template `plll::arithmetic::NIntContext<>` as described
   in \ref arithctxts.
   
   Concrete integers are represented by the `plll::arithmetic::NInt<>` template.
*/

#include <sstream>
#include <limits>
#include <cmath>
#include <plll/arithmetic.hpp>

namespace plll
{
    namespace arithmetic
    {
        template<typename Type>
        // Type should be long, long long, int, ....
        class NInt;
        // NInt stands for "Native Integers" (i.e. something native to the CPU/FPU)
    }
}

namespace plll
{
    namespace arithmetic
    {
        /**
           \brief Represents an arithmetic context for native integers.
           
           \tparam IType A native signed integer type, such as `int`, `long int` and `long long`.
           
           \sa \ref arithctxts_contexts.
         */
        template<class IType>
        class NIntContext
        // Emulates IntegerContext
        {
        public:
            /**
               \brief The integer type.
             */
            typedef arithmetic::NInt<IType> Integer;
            /**
               \brief The integer type.
             */
            typedef arithmetic::NInt<IType> Type;
            
            /**
               \brief The properties of native CPU integers.
             */
            enum { is_cputype = true, is_realtype = false, is_inttype = true, is_exact = true,
                   is_modulo = true, has_infinity = false, has_uniform_rng = true };
            
            /**
               \brief A uniform random number generator frontend.
               
               \sa \ref arithctxts_urng
             */
            class UniformRNG
            {
            private:
                RandomNumberGenerator & d_rng;
                
            public:
                /**
                   \brief Creates a new uniform random number generator based on the given random
                          number generator.
                */
                UniformRNG(RandomNumberGenerator & rng) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    : d_rng(rng)
                {
                }
                
                /**
                   \brief Creates a random native CPU integer in the range \f$[0, bound)\f$.
                   
                   \param res Where to store the result.
                   \param bound A bound on the maximal integer returned.
                */
                inline void random(NInt<IType> & res, const NInt<IType> & bound);
                
                /**
                   \brief Creates and returns a random native CPU integer in the range \f$[0, bound)\f$.
                   
                   \param bound A bound on the maximal integer returned.
                   \param ic An integer context.
                   \return The result.
                */
                inline NInt<IType> random(const NInt<IType> & bound, const NIntContext<IType> & ic)
                {
                    NInt<IType> res;
                    random(res, bound);
                    return res;
                }
                
                /**
                   \brief Creates random bits.
                   
                   \param res Will be filled with a random integer in range \f$[0, 2^{bits})\f$.
                   \param bits The number of bits.
                */
                inline void randomBits(NInt<IType> & res, unsigned long bits);
                
                /**
                   \brief Creates random bits.
                   
                   \param bits The number of bits.
                   \param ic An integer context.
                   \return A random integer in range \f$[0, 2^{bits})\f$.
                */
                inline NInt<IType> randomBits(unsigned long bits, const NIntContext<IType> & ic)
                {
                    NInt<IType> res;
                    randomBits(res, bits);
                    return res;
                }
                
                /**
                   \brief Creates a random integer of a fixed bit length.
                   
                   \param res Will be filled with a random integer in range \f$[2^{bits-1}, 2^{bits})\f$.
                   \param bits The number of bits.
                */
                inline void randomLen(NInt<IType> & res, unsigned long bits);
                
                /**
                   \brief Creates a random integer of a fixed bit length.
                   
                   \param bits The number of bits.
                   \param ic An integer context.
                   \return A random integer in range \f$[2^{bits-1}, 2^{bits})\f$.
                */
                inline NInt<IType> randomLen(unsigned long bits, const NIntContext<IType> & ic)
                {
                    NInt<IType> res;
                    randomLen(res, bits);
                    return res;
                }
            };
        };
        
        namespace traits
        {
            template<typename Type>
            struct type_traits<NInt<Type> >
            {
                enum
                {
                    is_number = true,
                    is_cputype = true,
                    is_realtype = false,
                    is_inttype = true,
                    is_string = false,
                    is_cpp_string = false,
                    is_c_string = false,
                    is_exact = true,
                    has_uniform_rng = true,
                    has_context = true,
                    is_native = false,
                    is_modulo = true,
                    has_infinity = false,
                    is_variable_precision = false,
                    has_squareroot = false,
                    has_full_power = false,
                    has_special_fns = false,
                    has_huge_exponent = false,
                    has_constants = false,
                    has_trigonometric = false
                };
                
                typedef NIntContext<Type> Context;
                typedef NInt<typename type_traits<Type>::PromoteType> PromoteType;
                typedef NInt<Type> ConstReferenceType;
            };
        }
        
        /**@{
           \name Comparisons.
        */
        /**
           \brief Compares the current integer with the given one for equality.
           \param a The first operand.
           \param b The second operand.
           \return `true` if `a` equals `b`.
        */
        template<typename Type>
        inline bool operator == (const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Compares the current integer with the given one for inequality.
           \param a The first operand.
           \param b The second operand.
           \return `true` if `a` does not equal `b`.
        */
        template<typename Type>
        inline bool operator != (const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Compares the current integer with the given one.
           \param a The first operand.
           \param b The second operand.
           \return `true` if `a` is less than or equal to `b`.
        */
        template<typename Type>
        inline bool operator <= (const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Compares the current integer with the given one.
           \param a The first operand.
           \param b The second operand.
           \return `true` if `a` is greater than or equal to `b`.
        */
        template<typename Type>
        inline bool operator >= (const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Compares the current integer with the given one.
           \param a The first operand.
           \param b The second operand.
           \return `true` if `a` is less than `b`.
        */
        template<typename Type>
        inline bool operator < (const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Compares the current integer with the given one.
           \param a The first operand.
           \param b The second operand.
           \return `true` if `a` is greater than `b`.
        */
        template<typename Type>
        inline bool operator > (const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        ///@}
        /**@{
           \name Operators.
        */
        /**
           \brief Negates the integer.
           
           \param a The integer.
           \return The result.
        */
        template<typename Type>
        inline NInt<Type> operator - (const NInt<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Adds the two integers and returns the result.
           
           \param a The first operand.
           \param b The second operand.
           \return The result.
        */
        template<typename Type>
        inline NInt<Type> operator + (const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Subtracts the second from the first integer and returns the result.
           
           \param a The first operand.
           \param b The second operand.
           \return The result.
        */
        template<typename Type>
        inline NInt<Type> operator - (const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Multiplies the two integers and returns the result.
           
           \param a The first operand.
           \param b The second operand.
           \return The result.
        */
        template<typename Type>
        inline NInt<Type> operator * (const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Divides the first by the second integer and returns the result.
           
           \param a The first operand.
           \param b The second operand.
           \return The result.
        */
        template<typename Type>
        inline NInt<Type> operator / (const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Divides the first by the second integer and returns the remainder.
           
           \param a The first operand.
           \param b The second operand.
           \return The result.
        */
        template<typename Type>
        inline NInt<Type> operator % (const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Computes the bitwise left shift of the first integer by the number of bits given
                  by the second integer, and returns the result.
           
           This effectively multiplies the first integer by 2 to the power of the second integer.
           
           \param a The first operand.
           \param b The second operand.
           \return The result.
        */
        template<typename Type>
        inline NInt<Type> operator << (const NInt<Type> & a, long b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Computes the bitwise right shift of the first integer by the number of bits given
                  by the second integer, and returns the result.
           
           This effectively divides the first integer by 2 to the power of the second integer, and
           rounds towards zero.
           
           \param a The first operand.
           \param b The second operand.
           \return The result.
        */
        template<typename Type>
        inline NInt<Type> operator >> (const NInt<Type> & a, long b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Computes the bitwise left shift of the first integer by the number of bits given
                  by the second integer, and returns the result.
           
           This effectively multiplies the first integer by 2 to the power of the second integer.
           
           \param a The first operand.
           \param b The second operand.
           \return The result.
        */
        template<typename Type>
        inline NInt<Type> operator << (const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Computes the bitwise right shift of the first integer by the number of bits given
                  by the second integer, and returns the result.
           
           This effectively divides the first integer by 2 to the power of the second integer, and
           rounds towards zero.
           
           \param a The first operand.
           \param b The second operand.
           \return The result.
        */
        template<typename Type>
        inline NInt<Type> operator >> (const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Increments the integer by one and returns the previous value.
           
           \param cur The integer to work on.
           \return The previous value.
        */
        template<typename Type>
        inline NInt<Type> operator ++ (NInt<Type> & cur, int) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Decrements the integer by one and returns the previous value.
           
           \param cur The integer to work on.
           \return The previous value.
        */
        template<typename Type>
        inline NInt<Type> operator -- (NInt<Type> & cur, int) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Increments the integer by one and returns the new value.
           
           \param cur The integer to work on.
           \return The new value (a reference to `cur`).
        */
        template<typename Type>
        inline NInt<Type> & operator ++ (NInt<Type> & cur) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Decreases the integer by one and returns the new value.
           
           \param cur The integer to work on.
           \return The new value (a reference to `cur`).
        */
        template<typename Type>
        inline NInt<Type> & operator -- (NInt<Type> & cur) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        ///@}
        
        /**@{
           \name Assignment-operation operators.
        */
        /**
           \brief Subtracts the integer `i` from `cur`.
           
           \param cur The integer to operate on.
           \param i The integer to subtract from `cur`.
           \return A reference to `cur` containing the result.
        */
        template<typename Type>
        inline NInt<Type> & operator -= (NInt<Type> & cur, const NInt<Type> & i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Adds the integer `i` to `cur`.
           
           \param cur The integer to operate on.
           \param i The integer to add to `cur`.
           \return A reference to `cur` containing the result.
        */
        template<typename Type>
        inline NInt<Type> & operator += (NInt<Type> & cur, const NInt<Type> & i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Multiplies the given integer `i` with `cur`.
               
           \param cur The integer to operate on.
           \param i The integer to multiply to `cur`.
           \return A reference to `cur` containing the result.
        */
        template<typename Type>
        inline NInt<Type> & operator *= (NInt<Type> & cur, const NInt<Type> & i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Divides `cur` by the given integer `i`.
           
           \param cur The integer to operate on.
           \param i The integer to divide by.
           \return A reference to `cur` containing the result.
        */
        template<typename Type>
        inline NInt<Type> & operator /= (NInt<Type> & cur, const NInt<Type> & i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Divides `cur` by the given integer `i` and
                  stores the remainder in `cur`.
           
           \param cur The integer to operate on.
           \param i The integer to divide by.
           \return A reference to `cur` containing the result.
        */
        template<typename Type>
        inline NInt<Type> & operator %= (NInt<Type> & cur, const NInt<Type> & i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Computes the bitwise left shift of `cur` by `i` bits, and stores the result in
                  `cur`.
           
           This effectively multiplies `cur` by 2 to the power of `i`.
           
           \param cur The integer to operate on.
           \param i The second operand.
           \return A reference to `cur` containing the result.
        */
        template<typename Type>
        inline NInt<Type> & operator <<= (NInt<Type> & cur, long i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Computes the bitwise right shift of `cur` by `i` bits, and stores the result in
                  `cur`.
           
           This effectively divides `cur` by 2 to the power of `i`, and rounds towards zero.
           
           \param cur The integer to operate on.
           \param i The second operand.
           \return A reference to `cur` containing the result.
        */
        template<typename Type>
        inline NInt<Type> & operator >>= (NInt<Type> & cur, long i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Computes the bitwise left shift of `cur` by `i` bits, and stores the result in
                  `cur`.
           
           This effectively multiplies `cur` by 2 to the power of `i`.
           
           \param cur The integer to operate on.
           \param i The second operand.
           \return A reference to `cur` containing the result.
        */
        template<typename Type>
        inline NInt<Type> & operator <<= (NInt<Type> & cur, const NInt<Type> & i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Computes the bitwise right shift of `cur` by `i` bits, and stores the result in
                  `cur`.
           
           This effectively divides `cur` by 2 to the power of `i`, and rounds towards zero.
           
           \param cur The integer to operate on.
           \param i The second operand.
           \return A reference to `cur` containing the result.
        */
        template<typename Type>
        inline NInt<Type> & operator >>= (NInt<Type> & cur, const NInt<Type> & i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        ///@}
        
        /**@{
           \name Predicates.
        */
        /**
           \brief Tests the given `plll::arithmetic::Integer` object for being zero.
           
           \return Returns `true` if and only if the argument is zero.
           \sa \ref arithctxts_preds
        */
        template<typename Type>
        inline bool isZero(const NInt<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        
        /**
           \brief Tests the given `plll::arithmetic::Integer` object for being one.
           
           \return Returns `true` if and only if the argument is one.
           \sa \ref arithctxts_preds
        */
        template<typename Type>
        inline bool isOne(const NInt<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        
        /**
           \brief Tests the given `plll::arithmetic::Integer` object for being one or minus one.
           
           \return Returns `true` if and only if the argument is \f$\pm 1\f$.
           \sa \ref arithctxts_preds
        */
        template<typename Type>
        inline bool isPMOne(const NInt<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        
        /**
           \brief Tests the given `plll::arithmetic::Integer` object for being two or minus two.
           
           \return Returns `true` if and only if the argument is \f$\pm 2\f$.
           \sa \ref arithctxts_preds
        */
        template<typename Type>
        inline bool isPMTwo(const NInt<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        
        /**
           \brief Tests the given `plll::arithmetic::Integer` object for being strictly positive.
           
           \return Returns `true` if and only if the argument is strictly positive.
           \sa \ref arithctxts_preds
        */
        template<typename Type>
        inline bool isPositive(const NInt<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        
        /**
           \brief Tests the given `plll::arithmetic::Integer` object for being positive or zero.
           
           \return Returns `true` if and only if the argument is positive or zero.
           \sa \ref arithctxts_preds
        */
        template<typename Type>
        inline bool isNonNegative(const NInt<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        
        /**
           \brief Tests the given `plll::arithmetic::Integer` object for being strictly negative.
           
           \return Returns `true` if and only if the argument is strictly negative.
           \sa \ref arithctxts_preds
        */
        template<typename Type>
        inline bool isNegative(const NInt<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        
        /**
           \brief Tests the given `plll::arithmetic::Integer` object for being negative or zero.
           
           \return Returns `true` if and only if the argument is negative or zero.
           \sa \ref arithctxts_preds
        */
        template<typename Type>
        inline bool isNonPositive(const NInt<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
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
        template<typename Type>
        inline void add(NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Subtracts `b` from `a` and stores the result in `r`.
           
           \param r The result.
           \param a The first operand.
           \param b The second operand.
         */
        template<typename Type>
        inline void sub(NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Multiplies `a` and `b` and adds the result to `r`.
           
           \param r The accumulator.
           \param a The first operand.
           \param b The second operand.
         */
        template<typename Type>
        inline void addmul(NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Multiplies `a` and `b` and subtracts the result from `r`.
           
           \param r The accumulator.
           \param a The first operand.
           \param b The second operand.
         */
        template<typename Type>
        inline void submul(NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Multiplies `a` with `b` and stores the result in `r`.
           
           \param r The result.
           \param a The first operand.
           \param b The second operand.
         */
        template<typename Type>
        inline void mul(NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Divides `a` by `b` and stores the result in `r`.
           
           \param r The result.
           \param a The first operand.
           \param b The second operand.
         */
        template<typename Type>
        inline void div(NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Stores quotient and remainder of the division of `a` by `b` in `q` respectively
                  `r`.
           
           \param q The quotient.
           \param r The remainder.
           \param a The first operand.
           \param b The second operand.
         */
        template<typename Type>
        inline void divmod(NInt<Type> & q, NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Takes the remainder of the division of `a` by `b` and stores it in `r`.
           
           \param r The result.
           \param a The first operand.
           \param b The second operand.
         */
        template<typename Type>
        inline void mod(NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Shifts `a` by `b` bits to the left and stores the result in `r`.
           
           \param r The result.
           \param a The first operand.
           \param b The second operand.
         */
        template<typename Type>
        inline void shl(NInt<Type> & r, const NInt<Type> & a, long b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Shifts `a` by `b` bits to the left and stores the result in `r`.
           
           \param r The result.
           \param a The first operand.
           \param b The second operand.
         */
        template<typename Type>
        inline void shr(NInt<Type> & r, const NInt<Type> & a, long b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Shifts `a` by `b` bits to the left and stores the result in `r`.
           
           \param r The result.
           \param a The first operand.
           \param b The second operand.
         */
        template<typename Type>
        inline void shl(NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Shifts `a` by `b` bits to the left and stores the result in `r`.
           
           \param r The result.
           \param a The first operand.
           \param b The second operand.
         */
        template<typename Type>
        inline void shr(NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Increments `a` by one and stores the result in `r`.
           
           \param r The result.
           \param a The operand.
           \sa \ref arithctxts_ints
         */
        template<typename Type>
        inline void increment(NInt<Type> & r, const NInt<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Decrements `a` by one and stores the result in `r`.
           
           \param r The result.
           \param a The operand.
         */
        template<typename Type>
        inline void decrement(NInt<Type> & r, const NInt<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Negates `a` and stores the result in `r`.
           
           \param r The result.
           \param a The operand.
         */
        template<typename Type>
        inline void neg(NInt<Type> & r, const NInt<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Takes the absolute value of `a` and stores the result in `r`.
           
           \param r The result.
           \param a The operand.
         */
        template<typename Type>
        inline void abs(NInt<Type> & r, const NInt<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Computes the square of `a` and stores the result in `r`.
           
           \param r The result.
           \param a The operand.
         */
        template<typename Type>
        inline void square(NInt<Type> & r, const NInt<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        ///@}
        
        /**@{
           \name Operators.
        */
        /**
           \brief Computes and returns the absolute value of `i`.
           \param i The operand.
           \return The result.
        */
        template<typename Type>
        inline NInt<Type> abs(const NInt<Type> & i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Computes and returns the square of `i`.
           \param i The operand.
           \return The result.
        */
        template<typename Type>
        inline NInt<Type> square(const NInt<Type> & i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        ///*}
        
        /**@{
           \name Stream input/output.
        */
        /**
           \brief Outputs the integer on the given output stream.
         */
        template<typename Type>
        std::ostream & operator << (std::ostream &, const NInt<Type> &);
        
        /**
           \brief Reads the integer from the given input stream.
         */
        template<typename Type>
        std::istream & operator >> (std::istream &, NInt<Type> &);
        ///@}
        
        /**@{
           \name Setting to specific constants.
        */
        /**
           \brief Sets the given integer to zero.
         */
        template<typename Type>
        inline void setZero(NInt<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Sets the given integer to one.
         */
        template<typename Type>
        inline void setOne(NInt<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
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
        template<typename Type>
        inline int compare(const NInt<Type> &, const NInt<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Compares the two integers in absolute value.

           \param a The first operand.
           \param b The second operand.
           \return Returns a negative number if \f$|a| < |b|\f$, zero if \f$|a| = |b|\f$ and a
                   positive number if \f$|a| > |b|\f$.
         */
        template<typename Type>
        inline int compareAbsValues(const NInt<Type> &, const NInt<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        ///@}
        
        /**@{
           \name Sign querying/modification.
        */
        /**
           \brief Returns the sign of the given integer.
           
           \return Returns a negative value if the value is negative, 0 if it is zero, and a
                   positive value if it is positive.
         */
        template<typename Type>
        inline int sign(const NInt<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Makes the operand non-negative.
           
           \param a The operand.
         */
        template<typename Type>
        inline void makeAbs(NInt<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        ///@}
        /**@{
           \name Bit manipulation.
        */
        /**
           \brief Returns the `n` bit of \f$|x|\f$ in the usual binary representation.
           \param x The integer whose bit to query.
           \param n The index of the bit to query.
         */
        template<typename Type>
        inline int bit(const NInt<Type> & x, long n) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Sets the `n`-th bit of \f$|x|\f$ to `value`.
           \param x The integer whose bits to modify.
           \param n The index of the bit to set or clear.
           \param value The new value of the `n`-th bit. The default value is `true`.
           \sa \ref arithctxts_ints
         */
        template<typename Type>
        inline void setbit(NInt<Type> & x, long n, bool value = true) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        ///@}
        
        /**@{
           \name Exponentiation.
        */
        /**
           \brief Computes and returns `a` raised to the power of `b`.
           \param a The first operand.
           \param b The second operand.
           \return The result.
        */
        template<typename Type>
        inline NInt<Type> power(const NInt<Type> & a, long b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Raises `a` to the power `b` and stores the result in `r`.
           
           \param r The result.
           \param a The base.
           \param b The exponent.
         */
        template<typename Type>
        inline void power(NInt<Type> & r, const NInt<Type> & a, long b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Computes and returns `a` raised to the power of `b`.
           \param a The first operand.
           \param b The second operand.
           \return The result.
        */
        template<typename Type>
        inline NInt<Type> power(const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Raises `a` to the power `b` and stores the result in `r`.
           
           \param r The result.
           \param a The base.
           \param b The exponent.
         */
        template<typename Type>
        inline void power(NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        ///@}
        
        /**@{
           \name Integer approximation.
        */
        /**
           \brief Computes \f$\lceil\sqrt{a}\rceil\f$ and stores the result in `r`.
           
           \param r The result.
           \param a The operand.
         */
        template<typename Type>
        inline void sqrtCeil(NInt<Type> & r, const NInt<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Computes and returns \f$\lceil\sqrt{a}\rceil\f$.
           
           \return The result.
           \param a The operand.
         */
        template<typename Type>
        inline NInt<Type> sqrtCeil(const NInt<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Computes \f$\lfloor\sqrt{a}\rfloor\f$ and stores the result in `r`.
           
           \param r The result.
           \param a The operand.
         */
        template<typename Type>
        inline void sqrtFloor(NInt<Type> & r, const NInt<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Computes and returns \f$\lfloor\sqrt{a}\rfloor\f$.
           
           \return The result.
           \param a The operand.
         */
        template<typename Type>
        inline NInt<Type> sqrtFloor(const NInt<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Computes and returns \f$\lceil \log_2 |x| \rceil\f$.
           
           \param x A non-zero integer.
           \return \f$n \in \mathbb{N}\f$ such that `\f$2^{n-1} < |x| \le 2^n\f$.
           \sa \ref arithctxts_ints
         */
        template<typename Type>
        long ceilOfLog2(const NInt<Type> & x) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE; // Returns   ceil(log_2(|x|)).
        /**
           \brief Computes and returns \f$\lfloor \log_2 |x| \rfloor\f$.
           
           \param x A non-zero integer.
           \return \f$n \in \mathbb{N}\f$ such that `\f$2^n \le |x| < 2^{n+1}\f$.
           \sa \ref arithctxts_ints
         */
        template<typename Type>
        long floorOfLog2(const NInt<Type> & x) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE; // Returns   floor(log_2(|x|)).
        /**
           \brief Quickly approximates \f$\log_2 |x|\f$ and returns the approximation.
           
           \param x A non-zero integer.
           \return \f$\approx \log_2 |x|\f$.
         */
        template<typename Type>
        long approxLog2(const NInt<Type> & x) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE; // Returns   ~log_2(|x|).
        /**
           \brief Computes and returns `n` such that \f$2^{n-1} \le |x| < 2^n\f$.
           
           \param x A non-zero integer.
           \return \f$n \in \mathbb{N}\f$ such that `\f$2^{n-1} \le |x| < 2^n\f$.
           \sa \ref arithctxts_ints
         */
        template<typename Type>
        inline long bitLength(const NInt<Type> & x) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE; // Returns n such that 2^{n-1} <= |x| < 2^n
        /**
           \brief Computes \f$\lfloor \tfrac{a}{b} \rfloor\f$ and stores the result in `r`.
           
           \param r The result.
           \param a The divident.
           \param b The divisor.
           \sa \ref arithctxts_ints
         */
        template<typename Type>
        void floorDiv(NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE; // Computes r = floor(a/b), where b \neq 0.
        /**
           \brief Computes and returns \f$\lfloor \tfrac{a}{b} \rfloor\f$.
           
           \return The result.
           \param a The divident.
           \param b The divisor.
           \sa \ref arithctxts_ints
         */
        template<typename Type>
        NInt<Type> floorDiv(const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE; // Computes r = floor(a/b), where b \neq 0.
        /**
           \brief Computes \f$\lceil \tfrac{a}{b} \rceil\f$ and stores the result in `r`.
           
           \param r The result.
           \param a The divident.
           \param b The divisor.
           \sa \ref arithctxts_ints
         */
        template<typename Type>
        void ceilDiv(NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE; // Computes r = ceil(a/b), where b \neq 0.
        /**
           \brief Computes and returns \f$\lceil \tfrac{a}{b} \rceil\f$.
           
           \return The result.
           \param a The divident.
           \param b The divisor.
           \sa \ref arithctxts_ints
         */
        template<typename Type>
        NInt<Type> ceilDiv(const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE; // Computes r = ceil(a/b), where b \neq 0.
        /**
           \brief Computes \f$\lfloor \tfrac{a}{b} \rceil\f$ (rounding to the next integer) and
                  stores the result in `r`.
           
           \param r The result.
           \param a The divident.
           \param b The divisor.
           \sa \ref arithctxts_ints
         */
        template<typename Type>
        void roundDiv(NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE; // Computes r = round(a/b), where b \neq 0.
        /**
           \brief Computes and returns \f$\lfloor \tfrac{a}{b} \rceil\f$ (rounding to the next
                  integer).
           
           \return The result.
           \param a The divident.
           \param b The divisor.
           \sa \ref arithctxts_ints
         */
        template<typename Type>
        NInt<Type> roundDiv(const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE; // Computes r = round(a/b), where b \neq 0.
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
        template<typename Type>
        inline void euclideanDivision(NInt<Type> & q, NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        
        /**
           \brief Computes an Euclidean Division of `a` by `b`.
           
           \param q The quotient of `a` divided by `b`.
           \param r The remainder of `a` divided by `b`, i.e. `a` modulo `b`.
           \param a The first operand.
           \param b The second operand. Must be non-zero.
           
           Computes `q` and `r` such that `a == q * b + r` and that \f$0 \le r \le |b|\f$.
         */
        template<typename Type>
        inline void euclideanDivisionPos(NInt<Type> & q, NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        
        /**
           \brief Computes the non-negative Greatest Common Divisior `r` of `x` and `y`.
           
           \param r The result is stored in here.
           \param x The first operand.
           \param y The first operand.
           \sa \ref arithctxts_ints
         */
        template<typename Type>
        inline void GCD(NInt<Type> & r, const NInt<Type> & x, const NInt<Type> & y) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Computes and returns the non-negative Greatest Common Divisior of `x` and `y`.
           
           \return The GCD of `x` and `y`.
           \param x The first operand.
           \param y The first operand.
           \sa \ref arithctxts_ints
         */
        template<typename Type>
        inline NInt<Type> GCD(const NInt<Type> & x, const NInt<Type> & y) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        
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
        template<typename Type>
        inline void XGCD(NInt<Type> & r, NInt<Type> & a, NInt<Type> & b, const NInt<Type> & x, const NInt<Type> & y) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        
        /**
           \brief Computes the non-negative Least Common Multiple `r` of `x` and `y`.
           
           \param r The result is stored in here.
           \param x The first operand.
           \param y The first operand.
           \sa \ref arithctxts_ints
         */
        template<typename Type>
        inline void LCM(NInt<Type> & r, const NInt<Type> & x, const NInt<Type> & y) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        /**
           \brief Computes and returns the non-negative Least Common Multiple of `x` and `y`.
           
           \return The LCM of `x` and `y`.
           \param x The first operand.
           \param y The first operand.
           \sa \ref arithctxts_ints
         */
        template<typename Type>
        inline NInt<Type> LCM(const NInt<Type> & x, const NInt<Type> & y) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        ///@}
        
        /**@{
           \name Swap functions.
        */
        /**
           \brief Swaps two `plll::arithmetic::NInt<>` objects.
         */
        template<typename Type>
        inline void swap(NInt<Type> &, NInt<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        ///@}
        
        /**
           \brief Represents a native integer.
           
           \tparam Type A native signed integer type, such as `int`, `long int` and `long long`.
         */
        template<typename Type>
        class NInt
        {
#ifdef PLLL_INTERNAL_NO_TEMPLATE_FRIENDS
        public:
#else
        private:
            friend class NIntContext<Type>;
            
            template<class X, class Y>
            friend class implementation::conversion_impl;
            
            friend class implementation::nativeconversion_impl<NInt<Type> >;
#endif
            
            Type d_value;
            
            inline NInt(bool, Type v) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_value(v)
            {
            }
            
            inline void assignType(Type v) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                d_value = v;
            }
            
        public:
            /**
               \brief The context type.
             */
            typedef NIntContext<Type> Context;
            
            /**
               \brief Creates a new integer. Default value is zero.
             */
            inline NInt() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_value(0)
            {
            }
            
            /**
               \brief Creates a new integer. Default value is zero.
               
               \param c A native CPU integer context.
             */
            inline NInt(const NIntContext<Type> & c) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_value(0)
            {
            }
            
            /**
               \brief Creates a copy of the given integer.
               
               \param i The integer to be copied.
             */
            inline NInt(const NInt & i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_value(i.d_value)
            {
            }
            
            /**
               \brief Creates a copy of the given integer (of another native type).
               
               \param i The integer to be copied.
             */
            template<typename T>
            inline NInt(const NInt<T> & i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_value(i.d_value)
            {
            }
            
            /**
               \brief Creates a copy of the given integer.
               
               \param i The integer to be copied.
               \param ic An integer context.
             */
            inline NInt(const NInt & i, const NIntContext<Type> & ic) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_value(i.d_value)
            {
            }
            
            /**
               \brief Creates a copy of the given integer (of another native type).
               
               \param i The integer to be copied.
               \param ic An integer context.
             */
            template<typename T>
            inline NInt(const NInt<T> & i, const NIntContext<Type> & ic) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_value(i.d_value)
            {
            }
            
            /**
               \brief Sets the integer context `c`.
               
               \param c The integer context.
             */
            static inline void setContext(const NIntContext<Type> & c) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
            }
            
            /**
               \brief Creates a native CPU integer from the given native floating point number.
               
               \param d The native floating point number.
             */
            inline explicit NInt(double d) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_value(d)
            {
            }
        
            /**
               \brief Creates a native CPU integer from the given native floating point number.
               
               \param d The native floating point number.
             */
            inline explicit NInt(long double d) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_value(d)
            {
            }
            
            /**
               \brief Creates a native CPU integer from the given native integer.
               
               \param i The native integer.
             */
            inline explicit NInt(long i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_value(i)
            {
            }
            
            /**
               \brief Creates a native CPU integer from the given native integer.
               
               \param i The native integer.
             */
            inline explicit NInt(unsigned long i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_value(i)
            {
            }
            
            /**
               \brief Creates a native CPU integer from the given arbitrary precision integer.
               
               \param i The native integer.
             */
            inline explicit NInt(const Integer & i)
                : d_value(convert<Type>(i))
            {
            }
            
            /**
               \brief Copies the given integer into this one.
               
               \param r A integer to be moved into this one.
             */
            inline NInt & operator = (const NInt & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                d_value = r.d_value;
                return *this;
            }
            
#ifndef PLLL_INTERNAL_NO_TEMPLATE_FRIENDS
            friend NInt<Type> operator - <>(const NInt<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NInt<Type> operator + <>(const NInt<Type> &, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NInt<Type> operator - <>(const NInt<Type> &, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NInt<Type> operator * <>(const NInt<Type> &, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NInt<Type> operator / <>(const NInt<Type> &, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NInt<Type> operator % <>(const NInt<Type> &, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NInt<Type> operator << <>(const NInt<Type> & a, long b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NInt<Type> operator >> <>(const NInt<Type> & a, long b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NInt<Type> operator << <>(const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NInt<Type> operator >> <>(const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NInt<Type> & operator -= <>(NInt<Type> &, const NInt<Type> & i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NInt<Type> & operator += <>(NInt<Type> &, const NInt<Type> & i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NInt<Type> & operator *= <>(NInt<Type> &, const NInt<Type> & i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NInt<Type> & operator /= <>(NInt<Type> &, const NInt<Type> & i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NInt<Type> & operator %= <>(NInt<Type> &, const NInt<Type> & i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NInt<Type> operator ++ <>(NInt<Type> &, int) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NInt<Type> operator -- <>(NInt<Type> &, int) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NInt<Type> & operator ++ <>(NInt<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NInt<Type> & operator -- <>(NInt<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend bool operator == <>(const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend bool operator != <>(const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend bool operator <= <>(const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend bool operator >= <>(const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend bool operator < <>(const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend bool operator > <>(const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend bool isZero<>(const NInt<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend bool isOne<>(const NInt<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend bool isPMOne<>(const NInt<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend bool isPMTwo<>(const NInt<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend bool isPositive<>(const NInt<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend bool isNonNegative<>(const NInt<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend bool isNegative<>(const NInt<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend bool isNonPositive<>(const NInt<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void add<>(NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void sub<>(NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void mul<>(NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void div<>(NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void mod<>(NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void shl<>(NInt<Type> & r, const NInt<Type> & a, long b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void shr<>(NInt<Type> & r, const NInt<Type> & a, long b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void shl<>(NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void shr<>(NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void increment<>(NInt<Type> & r, const NInt<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void decrement<>(NInt<Type> & r, const NInt<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void neg<>(NInt<Type> & r, const NInt<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void abs<>(NInt<Type> & r, const NInt<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NInt<Type> abs<>(const NInt<Type> & i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void makeAbs<>(NInt<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend std::ostream & operator << <>(std::ostream & s, const NInt<Type> & r);
            // Output to stream
            friend std::istream & operator >> <>(std::istream & s, NInt<Type> & r);
            // Input from stream
            friend void setZero<>(NInt<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            // Set to zero
            friend void setOne<>(NInt<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            // Set to one
            friend int compare<>(const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            // Tests whether the first real is < (-1), = (0) or > (1) than the second.
            friend int compareAbsValues<>(const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            // Tests whether the absolute value of the first NInt<Type> is < (-1), = (0) or > (1) than the
            // absolute value of the second NInt.
            friend int sign<>(const NInt<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            // Returns the sign (i.e. -1, 0, 1).
            friend int bit<>(const NInt<Type> &, long n) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NInt<Type> square<>(const NInt<Type> & i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void square<>(NInt<Type> & r, const NInt<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NInt<Type> power<>(const NInt<Type> &, long) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void power<>(NInt<Type> &, const NInt<Type> &, long) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NInt<Type> power<>(const NInt<Type> &, const NInt<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void power<>(NInt<Type> &, const NInt<Type> &, const NInt<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void sqrtCeil<>(NInt<Type> &, const NInt<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void sqrtFloor<>(NInt<Type> &, const NInt<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend long ceilOfLog2<>(const NInt<Type> & x) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend long floorOfLog2<>(const NInt<Type> & x) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend long approxLog2<>(const NInt<Type> & x) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend long bitLength<>(const NInt<Type> & x) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void floorDiv<>(NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NInt<Type> floorDiv<>(const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void ceilDiv<>(NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NInt<Type> ceilDiv<>(const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void roundDiv<>(NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NInt<Type> roundDiv<>(const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void euclideanDivision<>(NInt<Type> & q, NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void euclideanDivisionPos<>(NInt<Type> & q, NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void GCD<>(NInt<Type> & r, const NInt<Type> & x, const NInt<Type> & y) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NInt<Type> GCD<>(const NInt<Type> & x, const NInt<Type> & y) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void XGCD<>(NInt<Type> & r, NInt<Type> & a, NInt<Type> & b, const NInt<Type> & x, const NInt<Type> & y) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void LCM<>(NInt<Type> & r, const NInt<Type> & x, const NInt<Type> & y) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NInt<Type> LCM<>(const NInt<Type> & x, const NInt<Type> & y) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void swap<>(NInt<Type> & a, NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
#endif
        };
        
        template<typename Type>
        inline NInt<Type> operator - (const NInt<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NInt<Type>(true, -a.d_value);
        }
        
        template<typename Type>
        inline NInt<Type> operator + (const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NInt<Type>(true, a.d_value + b.d_value);
        }
        
        template<typename Type>
        inline NInt<Type> operator - (const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NInt<Type>(true, a.d_value - b.d_value);
        }
        
        template<typename Type>
        inline NInt<Type> operator * (const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NInt<Type>(true, a.d_value * b.d_value);
        }
        
        template<typename Type>
        inline NInt<Type> operator / (const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NInt<Type>(true, a.d_value / b.d_value);
        }
        
        template<typename Type>
        inline NInt<Type> operator % (const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NInt<Type>(true, a.d_value % b.d_value);
        }
        
        template<typename Type>
        inline NInt<Type> & operator -= (NInt<Type> & a, const NInt<Type> & i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            a.d_value -= i.d_value;
            return a;
        }
        
        template<typename Type>
        inline NInt<Type> & operator += (NInt<Type> & a, const NInt<Type> & i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            a.d_value += i.d_value;
            return a;
        }
        
        template<typename Type>
        inline NInt<Type> & operator *= (NInt<Type> & a, const NInt<Type> & i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            a.d_value *= i.d_value;
            return a;
        }
        
        template<typename Type>
        inline NInt<Type> & operator /= (NInt<Type> & a, const NInt<Type> & i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            a.d_value /= i.d_value;
            return a;
        }
        
        template<typename Type>
        inline NInt<Type> & operator %= (NInt<Type> & a, const NInt<Type> & i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            a.d_value %= i.d_value;
            return a;
        }
        
        template<typename Type>
        inline NInt<Type> & operator <<= (NInt<Type> & a, long r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            shl(a, a, r);
            return a;
        }
        
        template<typename Type>
        inline NInt<Type> & operator >>= (NInt<Type> & a, long r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            shr(a, a, r);
            return a;
        }
        
        template<typename Type>
        inline NInt<Type> & operator <<= (NInt<Type> & a, const NInt<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            shl(a, a, r);
            return a;
        }
        
        template<typename Type>
        inline NInt<Type> & operator >>= (NInt<Type> & a, const NInt<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            shr(a, a, r);
            return a;
        }
        
        template<typename Type>
        inline NInt<Type> operator ++ (NInt<Type> & a, int) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NInt<Type>(true, a.d_value++);
        }
        
        template<typename Type>
        inline NInt<Type> operator -- (NInt<Type> & a, int) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NInt<Type>(true, a.d_value--);
        }
        
        template<typename Type>
        inline NInt<Type> & operator ++ (NInt<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            ++a.d_value;
            return a;
        }
        
        template<typename Type>
        inline NInt<Type> & operator -- (NInt<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            --a.d_value;
            return a;
        }
        
        template<typename Type>
        inline NInt<Type> operator << (const NInt<Type> & a, long b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            NInt<Type> r;
            shl(r, a, b);
            return r;
        }
        
        template<typename Type>
        inline NInt<Type> operator >> (const NInt<Type> & a, long b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            NInt<Type> r;
            shr(r, a, b);
            return r;
        }
        
        template<typename Type>
        inline NInt<Type> operator << (const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            NInt<Type> r;
            shl(r, a, b);
            return r;
        }
        
        template<typename Type>
        inline NInt<Type> operator >> (const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            NInt<Type> r;
            shr(r, a, b);
            return r;
        }
        
        template<typename Type>
        inline bool operator == (const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return a.d_value == b.d_value;
        }
        
        template<typename Type>
        inline bool operator != (const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return a.d_value != b.d_value;
        }
        
        template<typename Type>
        inline bool operator <= (const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return a.d_value <= b.d_value;
        }
        
        template<typename Type>
        inline bool operator >= (const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return a.d_value >= b.d_value;
        }
        
        template<typename Type>
        inline bool operator < (const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return a.d_value < b.d_value;
        }
        
        template<typename Type>
        inline bool operator > (const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return a.d_value > b.d_value;
        }
        
        template<typename Type>
        inline bool isZero(const NInt<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Tests whether the given value is = 0.
        {
            return r.d_value == (Type)0;
        }
        
        template<typename Type>
        inline bool isOne(const NInt<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Tests whether the given value is = 1.
        {
            return r.d_value == (Type)1;
        }
        
        template<typename Type>
        inline bool isPMOne(const NInt<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Tests whether the given value is = +-1.
        {
            return (r.d_value == (Type)1) || (r.d_value == (Type)(-1));
        }
        
        template<typename Type>
        inline bool isPMTwo(const NInt<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Tests whether the given value is = +-2.
        {
            return (r.d_value == (Type)2) || (r.d_value == (Type)(-2));
        }
        
        template<typename Type>
        inline bool isPositive(const NInt<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Tests whether the given value is > 0.
        {
            return r.d_value > (Type)0;
        }
        
        template<typename Type>
        inline bool isNonNegative(const NInt<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Tests whether the given value is >= 0.
        {
            return r.d_value >= (Type)0;
        }
        
        template<typename Type>
        inline bool isNegative(const NInt<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Tests whether the given value is < 0.
        {
            return r.d_value < (Type)0;
        }
        
        template<typename Type>
        inline bool isNonPositive(const NInt<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Tests whether the given value is <= 0.
        {
            return r.d_value <= (Type)0;
        }
        
        template<typename Type>
        inline void add(NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = a.d_value + b.d_value;
        }
        
        template<typename Type>
        inline void sub(NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = a.d_value - b.d_value;
        }
        
        template<typename Type>
        inline void addmul(NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value += a.d_value * b.d_value;
        }
        
        template<typename Type>
        inline void submul(NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value -= a.d_value * b.d_value;
        }
        
        template<typename Type>
        inline void mul(NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = a.d_value * b.d_value;
        }
        
        template<typename Type>
        inline void div(NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = a.d_value / b.d_value;
        }
        
        template<typename Type>
        inline void mod(NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = a.d_value % b.d_value;
        }
        
        template<typename Type>
        inline void divmod(NInt<Type> & q, NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            Type qq = a.d_value / b.d_value, rr = a.d_value % b.d_value;
            q.d_value = qq;
            r.d_value = rr;
        }
        
        template<typename Type>
        inline void shl(NInt<Type> & r, const NInt<Type> & a, long b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            if (b < 0)
            {
                if (static_cast<unsigned long>(-b) >= static_cast<unsigned long>(std::numeric_limits<Type>::digits))
                    r.d_value = 0;
                else
                    r.d_value = a.d_value >> static_cast<unsigned long>(-b);
            }
            else
            {
                if (b >= std::numeric_limits<Type>::digits)
                    r.d_value = 0;
                else
                    r.d_value = a.d_value << b;
            }
        }
        
        template<typename Type>
        inline void shr(NInt<Type> & r, const NInt<Type> & a, long b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            if (b < 0)
            {
                if (static_cast<unsigned long>(-b) >= static_cast<unsigned long>(std::numeric_limits<Type>::digits))
                    r.d_value = 0;
                else
                    r.d_value = a.d_value << static_cast<unsigned long>(-b);
            }
            else
            {
                if (b >= std::numeric_limits<Type>::digits)
                    r.d_value = 0;
                else
                    r.d_value = a.d_value >> b;
            }
        }
        
        template<typename Type>
        inline void shl(NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = a.d_value << b.d_value;
        }
        
        template<typename Type>
        inline void shr(NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            if (b.d_value >= std::numeric_limits<Type>::digits)
                r.d_value = 0;
            else
                r.d_value = a.d_value >> b.d_value;
        }
        
        template<typename Type>
        inline void increment(NInt<Type> & r, const NInt<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = a.d_value + (Type)1;
        }
        
        template<typename Type>
        inline void decrement(NInt<Type> & r, const NInt<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = a.d_value - (Type)1;
        }
        
        template<typename Type>
        inline void neg(NInt<Type> & r, const NInt<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = -a.d_value;
        }
        
        template<typename Type>
        inline void abs(NInt<Type> & r, const NInt<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = a.d_value < 0 ? -a.d_value : a.d_value;
        }
        
        template<typename Type>
        inline NInt<Type> abs(const NInt<Type> & i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NInt<Type>(true, i.d_value < 0 ? -i.d_value : i.d_value);
        }
        
        template<typename Type>
        inline void makeAbs(NInt<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            if (a.d_value < 0)
                a.d_value = -a.d_value;
        }
        
        template<typename Type>
        std::ostream & operator << (std::ostream & s, const NInt<Type> & r)
        // Output to stream
        {
            return s << r.d_value;
        }
        
        template<typename Type>
        std::istream & operator >> (std::istream & s, NInt<Type> & r)
        // Input from stream
        {
            return s >> r.d_value;
        }
        
        template<typename Type>
        inline void setZero(NInt<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Set to zero
        {
            r.d_value = (Type)0;
        }
        
        template<typename Type>
        inline void setOne(NInt<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Set to one
        {
            r.d_value = (Type)1;
        }
        
        template<typename Type>
        inline int compare(const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Tests whether the first real is < (-1), = (0) or > (1) than the second.
        {
            return a.d_value < b.d_value ? -1 : (a.d_value > b.d_value ? 1 : 0);
        }
        
        template<typename Type>
        inline int compareAbsValues(const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Tests whether the absolute value of the first NInt is < (-1), = (0) or > (1) than the
        // absolute value of the second NInt.
        {
            Type aabs = a.d_value < 0 ? -a.d_value : a.d_value, babs = b.d_value < 0 ? -b.d_value : b.d_value;
            return aabs < babs ? -1 : (aabs > babs ? 1 : 0);
        }
        
        template<typename Type>
        inline int sign(const NInt<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Returns the sign (i.e. -1, 0, 1).
        {
            return r.d_value < 0 ? -1 : (r.d_value > 0 ? 1 : 0);
        }
        
        namespace implementation
        {
            template<typename Type>
            struct get_unsigned;
            
            template<>
            struct get_unsigned<int>
            {
                typedef unsigned int result;
            };
            
            template<>
            struct get_unsigned<long>
            {
                typedef unsigned long result;
            };
            
            template<>
            struct get_unsigned<long long>
            {
                typedef unsigned long long result;
            };
        }
        
        template<typename Type>
        inline int bit(const NInt<Type> & v, long n) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Returns the n-th bit of the given integer.
        {
            return v.d_value < 0
                               ? ((static_cast<typename implementation::get_unsigned<Type>::result>(-v.d_value) >> (Type)n) & 1)
                                    : ((v.d_value >> (Type)n) & 1);
        }
        
        template<typename Type>
        inline void setbit(NInt<Type> & v, long n, bool value) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Sets the n-th bit of the given integer.
        {
            if (value)
                v.d_value |= (static_cast<Type>(1) << n);
            else
                v.d_value &= ~(static_cast<Type>(1) << n);
        }
        
        template<typename Type>
        inline NInt<Type> square(const NInt<Type> & i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NInt<Type>(true, i.d_value * i.d_value);
        }
        
        template<typename Type>
        inline void square(NInt<Type> & r, const NInt<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = a.d_value * a.d_value;
        }
        
        template<typename Type>
        inline Type powerimpl(Type base, Type exp) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            assert(exp >= 0);
            Type result = (exp & 1) ? base : (Type)1;
            exp >>= 1;
            while (exp)
            {
                base = base * base;
                if (exp & 1)
                    result *= base;
                exp >>= 1;
            }
            return result;
        }
        
        template<typename Type>
        inline void power(NInt<Type> & r, const NInt<Type> & a, long e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = powerimpl(a.d_value, (Type)e);
        }
        
        template<typename Type>
        inline NInt<Type> power(const NInt<Type> & a, long e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NInt<Type>(true, powerimpl(a.d_value, (Type)e));
        }
        
        template<typename Type>
        inline void power(NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = powerimpl(a.d_value, e.d_value);
        }
        
        template<typename Type>
        inline NInt<Type> power(const NInt<Type> & a, const NInt<Type> & e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NInt<Type>(true, powerimpl(a.d_value, e.d_value));
        }
        
        template<typename Type>
        inline void sqrtCeil(NInt<Type> & v, const NInt<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            Type result = 0, rsq = 0, bm = 1;
            int bits = 0;
            while (bm + (bm - 1) < a.d_value)
            {
                bm <<= 2;
                ++bits;
            }
            while (bm)
            {
                Type nsq = rsq + bm + (result << (bits + 1));
                if (nsq <= a.d_value)
                {
                    result += (((Type)1) << bits);
                    rsq = nsq;
                }
                bm >>= 2;
                --bits;
            }
            v.d_value = (rsq < a.d_value) ? result + 1 : result;
        }
    
        template<typename Type>
        inline void sqrtFloor(NInt<Type> & v, const NInt<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            Type result = 0, rsq = 0, bm = 1;
            int bits = 0;
            while (bm + (bm - 1) < a.d_value)
            {
                bm <<= 2;
                ++bits;
            }
            while (bm)
            {
                Type nsq = rsq + bm + (result << (bits + 1));
                if (nsq <= a.d_value)
                {
                    result += (((Type)1) << bits);
                    rsq = nsq;
                }
                bm >>= 2;
                --bits;
            }
            v.d_value = result;
        }
    
        template<typename Type>
        long ceilOfLog2(const NInt<Type> & x) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Returns   ceil(log_2(|x|)).
        {
            Type v = x.d_value < 0 ? -x.d_value : x.d_value;
            bool pot = (v & (-v)) == v;
            unsigned res = 0;
            while (v)
            {
                --res;
                v >>= 1;
            }
            if (!pot)
                ++res;
            return res;
        }
        
        template<typename Type>
        long floorOfLog2(const NInt<Type> & x) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Returns   floor(log_2(|x|)).
        {
            Type v = x.d_value < 0 ? -x.d_value : x.d_value;
            unsigned res = 0;
            while (v)
            {
                --res;
                v >>= 1;
            }
            return res;
        }
        
        namespace implementation
        {
            template<class X> class GetUnsignedType;
            
            template<> class GetUnsignedType<signed short int>
            {
            public:
                typedef signed short int signed_type;
                typedef unsigned short int unsigned_type;
                enum { is_signed = true };
            };
            
            template<> class GetUnsignedType<signed int>
            {
            public:
                typedef signed int signed_type;
                typedef unsigned int unsigned_type;
                enum { is_signed = true };
            };
            
            template<> class GetUnsignedType<signed long>
            {
            public:
                typedef signed long signed_type;
                typedef unsigned long unsigned_type;
                enum { is_signed = true };
            };
            
            template<> class GetUnsignedType<signed long long>
            {
            public:
                typedef signed long long signed_type;
                typedef unsigned long long unsigned_type;
                enum { is_signed = true };
            };
            
            template<> class GetUnsignedType<unsigned short int>
            {
            public:
                typedef signed short int signed_type;
                typedef unsigned short int unsigned_type;
                enum { is_signed = false };
            };
            
            template<> class GetUnsignedType<unsigned int>
            {
            public:
                typedef signed int signed_type;
                typedef unsigned int unsigned_type;
                enum { is_signed = false };
            };
            
            template<> class GetUnsignedType<unsigned long>
            {
            public:
                typedef signed long signed_type;
                typedef unsigned long unsigned_type;
                enum { is_signed = false };
            };
            
            template<> class GetUnsignedType<unsigned long long>
            {
            public:
                typedef signed long long signed_type;
                typedef unsigned long long unsigned_type;
                enum { is_signed = false };
            };
            
            template<class UType>
            inline long approxLog2Impl(UType xx) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE // for unsigned types
            // Returns   ~log_2(xx).
            {
                // Uses the technique described by Sean E. Anderson in
                // http://graphics.stanford.edu/~seander/bithacks.html#IntegerLog. The original version is
                // from Eric Cole and was modified by Andrew Shapira.
                unsigned int res, shift;
                PLLL_INTERNAL_STATIC_CHECK((std::numeric_limits<UType>::digits == 7) ||
                                           (std::numeric_limits<UType>::digits == 8) ||
                                           (std::numeric_limits<UType>::digits == 15) ||
                                           (std::numeric_limits<UType>::digits == 16) ||
                                           (std::numeric_limits<UType>::digits == 31) ||
                                           (std::numeric_limits<UType>::digits == 32) ||
                                           (std::numeric_limits<UType>::digits == 63) ||
                                           (std::numeric_limits<UType>::digits == 64) ||
                                           (std::numeric_limits<UType>::digits == 127) ||
                                           (std::numeric_limits<UType>::digits == 128), RequiresTypeOf_8_16_32_64_128_Bits)
                    switch (std::numeric_limits<UType>::digits) // Note that CLANG seems to return bitsize-1
                        // instead of bitsize.
                    {
                    case 7:
                    case 8:
                    {
                        res =   (xx > (UType)0xF) << 2; xx >>= res;
                        shift = (xx > (UType)0x3) << 1; xx >>= shift; res |= shift;
                        res |= (xx >> (UType)1);
                    }
                    return res;
                    case 15:
                    case 16:
                    {
                        res =   (xx > (UType)0xFF) << 3; xx >>= res;
                        shift = (xx > (UType)0xF ) << 2; xx >>= shift; res |= shift;
                        shift = (xx > (UType)0x3 ) << 1; xx >>= shift; res |= shift;
                        res |= (xx >> (UType)1);
                    }
                    return res;
                    case 31:
                    case 32:
                    {
                        res =   (xx > (UType)0xFFFF) << 4; xx >>= res;
                        shift = (xx > (UType)0xFF  ) << 3; xx >>= shift; res |= shift;
                        shift = (xx > (UType)0xF   ) << 2; xx >>= shift; res |= shift;
                        shift = (xx > (UType)0x3   ) << 1; xx >>= shift; res |= shift;
                        res |= (xx >> (UType)1);
                    }
                    return res;
                    case 63:
                    case 64:
                    {
                        res =   (xx > (UType)0xFFFFFFFFul) << 5; xx >>= res;
                        shift = (xx > (UType)0xFFFF      ) << 4; xx >>= shift; res |= shift;
                        shift = (xx > (UType)0xFF        ) << 3; xx >>= shift; res |= shift;
                        shift = (xx > (UType)0xF         ) << 2; xx >>= shift; res |= shift;
                        shift = (xx > (UType)0x3         ) << 1; xx >>= shift; res |= shift;
                        res |= (xx >> (UType)1);
                    }
                    return res;
                    case 127:
                    case 128:
                    {
                        res =   (xx > (UType)0xFFFFFFFFFFFFFFFFul) << 6; xx >>= res;
                        shift = (xx > (UType)0xFFFFFFFFul        ) << 5; xx >>= shift; res |= shift;
                        shift = (xx > (UType)0xFFFF              ) << 4; xx >>= shift; res |= shift;
                        shift = (xx > (UType)0xFF                ) << 3; xx >>= shift; res |= shift;
                        shift = (xx > (UType)0xF                 ) << 2; xx >>= shift; res |= shift;
                        shift = (xx > (UType)0x3                 ) << 1; xx >>= shift; res |= shift;
                        res |= (xx >> (UType)1);
                    }
                    return res;
                    }
            }
            
            template<typename Type>
            inline long approxLog2_impl(Type x, helper::BoolToType<true>) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE // for signed types
            // Returns   ~log_2(|x|).
            {
                typedef typename GetUnsignedType<Type>::unsigned_type UType;
                return approxLog2Impl(x < 0 ? -(UType)x : (UType)x);
            }
            
            template<typename Type>
            inline long approxLog2_impl(Type x, helper::BoolToType<false>) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE // for unsigned types
            // Returns   ~log_2(|x|).
            {
                return approxLog2Impl(x);
            }
        }
        
        template<typename Type>
        inline long approxLog2(Type x) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE // for signed & unsigned types
        // Returns   ~log_2(|x|).
        {
            return implementation::approxLog2_impl(x, helper::BoolToType<implementation::GetUnsignedType<Type>::is_signed>());
        }
        
        template<typename Type>
        inline long approxLog2(const NInt<Type> & x) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Returns   ~log_2(|x|).
        {
            return approxLog2(x.d_value);
        }
        
        template<typename Type>
        inline long bitLength(const NInt<Type> & x) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Returns n such that 2^{n-1} <= |x| < 2^n
        {
            return isZero(x) ? 0 : floorOfLog2(x) + 1;
        }
        
        template<typename Type>
        void floorDiv(NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Computes r = floor(a/b), where b \neq 0.
        {
            Type aa = a.d_value < 0 ? -a.d_value : a.d_value;
            Type ba = b.d_value < 0 ? -b.d_value : b.d_value;
            bool neg = (b.d_value < 0) ^ (a.d_value < 0);
            if (neg)
                r.d_value = -(aa + (ba - 1)) / ba;
            else
                r.d_value = aa / ba;
        }
        
        template<typename Type>
        NInt<Type> floorDiv(const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            NInt<Type> result;
            floorDiv(result, a, b);
            return result;
        }
        
        template<typename Type>
        void ceilDiv(NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Computes r = ceil(a/b), where b \neq 0.
        {
            Type aa = a.d_value < 0 ? -a.d_value : a.d_value;
            Type ba = b.d_value < 0 ? -b.d_value : b.d_value;
            bool neg = (b.d_value < 0) ^ (a.d_value < 0);
            if (neg)
                r.d_value = -aa / ba;
            else
                r.d_value = (aa + (ba - 1)) / ba;
        }
        
        template<typename Type>
        NInt<Type> ceilDiv(const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            NInt<Type> result;
            ceilDiv(result, a, b);
            return result;
        }
        
        template<typename Type>
        void roundDiv(NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Computes r = round(a/b), where b \neq 0.
        {
            Type aa = a.d_value < 0 ? -a.d_value : a.d_value;
            Type ba = b.d_value < 0 ? -b.d_value : b.d_value;
            bool neg = (b.d_value < 0) ^ (a.d_value < 0);
            if (neg)
                r.d_value = -(aa + ba/2) / ba;
            else
                r.d_value = (aa + ba/2) / ba;
        }
        
        template<typename Type>
        NInt<Type> roundDiv(const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            NInt<Type> result;
            roundDiv(result, a, b);
            return result;
        }
        
        template<typename Type>
        inline void euclideanDivision(NInt<Type> & q, NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Given b \neq 0, computes q, r such that a = q b + r, 0 <= |r| < |b| and b*r >= 0.
        {
            q.d_value = a.d_value / b.d_value;
            r.d_value = a.d_value % b.d_value;
        }
        
        template<typename Type>
        inline void euclideanDivisionPos(NInt<Type> & q, NInt<Type> & r, const NInt<Type> & a, const NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Given b \neq 0, computes q, r such that a = q b + r, 0 <= r < |b|
        {
            Type aa = a.d_value < 0 ? -a.d_value : a.d_value;
            Type ba = b.d_value < 0 ? -b.d_value : b.d_value;
            bool neg = (b.d_value < 0) ^ (a.d_value < 0);
            if (neg)
                q.d_value = -(aa + (ba - 1)) / ba;
            else
                q.d_value = aa / ba;
            r.d_value = a.d_value - q.d_value * b.d_value;
        }
        
        template<typename Type>
        inline void GCD(NInt<Type> & r, const NInt<Type> & x, const NInt<Type> & y) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Computes the gcd of x and y.
        {
            Type a = x.d_value, b = y.d_value;
            if (a < 0)
                a = -a;
            if (b < 0)
                b = -b;
            while (a)
            {
                Type c = b % a;
                b = a;
                a = c;
            }
            r.d_value = b;
        }
        
        template<typename Type>
        NInt<Type> GCD(const NInt<Type> & x, const NInt<Type> & y) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            NInt<Type> result;
            GCD(result, x, y);
            return result;
        }
        
        template<typename Type>
        inline void XGCD(NInt<Type> & r, NInt<Type> & a, NInt<Type> & b, const NInt<Type> & x, const NInt<Type> & y) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Computes the gcd of x and y, as well as a, b such that   gcd = a*x + b*y.
        {
            Type aa = x.d_value, bb = y.d_value;
            Type x1, x2 = 0, y1 = 0, y2; // beginning: |a| = x1 * a + y1 * b, |b| = x2 * a + y2 * b
            if (aa < 0)
            {
                aa = -aa;
                x1 = -1;
            }
            else
                x1 = 1;
            if (bb < 0)
            {
                bb = -bb;
                y2 = -1;
            }
            else
                y2 = 1;
            while (aa)
            {
                Type r = bb % aa, q = bb / aa; // b = q * a + r
                // Have: a = x1  * a_orig + y1  * b_orig, b = x2  * a_orig + y2  * b_orig
                // Want: r = x1' * a_orig + y1' * b_orig, a = x2' * a_orig + y2' * b_orig
                // Use:  r = b - q * a
                // Thus: x1' = q * x2 - x1
                //       y1' = q * y2 - y1
                //       x2' = x1
                //       y2' = y1
                bb = aa;
                aa = r;
                Type x = x2 - q * x1;
                x2 = x1;
                x1 = x;
                Type y = y2 - q * y1;
                y2 = y1;
                y1 = y;
            }
            a.d_value = x2;
            b.d_value = y2;
            r.d_value = bb;
        }
        
        template<typename Type>
        inline void LCM(NInt<Type> & r, const NInt<Type> & x, const NInt<Type> & y) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE // Computes the lcm of x and y.
        {
            GCD(r, x, y);
            r.d_value = x.d_value * (y.d_value / r.d_value);
        }
        
        template<typename Type>
        NInt<Type> LCM(const NInt<Type> & x, const NInt<Type> & y) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            NInt<Type> result;
            LCM(result, x, y);
            return result;
        }
        
        namespace implementation
        {
            void randomInt(int & l, long b, RandomNumberGenerator & rng);
            void randomInt(long & l, long b, RandomNumberGenerator & rng);
            void randomInt(long long & ll, long long b, RandomNumberGenerator & rng);
            void randomBits(int & l, unsigned long b, RandomNumberGenerator & rng);
            void randomBits(long & l, unsigned long b, RandomNumberGenerator & rng);
            void randomBits(long long & ll, unsigned long b, RandomNumberGenerator & rng);
        }
        
        template<class IType>
        inline void NIntContext<IType>::UniformRNG::random(NInt<IType> & res, const NInt<IType> & bound)
        {
            implementation::randomInt(res.d_value, bound.d_value, d_rng);
        }
        
        template<class IType>
        inline void NIntContext<IType>::UniformRNG::randomBits(NInt<IType> & res, unsigned long bits)
        {
            implementation::randomBits(res.d_value, bits, d_rng);
        }
        
        template<class IType>
        inline void NIntContext<IType>::UniformRNG::randomLen(NInt<IType> & res, unsigned long bits)
        {
            if (bits)
            {
                randomBits(res, bits - 1);
                IType b = 1;
                b <<= bits;
                res.d_value |= b;
            }
            else
                setZero(res);
        }
        
        template<typename Type>
        inline void swap(NInt<Type> & a, NInt<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            std::swap(a.d_value, b.d_value);
        }
    }
}

#include "arithmetic-nint-conv.hpp"

#endif

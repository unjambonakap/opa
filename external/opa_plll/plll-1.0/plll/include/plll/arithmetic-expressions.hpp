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

#ifndef PLLL_INCLUDE_GUARD__ARITHMETIC_EXPRESSIONS_HPP
#define PLLL_INCLUDE_GUARD__ARITHMETIC_EXPRESSIONS_HPP

#include <plll/arithmetic.hpp>

/**
   \file
   \brief Main header for arithmetic expressions.
   
   This header provides expression templates for arbitrary precision integers, floating point
   numbers, for rational numbers and for conversions.
 */

namespace plll
{
    namespace arithmetic
    {
        /**
           \brief Expression templates.
           
           This namespace contains all facilities required to implement expression templates for
           arithmetic. It is used to provide arithmetic for the `arithmetic::Integer`,
           `arithmetic::Real` and `arithmetic::Rational` types.
           
           An expression `arithmetic::expressions::Expression<Context, Data, Op>` consists of three
           parts:
           
           - an arithmetic context `Context` describing the result;
           
           - arbitrary data of type `Data`, which describes the operands of the operation and can be
             other expressions, wrappers (see the `expressions::Wrapper<>` template)
           
           - a operation, modelled by a template template parameter `Op`, which is instanciated with
             the context and the data type.
           
           Types using expression templates should provide constructors and assignment operators
           which call the `evaluateTo()` method of the expression; otherwise, it will be evaluated
           by calling the cast operator which returns an object of type `Context::Type`.
         */
        namespace expressions
        {
            ////////////////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////
            ////// WRAPPER / FORWARDER /////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////
            
            template<typename Context_>
            class Wrapper
            /**
               \brief Wraps a variable to behave similarly to an `Expression<>` object.
               
               The wrapper object stores a reference to the original variable.
             */
            {
            private:
                const typename Context_::Type & d_value;
                
            public:
                enum { has_evaluate = true, has_context = false };
                typedef const typename Context_::Type & evaluate_type;
                typedef Context_ Context;
                
                inline Wrapper(const typename Context_::Type & value)
                /**
                   \brief Initializes the wrapper object with the given variable.
                 */
                    : d_value(value)
                {
                }
                
                inline unsigned long precision() const
                /**
                   \brief Returns the precision of the variable.
                   
                   Only available for real types.
                 */
                {
                    return d_value.precision();
                }
                
                inline operator const typename Context_::Type & () const
                /**
                   \brief Cast operator.
                   
                   \return a const reference to the the stored value.
                 */
                {
                    return d_value;
                }
                
                inline void assignTo(typename Context_::Type & x) const
                /**
                   \brief Assigns the value to the destination variable.
                   
                   \param x The varible where to store the value.
                 */
                {
                    x = d_value;
                }
                
                inline const typename Context_::Type & evaluate() const
                /**
                   \brief "Evaluates" (returns) the stored value.
                   
                   \return a const reference to the the stored value.
                 */
                {
                    return d_value;
                }
            };
            
            template<class Context, class Data>
            class NoneOp
            /**
               \brief Simply provides the data without any modifications.
             */
            {
            public:
                enum { has_evaluate = Data::has_evaluate, has_context = Data::has_context };
                typedef typename Data::evaluate_type evaluate_type;
                
                static inline unsigned long precision(const Data & a)
                {
                    return a.precision();
                }
                
                static inline const Context & context(const Data & a)
                {
                    return a.context();
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & a)
                {
                    a.assignTo(x);
                }
                
                static inline evaluate_type evaluate(const Data & a)
                {
                    return a.evaluate();
                }
            };
            
            struct NoData
            /**
               \brief A simple "no data" indicator.
             */
            {
            };
            
            ////////////////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////
            ////// BASIC ARITHMETIC ////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////
            
            template<typename Context, typename Data>
            class NegOp
            /**
               \brief Returns the negative of the operand.
             */
            {
            public:
                enum { has_evaluate = false, has_context = false };
                
                static inline unsigned long precision(const Data & a)
                {
                    return a.precision();
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & a)
                {
                    neg(x, a.evaluate());
                }
            };
            
            template<typename Context, typename Data>
            class AddOp
            /**
               \brief Returns the sum of the operands.
             */
            {
            public:
                enum { has_evaluate = false, has_context = false };
                
                static inline unsigned long precision(const Data & d)
                {
                    return std::max(d.first.precision(), d.second.precision());
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & d)
                {
                    add(x, d.first.evaluate(), d.second.evaluate());
                }
            };
            
            template<typename Context, typename Data>
            class SubOp
            /**
               \brief Returns the difference of the operands.
             */
            {
            public:
                enum { has_evaluate = false, has_context = false };
                
                static inline unsigned long precision(const Data & d)
                {
                    return std::max(d.first.precision(), d.second.precision());
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & d)
                {
                    sub(x, d.first.evaluate(), d.second.evaluate());
                }
            };
            
            template<typename Context, typename Data>
            class MulOp
            /**
               \brief Returns the product of the operands.
             */
            {
            public:
                enum { has_evaluate = false, has_context = false };
                
                static inline unsigned long precision(const Data & d)
                {
                    return std::max(d.first.precision(), d.second.precision());
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & d)
                {
                    mul(x, d.first.evaluate(), d.second.evaluate());
                }
            };
            
            template<typename Context, typename Data>
            class DivOp
            /**
               \brief Returns the quotient of the operands.
             */
            {
            public:
                enum { has_evaluate = false, has_context = false };
                
                static inline unsigned long precision(const Data & d)
                {
                    return std::max(d.first.precision(), d.second.precision());
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & d)
                {
                    div(x, d.first.evaluate(), d.second.evaluate());
                }
            };
            
            template<typename Context, typename Data>
            class ModOp
            /**
               \brief Returns the remainder of the operands.
             */
            {
            public:
                enum { has_evaluate = false, has_context = false };
                
                static inline unsigned long precision(const Data & d)
                {
                    return std::max(d.first.precision(), d.second.precision());
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & d)
                {
                    mod(x, d.first.evaluate(), d.second.evaluate());
                }
            };
            
            template<typename Context, typename Data>
            class ShLOp
            /**
               \brief Returns the first operand multiplied by 2 to the power of the second operand.
             */
            {
            public:
                enum { has_evaluate = false, has_context = false };
                
                static inline unsigned long precision(const Data & d)
                {
                    return std::max(d.first.precision(), d.second.precision());
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & d)
                {
                    shl(x, d.first.evaluate(), d.second.evaluate());
                }
            };
            
            template<typename Context, typename Data>
            class ShROp
            /**
               \brief Returns the first operand divided by 2 to the power of the second operand.
             */
            {
            public:
                enum { has_evaluate = false, has_context = false };
                
                static inline unsigned long precision(const Data & d)
                {
                    return std::max(d.first.precision(), d.second.precision());
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & d)
                {
                    shr(x, d.first.evaluate(), d.second.evaluate());
                }
            };
            
            template<typename Context, typename Data>
            class ShiftCOp
            /**
               \brief Returns the operand multiplied by 2 to the power of a constant.
               
               The constant is stored in the operator object.
             */
            {
            private:
                signed long d_v;
                
            public:
                enum { has_evaluate = false, has_context = false };
                
                signed long shift() const
                {
                    return d_v;
                }
                
                inline ShiftCOp(signed long v)
                : d_v(v)
                {
                }
                
                static inline unsigned long precision(const Data & a)
                {
                    return a.precision();
                }
                
                inline void assignTo(typename Context::Type & x, const Data & a) const
                {
                    shl(x, a.evaluate(), d_v);
                }
            };
            
            template<typename Context, typename Data>
            class AndOp
            /**
               \brief Returns the bitwise AND of the operands.
             */
            {
            public:
                enum { has_evaluate = false, has_context = false };
                
                static inline unsigned long precision(const Data & d)
                {
                    return std::max(d.first.precision(), d.second.precision());
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & d)
                {
                    band(x, d.first.evaluate(), d.second.evaluate());
                }
            };
            
            template<typename Context, typename Data>
            class OrOp
            /**
               \brief Returns the bitwise OR of the operands.
             */
            {
            public:
                enum { has_evaluate = false, has_context = false };
                
                static inline unsigned long precision(const Data & d)
                {
                    return std::max(d.first.precision(), d.second.precision());
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & d)
                {
                    bor(x, d.first.evaluate(), d.second.evaluate());
                }
            };
            
            template<typename Context, typename Data>
            class XOROp
            /**
               \brief Returns the bitwise XOR (exclusive or) of the operands.
             */
            {
            public:
                enum { has_evaluate = false, has_context = false };
                
                static inline unsigned long precision(const Data & d)
                {
                    return std::max(d.first.precision(), d.second.precision());
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & d)
                {
                    bxor(x, d.first.evaluate(), d.second.evaluate());
                }
            };
            
            template<typename Context, typename Data>
            class BitInvOp
            /**
               \brief Returns the bitwise inversion of the operands.
             */
            {
            public:
                enum { has_evaluate = false, has_context = false };
                
                static inline unsigned long precision(const Data & d)
                {
                    return d.precision();
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & d)
                {
                    bneg(x, d.evaluate());
                }
            };
            
            ////////////////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////
            ////// FUNCTIONAL ARITHMETIC W/O CONTEXT ///////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////
            
            template<typename Context, typename Data>
            class SqrtFloorOp
            /**
               \brief Returns the floor of the square root of the operand.
             */
            {
            public:
                enum { has_evaluate = false, has_context = false };
                
                static inline unsigned long precision(const Data & d)
                {
                    return d.precision();
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & d)
                {
                    sqrtFloor(x, d.evaluate());
                }
            };
            
            template<typename Context, typename Data>
            class SqrtCeilOp
            /**
               \brief Returns the ceil of the square root of the operand.
             */
            {
            public:
                enum { has_evaluate = false, has_context = false };
                
                static inline unsigned long precision(const Data & d)
                {
                    return d.precision();
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & d)
                {
                    sqrtCeil(x, d.evaluate());
                }
            };
            
            template<class Context, class Data>
            class GCDOp
            /**
               \brief Returns the GCD of the two operands.
             */
            {
            public:
                enum { has_evaluate = false, has_context = false };
                
                static inline unsigned long precision(const Data & d)
                {
                    return std::max(d.first.precision(), d.second.precision());
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & d)
                {
                    GCD(x, d.first.evaluate(), d.second.evaluate());
                }
            };
            
            template<class Context, class Data>
            class LCMOp
            /**
               \brief Returns the LCM of the two operands.
             */
            {
            public:
                enum { has_evaluate = false, has_context = false };
                
                static inline unsigned long precision(const Data & d)
                {
                    return std::max(d.first.precision(), d.second.precision());
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & d)
                {
                    LCM(x, d.first.evaluate(), d.second.evaluate());
                }
            };
            
            template<class Context, class Data>
            class FloorDivOp
            /**
               \brief Returns the floor of the quotient of the two operands.
             */
            {
            public:
                enum { has_evaluate = false, has_context = false };
                
                static inline unsigned long precision(const Data & d)
                {
                    return std::max(d.first.precision(), d.second.precision());
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & d)
                {
                    floorDiv(x, d.first.evaluate(), d.second.evaluate());
                }
            };
            
            template<class Context, class Data>
            class CeilDivOp
            /**
               \brief Returns the ceil of the quotient of the two operands.
             */
            {
            public:
                enum { has_evaluate = false, has_context = false };
                
                static inline unsigned long precision(const Data & d)
                {
                    return std::max(d.first.precision(), d.second.precision());
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & d)
                {
                    ceilDiv(x, d.first.evaluate(), d.second.evaluate());
                }
            };
            
            template<class Context, class Data>
            class RoundDivOp
            /**
               \brief Returns the rounded quotient of the two operands.
             */
            {
            public:
                enum { has_evaluate = false, has_context = false };
                
                static inline unsigned long precision(const Data & d)
                {
                    return std::max(d.first.precision(), d.second.precision());
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & d)
                {
                    roundDiv(x, d.first.evaluate(), d.second.evaluate());
                }
            };
            
            template<typename Context, typename Data>
            class AbsOp
            /**
               \brief Returns the absolute value of the operand.
             */
            {
            public:
                enum { has_evaluate = false, has_context = false };
                
                static inline unsigned long precision(const Data & a)
                {
                    return a.precision();
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & a)
                {
                    abs(x, a.evaluate());
                }
            };
            
            template<typename Context, typename Data>
            class SquareOp
            /**
               \brief Returns the square of the operand.
             */
            {
            public:
                enum { has_evaluate = false, has_context = false };
                
                static inline unsigned long precision(const Data & a)
                {
                    return a.precision();
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & a)
                {
                    square(x, a.evaluate());
                }
            };
            
            template<class Context, class Data>
            class SinOp
            /**
               \brief Returns the sine of the operand.
             */
            {
            public:
                enum { has_evaluate = false, has_context = false };
                
                static inline unsigned long precision(const Data & a)
                {
                    return a.precision();
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & a)
                {
                    sin(x, a.evaluate());
                }
            };
            
            template<class Context, class Data>
            class CosOp
            /**
               \brief Returns the cosine of the operand.
             */
            {
            public:
                enum { has_evaluate = false, has_context = false };
                
                static inline unsigned long precision(const Data & a)
                {
                    return a.precision();
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & a)
                {
                    cos(x, a.evaluate());
                }
            };
            
            template<class Context, class Data>
            class TanOp
            /**
               \brief Returns the tangent of the operand.
             */
            {
            public:
                enum { has_evaluate = false, has_context = false };
                
                static inline unsigned long precision(const Data & a)
                {
                    return a.precision();
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & a)
                {
                    tan(x, a.evaluate());
                }
            };
            
            template<class Context, class Data>
            class ASinOp
            /**
               \brief Returns the arcus sine of the operand.
             */
            {
            public:
                enum { has_evaluate = false, has_context = false };
                
                static inline unsigned long precision(const Data & a)
                {
                    return a.precision();
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & a)
                {
                    asin(x, a.evaluate());
                }
            };
            
            template<class Context, class Data>
            class ACosOp
            /**
               \brief Returns the arcus cosine of the operand.
             */
            {
            public:
                enum { has_evaluate = false, has_context = false };
                
                static inline unsigned long precision(const Data & a)
                {
                    return a.precision();
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & a)
                {
                    acos(x, a.evaluate());
                }
            };
            
            template<class Context, class Data>
            class ATanOp
            /**
               \brief Returns the arcus tangent of the operand.
             */
            {
            public:
                enum { has_evaluate = false, has_context = false };
                
                static inline unsigned long precision(const Data & a)
                {
                    return a.precision();
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & a)
                {
                    atan(x, a.evaluate());
                }
            };
            
            template<class Context, class Data>
            class ATan2Op
            /**
               \brief Returns the arcus tangent of the quotient of the operands.
             */
            {
            public:
                enum { has_evaluate = false, has_context = false };
                
                static inline unsigned long precision(const Data & d)
                {
                    return std::max(d.first.precision(), d.second.precision());
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & d)
                {
                    atan2(x, d.first.evaluate(), d.second.evaluate());
                }
            };
            
            template<class Context, class Data>
            class ExpOp
            /**
               \brief Returns the exponential of the operand.
             */
            {
            public:
                enum { has_evaluate = false, has_context = false };
                
                static inline unsigned long precision(const Data & a)
                {
                    return a.precision();
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & a)
                {
                    exp(x, a.evaluate());
                }
            };
            
            template<class Context, class Data>
            class LogOp
            /**
               \brief Returns the natural logarithm of the operand.
             */
            {
            public:
                enum { has_evaluate = false, has_context = false };
                
                static inline unsigned long precision(const Data & a)
                {
                    return a.precision();
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & a)
                {
                    log(x, a.evaluate());
                }
            };
            
            template<class Context, class Data>
            class Log2Op
            /**
               \brief Returns the logarithm to base 2 of the operand.
             */
            {
            public:
                enum { has_evaluate = false, has_context = false };
                
                static inline unsigned long precision(const Data & a)
                {
                    return a.precision();
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & a)
                {
                    log2(x, a.evaluate());
                }
            };
            
            template<class Context, class Data>
            class Log10Op
            /**
               \brief Returns the logarithm to base 10 of the operand.
             */
            {
            public:
                enum { has_evaluate = false, has_context = false };
                
                static inline unsigned long precision(const Data & a)
                {
                    return a.precision();
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & a)
                {
                    log10(x, a.evaluate());
                }
            };
            
            template<class Context, class Data>
            class SqrtOp
            /**
               \brief Returns the square root of the operand.
             */
            {
            public:
                enum { has_evaluate = false, has_context = false };
                
                static inline unsigned long precision(const Data & a)
                {
                    return a.precision();
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & a)
                {
                    sqrt(x, a.evaluate());
                }
            };
            
            template<class Context, class Data>
            class GammaOp
            /**
               \brief Returns the Gamma function evaluated of the operand.
             */
            {
            public:
                enum { has_evaluate = false, has_context = false };
                
                static inline unsigned long precision(const Data & a)
                {
                    return a.precision();
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & a)
                {
                    gamma(x, a.evaluate());
                }
            };
            
            template<class Context, class Data>
            class LGammaOp
            /**
               \brief Returns the logarithm of the absolute value of the Gamma function evaluated of
                      the operand.
             */
            {
            public:
                enum { has_evaluate = false, has_context = false };
                
                static inline unsigned long precision(const Data & a)
                {
                    return a.precision();
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & a)
                {
                    lgamma(x, a.evaluate());
                }
            };
            
            template<class Context, class Data>
            class LGamma2Op
            /**
               \brief Returns the logarithm of the absolute value of the Gamma function evaluated of
                      the operand together with the sign of the Gamma function.
             */
            {
            private:
                int & d_sign;
                
            public:
                enum { has_evaluate = false, has_context = false };
                
                inline LGamma2Op(int & sign)
                : d_sign(sign)
                {
                }
                
                static inline unsigned long precision(const Data & a)
                {
                    return a.precision();
                }
                
                inline void assignTo(typename Context::Type & x, const Data & a) const
                {
                    lgamma(x, d_sign, a.evaluate());
                }
            };
            
            template<typename Type>
            struct PowerCOp
            /**
               \brief Returns the power of the operand to a constant stored in the operator object.
               
               To obtain the operator template, instantiate this template with the desired exponent
               type, and use its member `impl`.
             */
            {
                template<class Context, class Data>
                class impl
                /**
                   \brief Implementation of raising to a power by a constant exponent.
                 */
                {
                private:
                    Type d_v;
                    
                public:
                    enum { has_evaluate = false, has_context = false };
                    
                    inline impl(Type v)
                    : d_v(v)
                    {
                    }
                    
                    static inline unsigned long precision(const Data & a)
                    {
                        return a.precision();
                    }
                    
                    inline void assignTo(typename Context::Type & x, const Data & a) const
                    {
                        power(x, a.evaluate(), d_v);
                    }
                };
            };
            
            template<class Context, class Data>
            class PowerOp
            /**
               \brief Returns the power of the first operand to the second operand.
             */
            {
            public:
                enum { has_evaluate = false, has_context = false };
                
                static inline unsigned long precision(const Data & d)
                {
                    return d.first.precision();
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & d)
                {
                    power(x, d.first.evaluate(), d.second.evaluate());
                }
            };
            
            ////////////////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////
            ////// FUNCTIONAL ARITHMETIC W/ CONTEXT ////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////
            
            template<typename Context, typename Data>
            class AbsOp_Context
            /**
               \brief Returns the absolute value of the operand. Uses the context stored in the
                      operator.
             */
            {
            private:
                const Context & d_context;
                
            public:
                enum { has_evaluate = false, has_context = true };
                
                inline AbsOp_Context(const Context & context)
                    : d_context(context)
                {
                }
                
                inline unsigned long precision(const Data & a)
                {
                    return d_context.getRealPrecision();
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & a)
                {
                    abs(x, a.evaluate());
                }
            };
            
            template<typename Context, typename Data>
            class SquareOp_Context
            /**
               \brief Returns the square of the operand.  Uses the context stored in the operator.
             */
            {
            private:
                const Context & d_context;
                
            public:
                enum { has_evaluate = false, has_context = true };
                
                inline SquareOp_Context(const Context & context)
                    : d_context(context)
                {
                }
                
                inline unsigned long precision(const Data & a)
                {
                    return d_context.getRealPrecision();
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & a)
                {
                    square(x, a.evaluate());
                }
            };
            
            template<class Context, class Data>
            class SinOp_Context
            /**
               \brief Returns the sine of the operand.  Uses the context stored in the operator.
             */
            {
            private:
                const Context & d_context;
                
            public:
                enum { has_evaluate = false, has_context = true };
                
                inline SinOp_Context(const Context & context)
                    : d_context(context)
                {
                }
                
                inline unsigned long precision(const Data & a)
                {
                    return d_context.getRealPrecision();
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & a)
                {
                    sin(x, a.evaluate());
                }
            };
            
            template<class Context, class Data>
            class CosOp_Context
            /**
               \brief Returns the cosine of the operand. Uses the context stored in the operator.
             */
            {
            private:
                const Context & d_context;
                
            public:
                enum { has_evaluate = false, has_context = true };
                
                inline CosOp_Context(const Context & context)
                    : d_context(context)
                {
                }
                
                inline unsigned long precision(const Data & a)
                {
                    return d_context.getRealPrecision();
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & a)
                {
                    cos(x, a.evaluate());
                }
            };
            
            template<class Context, class Data>
            class TanOp_Context
            /**
               \brief Returns the tangent of the operand. Uses the context stored in the operator.
             */
            {
            private:
                const Context & d_context;
                
            public:
                enum { has_evaluate = false, has_context = true };
                
                inline TanOp_Context(const Context & context)
                    : d_context(context)
                {
                }
                
                inline unsigned long precision(const Data & a)
                {
                    return d_context.getRealPrecision();
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & a)
                {
                    tan(x, a.evaluate());
                }
            };
            
            template<class Context, class Data>
            class ASinOp_Context
            /**
               \brief Returns the arcus sine of the operand.  Uses the context stored in the
                      operator.
             */
            {
            private:
                const Context & d_context;
                
            public:
                enum { has_evaluate = false, has_context = true };
                
                inline ASinOp_Context(const Context & context)
                    : d_context(context)
                {
                }
                
                inline unsigned long precision(const Data & a)
                {
                    return d_context.getRealPrecision();
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & a)
                {
                    asin(x, a.evaluate());
                }
            };
            
            template<class Context, class Data>
            class ACosOp_Context
            /**
               \brief Returns the arcus cosine of the operand. Uses the context stored in the
                      operator.
             */
            {
            private:
                const Context & d_context;
                
            public:
                enum { has_evaluate = false, has_context = true };
                
                inline ACosOp_Context(const Context & context)
                    : d_context(context)
                {
                }
                
                inline unsigned long precision(const Data & a)
                {
                    return d_context.getRealPrecision();
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & a)
                {
                    acos(x, a.evaluate());
                }
            };
            
            template<class Context, class Data>
            class ATanOp_Context
            /**
               \brief Returns the arcus tangent of the operand. Uses the context stored in the
                      operator.
             */
            {
            private:
                const Context & d_context;
                
            public:
                enum { has_evaluate = false, has_context = true };
                
                inline ATanOp_Context(const Context & context)
                    : d_context(context)
                {
                }
                
                inline unsigned long precision(const Data & a)
                {
                    return d_context.getRealPrecision();
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & a)
                {
                    atan(x, a.evaluate());
                }
            };
            
            template<class Context, class Data>
            class ATan2Op_Context
            /**
               \brief Returns the arcus tangent of the quotient of the operands. Uses the context
                      stored in the operator.
             */
            {
            private:
                const Context & d_context;
                
            public:
                enum { has_evaluate = false, has_context = true };
                
                inline ATan2Op_Context(const Context & context)
                    : d_context(context)
                {
                }
                
                inline unsigned long precision(const Data & d)
                {
                    return d_context.getRealPrecision();
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & d)
                {
                    atan2(x, d.first.evaluate(), d.second.evaluate());
                }
            };
            
            template<class Context, class Data>
            class ExpOp_Context
            /**
               \brief Returns the exponential of the operand. Uses the context stored in the
                      operator.
             */
            {
            private:
                const Context & d_context;
                
            public:
                enum { has_evaluate = false, has_context = true };
                
                inline ExpOp_Context(const Context & context)
                    : d_context(context)
                {
                }
                
                inline unsigned long precision(const Data & a)
                {
                    return d_context.getRealPrecision();
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & a)
                {
                    exp(x, a.evaluate());
                }
            };
            
            template<class Context, class Data>
            class LogOp_Context
            /**
               \brief Returns the natural logarithm of the operand. Uses the context stored in the
                      operator.
             */
            {
            private:
                const Context & d_context;
                
            public:
                enum { has_evaluate = false, has_context = true };
                
                inline LogOp_Context(const Context & context)
                    : d_context(context)
                {
                }
                
                inline unsigned long precision(const Data & a)
                {
                    return d_context.getRealPrecision();
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & a)
                {
                    log(x, a.evaluate());
                }
            };
            
            template<class Context, class Data>
            class Log2Op_Context
            /**
               \brief Returns the logarithm to base 2 of the operand. Uses the context stored in the
                      operator.
             */
            {
            private:
                const Context & d_context;
                
            public:
                enum { has_evaluate = false, has_context = true };
                
                inline Log2Op_Context(const Context & context)
                    : d_context(context)
                {
                }
                
                inline unsigned long precision(const Data & a)
                {
                    return d_context.getRealPrecision();
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & a)
                {
                    log2(x, a.evaluate());
                }
            };
            
            template<class Context, class Data>
            class Log10Op_Context
            /**
               \brief Returns the logarithm to base 10 of the operand. Uses the context stored in
                      the operator.
             */
            {
            private:
                const Context & d_context;
                
            public:
                enum { has_evaluate = false, has_context = true };
                
                inline Log10Op_Context(const Context & context)
                    : d_context(context)
                {
                }
                
                inline unsigned long precision(const Data & a)
                {
                    return d_context.getRealPrecision();
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & a)
                {
                    log10(x, a.evaluate());
                }
            };
            
            template<class Context, class Data>
            class SqrtOp_Context
            /**
               \brief Returns the square root of the operand. Uses the context stored in the
                      operator.
             */
            {
            private:
                const Context & d_context;
                
            public:
                enum { has_evaluate = false, has_context = true };
                
                inline SqrtOp_Context(const Context & context)
                    : d_context(context)
                {
                }
                
                inline unsigned long precision(const Data & a)
                {
                    return d_context.getRealPrecision();
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & a)
                {
                    sqrt(x, a.evaluate());
                }
            };
            
            template<class Context, class Data>
            class GammaOp_Context
            /**
               \brief Returns the Gamma function evaluated of the operand. Uses the context stored
                      in the operator.
             */
            {
            private:
                const Context & d_context;
                
            public:
                enum { has_evaluate = false, has_context = true };
                
                inline GammaOp_Context(const Context & context)
                    : d_context(context)
                {
                }
                
                inline unsigned long precision(const Data & a)
                {
                    return d_context.getRealPrecision();
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & a)
                {
                    gamma(x, a.evaluate());
                }
            };
            
            template<class Context, class Data>
            class LGammaOp_Context
            /**
               \brief Returns the logarithm of the absolute value of the Gamma function evaluated of
                      the operand. Uses the context stored in the operator.
             */
            {
            private:
                const Context & d_context;
                
            public:
                enum { has_evaluate = false, has_context = true };
                
                inline LGammaOp_Context(const Context & context)
                    : d_context(context)
                {
                }
                
                inline unsigned long precision(const Data & a)
                {
                    return d_context.getRealPrecision();
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & a)
                {
                    lgamma(x, a.evaluate());
                }
            };
            
            template<class Context, class Data>
            class LGamma2Op_Context
            /**
               \brief Returns the logarithm of the absolute value of the Gamma function evaluated of
                      the operand together with the sign of the Gamma function. Uses the context
                      stored in the operator.
             */
            {
            private:
                const Context & d_context;
                int & d_sign;
                
            public:
                enum { has_evaluate = false, has_context = true };
                
                inline LGamma2Op_Context(int & sign, const Context & context)
                    : d_context(context), d_sign(sign)
                {
                }
                
                inline unsigned long precision(const Data & a)
                {
                    return d_context.getRealPrecision();
                }
                
                inline void assignTo(typename Context::Type & x, const Data & a) const
                {
                    lgamma(x, d_sign, a.evaluate());
                }
            };
            
            template<typename Type>
            struct PowerCOp_Context
            /**
               \brief Returns the power of the operand to a constant stored in the operator
                      object. Uses the context stored in the operator.
               
               To obtain the operator template, instantiate this template with the desired exponent
               type, and use its member `impl`.
             */
            {
                template<class Context, class Data>
                class impl
                /**
                   \brief Implementation of raising to a power by a constant exponent.
                 */
                {
                private:
                    const Context & d_context;
                    Type d_v;
                    
                public:
                    enum { has_evaluate = false, has_context = true };
                    
                    inline impl(Type v, const Context & context)
                        : d_context(context), d_v(v)
                    {
                    }
                    
                    inline unsigned long precision(const Data & a)
                    {
                        return d_context.getRealPrecision();
                    }
                    
                    inline void assignTo(typename Context::Type & x, const Data & a) const
                    {
                        power(x, a.evaluate(), d_v);
                    }
                };
            };
            
            template<class Context, class Data>
            class PowerOp_Context
            /**
               \brief Returns the power of the first operand to the second operand. Uses the context
                      stored in the operator.
             */
            {
            private:
                const Context & d_context;
                
            public:
                enum { has_evaluate = false, has_context = true };
                
                inline PowerOp_Context(const Context & context)
                    : d_context(context)
                {
                }
                
                inline unsigned long precision(const Data & d)
                {
                    return d.first.precision();
                }
                
                static inline void assignTo(typename Context::Type & x, const Data & d)
                {
                    power(x, d.first.evaluate(), d.second.evaluate());
                }
            };
            
            ////////////////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////
            ////// CONVERSION OPERATORS ////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////
            
            template<typename Type>
            class ConversionWrapper
            /**
               \brief Wraps a native type for a conversion.
               
               \warning A copy of the type is stored in the `ConversionWrapper<>` object.
             */
            {
            private:
                Type d_value;
                
            public:
                inline ConversionWrapper(const Type & value)
                    : d_value(value)
                {
                }
                
                inline const Type & evaluate() const
                {
                    return d_value;
                }
            };
            
            template<typename Context, typename Data>
            class ConvertOp_Context
            /**
               \brief Converts the argument to the destination type. Uses the context stored in the
                      operator.
             */
            {
            private:
                const Context & d_context;
                
            public:
                enum { has_evaluate = false, has_context = true };
                
                inline ConvertOp_Context(const Context & context)
                    : d_context(context)
                {
                }
                
                const Context & context(const Data &) const
                {
                    return d_context;
                }
                
                inline unsigned long precision(const Data & a) const
                {
                    return d_context.getRealPrecision();
                }
                
                inline void assignTo(typename Context::Type & x, const Data & a) const
                {
                    convert(x, a.evaluate(), d_context);
                }
            };
            
            template<typename Context, typename Data>
            class ConvertFloorOp_Context
            /**
               \brief Converts the argument to the destination type by rounding down (floor). Uses
                      the context stored in the operator.
             */
            {
            private:
                const Context & d_context;
                
            public:
                enum { has_evaluate = false, has_context = true };
                
                inline ConvertFloorOp_Context(const Context & context)
                    : d_context(context)
                {
                }
                
                const Context & context(const Data &) const
                {
                    return d_context;
                }
                
                inline unsigned long precision(const Data & a) const
                {
                    return d_context.getRealPrecision();
                }
                
                inline void assignTo(typename Context::Type & x, const Data & a) const
                {
                    convert_floor(x, a.evaluate(), d_context);
                }
            };
            
            template<typename Context, typename Data>
            class ConvertRoundOp_Context
            /**
               \brief Converts the argument to the destination type by rounding. Uses the context
                      stored in the operator.
             */
            {
            private:
                const Context & d_context;
                
            public:
                enum { has_evaluate = false, has_context = true };
                
                inline ConvertRoundOp_Context(const Context & context)
                    : d_context(context)
                {
                }
                
                const Context & context(const Data &) const
                {
                    return d_context;
                }
                
                inline unsigned long precision(const Data & a) const
                {
                    return d_context.getRealPrecision();
                }
                
                inline void assignTo(typename Context::Type & x, const Data & a) const
                {
                    convert_round(x, a.evaluate(), d_context);
                }
            };
            
            template<typename Context, typename Data>
            class ConvertRound2Op_Context
            /**
               \brief Converts the argument to the destination type by rounding. Uses the context
                      stored in the operator and stores the rounding direction in the provided
                      `bool` reference.
             */
            {
            private:
                const Context & d_context;
                bool & d_rounded_up;
                
            public:
                enum { has_evaluate = false, has_context = true };
                
                inline ConvertRound2Op_Context(bool & rounded_up, const Context & context)
                    : d_context(context), d_rounded_up(rounded_up)
                {
                }
                
                const Context & context(const Data &) const
                {
                    return d_context;
                }
                
                inline unsigned long precision(const Data & a) const
                {
                    return d_context.getRealPrecision();
                }
                
                inline void assignTo(typename Context::Type & x, const Data & a) const
                {
                    convert_round(x, a.evaluate(), d_rounded_up, d_context);
                }
            };
            
            template<typename Context, typename Data>
            class ConvertCeilOp_Context
            /**
               \brief Converts the argument to the destination type by rounding up (ceil). Uses the
                      context stored in the operator.
             */
            {
            private:
                const Context & d_context;
                
            public:
                enum { has_evaluate = false, has_context = true };
                
                inline ConvertCeilOp_Context(const Context & context)
                    : d_context(context)
                {
                }
                
                const Context & context(const Data &) const
                {
                    return d_context;
                }
                
                inline unsigned long precision(const Data & a) const
                {
                    return d_context.getRealPrecision();
                }
                
                inline void assignTo(typename Context::Type & x, const Data & a) const
                {
                    convert_ceil(x, a.evaluate(), d_context);
                }
            };
            
            ////////////////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////
            ////// EXPRESSION //////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////
            
            template<typename Context_, class Data, template<typename, typename> class Op>
            class Expression : public Op<Context_, Data>
            /**
               \brief Represents an expression.
               
               \tparam Context_ The arithmetic context describing the result of the expression.
               \tparam Data Data for the operation.
               \tparam Op The operand template template class. Describes how the data should be
                          turned into the result.
             */
            {
            private:
                Data d_data;
                
                template<bool, class A, class B>
                struct SelectET
                {
                    typedef typename A::evaluate_type result;
                };
                
                template<class A, class B>
                struct SelectET<false, A, B>
                {
                    typedef B result;
                };
                
            public:
                enum { has_evaluate = true, has_context = Op<Context_, Data>::has_context };
                typedef typename SelectET<Op<Context_, Data>::has_evaluate,
                                          Op<Context_, Data>, typename Context_::Type>::result evaluate_type;
                typedef Context_ Context;
                
                inline const Op<Context_, Data> & op() const
                /**
                   \brief Allows to access the operator object.
                 */
                {
                    return *this;
                }
                
                inline const Data & data() const
                /**
                   \brief Allows to access the expression's data.
                 */
                {
                    return d_data;
                }
                
                inline const Context & context() const
                /**
                   \brief Allows to access the expression's context. This function only compiles if
                          a context is stored in the expression!
                 */
                {
                    return Op<Context_, Data>::context(d_data);
                }
                
                inline Expression(const Data & data, const Op<Context_, Data> & op)
                /**
                   \brief Creates an expression from the given data and operator object.
                 */
                    : Op<Context_, Data>(op), d_data(data)
                {
                }
                
                inline Expression(const Data & data)
                /**
                   \brief Creates an expression from the given data. The operator object is default
                          initialized.
                 */
                    : Op<Context_, Data>(), d_data(data)
                {
                }
                
                inline unsigned long precision() const
                /**
                   \brief Returns the precision of the expression. Only compiles if the expression
                          supports this.
                 */
                {
                    return Op<Context_, Data>::precision(d_data);
                }
                
                /**
                   \brief Evaluates the expression into the given object.
                 */
                inline void assignTo(typename Context::Type & x) const;
                
            private:
                template<bool b>
                inline evaluate_type do_cast(helper::BoolToType<true>, helper::BoolToType<b>) const
                {
                    return Op<Context_, Data>::evaluate(d_data);
                }
                
                inline evaluate_type do_cast(helper::BoolToType<false>, helper::BoolToType<true>) const
                {
                    return typename Context_::Type(*this, Op<Context_, Data>::context(d_data));
                }
                
                inline evaluate_type do_cast(helper::BoolToType<false>, helper::BoolToType<false>) const
                {
                    return typename Context_::Type(*this);
                }
                
            public:
                inline operator evaluate_type () const
                /**
                   \brief Returns the value of the expression.
                 */
                {
                    return do_cast(helper::BoolToType<Op<Context_, Data>::has_evaluate>(), helper::BoolToType<has_context>());
                }
                
                inline evaluate_type evaluate() const
                /**
                   \brief Returns the value of the expression.
                 */
                {
                    return do_cast(helper::BoolToType<Op<Context_, Data>::has_evaluate>(), helper::BoolToType<has_context>());
                }
            };
            
            template<typename Context, class Data, class Op>
            void do_assign(typename Context::Type & x, const Op & op, const Data & data)
            /**
               \brief Performs the assignment of the expression to `x`.
               
               This function is provided to allow overloading for specialized situations. For
               example, if the underlying type provides an `addmul()` function, overloading this
               function for certain expression constructs allows to use this function.
               
               By default, the operator object's `assignTo()` function is called to do this
               assignment.
               
               \param x Where to store the result in.
               \param op The operator object.
               \param data The provided data.
            */
            {
                op.assignTo(x, data);
            }
            
            template<typename Context, class Data, template<typename, typename> class Op>
            inline void Expression<Context, Data, Op>::assignTo(typename Context::Type & x) const
            {
                do_assign<Context, Data, Op<Context, Data> >(x, static_cast<const Op<Context, Data> &>(*this), d_data);
            }
            
            template<typename Context>
            inline Expression<Context, Wrapper<Context>, NoneOp> make_expression(const typename Context::Type & a)
            /**
               \brief Creates an expression from a type.
             */
            {
                return Expression<Context, Wrapper<Context>, NoneOp>(Wrapper<Context>(a));
            }
        }
        
        // generic conversion from expressions:
        namespace implementation
        {
            template<class SourceContext, class Data, template<typename, typename> class Op, class Dest>
            class conversion_impl<expressions::Expression<SourceContext, Data, Op>, Dest>
            /**
               \brief Implements a conversion from an expression to a destination context.
               
               Uses the conversion provided by the expression's result context to the destination
               context.
             */
            {
            public:
                typedef typename conversion_impl<typename SourceContext::Type, Dest>::RetVal RetVal;
                typedef typename conversion_impl<typename SourceContext::Type, Dest>::RetVal_Frac RetVal_Frac;
                typedef typename conversion_impl<typename SourceContext::Type, Dest>::RetVal_Floor RetVal_Floor;
                typedef typename conversion_impl<typename SourceContext::Type, Dest>::RetVal_Ceil RetVal_Ceil;
                typedef typename conversion_impl<typename SourceContext::Type, Dest>::RetVal_Round RetVal_Round;
                typedef typename conversion_impl<typename SourceContext::Type, Dest>::RetVal_Round2 RetVal_Round2;
                
                static void convert(typename Dest::Type & d, const expressions::Expression<SourceContext, Data, Op> & v, const Dest & c)
                {
                    conversion_impl<typename SourceContext::Type, Dest>::convert(d, v.evaluate(), c);
                }
                
                static RetVal convert(const expressions::Expression<SourceContext, Data, Op> & v, const Dest & c)
                {
                    return conversion_impl<typename SourceContext::Type, Dest>::convert(v.evaluate(), c);
                }
                
                static void convert_frac(typename Dest::Type & d, const expressions::Expression<SourceContext, Data, Op> & v1,
                                         const expressions::Expression<SourceContext, Data, Op> & v2, const Dest & c)
                {
                    conversion_impl<typename SourceContext::Type, Dest>::convert(d, v1.evaluate(), v2.evaluate(), c);
                }
                
                static RetVal_Frac convert_frac(const expressions::Expression<SourceContext, Data, Op> & v1,
                                                const expressions::Expression<SourceContext, Data, Op> & v2, const Dest & c)
                {
                    return conversion_impl<typename SourceContext::Type, Dest>::convert(v1.evaluate(), v2.evaluate(), c);
                }
                
                static void floor(typename Dest::Type & d, const expressions::Expression<SourceContext, Data, Op> & v, const Dest & c)
                {
                    conversion_impl<typename SourceContext::Type, Dest>::floor(d, v.evaluate(), c);
                }
                
                static RetVal_Floor floor(const expressions::Expression<SourceContext, Data, Op> & v, const Dest & c)
                {
                    return conversion_impl<typename SourceContext::Type, Dest>::floor(v.evaluate(), c);
                }
                
                static void round(typename Dest::Type & d, const expressions::Expression<SourceContext, Data, Op> & v, const Dest & c)
                {
                    conversion_impl<typename SourceContext::Type, Dest>::round(d, v.evaluate(), c);
                }
                
                static RetVal_Round round(const expressions::Expression<SourceContext, Data, Op> & v, const Dest & c)
                {
                    return conversion_impl<typename SourceContext::Type, Dest>::round(v.evaluate(), c);
                }
                
                static void round(typename Dest::Type & d, const expressions::Expression<SourceContext, Data, Op> & v, bool & up, const Dest & c)
                {
                    conversion_impl<typename SourceContext::Type, Dest>::round(d, v.evaluate(), up, c);
                }
                
                static RetVal_Round2 round(const expressions::Expression<SourceContext, Data, Op> & v, bool & up, const Dest & c)
                {
                    return conversion_impl<typename SourceContext::Type, Dest>::round(v.evaluate(), up, c);
                }
                
                static void ceil(typename Dest::Type & d, const expressions::Expression<SourceContext, Data, Op> & v, const Dest & c)
                {
                    conversion_impl<typename SourceContext::Type, Dest>::ceil(d, v.evaluate(), c);
                }
                
                static RetVal_Ceil ceil(const expressions::Expression<SourceContext, Data, Op> & v, const Dest & c)
                {
                    return conversion_impl<typename SourceContext::Type, Dest>::ceil(v.evaluate(), c);
                }
            };
            
            template<class SourceContext, class Data, template<typename, typename> class Op>
            class nativeconversion_impl<expressions::Expression<SourceContext, Data, Op> >
            /**
               \brief Implements native conversions from an expression.
               
               Uses the native conversions provided by the expression's result context.
             */
            {
            public:
                static int toInt(const expressions::Expression<SourceContext, Data, Op> & v)
                {
                    return nativeconversion_impl<typename SourceContext::Type>::toInt(v.evaluate());
                }
        
                static int toInt_Floor(const expressions::Expression<SourceContext, Data, Op> & v)
                {
                    return nativeconversion_impl<typename SourceContext::Type>::toInt_Floor(v.evaluate());
                }
        
                static int toInt_Round(const expressions::Expression<SourceContext, Data, Op> & v)
                {
                    return nativeconversion_impl<typename SourceContext::Type>::toInt_Round(v.evaluate());
                }
        
                static int toInt_Round(const expressions::Expression<SourceContext, Data, Op> & v, bool & up)
                {
                    return nativeconversion_impl<typename SourceContext::Type>::toInt_Round(v.evaluate(), up);
                }
        
                static int toInt_Ceil(const expressions::Expression<SourceContext, Data, Op> & v)
                {
                    return nativeconversion_impl<typename SourceContext::Type>::toInt_Ceil(v.evaluate());
                }
            
                static unsigned int toUInt(const expressions::Expression<SourceContext, Data, Op> & v)
                {
                    return nativeconversion_impl<typename SourceContext::Type>::toUInt(v.evaluate());
                }
        
                static unsigned int toUInt_Floor(const expressions::Expression<SourceContext, Data, Op> & v)
                {
                    return nativeconversion_impl<typename SourceContext::Type>::toUInt_Floor(v.evaluate());
                }
        
                static unsigned int toUInt_Round(const expressions::Expression<SourceContext, Data, Op> & v)
                {
                    return nativeconversion_impl<typename SourceContext::Type>::toUInt_Round(v.evaluate());
                }
        
                static unsigned int toUInt_Round(const expressions::Expression<SourceContext, Data, Op> & v, bool & up)
                {
                    return nativeconversion_impl<typename SourceContext::Type>::toUInt_Round(v.evaluate(), up);
                }
        
                static unsigned int toUInt_Ceil(const expressions::Expression<SourceContext, Data, Op> & v)
                {
                    return nativeconversion_impl<typename SourceContext::Type>::toUInt_Ceil(v.evaluate());
                }
        
                static long toLong(const expressions::Expression<SourceContext, Data, Op> & v)
                {
                    return nativeconversion_impl<typename SourceContext::Type>::toLong(v.evaluate());
                }
        
                static long toLong_Floor(const expressions::Expression<SourceContext, Data, Op> & v)
                {
                    return nativeconversion_impl<typename SourceContext::Type>::toLong_Floor(v.evaluate());
                }
        
                static long toLong_Round(const expressions::Expression<SourceContext, Data, Op> & v)
                {
                    return nativeconversion_impl<typename SourceContext::Type>::toLong_Round(v.evaluate());
                }
        
                static long toLong_Round(const expressions::Expression<SourceContext, Data, Op> & v, bool & up)
                {
                    return nativeconversion_impl<typename SourceContext::Type>::toLong_Round(v.evaluate(), up);
                }
        
                static long toLong_Ceil(const expressions::Expression<SourceContext, Data, Op> & v)
                {
                    return nativeconversion_impl<typename SourceContext::Type>::toLong_Ceil(v.evaluate());
                }
        
                static unsigned long toULong(const expressions::Expression<SourceContext, Data, Op> & v)
                {
                    return nativeconversion_impl<typename SourceContext::Type>::toULong(v.evaluate());
                }
        
                static unsigned long toULong_Floor(const expressions::Expression<SourceContext, Data, Op> & v)
                {
                    return nativeconversion_impl<typename SourceContext::Type>::toULong_Floor(v.evaluate());
                }
        
                static unsigned long toULong_Round(const expressions::Expression<SourceContext, Data, Op> & v)
                {
                    return nativeconversion_impl<typename SourceContext::Type>::toULong_Round(v.evaluate());
                }
        
                static unsigned long toULong_Round(const expressions::Expression<SourceContext, Data, Op> & v, bool & up)
                {
                    return nativeconversion_impl<typename SourceContext::Type>::toULong_Round(v.evaluate(), up);
                }
        
                static unsigned long toULong_Ceil(const expressions::Expression<SourceContext, Data, Op> & v)
                {
                    return nativeconversion_impl<typename SourceContext::Type>::toULong_Ceil(v.evaluate());
                }
        
                static long long toLongLong(const expressions::Expression<SourceContext, Data, Op> & v)
                {
                    return nativeconversion_impl<typename SourceContext::Type>::toLongLong(v.evaluate());
                }
        
                static long long toLongLong_Floor(const expressions::Expression<SourceContext, Data, Op> & v)
                {
                    return nativeconversion_impl<typename SourceContext::Type>::toLongLong_Floor(v.evaluate());
                }
        
                static long long toLongLong_Round(const expressions::Expression<SourceContext, Data, Op> & v)
                {
                    return nativeconversion_impl<typename SourceContext::Type>::toLongLong_Round(v.evaluate());
                }
        
                static long long toLongLong_Round(const expressions::Expression<SourceContext, Data, Op> & v, bool & up)
                {
                    return nativeconversion_impl<typename SourceContext::Type>::toLongLong_Round(v.evaluate(), up);
                }
        
                static long long toLongLong_Ceil(const expressions::Expression<SourceContext, Data, Op> & v)
                {
                    return nativeconversion_impl<typename SourceContext::Type>::toLongLong_Ceil(v.evaluate());
                }
        
                static float toFloat(const expressions::Expression<SourceContext, Data, Op> & v)
                {
                    return nativeconversion_impl<typename SourceContext::Type>::toFloat(v.evaluate());
                }
        
                static double toDouble(const expressions::Expression<SourceContext, Data, Op> & v)
                {
                    return nativeconversion_impl<typename SourceContext::Type>::toDouble(v.evaluate());
                }
        
                static long double toLongDouble(const expressions::Expression<SourceContext, Data, Op> & v)
                {
                    return nativeconversion_impl<typename SourceContext::Type>::toLongDouble(v.evaluate());
                }
            };
        }
    }
}

#endif

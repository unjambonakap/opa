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

#ifndef PLLL_INCLUDE_GUARD__RATIONAL_OPS_HPP
#define PLLL_INCLUDE_GUARD__RATIONAL_OPS_HPP

/**
   \file
   \brief Operator definitions for rational numbers.
   
   This header contains templates and instantiations to implement all operations on rational numbers.
*/
namespace plll
{
    namespace arithmetic
    {
        namespace expressions
        {
            // add
            template<class D1, class D2, template<typename, typename> class O2, class E>
            void do_assign(arithmetic::Rational & x,
                           const AddOp<RationalContext, std::pair<D1, Expression<RationalContext, Expression<IntegerContext, D2, O2>, ConvertOp_Context> > > &,
                           const std::pair<D1, Expression<RationalContext, Expression<IntegerContext, D2, O2>, ConvertOp_Context> > & data)
            {
                add_z(x, data.first.evaluate(), data.second.data().evaluate());
            }
            
            template<class D1, class D2, template<typename, typename> class O2, class E>
            void do_assign(arithmetic::Rational & x,
                           const AddOp<RationalContext, std::pair<Expression<RationalContext, Expression<IntegerContext, D2, O2>, ConvertOp_Context>, D1> > &,
                           const std::pair<Expression<RationalContext, Expression<IntegerContext, D2, O2>, ConvertOp_Context>, D1> & data)
            {
                add_z(x, data.second.evaluate(), data.first.data().evaluate());
            }
            
            // subtract
            template<class D1, class D2, template<typename, typename> class O2, class E>
            void do_assign(arithmetic::Rational & x,
                           const SubOp<RationalContext, std::pair<D1, Expression<RationalContext, Expression<IntegerContext, D2, O2>, ConvertOp_Context> > > &,
                           const std::pair<D1, Expression<RationalContext, Expression<IntegerContext, D2, O2>, ConvertOp_Context> > & data)
            {
                sub_z(x, data.first.evaluate(), data.second.data().evaluate());
            }
            
            template<class D1, class D2, template<typename, typename> class O2, class E>
            void do_assign(arithmetic::Rational & x,
                           const SubOp<RationalContext, std::pair<Expression<RationalContext, Expression<IntegerContext, D2, O2>, ConvertOp_Context>, D1> > &,
                           const std::pair<Expression<RationalContext, Expression<IntegerContext, D2, O2>, ConvertOp_Context>, D1> & data)
            {
                z_sub(x, data.first.data().evaluate(), data.second.evaluate());
            }
            
            // multiply
            template<class D1, class D2, template<typename, typename> class O2, class E>
            void do_assign(arithmetic::Rational & x,
                           const MulOp<RationalContext, std::pair<D1, Expression<RationalContext, Expression<IntegerContext, D2, O2>, ConvertOp_Context> > > &,
                           const std::pair<D1, Expression<RationalContext, Expression<IntegerContext, D2, O2>, ConvertOp_Context> > & data)
            {
                mul_z(x, data.first.evaluate(), data.second.data().evaluate());
            }
            
            template<class D1, class D2, template<typename, typename> class O2, class E>
            void do_assign(arithmetic::Rational & x,
                           const MulOp<RationalContext, std::pair<Expression<RationalContext, Expression<IntegerContext, D2, O2>, ConvertOp_Context>, D1> > &,
                           const std::pair<Expression<RationalContext, Expression<IntegerContext, D2, O2>, ConvertOp_Context>, D1> & data)
            {
                mul_z(x, data.second.evaluate(), data.first.data().evaluate());
            }
            
            // divide
            template<class D1, class D2, template<typename, typename> class O2, class E>
            void do_assign(arithmetic::Rational & x,
                           const DivOp<RationalContext, std::pair<D1, Expression<RationalContext, Expression<IntegerContext, D2, O2>, ConvertOp_Context> > > &,
                           const std::pair<D1, Expression<RationalContext, Expression<IntegerContext, D2, O2>, ConvertOp_Context> > & data)
            {
                div_z(x, data.first.evaluate(), data.second.data().evaluate());
            }
            
            template<class D1, class D2, template<typename, typename> class O2, class E>
            void do_assign(arithmetic::Rational & x,
                           const DivOp<RationalContext, std::pair<Expression<RationalContext, Expression<IntegerContext, D2, O2>, ConvertOp_Context>, D1> > &,
                           const std::pair<Expression<RationalContext, Expression<IntegerContext, D2, O2>, ConvertOp_Context>, D1> & data)
            {
                z_div(x, data.first.data().evaluate(), data.second.evaluate());
            }
            
            // modulo
            template<class D1, class D2, template<typename, typename> class O2, class E>
            void do_assign(arithmetic::Rational & x,
                           const ModOp<RationalContext, std::pair<D1, Expression<RationalContext, Expression<IntegerContext, D2, O2>, ConvertOp_Context> > > &,
                           const std::pair<D1, Expression<RationalContext, Expression<IntegerContext, D2, O2>, ConvertOp_Context> > & data)
            {
                mod_z(x, data.first.evaluate(), data.second.data().evaluate());
            }
            
            template<class D1, class D2, template<typename, typename> class O2, class E>
            void do_assign(arithmetic::Rational & x,
                           const ModOp<RationalContext, std::pair<Expression<RationalContext, Expression<IntegerContext, D2, O2>, ConvertOp_Context>, D1> > &,
                           const std::pair<Expression<RationalContext, Expression<IntegerContext, D2, O2>, ConvertOp_Context>, D1> & data)
            {
                z_mod(x, data.first.data().evaluate(), data.second.evaluate());
            }
            
            // shift
            template<class D, template<typename, typename> class O, class E>
            void do_assign(arithmetic::Rational & x,
                           const ShiftCOp<RationalContext, Expression<RationalContext, Expression<IntegerContext, D, O>, ConvertOp_Context> > & op,
                           Expression<RationalContext, Expression<IntegerContext, D, O>, ConvertOp_Context> & data)
            {
                shl_z(x, data.data().evaluate(), op.shift());
            }
            
            // negation, absolute value, square
            template<class D, template<typename, typename> class O>
            void do_assign(arithmetic::Rational & x,
                           const NegOp<RationalContext, Expression<RationalContext, Expression<IntegerContext, D, O>, ConvertOp_Context> > & op,
                           Expression<RationalContext, Expression<IntegerContext, D, O>, ConvertOp_Context> & data)
            {
                neg_z(x, data.op().value());
            }
            
            template<class D, template<typename, typename> class O>
            void do_assign(arithmetic::Rational & x,
                           const AbsOp<RationalContext, Expression<RationalContext, Expression<IntegerContext, D, O>, ConvertOp_Context> > & op,
                           Expression<RationalContext, Expression<IntegerContext, D, O>, ConvertOp_Context> & data)
            {
                abs_z(x, data.op().value());
            }
            
            template<class D, template<typename, typename> class O>
            void do_assign(arithmetic::Rational & x,
                           const AbsOp_Context<RationalContext, Expression<RationalContext, Expression<IntegerContext, D, O>, ConvertOp_Context> > & op,
                           Expression<RationalContext, Expression<IntegerContext, D, O>, ConvertOp_Context> & data)
            {
                abs_z(x, data.op().value());
            }
            
            template<class D, template<typename, typename> class O>
            void do_assign(arithmetic::Rational & x,
                           const SquareOp<RationalContext, Expression<RationalContext, Expression<IntegerContext, D, O>, ConvertOp_Context> > & op,
                           Expression<RationalContext, Expression<IntegerContext, D, O>, ConvertOp_Context> & data)
            {
                square_z(x, data.op().value());
            }
            
            template<class D, template<typename, typename> class O>
            void do_assign(arithmetic::Rational & x,
                           const SquareOp_Context<RationalContext, Expression<RationalContext, Expression<IntegerContext, D, O>, ConvertOp_Context> > & op,
                           Expression<RationalContext, Expression<IntegerContext, D, O>, ConvertOp_Context> & data)
            {
                square_z(x, data.op().value());
            }
            
            // power
            template<class D, template<typename, typename> class O>
            void do_assign(arithmetic::Rational & x,
                           const PowerCOp<signed long>::impl<RationalContext, Expression<RationalContext, Expression<IntegerContext, D, O>, ConvertOp_Context> > & op,
                           Expression<RationalContext, Expression<IntegerContext, D, O>, ConvertOp_Context> & data)
            {
                power_z(x, data.op().value(), op.value());
            }
            
            template<class D, template<typename, typename> class O>
            void do_assign(arithmetic::Rational & x,
                           const PowerCOp<unsigned long>::impl<RationalContext, Expression<RationalContext, Expression<IntegerContext, D, O>, ConvertOp_Context> > & op,
                           Expression<RationalContext, Expression<IntegerContext, D, O>, ConvertOp_Context> & data)
            {
                power_z(x, data.op().value(), op.value());
            }
            
            template<class D, template<typename, typename> class O, template<typename, typename> class IO, typename ID>
            void do_assign(arithmetic::Rational & x,
                           const PowerOp<RationalContext, std::pair<Expression<RationalContext, Expression<IntegerContext, D, O>, ConvertOp_Context>, Expression<IntegerContext, ID, IO> > > & op,
                           std::pair<Expression<RationalContext, Expression<IntegerContext, D, O>, ConvertOp_Context>, Expression<IntegerContext, ID, IO> > & data)
            {
                power_z(x, data.first.op().value(), data.second.evaluate());
            }
        }
        
        /**@{
           \name Assignment-operation operators.
        */
        
        /**
           \brief Subtracts the given rational number `r` from this rational number.
           
           \param cur The rational number to work on.
           \param r The rational number to subtract from this rational number.
           \return A reference to this rational number containing the result.
        */
        inline Rational & operator -= (Rational & cur, const Rational & r)
        {
            sub(cur, cur, r);
            return cur;
        }
        
        /**
           \brief Subtracts the given multiplication expression from this rational number.
           
           Note that the faster `submul()` command is used for this special case where the product
           of two rational numbers is subtracted from a third one.
           
           \param cur The rational number to work on.
           \param E The multiplication expression to subtract from this rational number.
           \return A reference to this rational number containing the result.
        */
        template<class D>
        inline Rational & operator -= (Rational & cur, const expressions::Expression<RationalContext, D, expressions::MulOp> & E)
        {
            submul(cur, E.data().first.evaluate(), E.data().second.evaluate());
            return cur;
        }
        
        /**
           \brief Adds the given rational number `r` to this rational number.
           
           \param cur The rational number to work on.
           \param r The rational number to add to this rational number.
           \return A reference to this rational number containing the result.
        */
        inline Rational & operator += (Rational & cur, const Rational & r)
        {
            add(cur, cur, r);
            return cur;
        }
        
        /**
           \brief Adds the given multiplication expression to this rational number.
           
           Note that the faster `addmul()` command is used for this special case where the product
           of two rational numbers is added to a third one.
           
           \param cur The rational number to work on.
           \param E The multiplication expression to add to this rational number.
           \return A reference to this rational number containing the result.
        */
        template<class D>
        inline Rational & operator += (Rational & cur, const expressions::Expression<RationalContext, D, expressions::MulOp> & E)
        {
            addmul(cur, E.data().first.evaluate(), E.data().second.evaluate());
            return cur;
        }
        
        /**
           \brief Multiplies the given rational number `r` with this rational number.
           
           \param cur The rational number to work on.
           \param r The rational number to multiply to this rational number.
           \return A reference to this rational number containing the result.
        */
        inline Rational & operator *= (Rational & cur, const Rational & r)
        {
            mul(cur, cur, r);
            return cur;
        }
        
        /**
           \brief Divides this rational number by the given rational number `r`.
           
           \param cur The rational number to work on.
           \param r The rational number to divide by.
           \return A reference to this rational number containing the result.
        */
        inline Rational & operator /= (Rational & cur, const Rational & r)
        {
            div(cur, cur, r);
            return cur;
        }
        
        /**
           \brief Divides this rational number by the given rational number `r` and stores the
                  remainder in this rational number.
                  
           \param cur The rational number to work on.
           \param r The rational number to divide by.
           \return A reference to this rational number containing the result.
        */
        inline Rational & operator %= (Rational & cur, const Rational & r)
        {
            mod(cur, cur, r);
            return cur;
        }
        
        /**
           \brief Multiplies this rational number by 2 to the power of `r`.
           
           \param cur The rational number to work on.
           \param r The exponent.
           \return A reference to this rational number containing the result.
        */
        inline Rational & operator <<= (Rational & cur, long r)
        {
            shl(cur, cur, r);
            return cur;
        }
        
        /**
           \brief Divides this rational number by 2 to the power of `r`.
           
           \param cur The rational number to work on.
           \param r The exponent.
           \return A reference to this rational number containing the result.
        */
        inline Rational & operator >>= (Rational & cur, long r)
        {
            shr(cur, cur, r);
            return cur;
        }
        ///@}
        
        template<class D, typename D2, template<typename, typename> class O2>
        inline Rational & operator -= (Rational & cur,
                                       const expressions::Expression<RationalContext,
                                                                     std::pair<D, expressions::Expression<RationalContext,
                                                                                                          expressions::Expression<IntegerContext, D2, O2>,
                                                                                                          expressions::ConvertOp_Context> >,
                                                                     expressions::MulOp> & E)
        {
            submul_z(cur, E.data().first.evaluate(), E.data().second.data().evaluate());
            return cur;
        }
        
        template<class D, typename D2, template<typename, typename> class O2>
        inline Rational & operator -= (Rational & cur,
                                       const expressions::Expression<RationalContext,
                                                                     std::pair<expressions::Expression<RationalContext,
                                                                                                       expressions::Expression<IntegerContext, D2, O2>,
                                                                                                       expressions::ConvertOp_Context>, D>,
                                                                     expressions::MulOp> & E)
        {
            submul_z(cur, E.data().second.evaluate(), E.data().first.data().evaluate());
            return cur;
        }
        
        template<class D, typename D2, template<typename, typename> class O2>
        inline Rational & operator += (Rational & cur,
                                       const expressions::Expression<RationalContext,
                                                                     std::pair<D, expressions::Expression<RationalContext,
                                                                                                          expressions::Expression<IntegerContext, D2, O2>,
                                                                                                          expressions::ConvertOp_Context> >,
                                                                     expressions::MulOp> & E)
        {
            addmul_z(cur, E.data().first.evaluate(), E.data().second.data().evaluate());
            return cur;
        }
        
        template<class D, typename D2, template<typename, typename> class O2>
        inline Rational & operator += (Rational & cur,
                                       const expressions::Expression<RationalContext,
                                                                     std::pair<expressions::Expression<RationalContext,
                                                                                                       expressions::Expression<IntegerContext, D2, O2>,
                                                                                                       expressions::ConvertOp_Context>, D>,
                                                                     expressions::MulOp> & E)
        {
            addmul_z(cur, E.data().second.evaluate(), E.data().first.data().evaluate());
            return cur;
        }
        
        /**@{
           \name Operators.
        */
        
        /**
           \brief Increments this rational number by one and returns the previous value.
           
           \param cur The rational number to work on.
           \return The previous value before incrementing.
        */
        inline Rational operator ++ (Rational & cur, int)
        {
            Rational r(cur);
            increment(cur, cur);
            return r;
        }
        
        /**
           \brief Decrements this rational number by one and returns the previous value.
           
           \param cur The rational number to work on.
           \return The previous value before decrementing.
        */
        inline Rational operator -- (Rational & cur, int)
        {
            Rational r(cur);
            decrement(cur, cur);
            return r;
        }
        
        /**
           \brief Increments this rational number by one and returns the new value.
           
           \param cur The rational number to work on.
           \return The new value.
        */
        inline Rational & operator ++ (Rational & cur)
        {
            increment(cur, cur);
            return cur;
        }
        
        /**
           \brief Decrements this rational number by one and returns the new value.
           
           \param cur The rational number to work on.
           \return The new value.
        */
        inline Rational & operator -- (Rational & cur)
        {
            decrement(cur, cur);
            return cur;
        }

        ///@}
        
        /**@{
           \name Comparisons.
        */
        
        /**
           \brief Compares the two rational numbers `a` and `b` for equality.
           
           \param a The first operand.
           \param b The second operand.
           \return `true` if and only if `a` equals `b`.
        */
        inline bool operator == (const Rational & a, const Rational & b)
        {
            return compare(a, b) == 0;
        }
        
        /**
           \brief Compares the two rational numbers `a` and `b` for inequality.
           
           \param a The first operand.
           \param b The second operand.
           \return `true` if and only if `a` does not equal `b`.
        */
        inline bool operator != (const Rational & a, const Rational & b)
        {
            return compare(a, b) != 0;
        }
        
        /**
           \brief Compares the two rational numbers `a` and `b`.
           
           \param a The first operand.
           \param b The second operand.
           \return `true` if and only if `a` is less than or equal to `b`.
        */
        inline bool operator <= (const Rational & a, const Rational & b)
        {
            return compare(a, b) <= 0;
        }
        
        /**
           \brief Compares the two rational numbers `a` and `b`.
               
           \param a The first operand.
           \param b The second operand.
           \return `true` if and only if `a` is greater than or equal to `b`.
        */
        inline bool operator >= (const Rational & a, const Rational & b)
        {
            return compare(a, b) >= 0;
        }
        
        /**
           \brief Compares the two rational numbers `a` and `b`.
               
           \param a The first operand.
           \param b The second operand.
           \return `true` if and only if `a` is less than `b`.
        */
        inline bool operator < (const Rational & a, const Rational & b)
        {
            return compare(a, b) < 0;
        }
        
        /**
           \brief Compares the two rational numbers `a` and `b`.
               
           \param a The first operand.
           \param b The second operand.
           \return `true` if and only if `a` is greater than `b`.
        */
        inline bool operator > (const Rational & a, const Rational & b)
        {
            return compare(a, b) > 0;
        }
        
        ///@}
        
        // Comparsions with second operand expression: make the first operand also an expression
        template<class A2, template<typename, typename> class O2>
        inline bool operator == (const Rational & a, const expressions::Expression<RationalContext, A2, O2> & b)
        {
            return expressions::make_expression<RationalContext>(a) == b;
        }
        
        template<class A2, template<typename, typename> class O2>
        inline bool operator != (const Rational & a, const expressions::Expression<RationalContext, A2, O2> & b)
        {
            return expressions::make_expression<RationalContext>(a) != b;
        }
        
        template<class A2, template<typename, typename> class O2>
        inline bool operator <= (const Rational & a, const expressions::Expression<RationalContext, A2, O2> & b)
        {
            return expressions::make_expression<RationalContext>(a) <= b;
        }
        
        template<class A2, template<typename, typename> class O2>
        inline bool operator >= (const Rational & a, const expressions::Expression<RationalContext, A2, O2> & b)
        {
            return expressions::make_expression<RationalContext>(a) >= b;
        }
            
        template<class A2, template<typename, typename> class O2>
        inline bool operator < (const Rational & a, const expressions::Expression<RationalContext, A2, O2> & b)
        {
            return expressions::make_expression<RationalContext>(a) < b;
        }
            
        template<class A2, template<typename, typename> class O2>
        inline bool operator > (const Rational & a, const expressions::Expression<RationalContext, A2, O2> & b)
        {
            return expressions::make_expression<RationalContext>(a) > b;
        }
        
        namespace expressions
        {
            // Comparsions with first operand expression: make the second operand also an expression
            template<class A1, template<typename, typename> class O1>
            inline bool operator == (const expressions::Expression<RationalContext, A1, O1> & a, const arithmetic::Rational & b)
            {
                return a == expressions::make_expression<RationalContext>(b);
            }
            
            template<class A1, template<typename, typename> class O1>
            inline bool operator != (const expressions::Expression<RationalContext, A1, O1> & a, const arithmetic::Rational & b)
            {
                return a != expressions::make_expression<RationalContext>(b);
            }
            
            template<class A1, template<typename, typename> class O1>
            inline bool operator <= (const expressions::Expression<RationalContext, A1, O1> & a, const arithmetic::Rational & b)
            {
                return a <= expressions::make_expression<RationalContext>(b);
            }
            
            template<class A1, template<typename, typename> class O1>
            inline bool operator >= (const expressions::Expression<RationalContext, A1, O1> & a, const arithmetic::Rational & b)
            {
                return a >= expressions::make_expression<RationalContext>(b);
            }
            
            template<class A1, template<typename, typename> class O1>
            inline bool operator < (const expressions::Expression<RationalContext, A1, O1> & a, const arithmetic::Rational & b)
            {
                return a < expressions::make_expression<RationalContext>(b);
            }
            
            template<class A1, template<typename, typename> class O1>
            inline bool operator > (const expressions::Expression<RationalContext, A1, O1> & a, const arithmetic::Rational & b)
            {
                return a > expressions::make_expression<RationalContext>(b);
            }
            
            // Comparsions with both operand expression
            template<class A1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
            inline bool operator == (const expressions::Expression<RationalContext, A1, O1> & a,
                                     const expressions::Expression<RationalContext, A2, O2> & b)
            {
                return compare(a, b) == 0;
            }
            
            template<class A1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
            inline bool operator != (const expressions::Expression<RationalContext, A1, O1> & a,
                                     const expressions::Expression<RationalContext, A2, O2> & b)
            {
                return compare(a, b) != 0;
            }
            
            template<class A1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
            inline bool operator <= (const expressions::Expression<RationalContext, A1, O1> & a,
                                     const expressions::Expression<RationalContext, A2, O2> & b)
            {
                return compare(a, b) <= 0;
            }
            
            template<class A1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
            inline bool operator >= (const expressions::Expression<RationalContext, A1, O1> & a,
                                     const expressions::Expression<RationalContext, A2, O2> & b)
            {
                return compare(a, b) >= 0;
            }
            
            template<class A1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
            inline bool operator < (const expressions::Expression<RationalContext, A1, O1> & a,
                                    const expressions::Expression<RationalContext, A2, O2> & b)
            {
                return compare(a, b) < 0;
            }
            
            template<class A1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
            inline bool operator > (const expressions::Expression<RationalContext, A1, O1> & a,
                                    const expressions::Expression<RationalContext, A2, O2> & b)
            {
                return compare(a, b) > 0;
            }
            
            // First operand is integer conversion expression, second operand arbitrary expression
            
            template<class A2, template<typename, typename> class O2, class A1, template<typename, typename> class O1>
            inline bool operator == (const expressions::Expression<RationalContext, expressions::Expression<IntegerContext, A1, O1>,
                                     expressions::ConvertOp_Context> & a,
                                     const expressions::Expression<RationalContext, A2, O2> & b)
            { return compare_z(b.evaluate(), a.op().value()) == 0; }
            template<class A2, template<typename, typename> class O2, class A1, template<typename, typename> class O1>
            inline bool operator != (const expressions::Expression<RationalContext, expressions::Expression<IntegerContext, A1, O1>,
                                     expressions::ConvertOp_Context> & a,
                                     const expressions::Expression<RationalContext, A2, O2> & b)
            { return compare_z(b.evaluate(), a.op().value()) != 0; }
            template<class A2, template<typename, typename> class O2, class A1, template<typename, typename> class O1>
            inline bool operator <= (const expressions::Expression<RationalContext, expressions::Expression<IntegerContext, A1, O1>,
                                     expressions::ConvertOp_Context> & a,
                                     const expressions::Expression<RationalContext, A2, O2> & b)
            { return compare_z(b.evaluate(), a.op().value()) >= 0; }
            template<class A2, template<typename, typename> class O2, class A1, template<typename, typename> class O1>
            inline bool operator >= (const expressions::Expression<RationalContext, expressions::Expression<IntegerContext, A1, O1>,
                                     expressions::ConvertOp_Context> & a,
                                     const expressions::Expression<RationalContext, A2, O2> & b)
            { return compare_z(b.evaluate(), a.op().value()) <= 0; }
            template<class A2, template<typename, typename> class O2, class A1, template<typename, typename> class O1>
            inline bool operator < (const expressions::Expression<RationalContext, expressions::Expression<IntegerContext, A1, O1>,
                                    expressions::ConvertOp_Context> & a,
                                    const expressions::Expression<RationalContext, A2, O2> & b)
            { return compare_z(b.evaluate(), a.op().value()) > 0; }
            template<class A2, template<typename, typename> class O2, class A1, template<typename, typename> class O1>
            inline bool operator > (const expressions::Expression<RationalContext, expressions::Expression<IntegerContext, A1, O1>,
                                    expressions::ConvertOp_Context> & a,
                                    const expressions::Expression<RationalContext, A2, O2> & b)
            { return compare_z(b.evaluate(), a.op().value()) < 0; }
            
            // First operand is arbitrary expression, Second operand is integer conversion expression
            template<class A1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
            inline bool operator == (const expressions::Expression<RationalContext, A1, O1> & a,
                                     const expressions::Expression<RationalContext, expressions::Expression<IntegerContext, A2, O2>,
                                     expressions::ConvertOp_Context> & b)
            { return compare_z(a.evaluate(), b.op().value()) == 0; }
            template<class A1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
            inline bool operator != (const expressions::Expression<RationalContext, A1, O1> & a,
                                     const expressions::Expression<RationalContext, expressions::Expression<IntegerContext, A2, O2>,
                                     expressions::ConvertOp_Context> & b)
            { return compare_z(a.evaluate(), b.op().value()) != 0; }
            template<class A1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
            inline bool operator <= (const expressions::Expression<RationalContext, A1, O1> & a,
                                     const expressions::Expression<RationalContext, expressions::Expression<IntegerContext, A2, O2>,
                                     expressions::ConvertOp_Context> & b)
            { return compare_z(a.evaluate(), b.op().value()) <= 0; }
            template<class A1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
            inline bool operator >= (const expressions::Expression<RationalContext, A1, O1> & a,
                                     const expressions::Expression<RationalContext, expressions::Expression<IntegerContext, A2, O2>,
                                     expressions::ConvertOp_Context> & b)
            { return compare_z(a.evaluate(), b.op().value()) >= 0; }
            template<class A1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
            inline bool operator < (const expressions::Expression<RationalContext, A1, O1> & a,
                                    const expressions::Expression<RationalContext, expressions::Expression<IntegerContext, A2, O2>,
                                    expressions::ConvertOp_Context> & b)
            { return compare_z(a.evaluate(), b.op().value()) < 0; }
            template<class A1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
            inline bool operator > (const expressions::Expression<RationalContext, A1, O1> & a,
                                    const expressions::Expression<RationalContext, expressions::Expression<IntegerContext, A2, O2>,
                                    expressions::ConvertOp_Context> & b)
            { return compare_z(a.evaluate(), b.op().value()) > 0; }
            
            // Predicates
            template<class A, template<typename, typename> class O>
            inline bool isZero(const expressions::Expression<RationalContext, A, O> & r)
            {
                return isZero(Rational(r));
            }
            
            template<class A, template<typename, typename> class O>
            inline bool isOne(const expressions::Expression<RationalContext, A, O> & r)
            {
                return isOne(Rational(r));
            }
            
            template<class A, template<typename, typename> class O>
            inline bool isPositive(const expressions::Expression<RationalContext, A, O> & r)
            {
                return isPositive(Rational(r));
            }
            
            template<class A, template<typename, typename> class O>
            inline bool isNonNegative(const expressions::Expression<RationalContext, A, O> & r)
            {
                return isNonNegative(Rational(r));
            }
            
            template<class A, template<typename, typename> class O>
            inline bool isNegative(const expressions::Expression<RationalContext, A, O> & r)
            {
                return isNegative(Rational(r));
            }
            
            template<class A, template<typename, typename> class O>
            inline bool isNonPositive(const expressions::Expression<RationalContext, A, O> & r)
            {
                return isZero(Rational(r));
            }
        }
        
        /**@{
           \name Operators.
        */
        
        /**
           \brief Negates this rational number.
           
           \param a The operand.
           \return The result.
        */
        inline expressions::Expression<RationalContext, expressions::Wrapper<RationalContext>,
                                       expressions::NegOp> operator - (const Rational & a)
        { return expressions::Expression<RationalContext, expressions::Wrapper<RationalContext>,
                                         expressions::NegOp>(expressions::Wrapper<RationalContext>(a)); }
        /**
           \brief Adds the argument to this rational number and returns the result.
               
           \param a The first operand.
           \param b The second operand.
           \return The result.
        */
        inline expressions::Expression<RationalContext, std::pair<expressions::Wrapper<RationalContext>,
                                                                  expressions::Wrapper<RationalContext> >,
                                       expressions::AddOp> operator + (const Rational & a, const Rational & b)
        { return expressions::Expression<RationalContext, std::pair<expressions::Wrapper<RationalContext>, expressions::Wrapper<RationalContext> >,
                                         expressions::AddOp>(std::make_pair(expressions::Wrapper<RationalContext>(a),
                                                                                             expressions::Wrapper<RationalContext>(b))); }
        /**
           \brief Subtracts the argument from this rational number and returns the result.
               
           \param a The first operand.
           \param b The second operand.
           \return The result.
        */
        inline expressions::Expression<RationalContext, std::pair<expressions::Wrapper<RationalContext>,
                                                                  expressions::Wrapper<RationalContext> >,
                                       expressions::SubOp> operator - (const Rational & a, const Rational & b)
        { return expressions::Expression<RationalContext, std::pair<expressions::Wrapper<RationalContext>,
                                                                    expressions::Wrapper<RationalContext> >,
                                         expressions::SubOp>(std::make_pair(expressions::Wrapper<RationalContext>(a),
                                                                            expressions::Wrapper<RationalContext>(b))); }
        /**
           \brief Multiplies the argument with this rational number and returns the result.
           
           \param a The first operand.
           \param b The second operand.
           \return The result.
        */
        inline expressions::Expression<RationalContext, std::pair<expressions::Wrapper<RationalContext>,
                                                                  expressions::Wrapper<RationalContext> >,
                                       expressions::MulOp> operator * (const Rational & a, const Rational & b)
        { return expressions::Expression<RationalContext, std::pair<expressions::Wrapper<RationalContext>,
                                                                    expressions::Wrapper<RationalContext> >,
                                         expressions::MulOp>(std::make_pair(expressions::Wrapper<RationalContext>(a),
                                                                            expressions::Wrapper<RationalContext>(b))); }
        /**
           \brief Divides this rational number by the argument and returns the result.
           
           \param a The first operand.
           \param b The second operand.
           \return The result.
        */
        inline expressions::Expression<RationalContext, std::pair<expressions::Wrapper<RationalContext>,
                                                                  expressions::Wrapper<RationalContext> >,
                                       expressions::DivOp> operator / (const Rational & a, const Rational & b)
        { return expressions::Expression<RationalContext, std::pair<expressions::Wrapper<RationalContext>,
                                                                    expressions::Wrapper<RationalContext> >,
                                         expressions::DivOp>(std::make_pair(expressions::Wrapper<RationalContext>(a),
                                                                            expressions::Wrapper<RationalContext>(b))); }
        /**
           \brief Returns the remainder of this rational number divided by the by the argument.
               
           \param a The first operand.
           \param b The second operand.
           \return The result.
        */
        inline expressions::Expression<RationalContext, std::pair<expressions::Wrapper<RationalContext>,
                                                                  expressions::Wrapper<RationalContext> >,
                                       expressions::ModOp> operator % (const Rational & a, const Rational & b)
        { return expressions::Expression<RationalContext, std::pair<expressions::Wrapper<RationalContext>,
                                                                    expressions::Wrapper<RationalContext> >,
                                         expressions::ModOp>(std::make_pair(expressions::Wrapper<RationalContext>(a),
                                                                            expressions::Wrapper<RationalContext>(b))); }
        /**
           \brief Multiplies this rational number with 2 to the power of `b` and returns the
           result.
               
           \param a The first operand.
           \param b The second operand.
           \return The result.
        */
        inline expressions::Expression<RationalContext, expressions::Wrapper<RationalContext>,
                                       expressions::ShiftCOp> operator << (const Rational & a, signed long b)
        { return expressions::Expression<RationalContext, expressions::Wrapper<RationalContext>,
                                         expressions::ShiftCOp>(expressions::Wrapper<RationalContext>(a), b); }
        
        /**
           \brief Divides this rational number by 2 to the power of `b` and returns the result.
                       
           \param a The first operand.
           \param b The second operand.
           \return The result.
        */
        inline expressions::Expression<RationalContext, expressions::Wrapper<RationalContext>,
                                       expressions::ShiftCOp> operator >> (const Rational & a, signed long b)
        { return expressions::Expression<RationalContext, expressions::Wrapper<RationalContext>,
                                         expressions::ShiftCOp>(expressions::Wrapper<RationalContext>(a), -b); }
        
        ///@}
        
        // Basic operations with second operand an expression
        template<class A, template<typename, typename> class O>
        inline expressions::Expression<RationalContext, std::pair<expressions::Wrapper<RationalContext>,
                                                                  expressions::Expression<RationalContext, A, O> >,
                                       expressions::AddOp> operator + (const Rational & a, const expressions::Expression<RationalContext, A, O> & b)
        { return expressions::Expression<RationalContext, std::pair<expressions::Wrapper<RationalContext>,
                                                                    expressions::Expression<RationalContext, A, O> >,
                                         expressions::AddOp>(std::make_pair(expressions::Wrapper<RationalContext>(a), b)); }
        template<class A, template<typename, typename> class O>
        inline expressions::Expression<RationalContext, std::pair<expressions::Wrapper<RationalContext>,
                                                                  expressions::Expression<RationalContext, A, O> >,
                                       expressions::SubOp> operator - (const Rational & a, const expressions::Expression<RationalContext, A, O> & b)
        { return expressions::Expression<RationalContext, std::pair<expressions::Wrapper<RationalContext>,
                                                                    expressions::Expression<RationalContext, A, O> >,
                                         expressions::SubOp>(std::make_pair(expressions::Wrapper<RationalContext>(a), b)); }
        template<class A, template<typename, typename> class O>
        inline expressions::Expression<RationalContext, std::pair<expressions::Wrapper<RationalContext>,
                                                                  expressions::Expression<RationalContext, A, O> >,
                                       expressions::MulOp> operator * (const Rational & a, const expressions::Expression<RationalContext, A, O> & b)
        { return expressions::Expression<RationalContext, std::pair<expressions::Wrapper<RationalContext>,
                                                                    expressions::Expression<RationalContext, A, O> >,
                                         expressions::MulOp>(std::make_pair(expressions::Wrapper<RationalContext>(a), b)); }
        template<class A, template<typename, typename> class O>
        inline expressions::Expression<RationalContext, std::pair<expressions::Wrapper<RationalContext>,
                                                                  expressions::Expression<RationalContext, A, O> >,
                                       expressions::DivOp> operator / (const Rational & a, const expressions::Expression<RationalContext, A, O> & b)
        { return expressions::Expression<RationalContext, std::pair<expressions::Wrapper<RationalContext>,
                                                                    expressions::Expression<RationalContext, A, O> >,
                                         expressions::DivOp>(std::make_pair(expressions::Wrapper<RationalContext>(a), b)); }
        template<class A, template<typename, typename> class O>
        inline expressions::Expression<RationalContext, std::pair<expressions::Wrapper<RationalContext>,
                                                                  expressions::Expression<RationalContext, A, O> >,
                                       expressions::ModOp> operator % (const Rational & a, const expressions::Expression<RationalContext, A, O> & b)
        { return expressions::Expression<RationalContext, std::pair<expressions::Wrapper<RationalContext>,
                                                                    expressions::Expression<RationalContext, A, O> >,
                                         expressions::ModOp>(std::make_pair(expressions::Wrapper<RationalContext>(a), b)); }
        
        namespace expressions
        {
            // Basic operations with first operand an expression
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<RationalContext, expressions::Expression<RationalContext, A, O>,
                                           expressions::NegOp> operator - (const expressions::Expression<RationalContext, A, O> & a)
            { return expressions::Expression<RationalContext, expressions::Expression<RationalContext, A, O>,
                                             expressions::NegOp>(a); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<RationalContext, std::pair<expressions::Expression<RationalContext, A, O>,
                                                                      expressions::Wrapper<RationalContext> >,
                                           expressions::AddOp> operator + (const expressions::Expression<RationalContext, A, O> & a,
                                                                           const arithmetic::Rational & b)
            { return expressions::Expression<RationalContext, std::pair<expressions::Expression<RationalContext, A, O>, expressions::Wrapper<RationalContext> >,
                                             expressions::AddOp>(std::make_pair(a, expressions::Wrapper<RationalContext>(b))); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<RationalContext, std::pair<expressions::Expression<RationalContext, A, O>,
                                                                      expressions::Wrapper<RationalContext> >,
                                           expressions::SubOp> operator - (const expressions::Expression<RationalContext, A, O> & a,
                                                                           const arithmetic::Rational & b)
            { return expressions::Expression<RationalContext, std::pair<expressions::Expression<RationalContext, A, O>,
                                                                        expressions::Wrapper<RationalContext> >,
                                             expressions::SubOp>(std::make_pair(a, expressions::Wrapper<RationalContext>(b))); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<RationalContext, std::pair<expressions::Expression<RationalContext, A, O>,
                                                                      expressions::Wrapper<RationalContext> >,
                                           expressions::MulOp> operator * (const expressions::Expression<RationalContext, A, O> & a,
                                                                           const arithmetic::Rational & b)
            { return expressions::Expression<RationalContext, std::pair<expressions::Expression<RationalContext, A, O>,
                                                                        expressions::Wrapper<RationalContext> >,
                                             expressions::MulOp>(std::make_pair(a, expressions::Wrapper<RationalContext>(b))); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<RationalContext, std::pair<expressions::Expression<RationalContext, A, O>,
                                                                      expressions::Wrapper<RationalContext> >,
                                           expressions::DivOp> operator / (const expressions::Expression<RationalContext, A, O> & a,
                                                                           const arithmetic::Rational & b)
            { return expressions::Expression<RationalContext, std::pair<expressions::Expression<RationalContext, A, O>,
                                                                        expressions::Wrapper<RationalContext> >,
                                             expressions::DivOp>(std::make_pair(a, expressions::Wrapper<RationalContext>(b))); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<RationalContext, std::pair<expressions::Expression<RationalContext, A, O>,
                                                                      expressions::Wrapper<RationalContext> >,
                                           expressions::ModOp> operator % (const expressions::Expression<RationalContext, A, O> & a,
                                                                           const arithmetic::Rational & b)
            { return expressions::Expression<RationalContext, std::pair<expressions::Expression<RationalContext, A, O>,
                                                                        expressions::Wrapper<RationalContext> >,
                                             expressions::ModOp>(std::make_pair(a, expressions::Wrapper<RationalContext>(b))); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<RationalContext, expressions::Expression<RationalContext, A, O>,
                                           expressions::ShiftCOp> operator << (const expressions::Expression<RationalContext, A, O> & a,
                                                                               signed long b)
            { return expressions::Expression<RationalContext, expressions::Expression<RationalContext, A, O>,
                                             expressions::ShiftCOp>(a, b); }
            
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<RationalContext, expressions::Expression<RationalContext, A, O>,
                                           expressions::ShiftCOp> operator >> (const expressions::Expression<RationalContext, A, O> & a,
                                                                               signed long b)
            { return expressions::Expression<RationalContext, expressions::Expression<RationalContext, A, O>,
                                             expressions::ShiftCOp>(a, -b); }
            
            // Basic operations with both operands an expression
            template<class A1, template<typename, typename> class O1, class A, template<typename, typename> class O>
            inline expressions::Expression<RationalContext, std::pair<expressions::Expression<RationalContext, A1, O1>,
                                                                      expressions::Expression<RationalContext, A, O> >,
                                           expressions::AddOp> operator + (const expressions::Expression<RationalContext, A1, O1> & a,
                                                                           const expressions::Expression<RationalContext, A, O> & b)
            { return expressions::Expression<RationalContext, std::pair<expressions::Expression<RationalContext, A1, O1>,
                                                                        expressions::Expression<RationalContext, A, O> >,
                                             expressions::AddOp>(std::make_pair(a, b)); }
            template<class A1, template<typename, typename> class O1, class A, template<typename, typename> class O>
            inline expressions::Expression<RationalContext, std::pair<expressions::Expression<RationalContext, A1, O1>,
                                                                      expressions::Expression<RationalContext, A, O> >,
                                           expressions::SubOp> operator - (const expressions::Expression<RationalContext, A1, O1> & a,
                                                                           const expressions::Expression<RationalContext, A, O> & b)
            { return expressions::Expression<RationalContext, std::pair<expressions::Expression<RationalContext, A1, O1>,
                                                                        expressions::Expression<RationalContext, A, O> >,
                                             expressions::SubOp>(std::make_pair(a, b)); }
            template<class A1, template<typename, typename> class O1, class A, template<typename, typename> class O>
            inline expressions::Expression<RationalContext, std::pair<expressions::Expression<RationalContext, A1, O1>,
                                                                      expressions::Expression<RationalContext, A, O> >,
                                           expressions::MulOp> operator * (const expressions::Expression<RationalContext, A1, O1> & a,
                                                                           const expressions::Expression<RationalContext, A, O> & b)
            { return expressions::Expression<RationalContext, std::pair<expressions::Expression<RationalContext, A1, O1>,
                                                                        expressions::Expression<RationalContext, A, O> >,
                                             expressions::MulOp>(std::make_pair(a, b)); }
            template<class A1, template<typename, typename> class O1, class A, template<typename, typename> class O>
            inline expressions::Expression<RationalContext, std::pair<expressions::Expression<RationalContext, A1, O1>,
                                                                      expressions::Expression<RationalContext, A, O> >,
                                           expressions::DivOp> operator / (const expressions::Expression<RationalContext, A1, O1> & a,
                                                                           const expressions::Expression<RationalContext, A, O> & b)
            { return expressions::Expression<RationalContext, std::pair<expressions::Expression<RationalContext, A1, O1>,
                                                                        expressions::Expression<RationalContext, A, O> >,
                                             expressions::DivOp>(std::make_pair(a, b)); }
            template<class A1, template<typename, typename> class O1, class A, template<typename, typename> class O>
            inline expressions::Expression<RationalContext, std::pair<expressions::Expression<RationalContext, A1, O1>,
                                                                      expressions::Expression<RationalContext, A, O> >,
                                           expressions::ModOp> operator % (const expressions::Expression<RationalContext, A1, O1> & a,
                                                                           const expressions::Expression<RationalContext, A, O> & b)
            { return expressions::Expression<RationalContext, std::pair<expressions::Expression<RationalContext, A1, O1>,
                                                                        expressions::Expression<RationalContext, A, O> >,
                                             expressions::ModOp>(std::make_pair(a, b)); }
        }
        
        /**@{
           \name Operators.
        */
        /**
           \brief Returns the absolute value of the given rational number.
           
           \param i The operand.
           \return The result.
        */
        inline expressions::Expression<RationalContext, expressions::Wrapper<RationalContext>,
                                       expressions::AbsOp> abs(const Rational & i)
        { return expressions::Expression<RationalContext, expressions::Wrapper<RationalContext>,
                                         expressions::AbsOp>(expressions::Wrapper<RationalContext>(i)); }
        /**
           \brief Returns the square of the given rational number.
           
           \param i The operand.
           \return The result.
        */
        inline expressions::Expression<RationalContext, expressions::Wrapper<RationalContext>,
                                       expressions::SquareOp> square(const Rational & i)
        { return expressions::Expression<RationalContext, expressions::Wrapper<RationalContext>,
                                         expressions::SquareOp>(expressions::Wrapper<RationalContext>(i)); }
        ///@}
        /**@{
           \name Exponentiation.
        */
        /**
           \brief Returns `a` raised to the power of `b`.
           
           \param a The base.
           \param b The exponent.
           \return The result.
        */
        inline expressions::Expression<RationalContext, expressions::Wrapper<RationalContext>,
                                       expressions::PowerCOp<signed long>::impl> power(const Rational & a, signed long b)
        { return expressions::Expression<RationalContext, expressions::Wrapper<RationalContext>,
                                         expressions::PowerCOp<signed long>::impl>(expressions::Wrapper<RationalContext>(a), b); }
        /**
           \brief Returns `a` raised to the power of `b`.
           
           \param a The base.
           \param b The exponent.
           \return The result.
        */
        inline expressions::Expression<RationalContext, expressions::Wrapper<RationalContext>,
                                       expressions::PowerCOp<unsigned long>::impl> power(const Rational & a, unsigned long b)
        { return expressions::Expression<RationalContext, expressions::Wrapper<RationalContext>,
                                         expressions::PowerCOp<unsigned long>::impl>(expressions::Wrapper<RationalContext>(a), b); }
        /**
           \brief Returns `a` raised to the power of `b`.
           
           \param a The base.
           \param b The exponent.
           \return The result.
        */
        inline expressions::Expression<RationalContext, std::pair<expressions::Wrapper<RationalContext>,
                                                                  expressions::Wrapper<IntegerContext> >,
                                       expressions::PowerOp> power(const Rational & a, const Integer & b)
        { return expressions::Expression<RationalContext, std::pair<expressions::Wrapper<RationalContext>,
                                                                    expressions::Wrapper<IntegerContext> >,
                                         expressions::PowerOp>(std::make_pair(expressions::Wrapper<RationalContext>(a),
                                                                              expressions::Wrapper<IntegerContext>(b))); }
        ///@}
        /**@{
           \name Operators.
        */
        /**
           \brief Returns the absolute value of the given rational number.
           
           \param i The operand.
           \param rc The context whose precision to use for the result.
           \return The result.
        */
        inline expressions::Expression<RationalContext, expressions::Wrapper<RationalContext>,
                                       expressions::AbsOp> abs(const Rational & i, const RationalContext & rc)
        { return expressions::Expression<RationalContext, expressions::Wrapper<RationalContext>,
                                         expressions::AbsOp>(expressions::Wrapper<RationalContext>(i)); }
        /**
           \brief Returns the square of the given rational number.
           
           \param i The operand.
           \param rc The context whose precision to use for the result.
           \return The result.
        */
        inline expressions::Expression<RationalContext, expressions::Wrapper<RationalContext>,
                                       expressions::SquareOp> square(const Rational & i, const RationalContext & rc)
        { return expressions::Expression<RationalContext, expressions::Wrapper<RationalContext>,
                                         expressions::SquareOp>(expressions::Wrapper<RationalContext>(i)); }
        ///@}
        /**@{
           \name Exponentiation.
        */
        /**
           \brief Returns `a` raised to the power of `b`.
           
           \param a The base.
           \param b The exponent.
           \param rc The context whose precision to use for the result.
           \return The result.
        */
        inline expressions::Expression<RationalContext, expressions::Wrapper<RationalContext>,
                                       expressions::PowerCOp<signed long>::impl> power(const Rational & a, signed long b, const RationalContext & rc)
        { return expressions::Expression<RationalContext, expressions::Wrapper<RationalContext>,
                                         expressions::PowerCOp<signed long>::impl>(expressions::Wrapper<RationalContext>(a), b); }
        /**
           \brief Returns `a` raised to the power of `b`.
           
           \param a The base.
           \param b The exponent.
           \param rc The context whose precision to use for the result.
           \return The result.
        */
        inline expressions::Expression<RationalContext, expressions::Wrapper<RationalContext>,
                                       expressions::PowerCOp<unsigned long>::impl> power(const Rational & a, unsigned long b, const RationalContext & rc)
        { return expressions::Expression<RationalContext, expressions::Wrapper<RationalContext>,
                                         expressions::PowerCOp<unsigned long>::impl>(expressions::Wrapper<RationalContext>(a), b); }
        /**
           \brief Returns `a` raised to the power of `b`.
           
           \param a The base.
           \param b The exponent.
           \param rc The context whose precision to use for the result.
           \return The result.
        */
        inline expressions::Expression<RationalContext, std::pair<expressions::Wrapper<RationalContext>,
                                                                  expressions::Wrapper<IntegerContext> >,
                                       expressions::PowerOp> power(const Rational & a, const Integer & b, const RationalContext & rc)
        { return expressions::Expression<RationalContext, std::pair<expressions::Wrapper<RationalContext>,
                                                                    expressions::Wrapper<IntegerContext> >,
                                         expressions::PowerOp>(std::make_pair(expressions::Wrapper<RationalContext>(a),
                                                                              expressions::Wrapper<IntegerContext>(b))); }
        ///@}
        
        // Functions applied to expressions
        template<class A, template<typename, typename> class O>
        inline expressions::Expression<RationalContext, expressions::Expression<RationalContext, A, O>,
                                       expressions::AbsOp> abs(const expressions::Expression<RationalContext, A, O> & b)
        { return expressions::Expression<RationalContext, expressions::Expression<RationalContext, A, O>,
                                         expressions::AbsOp>(b); }
        template<class A, template<typename, typename> class O>
        inline expressions::Expression<RationalContext, expressions::Expression<RationalContext, A, O>,
                                       expressions::SquareOp> square(const expressions::Expression<RationalContext, A, O> & b)
        { return expressions::Expression<RationalContext, expressions::Expression<RationalContext, A, O>,
                                         expressions::SquareOp>(b); }
        template<class A, template<typename, typename> class O>
        inline expressions::Expression<RationalContext, expressions::Expression<RationalContext, A, O>,
                                       expressions::PowerCOp<signed long>::impl> power(const expressions::Expression<RationalContext, A, O> & a,
                                                                                       signed long b)
        { return expressions::Expression<RationalContext, expressions::Expression<RationalContext, A, O>,
                                         expressions::PowerCOp<signed long>::impl>(a, b); }
        template<class A, template<typename, typename> class O>
        inline expressions::Expression<RationalContext, expressions::Expression<RationalContext, A, O>,
                                       expressions::PowerCOp<unsigned long>::impl> power(const expressions::Expression<RationalContext, A, O> & a,
                                                                                         unsigned long b)
        { return expressions::Expression<RationalContext, expressions::Expression<RationalContext, A, O>,
                                         expressions::PowerCOp<unsigned long>::impl>(a, b); }
        template<class A, template<typename, typename> class O>
        inline expressions::Expression<RationalContext, std::pair<expressions::Expression<RationalContext, A, O>,
                                                                  expressions::Wrapper<IntegerContext> >,
                                       expressions::PowerOp> power(const expressions::Expression<RationalContext, A, O> & a,
                                                                   const Integer & b)
        { return expressions::Expression<RationalContext, std::pair<expressions::Expression<RationalContext, A, O>,
                                                                    expressions::Wrapper<IntegerContext> >,
                                         expressions::PowerOp>(std::make_pair(a, expressions::Wrapper<IntegerContext>(b))); }
        template<class A, template<typename, typename> class O>
        inline expressions::Expression<RationalContext, expressions::Expression<RationalContext, A, O>,
                                       expressions::AbsOp> abs(const expressions::Expression<RationalContext, A, O> & b,
                                                               const RationalContext &)
        { return expressions::Expression<RationalContext, expressions::Expression<RationalContext, A, O>,
                                         expressions::AbsOp>(b); }
        template<class A, template<typename, typename> class O>
        inline expressions::Expression<RationalContext, expressions::Expression<RationalContext, A, O>,
                                       expressions::SquareOp> square(const expressions::Expression<RationalContext, A, O> & b,
                                                                     const RationalContext &)
        { return expressions::Expression<RationalContext, expressions::Expression<RationalContext, A, O>,
                                         expressions::SquareOp>(b); }
        template<class A, template<typename, typename> class O>
        inline expressions::Expression<RationalContext, expressions::Expression<RationalContext, A, O>,
                                       expressions::PowerCOp<signed long>::impl> power(const expressions::Expression<RationalContext, A, O> & a,
                                                                                       signed long b, const RationalContext &)
        { return expressions::Expression<RationalContext, expressions::Expression<RationalContext, A, O>,
                                         expressions::PowerCOp<signed long>::impl>(a, b); }
        template<class A, template<typename, typename> class O>
        inline expressions::Expression<RationalContext, expressions::Expression<RationalContext, A, O>,
                                       expressions::PowerCOp<unsigned long>::impl> power(const expressions::Expression<RationalContext, A, O> & a,
                                                                                         unsigned long b, const RationalContext &)
        { return expressions::Expression<RationalContext, expressions::Expression<RationalContext, A, O>,
                                         expressions::PowerCOp<unsigned long>::impl>(a, b); }
        template<class A, template<typename, typename> class O>
        inline expressions::Expression<RationalContext, std::pair<expressions::Expression<RationalContext, A, O>,
                                                                  expressions::Wrapper<IntegerContext> >,
                                       expressions::PowerOp> power(const expressions::Expression<RationalContext, A, O> & a,
                                                                   const Integer & b, const RationalContext &)
        { return expressions::Expression<RationalContext, std::pair<expressions::Expression<RationalContext, A, O>,
                                                                    expressions::Wrapper<IntegerContext> >,
                                         expressions::PowerOp>(std::make_pair(a, expressions::Wrapper<IntegerContext>(b))); }
        
        namespace expressions
        {
            // Functions applied to expressions
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<RationalContext, expressions::Expression<RationalContext, A, O>,
                                           expressions::AbsOp> abs(const expressions::Expression<RationalContext, A, O> & b)
            { return expressions::Expression<RationalContext, expressions::Expression<RationalContext, A, O>,
                                             expressions::AbsOp>(b); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<RationalContext, expressions::Expression<RationalContext, A, O>,
                                           expressions::SquareOp> square(const expressions::Expression<RationalContext, A, O> & b)
            { return expressions::Expression<RationalContext, expressions::Expression<RationalContext, A, O>,
                                             expressions::SquareOp>(b); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<RationalContext, expressions::Expression<RationalContext, A, O>,
                                           expressions::PowerCOp<signed long>::impl> power(const expressions::Expression<RationalContext, A, O> & a,
                                                                                           signed long b)
            { return expressions::Expression<RationalContext, expressions::Expression<RationalContext, A, O>,
                                             expressions::PowerCOp<signed long>::impl>(a, b); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<RationalContext, expressions::Expression<RationalContext, A, O>,
                                           expressions::PowerCOp<unsigned long>::impl> power(const expressions::Expression<RationalContext, A, O> & a,
                                                                                             unsigned long b)
            { return expressions::Expression<RationalContext, expressions::Expression<RationalContext, A, O>,
                                             expressions::PowerCOp<unsigned long>::impl>(a, b); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<RationalContext, std::pair<expressions::Expression<RationalContext, A, O>,
                                                                      expressions::Wrapper<IntegerContext> >,
                                           expressions::PowerOp> power(const expressions::Expression<RationalContext, A, O> & a,
                                                                       const Integer & b)
            { return expressions::Expression<RationalContext, std::pair<expressions::Expression<RationalContext, A, O>,
                                                                        expressions::Wrapper<IntegerContext> >,
                                             expressions::PowerOp>(std::make_pair(a, expressions::Wrapper<IntegerContext>(b))); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<RationalContext, expressions::Expression<RationalContext, A, O>,
                                           expressions::AbsOp> abs(const expressions::Expression<RationalContext, A, O> & b,
                                                                   const RationalContext &)
            { return expressions::Expression<RationalContext, expressions::Expression<RationalContext, A, O>,
                                             expressions::AbsOp>(b); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<RationalContext, expressions::Expression<RationalContext, A, O>,
                                           expressions::SquareOp> square(const expressions::Expression<RationalContext, A, O> & b,
                                                                         const RationalContext &)
            { return expressions::Expression<RationalContext, expressions::Expression<RationalContext, A, O>,
                                             expressions::SquareOp>(b); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<RationalContext, expressions::Expression<RationalContext, A, O>,
                                           expressions::PowerCOp<signed long>::impl> power(const expressions::Expression<RationalContext, A, O> & a,
                                                                                           signed long b, const RationalContext &)
            { return expressions::Expression<RationalContext, expressions::Expression<RationalContext, A, O>,
                                             expressions::PowerCOp<signed long>::impl>(a, b); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<RationalContext, expressions::Expression<RationalContext, A, O>,
                                           expressions::PowerCOp<unsigned long>::impl> power(const expressions::Expression<RationalContext, A, O> & a,
                                                                                             unsigned long b, const RationalContext &)
            { return expressions::Expression<RationalContext, expressions::Expression<RationalContext, A, O>,
                                             expressions::PowerCOp<unsigned long>::impl>(a, b); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<RationalContext, std::pair<expressions::Expression<RationalContext, A, O>,
                                                                      expressions::Wrapper<IntegerContext> >,
                                           expressions::PowerOp> power(const expressions::Expression<RationalContext, A, O> & a,
                                                                       const Integer & b, const RationalContext &)
            { return expressions::Expression<RationalContext, std::pair<expressions::Expression<RationalContext, A, O>,
                                                                        expressions::Wrapper<IntegerContext> >,
                                             expressions::PowerOp>(std::make_pair(a, expressions::Wrapper<IntegerContext>(b))); }
            
            template<class A, template<typename, typename> class O>
            inline std::ostream & operator << (std::ostream & s, const expressions::Expression<RationalContext, A, O> & E)
            {
                return s << Rational(E);
            }
        }
        
        inline void addmul(Rational & r, const Rational & a, const Rational & b)
        {
            Rational p = a * b;
            r.add_r_neq_1st(p, r);
        }
        
        inline void submul(Rational & r, const Rational & a, const Rational & b)
        {
            Rational p = a * b;
            r.sub_r_neq_1st(p, r, true);
        }
    }
}

#endif

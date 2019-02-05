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

#ifndef PLLL_INCLUDE_GUARD__ARITHMETIC_GMP_INTEGER_OPS_HPP
#define PLLL_INCLUDE_GUARD__ARITHMETIC_GMP_INTEGER_OPS_HPP

/**
   \file
   \brief Operator definitions for integers.
   
   This header contains templates and instantiations to implement all operations on integers.
*/
namespace plll
{
    namespace arithmetic
    {
        namespace expressions
        {
            // Adding integer and (un)signed long
            template<class D>
            void do_assign(arithmetic::Integer & x,
                           const AddOp<IntegerContext, std::pair<D, expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> > > &,
                           const std::pair<D, expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> > & data)
            {
                add_ui(x, data.first.evaluate(), data.second.data().evaluate());
            }
            
            template<class D>
            void do_assign(arithmetic::Integer & x,
                           const AddOp<IntegerContext, std::pair<D, expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> > > &,
                           const std::pair<D, expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> > & data)
            {
                add_si(x, data.first.evaluate(), data.second.data().evaluate());
            }
            
            template<class D>
            void do_assign(arithmetic::Integer & x,
                           const AddOp<IntegerContext, std::pair<expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context>, D> > &,
                           const std::pair< expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context>, D> & data)
            {
                add_ui(x, data.second.evaluate(), data.first.data().evaluate());
            }
            
            template<class D>
            void do_assign(arithmetic::Integer & x,
                           const AddOp<IntegerContext, std::pair<expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context>, D> > &,
                           const std::pair<expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context>, D> & data)
            {
                add_si(x, data.second.evaluate(), data.first.data().evaluate());
            }
            
            // Subtracting integer and (un)signed long
            template<class D>
            void do_assign(arithmetic::Integer & x,
                           const SubOp<IntegerContext, std::pair<D, expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> > > &,
                           const std::pair<D, expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> > & data)
            {
                sub_ui(x, data.first.evaluate(), data.second.data().evaluate());
            }
            
            template<class D>
            void do_assign(arithmetic::Integer & x,
                           const SubOp<IntegerContext, std::pair<D, expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> > > &,
                           const std::pair<D, expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> > & data)
            {
                sub_si(x, data.first.evaluate(), data.second.data().evaluate());
            }
            
            template<class D>
            void do_assign(arithmetic::Integer & x,
                           const SubOp<IntegerContext, std::pair<expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context>, D> > &,
                           const std::pair< expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context>, D> & data)
            {
                ui_sub(x, data.first.data().evaluate(), data.second.evaluate());
            }
            
            template<class D>
            void do_assign(arithmetic::Integer & x,
                           const SubOp<IntegerContext, std::pair<expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context>, D> > &,
                           const std::pair<expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context>, D> & data)
            {
                si_sub(x, data.first.data().evaluate(), data.second.evaluate());
            }
            
            // Multiplying integer and (un)signed long
            template<class D>
            void do_assign(arithmetic::Integer & x,
                           const MulOp<IntegerContext, std::pair<D, expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> > > &,
                           const std::pair<D, expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> > & data)
            {
                mul_ui(x, data.first.evaluate(), data.second.data().evaluate());
            }
            
            template<class D>
            void do_assign(arithmetic::Integer & x,
                           const MulOp<IntegerContext, std::pair<D, expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> > > &,
                           const std::pair<D, expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> > & data)
            {
                mul_si(x, data.first.evaluate(), data.second.data().evaluate());
            }
            
            template<class D>
            void do_assign(arithmetic::Integer & x,
                           const MulOp<IntegerContext, std::pair<expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context>, D> > &,
                           const std::pair< expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context>, D> & data)
            {
                mul_ui(x, data.second.evaluate(), data.first.data().evaluate());
            }
            
            template<class D>
            void do_assign(arithmetic::Integer & x,
                           const MulOp<IntegerContext, std::pair<expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context>, D> > &,
                           const std::pair<expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context>, D> & data)
            {
                mul_si(x, data.second.evaluate(), data.first.data().evaluate());
            }
            
            // Dividing integer by (un)signed long
            template<class D>
            void do_assign(arithmetic::Integer & x,
                           const DivOp<IntegerContext, std::pair<D, expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> > > &,
                           const std::pair<D, expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> > & data)
            {
                div_ui(x, data.first.evaluate(), data.second.data().evaluate());
            }
            
            template<class D>
            void do_assign(arithmetic::Integer & x,
                           const DivOp<IntegerContext, std::pair<D, expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> > > &,
                           const std::pair<D, expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> > & data)
            {
                div_si(x, data.first.evaluate(), data.second.data().evaluate());
            }
            
            // Modulo integer by (un)signed long
            template<class D>
            void do_assign(arithmetic::Integer & x,
                           const ModOp<IntegerContext, std::pair<D, expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> > > &,
                           const std::pair<D, expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> > & data)
            {
                mod_ui(x, data.first.evaluate(), data.second.data().evaluate());
            }
            
            template<class D>
            void do_assign(arithmetic::Integer & x,
                           const ModOp<IntegerContext, std::pair<D, expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> > > &,
                           const std::pair<D, expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> > & data)
            {
                mod_si(x, data.first.evaluate(), data.second.data().evaluate());
            }
            
            // GCD with (un)signed long
            template<class D>
            void do_assign(arithmetic::Integer & x,
                           const GCDOp<IntegerContext, std::pair<D, expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> > > &,
                           const std::pair<D, expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> > & data)
            {
                gcd_ui(x, data.first.evaluate(), data.second.data().evaluate());
            }
            
            template<class D>
            void do_assign(arithmetic::Integer & x,
                           const GCDOp<IntegerContext, std::pair<D, expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> > > &,
                           const std::pair<D, expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> > & data)
            {
                gcd_si(x, data.first.evaluate(), data.second.data().evaluate());
            }
            
            template<class D>
            void do_assign(arithmetic::Integer & x,
                           const GCDOp<IntegerContext, std::pair<expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context>, D> > &,
                           const std::pair< expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context>, D> & data)
            {
                gcd_ui(x, data.second.evaluate(), data.first.data().evaluate());
            }
            
            template<class D>
            void do_assign(arithmetic::Integer & x,
                           const GCDOp<IntegerContext, std::pair<expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context>, D> > &,
                           const std::pair<expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context>, D> & data)
            {
                gcd_si(x, data.second.evaluate(), data.first.data().evaluate());
            }
            
            // LCM with (un)signed long
            template<class D>
            void do_assign(arithmetic::Integer & x,
                           const LCMOp<IntegerContext, std::pair<D, expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> > > &,
                           const std::pair<D, expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> > & data)
            {
                lcm_ui(x, data.first.evaluate(), data.second.data().evaluate());
            }
            
            template<class D>
            void do_assign(arithmetic::Integer & x,
                           const LCMOp<IntegerContext, std::pair<D, expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> > > &,
                           const std::pair<D, expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> > & data)
            {
                lcm_si(x, data.first.evaluate(), data.second.data().evaluate());
            }
            
            template<class D>
            void do_assign(arithmetic::Integer & x,
                           const LCMOp<IntegerContext, std::pair<expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context>, D> > &,
                           const std::pair< expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context>, D> & data)
            {
                lcm_ui(x, data.second.evaluate(), data.first.data().evaluate());
            }
            
            template<class D>
            void do_assign(arithmetic::Integer & x,
                           const LCMOp<IntegerContext, std::pair<expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context>, D> > &,
                           const std::pair<expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context>, D> & data)
            {
                lcm_si(x, data.second.evaluate(), data.first.data().evaluate());
            }
        }
        
        /**@{
           \name Operators.
        */
        
        /**
           \brief Negates the integer.
           
           \param a The integer.
           \return The result.
        */
        inline expressions::Expression<IntegerContext, expressions::Wrapper<IntegerContext>,
                                                       expressions::NegOp> operator - (const Integer & a)
        { return expressions::Expression<IntegerContext, expressions::Wrapper<IntegerContext>,
                                                         expressions::NegOp>(expressions::Wrapper<IntegerContext>(a)); }
        /**
           \brief Bitwise negates the integer.
           
           \param a The integer.
           \return The result.
        */
        inline expressions::Expression<IntegerContext, expressions::Wrapper<IntegerContext>,
                                                        expressions::BitInvOp> operator ~ (const Integer & a)
        { return expressions::Expression<IntegerContext, expressions::Wrapper<IntegerContext>,
                                                          expressions::BitInvOp>(expressions::Wrapper<IntegerContext>(a)); }
        /**
           \brief Adds the two integers and returns the result.
           
           \param a The first operand.
           \param b The second operand.
           \return The result.
        */
        inline expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                 expressions::Wrapper<IntegerContext> >,
                                                       expressions::AddOp> operator + (const Integer & a, const Integer & b)
        { return expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                   expressions::Wrapper<IntegerContext> >,
                                                         expressions::AddOp>(std::make_pair(expressions::Wrapper<IntegerContext>(a),
                                                                                            expressions::Wrapper<IntegerContext>(b))); }
        /**
           \brief Subtracts the second from the first integer and returns the result.
           
           \param a The first operand.
           \param b The second operand.
           \return The result.
        */
        inline expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                 expressions::Wrapper<IntegerContext> >,
                                                       expressions::SubOp> operator - (const Integer & a, const Integer & b)
        { return expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                   expressions::Wrapper<IntegerContext> >,
                                                         expressions::SubOp>(std::make_pair(expressions::Wrapper<IntegerContext>(a),
                                                                                            expressions::Wrapper<IntegerContext>(b))); }
        /**
           \brief Multiplies the two integers and returns the result.
           
           \param a The first operand.
           \param b The second operand.
           \return The result.
        */
        inline expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                 expressions::Wrapper<IntegerContext> >,
                                                       expressions::MulOp> operator * (const Integer & a, const Integer & b)
        { return expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                   expressions::Wrapper<IntegerContext> >,
                                                         expressions::MulOp>(std::make_pair(expressions::Wrapper<IntegerContext>(a),
                                                                                            expressions::Wrapper<IntegerContext>(b))); }
        /**
           \brief Divides the first by the second integer and returns the result.
           
           \param a The first operand.
           \param b The second operand.
           \return The result.
        */
        inline expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                 expressions::Wrapper<IntegerContext> >,
                                                       expressions::DivOp> operator / (const Integer & a, const Integer & b)
        { return expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                   expressions::Wrapper<IntegerContext> >,
                                                         expressions::DivOp>(std::make_pair(expressions::Wrapper<IntegerContext>(a),
                                                                                            expressions::Wrapper<IntegerContext>(b))); }
        /**
           \brief Divides the first by the second integer and returns the remainder.
           
           \param a The first operand.
           \param b The second operand.
           \return The result.
        */
        inline expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                 expressions::Wrapper<IntegerContext> >,
                                                       expressions::ModOp> operator % (const Integer & a, const Integer & b)
        { return expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                   expressions::Wrapper<IntegerContext> >,
                                                         expressions::ModOp>(std::make_pair(expressions::Wrapper<IntegerContext>(a),
                                                                                            expressions::Wrapper<IntegerContext>(b))); }
        /**
           \brief Computes the bitwise and of the two integers and returns the result.
           
           \param a The first operand.
           \param b The second operand.
           \return The result.
        */
        inline expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                 expressions::Wrapper<IntegerContext> >,
                                                       expressions::AndOp> operator & (const Integer & a, const Integer & b)
        { return expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                   expressions::Wrapper<IntegerContext> >,
                                                         expressions::AndOp>(std::make_pair(expressions::Wrapper<IntegerContext>(a),
                                                                                            expressions::Wrapper<IntegerContext>(b))); }
        /**
           \brief Computes the bitwise or of the two integers and returns the result.
           
           \param a The first operand.
           \param b The second operand.
           \return The result.
        */
        inline expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                 expressions::Wrapper<IntegerContext> >,
                                                       expressions::OrOp> operator | (const Integer & a, const Integer & b)
        { return expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                   expressions::Wrapper<IntegerContext> >,
                                                         expressions::OrOp>(std::make_pair(expressions::Wrapper<IntegerContext>(a),
                                                                                           expressions::Wrapper<IntegerContext>(b))); }
        /**
           \brief Computes the bitwise exclusive or of the two integers and returns the result.
           
           \param a The first operand.
           \param b The second operand.
           \return The result.
        */
        inline expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                 expressions::Wrapper<IntegerContext> >,
                                                       expressions::XOROp> operator ^ (const Integer & a, const Integer & b)
        { return expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                   expressions::Wrapper<IntegerContext> >,
                                                         expressions::XOROp>(std::make_pair(expressions::Wrapper<IntegerContext>(a),
                                                                                            expressions::Wrapper<IntegerContext>(b))); }
        /**
           \brief Computes the bitwise left shift of the first integer by the number of bits given
                  by the second integer, and returns the result.
           
           This effectively multiplies the first integer by 2 to the power of the second integer.
           
           \param a The first operand.
           \param b The second operand.
           \return The result.
        */
        inline expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                 expressions::Wrapper<IntegerContext> >,
                                                       expressions::ShLOp> operator << (const Integer & a, const Integer & b)
        { return expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                   expressions::Wrapper<IntegerContext> >,
                                                         expressions::ShLOp>(std::make_pair(expressions::Wrapper<IntegerContext>(a),
                                                                                            expressions::Wrapper<IntegerContext>(b))); }
        /**
           \brief Computes the bitwise right shift of the first integer by the number of bits given
                  by the second integer, and returns the result.
           
           This effectively divides the first integer by 2 to the power of the second integer, and
           rounds towards zero.
           
           \param a The first operand.
           \param b The second operand.
           \return The result.
        */
        inline expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                 expressions::Wrapper<IntegerContext> >,
                                                       expressions::ShROp> operator >> (const Integer & a, const Integer & b)
        { return expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                   expressions::Wrapper<IntegerContext> >,
                                                         expressions::ShROp>(std::make_pair(expressions::Wrapper<IntegerContext>(a),
                                                                                            expressions::Wrapper<IntegerContext>(b))); }
        /**
           \brief Computes the bitwise left shift of the first integer by the number of bits given
                  by the second integer, and returns the result.
           
           This effectively multiplies the first integer by 2 to the power of the second integer.
           
           \param a The first operand.
           \param b The second operand.
           \return The result.
        */
        inline expressions::Expression<IntegerContext, expressions::Wrapper<IntegerContext>,
                                                       expressions::ShiftCOp> operator << (const Integer & a, signed long b)
        { return expressions::Expression<IntegerContext, expressions::Wrapper<IntegerContext>,
                                                         expressions::ShiftCOp>(expressions::Wrapper<IntegerContext>(a), b); }
        /**
           \brief Computes the bitwise right shift of the first integer by the number of bits given
                  by the second integer, and returns the result.
           
           This effectively divides the first integer by 2 to the power of the second integer, and
           rounds towards zero.
           
           \param a The first operand.
           \param b The second operand.
           \return The result.
        */
        inline expressions::Expression<IntegerContext, expressions::Wrapper<IntegerContext>,
                                                       expressions::ShiftCOp> operator >> (const Integer & a, signed long b)
        { return expressions::Expression<IntegerContext, expressions::Wrapper<IntegerContext>,
                                                         expressions::ShiftCOp>(expressions::Wrapper<IntegerContext>(a), -b); }
        ///@}
        
        // second operand is expression
        template<class A2, template<typename, typename> class O2>
        inline expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                 expressions::Expression<IntegerContext, A2, O2> >,
                                       expressions::AddOp> operator + (const Integer & a, const expressions::Expression<IntegerContext, A2, O2> & b)
        { return expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                   expressions::Expression<IntegerContext, A2, O2> >,
                                         expressions::AddOp>(std::make_pair(expressions::Wrapper<IntegerContext>(a), b)); }
        template<class A2, template<typename, typename> class O2>
        inline expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                 expressions::Expression<IntegerContext, A2, O2> >,
                                       expressions::SubOp> operator - (const Integer & a, const expressions::Expression<IntegerContext, A2, O2> & b)
        { return expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                   expressions::Expression<IntegerContext, A2, O2> >,
                                         expressions::SubOp>(std::make_pair(expressions::Wrapper<IntegerContext>(a), b)); }
        template<class A2, template<typename, typename> class O2>
        inline expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                 expressions::Expression<IntegerContext, A2, O2> >,
                                       expressions::MulOp> operator * (const Integer & a, const expressions::Expression<IntegerContext, A2, O2> & b)
        { return expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                   expressions::Expression<IntegerContext, A2, O2> >,
                                         expressions::MulOp>(std::make_pair(expressions::Wrapper<IntegerContext>(a), b)); }
        template<class A2, template<typename, typename> class O2>
        inline expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                 expressions::Expression<IntegerContext, A2, O2> >,
                                       expressions::DivOp> operator / (const Integer & a, const expressions::Expression<IntegerContext, A2, O2> & b)
        { return expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                   expressions::Expression<IntegerContext, A2, O2> >,
                                         expressions::DivOp>(std::make_pair(expressions::Wrapper<IntegerContext>(a), b)); }
        template<class A2, template<typename, typename> class O2>
        inline expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                 expressions::Expression<IntegerContext, A2, O2> >,
                                       expressions::ModOp> operator % (const Integer & a, const expressions::Expression<IntegerContext, A2, O2> & b)
        { return expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                   expressions::Expression<IntegerContext, A2, O2> >,
                                         expressions::ModOp>(std::make_pair(expressions::Wrapper<IntegerContext>(a), b)); }
        template<class A2, template<typename, typename> class O2>
        inline expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                 expressions::Expression<IntegerContext, A2, O2> >,
                                       expressions::AndOp> operator & (const Integer & a, const expressions::Expression<IntegerContext, A2, O2> & b)
        { return expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                   expressions::Expression<IntegerContext, A2, O2> >,
                                         expressions::AndOp>(std::make_pair(expressions::Wrapper<IntegerContext>(a), b)); }
        template<class A2, template<typename, typename> class O2>
        inline expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                 expressions::Expression<IntegerContext, A2, O2> >,
                                       expressions::OrOp> operator | (const Integer & a, const expressions::Expression<IntegerContext, A2, O2> & b)
        { return expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                   expressions::Expression<IntegerContext, A2, O2> >,
                                         expressions::OrOp>(std::make_pair(expressions::Wrapper<IntegerContext>(a), b)); }
        template<class A2, template<typename, typename> class O2>
        inline expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                 expressions::Expression<IntegerContext, A2, O2> >,
                                       expressions::XOROp> operator ^ (const Integer & a, const expressions::Expression<IntegerContext, A2, O2> & b)
        { return expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                   expressions::Expression<IntegerContext, A2, O2> >,
                                         expressions::XOROp>(std::make_pair(expressions::Wrapper<IntegerContext>(a), b)); }
        template<class A2, template<typename, typename> class O2>
        inline expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                 expressions::Expression<IntegerContext, A2, O2> >,
                                       expressions::ShLOp> operator << (const Integer & a, const expressions::Expression<IntegerContext, A2, O2> & b)
        { return expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                   expressions::Expression<IntegerContext, A2, O2> >,
                                         expressions::ShLOp>(std::make_pair(expressions::Wrapper<IntegerContext>(a), b)); }
        template<class A2, template<typename, typename> class O2>
        inline expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                 expressions::Expression<IntegerContext, A2, O2> >,
                                       expressions::ShROp> operator >> (const Integer & a, const expressions::Expression<IntegerContext, A2, O2> & b)
        { return expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                   expressions::Expression<IntegerContext, A2, O2> >,
                                         expressions::ShROp>(std::make_pair(expressions::Wrapper<IntegerContext>(a), b)); }
        
        namespace expressions
        {
            // first operand is expression
            template<class A1, template<typename, typename> class O1>
            inline expressions::Expression<IntegerContext, expressions::Expression<IntegerContext, A1, O1>,
                                           expressions::NegOp> operator - (const expressions::Expression<IntegerContext, A1, O1> & a)
            { return expressions::Expression<IntegerContext, expressions::Expression<IntegerContext, A1, O1>,
                                             expressions::NegOp>(a); }
            template<class A1, template<typename, typename> class O1>
            inline expressions::Expression<IntegerContext, expressions::Expression<IntegerContext, A1, O1>,
                                           expressions::BitInvOp> operator ~ (const expressions::Expression<IntegerContext, A1, O1> & a)
            { return expressions::Expression<IntegerContext, expressions::Expression<IntegerContext, A1, O1>,
                                             expressions::BitInvOp>(a); }
            template<class A1, template<typename, typename> class O1>
            inline expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A1, O1>,
                                                                     expressions::Wrapper<IntegerContext> >,
                                           expressions::AddOp> operator + (const expressions::Expression<IntegerContext, A1, O1> & a, const Integer & b)
            { return expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A1, O1>,
                                                                       expressions::Wrapper<IntegerContext> >,
                                             expressions::AddOp>(std::make_pair(a, expressions::Wrapper<IntegerContext>(b))); }
            template<class A1, template<typename, typename> class O1>
            inline expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A1, O1>,
                                                                     expressions::Wrapper<IntegerContext> >,
                                           expressions::SubOp> operator - (const expressions::Expression<IntegerContext, A1, O1> & a, const Integer & b)
            { return expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A1, O1>,
                                                                       expressions::Wrapper<IntegerContext> >,
                                             expressions::SubOp>(std::make_pair(a, expressions::Wrapper<IntegerContext>(b))); }
            template<class A1, template<typename, typename> class O1>
            inline expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A1, O1>,
                                                                     expressions::Wrapper<IntegerContext> >,
                                           expressions::MulOp> operator * (const expressions::Expression<IntegerContext, A1, O1> & a, const Integer & b)
            { return expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A1, O1>,
                                                                       expressions::Wrapper<IntegerContext> >,
                                             expressions::MulOp>(std::make_pair(a, expressions::Wrapper<IntegerContext>(b))); }
            template<class A1, template<typename, typename> class O1>
            inline expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A1, O1>,
                                                                     expressions::Wrapper<IntegerContext> >,
                                           expressions::DivOp> operator / (const expressions::Expression<IntegerContext, A1, O1> & a, const Integer & b)
            { return expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A1, O1>,
                                                                       expressions::Wrapper<IntegerContext> >,
                                             expressions::DivOp>(std::make_pair(a, expressions::Wrapper<IntegerContext>(b))); }
            template<class A1, template<typename, typename> class O1>
            inline expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A1, O1>,
                                                                     expressions::Wrapper<IntegerContext> >,
                                           expressions::ModOp> operator % (const expressions::Expression<IntegerContext, A1, O1> & a, const Integer & b)
            { return expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A1, O1>,
                                                                       expressions::Wrapper<IntegerContext> >,
                                             expressions::ModOp>(std::make_pair(a, expressions::Wrapper<IntegerContext>(b))); }
            template<class A1, template<typename, typename> class O1>
            inline expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A1, O1>,
                                                                     expressions::Wrapper<IntegerContext> >,
                                           expressions::AndOp> operator & (const expressions::Expression<IntegerContext, A1, O1> & a, const Integer & b)
            { return expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A1, O1>,
                                                                       expressions::Wrapper<IntegerContext> >,
                                             expressions::AndOp>(std::make_pair(a, expressions::Wrapper<IntegerContext>(b))); }
            template<class A1, template<typename, typename> class O1>
            inline expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A1, O1>,
                                                                     expressions::Wrapper<IntegerContext> >,
                                           expressions::OrOp> operator | (const expressions::Expression<IntegerContext, A1, O1> & a, const Integer & b)
            { return expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A1, O1>,
                                                                       expressions::Wrapper<IntegerContext> >,
                                             expressions::OrOp>(std::make_pair(a, expressions::Wrapper<IntegerContext>(b))); }
            template<class A1, template<typename, typename> class O1>
            inline expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A1, O1>,
                                                                     expressions::Wrapper<IntegerContext> >,
                                           expressions::XOROp> operator ^ (const expressions::Expression<IntegerContext, A1, O1> & a, const Integer & b)
            { return expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A1, O1>,
                                                                       expressions::Wrapper<IntegerContext> >,
                                             expressions::XOROp>(std::make_pair(a, expressions::Wrapper<IntegerContext>(b))); }
            template<class A1, template<typename, typename> class O1>
            inline expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A1, O1>,
                                                                     expressions::Wrapper<IntegerContext> >,
                                           expressions::ShLOp> operator << (const expressions::Expression<IntegerContext, A1, O1> & a, const Integer & b)
            { return expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A1, O1>,
                                                                       expressions::Wrapper<IntegerContext> >,
                                             expressions::ShLOp>(std::make_pair(a, expressions::Wrapper<IntegerContext>(b))); }
            template<class A1, template<typename, typename> class O1>
            inline expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A1, O1>,
                                                                     expressions::Wrapper<IntegerContext> >,
                                           expressions::ShROp> operator >> (const expressions::Expression<IntegerContext, A1, O1> & a, const Integer & b)
            { return expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A1, O1>,
                                                                       expressions::Wrapper<IntegerContext> >,
                                             expressions::ShROp>(std::make_pair(a, expressions::Wrapper<IntegerContext>(b))); }
            template<class A1, template<typename, typename> class O1>
            inline expressions::Expression<IntegerContext, expressions::Expression<IntegerContext, A1, O1>,
                                           expressions::ShiftCOp> operator << (const expressions::Expression<IntegerContext, A1, O1> & a, signed long b)
            { return expressions::Expression<IntegerContext, expressions::Expression<IntegerContext, A1, O1>,
                                             expressions::ShiftCOp>(a, b); }
            template<class A1, template<typename, typename> class O1>
            inline expressions::Expression<IntegerContext, expressions::Expression<IntegerContext, A1, O1>,
                                           expressions::ShiftCOp> operator >> (const expressions::Expression<IntegerContext, A1, O1> & a, signed long b)
            { return expressions::Expression<IntegerContext, expressions::Expression<IntegerContext, A1, O1>,
                                             expressions::ShiftCOp>(a, -b); }
            
            // both operands are expressions
            template<class A1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
            inline expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A1, O1>,
                                                                     expressions::Expression<IntegerContext, A2, O2> >,
                                           expressions::AddOp> operator + (const expressions::Expression<IntegerContext, A1, O1> & a,
                                                                           const expressions::Expression<IntegerContext, A2, O2> & b)
            { return expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A1, O1>,
                                                                       expressions::Expression<IntegerContext, A2, O2> >,
                                             expressions::AddOp>(std::make_pair(a, b)); }
            template<class A1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
            inline expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A1, O1>,
                                                                     expressions::Expression<IntegerContext, A2, O2> >,
                                           expressions::SubOp> operator - (const expressions::Expression<IntegerContext, A1, O1> & a,
                                                                           const expressions::Expression<IntegerContext, A2, O2> & b)
            { return expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A1, O1>,
                                                                       expressions::Expression<IntegerContext, A2, O2> >,
                                             expressions::SubOp>(std::make_pair(a, b)); }
            template<class A1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
            inline expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A1, O1>,
                                                                     expressions::Expression<IntegerContext, A2, O2> >,
                                           expressions::MulOp> operator * (const expressions::Expression<IntegerContext, A1, O1> & a,
                                                                           const expressions::Expression<IntegerContext, A2, O2> & b)
            { return expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A1, O1>,
                                                                       expressions::Expression<IntegerContext, A2, O2> >,
                                             expressions::MulOp>(std::make_pair(a, b)); }
            template<class A1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
            inline expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A1, O1>,
                                                                     expressions::Expression<IntegerContext, A2, O2> >,
                                           expressions::DivOp> operator / (const expressions::Expression<IntegerContext, A1, O1> & a,
                                                                           const expressions::Expression<IntegerContext, A2, O2> & b)
            { return expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A1, O1>,
                                                                       expressions::Expression<IntegerContext, A2, O2> >,
                                             expressions::DivOp>(std::make_pair(a, b)); }
            template<class A1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
            inline expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A1, O1>,
                                                                     expressions::Expression<IntegerContext, A2, O2> >,
                                           expressions::ModOp> operator % (const expressions::Expression<IntegerContext, A1, O1> & a,
                                                                           const expressions::Expression<IntegerContext, A2, O2> & b)
            { return expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A1, O1>,
                                                                       expressions::Expression<IntegerContext, A2, O2> >,
                                             expressions::ModOp>(std::make_pair(a, b)); }
            template<class A1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
            inline expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A1, O1>,
                                                                     expressions::Expression<IntegerContext, A2, O2> >,
                                           expressions::AndOp> operator & (const expressions::Expression<IntegerContext, A1, O1> & a,
                                                                           const expressions::Expression<IntegerContext, A2, O2> & b)
            { return expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A1, O1>,
                                                                       expressions::Expression<IntegerContext, A2, O2> >,
                                             expressions::AndOp>(std::make_pair(a, b)); }
            template<class A1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
            inline expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A1, O1>,
                                                                     expressions::Expression<IntegerContext, A2, O2> >,
                                           expressions::OrOp> operator | (const expressions::Expression<IntegerContext, A1, O1> & a,
                                                                          const expressions::Expression<IntegerContext, A2, O2> & b)
            { return expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A1, O1>,
                                                                       expressions::Expression<IntegerContext, A2, O2> >,
                                             expressions::OrOp>(std::make_pair(a, b)); }
            template<class A1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
            inline expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A1, O1>,
                                                                     expressions::Expression<IntegerContext, A2, O2> >,
                                           expressions::XOROp> operator ^ (const expressions::Expression<IntegerContext, A1, O1> & a,
                                                                           const expressions::Expression<IntegerContext, A2, O2> & b)
            { return expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A1, O1>,
                                                                       expressions::Expression<IntegerContext, A2, O2> >,
                                             expressions::XOROp>(std::make_pair(a, b)); }
            template<class A1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
            inline expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A1, O1>,
                                                                     expressions::Expression<IntegerContext, A2, O2> >,
                                           expressions::ShLOp> operator << (const expressions::Expression<IntegerContext, A1, O1> & a,
                                                                            const expressions::Expression<IntegerContext, A2, O2> & b)
            { return expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A1, O1>,
                                                                       expressions::Expression<IntegerContext, A2, O2> >,
                                             expressions::ShLOp>(std::make_pair(a, b)); }
            template<class A1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
            inline expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A1, O1>,
                                                                     expressions::Expression<IntegerContext, A2, O2> >,
                                           expressions::ShROp> operator >> (const expressions::Expression<IntegerContext, A1, O1> & a,
                                                                            const expressions::Expression<IntegerContext, A2, O2> & b)
            { return expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A1, O1>,
                                                                       expressions::Expression<IntegerContext, A2, O2> >,
                                             expressions::ShROp>(std::make_pair(a, b)); }
            
            template<class A, template<typename, typename> class O>
            inline bool isZero(const expressions::Expression<IntegerContext, A, O> & i)
            {
                return isZero(i.evaluate());
            }
            
            template<class A, template<typename, typename> class O>
            inline bool isOne(const expressions::Expression<IntegerContext, A, O> & i)
            {
                return isOne(i.evaluate());
            }
            
            template<class A, template<typename, typename> class O>
            inline bool isPMOne(const expressions::Expression<IntegerContext, A, O> & i)
            {
                return isPMOne(i.evaluate());
            }
            
            template<class A, template<typename, typename> class O>
            inline bool isPMTwo(const expressions::Expression<IntegerContext, A, O> & i)
            {
                return isPMTwo(i.evaluate());
            }
            
            template<class A, template<typename, typename> class O>
            inline bool isPositive(const expressions::Expression<IntegerContext, A, O> & i)
            {
                return isPositive(i.evaluate());
            }
            
            template<class A, template<typename, typename> class O>
            inline bool isNonNegative(const expressions::Expression<IntegerContext, A, O> & i)
            {
                return isNonNegative(i.evaluate());
            }
            
            template<class A, template<typename, typename> class O>
            inline bool isNegative(const expressions::Expression<IntegerContext, A, O> & i)
            {
                return isNegative(i.evaluate());
            }
            
            template<class A, template<typename, typename> class O>
            inline bool isNonPositive(const expressions::Expression<IntegerContext, A, O> & i)
            {
                return isNonPositive(i.evaluate());
            }
            
            template<class A, template<typename, typename> class O>
            inline std::ostream & operator << (std::ostream & s, const expressions::Expression<IntegerContext, A, O> & i)
            {
                return s << i.evaluate();
            }
        }
        
        /**@{
           \name Operators.
        */
        
        /**
           \brief Increments the integer by one and returns the previous value.
           
           \param cur The integer to work on.
           \return The previous value.
        */
        inline Integer operator ++ (Integer & cur, int)
        {
            Integer r(cur);
            increment(cur, cur);
            return r;
        }
        
        /**
           \brief Decrements the integer by one and returns the previous value.
           
           \param cur The integer to work on.
           \return The previous value.
        */
        inline Integer operator -- (Integer & cur, int)
        {
            Integer r(cur);
            decrement(cur, cur);
            return r;
        }
        
        /**
           \brief Increments the integer by one and returns the new value.
           
           \param cur The integer to work on.
           \return The new value (a reference to `cur`).
        */
        inline Integer & operator ++ (Integer & cur)
        {
            increment(cur, cur);
            return cur;
        }
        
        /**
           \brief Decreases the integer by one and returns the new value.
           
           \param cur The integer to work on.
           \return The new value (a reference to `cur`).
        */
        inline Integer & operator -- (Integer & cur)
        {
            decrement(cur, cur);
            return cur;
        }
        
        ///@}
        
        /**@{
           \name Assignment-operation operators.
        */
        
        /**
           \brief Subtracts the multiplication expression `E` from `cur`.
           
           Note that the faster `submul()` command is used for this special case where the product
           of two integers is subtracted from a third one.
           
           \param cur The integer to operate on.
           \param E The multiplication expression to subtract from `cur`.
           \return A reference to `cur` containing the result.
        */
        template<class D>
        inline Integer & operator -= (Integer & cur, const expressions::Expression<IntegerContext, D, expressions::MulOp> & E)
        {
            submul(cur, E.data().first.evaluate(), E.data().second.evaluate());
            return cur;
        }
        
        /**
           \brief Subtracts the integer `i` from `cur`.
           
           \param cur The integer to operate on.
           \param i The integer to subtract from `cur`.
           \return A reference to `cur` containing the result.
        */
        inline Integer & operator -= (Integer & cur, const Integer & i)
        {
            sub(cur, cur, i);
            return cur;
        }
        
        /**
           \brief Adds the given multiplication expression to `cur`.
           
           Note that the faster `addmul()` command is used for this special case where the product
           of two integers is added to a third one.
           
           \param cur The integer to operate on.
           \param E The multiplication expression to add to `cur`.
           \return A reference to `cur` containing the result.
        */
        template<class D>
        inline Integer & operator += (Integer & cur, const expressions::Expression<IntegerContext, D, expressions::MulOp> & E)
        {
            addmul(cur, E.data().first.evaluate(), E.data().second.evaluate());
            return cur;
        }
        
        /**
           \brief Adds the integer `i` to `cur`.
           
           \param cur The integer to operate on.
           \param i The integer to add to `cur`.
           \return A reference to `cur` containing the result.
        */
        inline Integer & operator += (Integer & cur, const Integer & i)
        {
            add(cur, cur, i);
            return cur;
        }
        
        /**
           \brief Multiplies the given integer `i` with `cur`.
               
           \param cur The integer to operate on.
           \param i The integer to multiply to `cur`.
           \return A reference to `cur` containing the result.
        */
        inline Integer & operator *= (Integer & cur, const Integer & i)
        {
            mul(cur, cur, i);
            return cur;
        }
        
        /**
           \brief Divides `cur` by the given integer `i`.
           
           \param cur The integer to operate on.
           \param i The integer to divide by.
           \return A reference to `cur` containing the result.
        */
        inline Integer & operator /= (Integer & cur, const Integer & i)
        {
            div(cur, cur, i);
            return cur;
        }
        
        /**
           \brief Divides `cur` by the given integer `i` and
                  stores the remainder in `cur`.
           
           \param cur The integer to operate on.
           \param i The integer to divide by.
           \return A reference to `cur` containing the result.
        */
        inline Integer & operator %= (Integer & cur, const Integer & i)
        {
            mod(cur, cur, i);
            return cur;
        }
        
        /**
           \brief Computes the bitwise or of `cur` with `i` and stores the result in `cur`.
           
           \param cur The integer to operate on.
           \param i The second operand.
           \return A reference to `cur` containing the result.
        */
        inline Integer & operator |= (Integer & cur, const Integer & i)
        {
            bor(cur, cur, i);
            return cur;
        }
        
        /**
           \brief Computes the bitwise and of `cur` with `i` and stores the result in `cur`.
           
           \param cur The integer to operate on.
           \param i The second operand.
           \return A reference to `cur` containing the result.
        */
        inline Integer & operator &= (Integer & cur, const Integer & i)
        {
            band(cur, cur, i);
            return cur;
        }
        
        /**
           \brief Computes the bitwise exclusive or of `cur` with `i` and stores the result in
                  `cur`.
           
           \param cur The integer to operate on.
           \param i The second operand.
           \return A reference to `cur` containing the result.
        */
        inline Integer & operator ^= (Integer & cur, const Integer & i)
        {
            bxor(cur, cur, i);
            return cur;
        }
        
        /**
           \brief Computes the bitwise left shift of `cur` by `i` bits, and stores the result in
                  `cur`.
           
           This effectively multiplies `cur` by 2 to the power of `i`.
           
           \param cur The integer to operate on.
           \param i The second operand.
           \return A reference to `cur` containing the result.
        */
        inline Integer & operator <<= (Integer & cur, const Integer & i)
        {
            shl(cur, cur, i);
            return cur;
        }
        
        /**
           \brief Computes the bitwise right shift of `cur` by `i` bits, and stores the result in
                  `cur`.
           
           This effectively divides `cur` by 2 to the power of `i`, and rounds towards zero.
           
           \param cur The integer to operate on.
           \param i The second operand.
           \return A reference to `cur` containing the result.
        */
        inline Integer & operator >>= (Integer & cur, const Integer & i)
        {
            shr(cur, cur, i);
            return cur;
        }
        
        /**
           \brief Computes the bitwise left shift of `cur` by `i` bits, and stores the result in
                  `cur`.
           
           This effectively multiplies `cur` by 2 to the power of `i`.
           
           \param cur The integer to operate on.
           \param i The second operand.
           \return A reference to `cur` containing the result.
        */
        inline Integer & operator <<= (Integer & cur, long i)
        {
            shl(cur, cur, i);
            return cur;
        }
        
        /**
           \brief Computes the bitwise right shift of `cur` by `i` bits, and stores the result in
                  `cur`.
           
           This effectively divides `cur` by 2 to the power of `i`, and rounds towards zero.
           
           \param cur The integer to operate on.
           \param i The second operand.
           \return A reference to `cur` containing the result.
        */
        inline Integer & operator >>= (Integer & cur, long i)
        {
            shr(cur, cur, i);
            return cur;
        }
        
        ///@}
        
        inline Integer & operator -= (Integer & cur, const expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> & E)
        {
            sub_si(cur, cur, E.data().evaluate());
            return cur;
        }
        
        inline Integer & operator -= (Integer & cur, const expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> & E)
        {
            sub_ui(cur, cur, E.data().evaluate());
            return cur;
        }
        
        template<class D>
        inline Integer & operator -= (Integer & cur, const expressions::Expression<IntegerContext, std::pair<D, expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> >, expressions::MulOp> & E)
        {
            submul_ui(cur, E.data().first.evaluate(), E.data().second.data().evaluate());
            return cur;
        }
        
        template<class D>
        inline Integer & operator -= (Integer & cur, const expressions::Expression<IntegerContext, std::pair<D, expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> >, expressions::MulOp> & E)
        {
            submul_si(cur, E.data().first.evaluate(), E.data().second.data().evaluate());
            return cur;
        }
        
        template<class D>
        inline Integer & operator -= (Integer & cur, const expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context>, D>, expressions::MulOp> & E)
        {
            submul_ui(cur, E.data().second.evaluate(), E.data().first.data().evaluate());
            return cur;
        }
        
        template<class D>
        inline Integer & operator -= (Integer & cur, const expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context>, D>, expressions::MulOp> & E)
        {
            mpz_submul_ui(cur, E.data().second.evaluate(), E.data().first.data().evaluate());
            return cur;
        }
        
        inline Integer & operator += (Integer & cur, const expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> & E)
        {
            add_si(cur, cur, E.data().evaluate());
            return cur;
        }
        
        inline Integer & operator += (Integer & cur, const expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> & E)
        {
            add_ui(cur, cur, E.data().evaluate());
            return cur;
        }
        
        template<class D>
        inline Integer & operator += (Integer & cur, const expressions::Expression<IntegerContext, std::pair<D, expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> >, expressions::MulOp> & E)
        {
            addmul_ui(cur, E.data().first.evaluate(), E.data().second.data().evaluate());
            return cur;
        }
        
        template<class D>
        inline Integer & operator += (Integer & cur, const expressions::Expression<IntegerContext, std::pair<D, expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> >, expressions::MulOp> & E)
        {
            addmul_si(cur, E.data().first.evaluate(), E.data().second.data().evaluate());
            return cur;
        }
        
        template<class D>
        inline Integer & operator += (Integer & cur, const expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context>, D>, expressions::MulOp> & E)
        {
            addmul_ui(cur, E.data().second.evaluate(), E.data().first.data().evaluate());
            return cur;
        }
        
        template<class D>
        inline Integer & operator += (Integer & cur, const expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context>, D>, expressions::MulOp> & E)
        {
            addmul_si(cur, E.data().second.evaluate(), E.data().first.data().evaluate());
            return cur;
        }
        
        inline Integer & operator *= (Integer & cur, const expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> & E)
        {
            mpz_mul_si(cur.getInternal(), cur.getInternal(), E.data().evaluate());
            return cur;
        }
        
        inline Integer & operator *= (Integer & cur, const expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> & E)
        {
            mpz_mul_ui(cur.getInternal(), cur.getInternal(), E.data().evaluate());
            return cur;
        }
        
        inline Integer & operator /= (Integer & cur, const expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> & E)
        {
            if (E.data().evaluate() < 0)
            {
                mpz_tdiv_q_ui(cur.getInternal(), cur.getInternal(), -E.data().evaluate());
                neg(cur, cur);
            }
            else
                mpz_tdiv_q_ui(cur.getInternal(), cur.getInternal(), E.data().evaluate());
            return cur;
        }
        
        inline Integer & operator /= (Integer & cur, const expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> & E)
        {
            mpz_tdiv_q_ui(cur.getInternal(), cur.getInternal(), E.data().evaluate());
            return cur;
        }
        
        inline Integer & operator %= (Integer & cur, const expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> & E)
        {
            if (E.data().evaluate() < 0)
                mpz_tdiv_r_ui(cur.getInternal(), cur.getInternal(), -E.data().evaluate());
            else
                mpz_tdiv_r_ui(cur.getInternal(), cur.getInternal(), E.data().evaluate());
            return cur;
        }
        
        inline Integer & operator %= (Integer & cur, const expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> & E)
        {
            mpz_tdiv_r_ui(cur.getInternal(), cur.getInternal(), E.data().evaluate());
            return cur;
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
        inline bool operator == (const Integer & a, const Integer & b)
        {
            return mpz_cmp(a.d_value, b.d_value) == 0;
        }
        
        /**
           \brief Compares the current integer with the given one for inequality.
           \param a The first operand.
           \param b The second operand.
           \return `true` if `a` does not equal `b`.
        */
        inline bool operator != (const Integer & a, const Integer & b)
        {
            return mpz_cmp(a.d_value, b.d_value) != 0;
        }
        
        /**
           \brief Compares the current integer with the given one.
           \param a The first operand.
           \param b The second operand.
           \return `true` if `a` is less than or equal to `b`.
        */
        inline bool operator <= (const Integer & a, const Integer & b)
        {
            return mpz_cmp(a.d_value, b.d_value) <= 0;
        }
        
        /**
           \brief Compares the current integer with the given one.
           \param a The first operand.
           \param b The second operand.
           \return `true` if `a` is greater than or equal to `b`.
        */
        inline bool operator >= (const Integer & a, const Integer & b)
        {
            return mpz_cmp(a.d_value, b.d_value) >= 0;
        }
        
        /**
           \brief Compares the current integer with the given one.
           \param a The first operand.
           \param b The second operand.
           \return `true` if `a` is less than `b`.
        */
        inline bool operator < (const Integer & a, const Integer & b)
        {
            return mpz_cmp(a.d_value, b.d_value) < 0;
        }
        
        /**
           \brief Compares the current integer with the given one.
           \param a The first operand.
           \param b The second operand.
           \return `true` if `a` is greater than `b`.
        */
        inline bool operator > (const Integer & a, const Integer & b)
        {
            return mpz_cmp(a.d_value, b.d_value) > 0;
        }
        
        ///@}
        
        // Second operand is expression: make both operands expressions
        template<class A2, template<typename, typename> class O2>
        inline bool operator == (const Integer & a, const expressions::Expression<IntegerContext, A2, O2> & b)
        {
            return expressions::make_expression<IntegerContext>(a) == b;
        }
        
        template<class A2, template<typename, typename> class O2>
        inline bool operator != (const Integer & a, const expressions::Expression<IntegerContext, A2, O2> & b)
        {
            return expressions::make_expression<IntegerContext>(a) != b;
        }
        
        template<class A2, template<typename, typename> class O2>
        inline bool operator <= (const Integer & a, const expressions::Expression<IntegerContext, A2, O2> & b)
        {
            return expressions::make_expression<IntegerContext>(a) <= b;
        }
        
        template<class A2, template<typename, typename> class O2>
        inline bool operator >= (const Integer & a, const expressions::Expression<IntegerContext, A2, O2> & b)
        {
            return expressions::make_expression<IntegerContext>(a) >= b;
        }
        
        template<class A2, template<typename, typename> class O2>
        inline bool operator < (const Integer & a, const expressions::Expression<IntegerContext, A2, O2> & b)
        {
            return expressions::make_expression<IntegerContext>(a) < b;
        }
        
        template<class A2, template<typename, typename> class O2>
        inline bool operator > (const Integer & a, const expressions::Expression<IntegerContext, A2, O2> & b)
        {
            return expressions::make_expression<IntegerContext>(a) > b;
        }
        
        namespace implementation
        {
            namespace Integer_impl
            {
                // First operand is expression: make both operands expressions
                template<class A1, template<typename, typename> class O1>
                inline bool operator == (const expressions::Expression<IntegerContext, A1, O1> & a, const arithmetic::Integer & b)
                {
                    return a == expressions::make_expression<IntegerContext>(b);
                }
                
                template<class A1, template<typename, typename> class O1>
                inline bool operator != (const expressions::Expression<IntegerContext, A1, O1> & a, const arithmetic::Integer & b)
                {
                    return a != expressions::make_expression<IntegerContext>(b);
                }
                
                template<class A1, template<typename, typename> class O1>
                inline bool operator <= (const expressions::Expression<IntegerContext, A1, O1> & a, const arithmetic::Integer & b)
                {
                    return a <= expressions::make_expression<IntegerContext>(b);
                }
                
                template<class A1, template<typename, typename> class O1>
                inline bool operator >= (const expressions::Expression<IntegerContext, A1, O1> & a, const arithmetic::Integer & b)
                {
                    return a >= expressions::make_expression<IntegerContext>(b);
                }
                
                template<class A1, template<typename, typename> class O1>
                inline bool operator < (const expressions::Expression<IntegerContext, A1, O1> & a, const arithmetic::Integer & b)
                {
                    return a < expressions::make_expression<IntegerContext>(b);
                }
                
                template<class A1, template<typename, typename> class O1>
                inline bool operator > (const expressions::Expression<IntegerContext, A1, O1> & a, const arithmetic::Integer & b)
                {
                    return a > expressions::make_expression<IntegerContext>(b);
                }
                
                // Both operands are expressions
                template<class A1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
                inline bool operator == (const expressions::Expression<IntegerContext, A1, O1> & a,
                                         const expressions::Expression<IntegerContext, A2, O2> & b)
                {
                    return a.evaluate() == b.evaluate();
                }
                
                template<class A1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
                inline bool operator != (const expressions::Expression<IntegerContext, A1, O1> & a,
                                         const expressions::Expression<IntegerContext, A2, O2> & b)
                {
                    return a.evaluate() != b.evaluate();
                }
                
                template<class A1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
                inline bool operator <= (const expressions::Expression<IntegerContext, A1, O1> & a,
                                         const expressions::Expression<IntegerContext, A2, O2> & b)
                {
                    return a.evaluate() <= b.evaluate();
                }
                
                template<class A1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
                inline bool operator >= (const expressions::Expression<IntegerContext, A1, O1> & a,
                                         const expressions::Expression<IntegerContext, A2, O2> & b)
                {
                    return a.evaluate() >= b.evaluate();
                }
                
                template<class A1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
                inline bool operator < (const expressions::Expression<IntegerContext, A1, O1> & a,
                                        const expressions::Expression<IntegerContext, A2, O2> & b)
                {
                    return a.evaluate() < b.evaluate();
                }
                
                template<class A1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
                inline bool operator > (const expressions::Expression<IntegerContext, A1, O1> & a,
                                        const expressions::Expression<IntegerContext, A2, O2> & b)
                {
                    return a.evaluate() > b.evaluate();
                }
                
                // First operand is integer conversion expression, second operand arbitrary expression
                template<class A2, template<typename, typename> class O2>
                inline bool operator == (const expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> & a,
                                         const expressions::Expression<IntegerContext, A2, O2> & b)
                { return compare_ui(b.evaluate(), a.data().evaluate()) == 0; }
                template<class A2, template<typename, typename> class O2>
                inline bool operator != (const expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> & a,
                                         const expressions::Expression<IntegerContext, A2, O2> & b)
                { return compare_ui(b.evaluate(), a.data().evaluate()) != 0; }
                template<class A2, template<typename, typename> class O2>
                inline bool operator <= (const expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> & a,
                                         const expressions::Expression<IntegerContext, A2, O2> & b)
                { return compare_ui(b.evaluate(), a.data().evaluate()) >= 0; }
                template<class A2, template<typename, typename> class O2>
                inline bool operator >= (const expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> & a,
                                         const expressions::Expression<IntegerContext, A2, O2> & b)
                { return compare_ui(b.evaluate(), a.data().evaluate()) <= 0; }
                template<class A2, template<typename, typename> class O2>
                inline bool operator < (const expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> & a,
                                        const expressions::Expression<IntegerContext, A2, O2> & b)
                { return compare_ui(b.evaluate(), a.data().evaluate()) > 0; }
                template<class A2, template<typename, typename> class O2>
                inline bool operator > (const expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> & a,
                                        const expressions::Expression<IntegerContext, A2, O2> & b)
                { return compare_ui(b.evaluate(), a.data().evaluate()) < 0; }
                
                template<class A2, template<typename, typename> class O2>
                inline bool operator == (const expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> & a,
                                         const expressions::Expression<IntegerContext, A2, O2> & b)
                { return compare_si(b.evaluate(), a.data().evaluate()) == 0; }
                template<class A2, template<typename, typename> class O2>
                inline bool operator != (const expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> & a,
                                         const expressions::Expression<IntegerContext, A2, O2> & b)
                { return compare_si(b.evaluate(), a.data().evaluate()) != 0; }
                template<class A2, template<typename, typename> class O2>
                inline bool operator <= (const expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> & a,
                                         const expressions::Expression<IntegerContext, A2, O2> & b)
                { return compare_si(b.evaluate(), a.data().evaluate()) >= 0; }
                template<class A2, template<typename, typename> class O2>
                inline bool operator >= (const expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> & a,
                                         const expressions::Expression<IntegerContext, A2, O2> & b)
                { return compare_si(b.evaluate(), a.data().evaluate()) <= 0; }
                template<class A2, template<typename, typename> class O2>
                inline bool operator < (const expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> & a,
                                        const expressions::Expression<IntegerContext, A2, O2> & b)
                { return compare_si(b.evaluate(), a.data().evaluate()) > 0; }
                template<class A2, template<typename, typename> class O2>
                inline bool operator > (const expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> & a,
                                        const expressions::Expression<IntegerContext, A2, O2> & b)
                { return compare_si(b.evaluate(), a.data().evaluate()) < 0; }
                
                template<class A2, template<typename, typename> class O2>
                inline bool operator == (const expressions::Expression<IntegerContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> & a,
                                         const expressions::Expression<IntegerContext, A2, O2> & b)
                { return compare_d(b.evaluate(), a.data().evaluate()) == 0; }
                template<class A2, template<typename, typename> class O2>
                inline bool operator != (const expressions::Expression<IntegerContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> & a,
                                         const expressions::Expression<IntegerContext, A2, O2> & b)
                { return compare_d(b.evaluate(), a.data().evaluate()) != 0; }
                template<class A2, template<typename, typename> class O2>
                inline bool operator <= (const expressions::Expression<IntegerContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> & a,
                                         const expressions::Expression<IntegerContext, A2, O2> & b)
                { return compare_d(b.evaluate(), a.data().evaluate()) >= 0; }
                template<class A2, template<typename, typename> class O2>
                inline bool operator >= (const expressions::Expression<IntegerContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> & a,
                                         const expressions::Expression<IntegerContext, A2, O2> & b)
                { return compare_d(b.evaluate(), a.data().evaluate()) <= 0; }
                template<class A2, template<typename, typename> class O2>
                inline bool operator < (const expressions::Expression<IntegerContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> & a,
                                        const expressions::Expression<IntegerContext, A2, O2> & b)
                { return compare_d(b.evaluate(), a.data().evaluate()) > 0; }
                template<class A2, template<typename, typename> class O2>
                inline bool operator > (const expressions::Expression<IntegerContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> & a,
                                        const expressions::Expression<IntegerContext, A2, O2> & b)
                { return compare_d(b.evaluate(), a.data().evaluate()) < 0; }
                
                // First operand is arbitrary expression, Second operand is integer conversion expression
                template<class A1, template<typename, typename> class O1>
                inline bool operator == (const expressions::Expression<IntegerContext, A1, O1> & a,
                                         const expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> & b)
                { return compare_ui(a.evaluate(), b.data().evaluate()) == 0; }
                template<class A1, template<typename, typename> class O1>
                inline bool operator != (const expressions::Expression<IntegerContext, A1, O1> & a,
                                         const expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> & b)
                { return compare_ui(a.evaluate(), b.data().evaluate()) != 0; }
                template<class A1, template<typename, typename> class O1>
                inline bool operator <= (const expressions::Expression<IntegerContext, A1, O1> & a,
                                         const expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> & b)
                { return compare_ui(a.evaluate(), b.data().evaluate()) <= 0; }
                template<class A1, template<typename, typename> class O1>
                inline bool operator >= (const expressions::Expression<IntegerContext, A1, O1> & a,
                                         const expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> & b)
                { return compare_ui(a.evaluate(), b.data().evaluate()) >= 0; }
                template<class A1, template<typename, typename> class O1>
                inline bool operator < (const expressions::Expression<IntegerContext, A1, O1> & a,
                                        const expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> & b)
                { return compare_ui(a.evaluate(), b.data().evaluate()) < 0; }
                template<class A1, template<typename, typename> class O1>
                inline bool operator > (const expressions::Expression<IntegerContext, A1, O1> & a,
                                        const expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> & b)
                { return compare_ui(a.evaluate(), b.data().evaluate()) > 0; }
                
                template<class A1, template<typename, typename> class O1>
                inline bool operator == (const expressions::Expression<IntegerContext, A1, O1> & a,
                                         const expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> & b)
                { return compare_si(a.evaluate(), b.data().evaluate()) == 0; }
                template<class A1, template<typename, typename> class O1>
                inline bool operator != (const expressions::Expression<IntegerContext, A1, O1> & a,
                                         const expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> & b)
                { return compare_si(a.evaluate(), b.data().evaluate()) != 0; }
                template<class A1, template<typename, typename> class O1>
                inline bool operator <= (const expressions::Expression<IntegerContext, A1, O1> & a,
                                         const expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> & b)
                { return compare_si(a.evaluate(), b.data().evaluate()) <= 0; }
                template<class A1, template<typename, typename> class O1>
                inline bool operator >= (const expressions::Expression<IntegerContext, A1, O1> & a,
                                         const expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> & b)
                { return compare_si(a.evaluate(), b.data().evaluate()) >= 0; }
                template<class A1, template<typename, typename> class O1>
                inline bool operator < (const expressions::Expression<IntegerContext, A1, O1> & a,
                                        const expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> & b)
                { return compare_si(a.evaluate(), b.data().evaluate()) < 0; }
                template<class A1, template<typename, typename> class O1>
                inline bool operator > (const expressions::Expression<IntegerContext, A1, O1> & a,
                                        const expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> & b)
                { return compare_si(a.evaluate(), b.data().evaluate()) > 0; }
                
                template<class A1, template<typename, typename> class O1>
                inline bool operator == (const expressions::Expression<IntegerContext, A1, O1> & a,
                                         const expressions::Expression<IntegerContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> & b)
                { return compare_d(a.evaluate(), b.data().evaluate()) == 0; }
                template<class A1, template<typename, typename> class O1>
                inline bool operator != (const expressions::Expression<IntegerContext, A1, O1> & a,
                                         const expressions::Expression<IntegerContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> & b)
                { return compare_d(a.evaluate(), b.data().evaluate()) != 0; }
                template<class A1, template<typename, typename> class O1>
                inline bool operator <= (const expressions::Expression<IntegerContext, A1, O1> & a,
                                         const expressions::Expression<IntegerContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> & b)
                { return compare_d(a.evaluate(), b.data().evaluate()) <= 0; }
                template<class A1, template<typename, typename> class O1>
                inline bool operator >= (const expressions::Expression<IntegerContext, A1, O1> & a,
                                         const expressions::Expression<IntegerContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> & b)
                { return compare_d(a.evaluate(), b.data().evaluate()) >= 0; }
                template<class A1, template<typename, typename> class O1>
                inline bool operator < (const expressions::Expression<IntegerContext, A1, O1> & a,
                                        const expressions::Expression<IntegerContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> & b)
                { return compare_d(a.evaluate(), b.data().evaluate()) < 0; }
                template<class A1, template<typename, typename> class O1>
                inline bool operator > (const expressions::Expression<IntegerContext, A1, O1> & a,
                                        const expressions::Expression<IntegerContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> & b)
                { return compare_d(a.evaluate(), b.data().evaluate()) > 0; }
            }
        }
        
        /**@{
           \name Operators.
        */
        /**
           \brief Computes and returns the absolute value of `i`.
           \param i The operand.
           \return The result.
        */
        inline expressions::Expression<IntegerContext, expressions::Wrapper<IntegerContext>,
                                                       expressions::AbsOp> abs(const Integer & i)
        { return expressions::Expression<IntegerContext, expressions::Wrapper<IntegerContext>,
                                                         expressions::AbsOp>(expressions::Wrapper<IntegerContext>(i)); }
        /**
           \brief Computes and returns the square of `i`.
           \param i The operand.
           \return The result.
        */
        inline expressions::Expression<IntegerContext, expressions::Wrapper<IntegerContext>,
                                                       expressions::SquareOp> square(const Integer & i)
        { return expressions::Expression<IntegerContext, expressions::Wrapper<IntegerContext>,
                                                         expressions::SquareOp>(expressions::Wrapper<IntegerContext>(i)); }
        ///@}
        /**@{
           \name Integer approximation.
        */
        /**
           \brief Computes and returns the ceil of the square root of `i`.
           \param i The operand.
           \return The result.
        */
        inline expressions::Expression<IntegerContext, expressions::Wrapper<IntegerContext>,
                                                       expressions::SqrtCeilOp> sqrtCeil(const Integer & i)
        { return expressions::Expression<IntegerContext, expressions::Wrapper<IntegerContext>,
                                                         expressions::SqrtCeilOp>(expressions::Wrapper<IntegerContext>(i)); }
        /**
           \brief Computes and returns the floor of the square root of `i`.
           \param i The operand.
           \return The result.
        */
        inline expressions::Expression<IntegerContext, expressions::Wrapper<IntegerContext>,
                                                       expressions::SqrtFloorOp> sqrtFloor(const Integer & i)
        { return expressions::Expression<IntegerContext, expressions::Wrapper<IntegerContext>,
                                                         expressions::SqrtFloorOp>(expressions::Wrapper<IntegerContext>(i)); }
        ///@}
        /**@{
           \name Euclidean ring functions.
        */
        /**
           \brief Computes and returns the non-negative Greatest Common Divisor of `a` and `b`.
           \param a The first operand.
           \param b The second operand.
           \return The result.
        */
        inline expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                 expressions::Wrapper<IntegerContext> >,
                                                       expressions::GCDOp> GCD(const Integer & a, const Integer & b)
        { return expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                   expressions::Wrapper<IntegerContext> >,
                                                         expressions::GCDOp>(std::make_pair(expressions::Wrapper<IntegerContext>(a),
                                                                                            expressions::Wrapper<IntegerContext>(b))); }
        /**
           \brief Computes and returns the non-negative Least Common Multiple of `a` and `b`.
           \param a The first operand.
           \param b The second operand.
           \return The result.
        */
        inline expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                 expressions::Wrapper<IntegerContext> >,
                                                       expressions::LCMOp> LCM(const Integer & a, const Integer & b)
        { return expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                   expressions::Wrapper<IntegerContext> >,
                                                         expressions::LCMOp>(std::make_pair(expressions::Wrapper<IntegerContext>(a),
                                                                                            expressions::Wrapper<IntegerContext>(b))); }
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
        inline expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                 expressions::Wrapper<IntegerContext> >,
                                                       expressions::PowerOp> power(const Integer & a, const Integer & b)
        { return expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                   expressions::Wrapper<IntegerContext> >,
                                                         expressions::PowerOp>(std::make_pair(expressions::Wrapper<IntegerContext>(a),
                                                                                              expressions::Wrapper<IntegerContext>(b))); }
        ///@}
        /**@{
           \name Integer approximation.
        */
        /**
           \brief Computes and returns \f$\lfloor \tfrac{a}{b} \rfloor\f$.
           \param a The base.
           \param b The exponent.
           \return The result.
        */
        inline expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                 expressions::Wrapper<IntegerContext> >,
                                                       expressions::FloorDivOp> floorDiv(const Integer & a, const Integer & b)
        { return expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                   expressions::Wrapper<IntegerContext> >,
                                                         expressions::FloorDivOp>(std::make_pair(expressions::Wrapper<IntegerContext>(a),
                                                                                                 expressions::Wrapper<IntegerContext>(b))); }
        /**
           \brief Computes and returns \f$\lceil \tfrac{a}{b} \rceil\f$.
           \param a The first operand.
           \param b The second operand.
           \return The result.
        */
        inline expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                 expressions::Wrapper<IntegerContext> >,
                                                       expressions::CeilDivOp> ceilDiv(const Integer & a, const Integer & b)
        { return expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                   expressions::Wrapper<IntegerContext> >,
                                                         expressions::CeilDivOp>(std::make_pair(expressions::Wrapper<IntegerContext>(a),
                                                                                                expressions::Wrapper<IntegerContext>(b))); }
        /**
           \brief Computes and returns \f$\lfloor \tfrac{a}{b} \rceil\f$.
           \param a The first operand.
           \param b The second operand.
           \return The result.
        */
        inline expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                 expressions::Wrapper<IntegerContext> >,
                                                       expressions::RoundDivOp> roundDiv(const Integer & a, const Integer & b)
        { return expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                   expressions::Wrapper<IntegerContext> >,
                                                         expressions::RoundDivOp>(std::make_pair(expressions::Wrapper<IntegerContext>(a),
                                                                                                 expressions::Wrapper<IntegerContext>(b))); }
        ///@}
        /**@{
           \name Exponentiation.
        */
        /**
           \brief Computes and returns `a` raised to the power of `b`.
           \param a The base.
           \param b The exponent.
           \return The result.
        */
        inline expressions::Expression<IntegerContext, expressions::Wrapper<IntegerContext>,
                                                       expressions::PowerCOp<signed long>::impl> power(const Integer & a, signed long b)
        { return expressions::Expression<IntegerContext, expressions::Wrapper<IntegerContext>,
                                                         expressions::PowerCOp<signed long>::impl>(expressions::Wrapper<IntegerContext>(a), b); }
        ///@}
        
        template<class A2, template<typename, typename> class O2>
        inline expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                 expressions::Expression<IntegerContext, A2, O2> >,
                                                       expressions::GCDOp> GCD(const Integer & a, const expressions::Expression<IntegerContext, A2, O2> & b)
        { return expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                   expressions::Expression<IntegerContext, A2, O2> >,
                                                         expressions::GCDOp>(std::make_pair(expressions::Wrapper<IntegerContext>(a), b)); }
        template<class A2, template<typename, typename> class O2>
        inline expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                 expressions::Expression<IntegerContext, A2, O2> >,
                                                       expressions::LCMOp> LCM(const Integer & a, const expressions::Expression<IntegerContext, A2, O2> & b)
        { return expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                   expressions::Expression<IntegerContext, A2, O2> >,
                                                         expressions::LCMOp>(std::make_pair(expressions::Wrapper<IntegerContext>(a), b)); }
        
        template<class A2, template<typename, typename> class O2>
        inline expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                 expressions::Expression<IntegerContext, A2, O2> >,
                                                       expressions::PowerOp> power(const Integer & a, const expressions::Expression<IntegerContext, A2, O2> & b)
        { return expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                   expressions::Expression<IntegerContext, A2, O2> >,
                                                         expressions::PowerOp>(std::make_pair(expressions::Wrapper<IntegerContext>(a), b)); }
        template<class A2, template<typename, typename> class O2>
        inline expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                 expressions::Expression<IntegerContext, A2, O2> >,
                                                       expressions::FloorDivOp> floorDiv(const Integer & a, const expressions::Expression<IntegerContext, A2, O2> & b)
        { return expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                   expressions::Expression<IntegerContext, A2, O2> >,
                                                         expressions::FloorDivOp>(std::make_pair(expressions::Wrapper<IntegerContext>(a), b)); }
        template<class A2, template<typename, typename> class O2>
        inline expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                 expressions::Expression<IntegerContext, A2, O2> >,
                                                       expressions::CeilDivOp> ceilDiv(const Integer & a, const expressions::Expression<IntegerContext, A2, O2> & b)
        { return expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                   expressions::Expression<IntegerContext, A2, O2> >,
                                                         expressions::CeilDivOp>(std::make_pair(expressions::Wrapper<IntegerContext>(a), b)); }
        template<class A2, template<typename, typename> class O2>
        inline expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                 expressions::Expression<IntegerContext, A2, O2> >,
                                                       expressions::RoundDivOp> roundDiv(const Integer & a, const expressions::Expression<IntegerContext, A2, O2> & b)
        { return expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,
                                                                   expressions::Expression<IntegerContext, A2, O2> >,
                                                         expressions::RoundDivOp>(std::make_pair(expressions::Wrapper<IntegerContext>(a), b)); }
        
        namespace implementation
        {
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<IntegerContext, expressions::Expression<IntegerContext, A, O>,
                                           expressions::AbsOp> abs(const expressions::Expression<IntegerContext, A, O> & b)
            { return expressions::Expression<IntegerContext, expressions::Expression<IntegerContext, A, O>,
                                             expressions::AbsOp>(b); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<IntegerContext, expressions::Expression<IntegerContext, A, O>,
                                           expressions::SquareOp> square(const expressions::Expression<IntegerContext, A, O> & b)
            { return expressions::Expression<IntegerContext, expressions::Expression<IntegerContext, A, O>,
                                             expressions::SquareOp>(b); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<IntegerContext, expressions::Expression<IntegerContext, A, O>,
                                           expressions::SqrtCeilOp> sqrtCeil(const expressions::Expression<IntegerContext, A, O> & b)
            { return expressions::Expression<IntegerContext, expressions::Expression<IntegerContext, A, O>,
                                             expressions::SqrtCeilOp>(b); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<IntegerContext, expressions::Expression<IntegerContext, A, O>,
                                           expressions::SqrtFloorOp> sqrtFloor(const expressions::Expression<IntegerContext, A, O> & b)
            { return expressions::Expression<IntegerContext, expressions::Expression<IntegerContext, A, O>,
                                             expressions::SqrtFloorOp>(b); }
            template<class A, template<typename, typename> class O, class A2, template<typename, typename> class O2>
            inline expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A, O>,
                                                                     expressions::Expression<IntegerContext, A2, O2> >,
                                           expressions::GCDOp> GCD(const expressions::Expression<IntegerContext, A, O> & a,
                                                                   const expressions::Expression<IntegerContext, A2, O2> & b)
            { return expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A, O>,
                                                                       expressions::Expression<IntegerContext, A2, O2> >,
                                             expressions::GCDOp>(std::make_pair(a, b)); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A, O>,
                                                                     expressions::Wrapper<IntegerContext> >,
                                           expressions::GCDOp> GCD(const expressions::Expression<IntegerContext, A, O> & a, const arithmetic::Integer & b)
            { return expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A, O>,
                                                                       expressions::Wrapper<IntegerContext> >,
                                             expressions::GCDOp>(std::make_pair(a, expressions::Wrapper<IntegerContext>(b))); }
            template<class A, template<typename, typename> class O, class A2, template<typename, typename> class O2>
            inline expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A, O>,
                                                                     expressions::Expression<IntegerContext, A2, O2> >,
                                           expressions::LCMOp> LCM(const expressions::Expression<IntegerContext, A, O> & a,
                                                                   const expressions::Expression<IntegerContext, A2, O2> & b)
            { return expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A, O>,
                                                                       expressions::Expression<IntegerContext, A2, O2> >,
                                             expressions::LCMOp>(std::make_pair(a, b)); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A, O>,
                                                                     expressions::Wrapper<IntegerContext> >,
                                           expressions::LCMOp> LCM(const expressions::Expression<IntegerContext, A, O> & a, const arithmetic::Integer & b)
            { return expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A, O>,
                                                                       expressions::Wrapper<IntegerContext> >,
                                             expressions::LCMOp>(std::make_pair(a, expressions::Wrapper<IntegerContext>(b))); }
            template<class A, template<typename, typename> class O, class A2, template<typename, typename> class O2>
            inline expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A, O>,
                                                                     expressions::Expression<IntegerContext, A2, O2> >,
                                           expressions::PowerOp> power(const expressions::Expression<IntegerContext, A, O> & a,
                                                                        const expressions::Expression<IntegerContext, A2, O2> & b)
            { return expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A, O>,
                                                                       expressions::Expression<IntegerContext, A2, O2> >,
                                             expressions::PowerOp>(std::make_pair(a, b)); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A, O>,
                                                                     expressions::Wrapper<IntegerContext> >,
                                           expressions::PowerOp> power(const expressions::Expression<IntegerContext, A, O> & a, const arithmetic::Integer & b)
            { return expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A, O>,
                                                                       expressions::Wrapper<IntegerContext> >,
                                             expressions::PowerOp>(std::make_pair(a, expressions::Wrapper<IntegerContext>(b))); }
            template<class A, template<typename, typename> class O, class A2, template<typename, typename> class O2>
            inline expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A, O>,
                                                                     expressions::Expression<IntegerContext, A2, O2> >,
                                           expressions::FloorDivOp> floorDiv(const expressions::Expression<IntegerContext, A, O> & a,
                                                                             const expressions::Expression<IntegerContext, A2, O2> & b)
            { return expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A, O>,
                                                                       expressions::Expression<IntegerContext, A2, O2> >,
                                             expressions::FloorDivOp>(std::make_pair(a, b)); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A, O>,
                                                                     expressions::Wrapper<IntegerContext> >,
                                           expressions::FloorDivOp> floorDiv(const expressions::Expression<IntegerContext, A, O> & a, const arithmetic::Integer & b)
            { return expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A, O>,
                                                                       expressions::Wrapper<IntegerContext> >,
                                             expressions::FloorDivOp>(std::make_pair(a, expressions::Wrapper<IntegerContext>(b))); }
            template<class A, template<typename, typename> class O, class A2, template<typename, typename> class O2>
            inline expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A, O>,
                                                                     expressions::Expression<IntegerContext, A2, O2> >,
                                           expressions::CeilDivOp> ceilDiv(const expressions::Expression<IntegerContext, A, O> & a,
                                                                           const expressions::Expression<IntegerContext, A2, O2> & b)
            { return expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A, O>,
                                                                       expressions::Expression<IntegerContext, A2, O2> >,
                                             expressions::CeilDivOp>(std::make_pair(a, b)); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A, O>,
                                                                     expressions::Wrapper<IntegerContext> >,
                                           expressions::CeilDivOp> ceilDiv(const expressions::Expression<IntegerContext, A, O> & a, const arithmetic::Integer & b)
            { return expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A, O>,
                                                                       expressions::Wrapper<IntegerContext> >,
                                             expressions::CeilDivOp>(std::make_pair(a, expressions::Wrapper<IntegerContext>(b))); }
            template<class A, template<typename, typename> class O, class A2, template<typename, typename> class O2>
            inline expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A, O>,
                                                                     expressions::Expression<IntegerContext, A2, O2> >,
                                           expressions::RoundDivOp> roundDiv(const expressions::Expression<IntegerContext, A, O> & a,
                                                                             const expressions::Expression<IntegerContext, A2, O2> & b)
            { return expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A, O>,
                                                                       expressions::Expression<IntegerContext, A2, O2> >,
                                             expressions::RoundDivOp>(std::make_pair(a, b)); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A, O>,
                                                                     expressions::Wrapper<IntegerContext> >,
                                           expressions::RoundDivOp> roundDiv(const expressions::Expression<IntegerContext, A, O> & a, const arithmetic::Integer & b)
            { return expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, A, O>,
                                                                       expressions::Wrapper<IntegerContext> >,
                                             expressions::RoundDivOp>(std::make_pair(a, expressions::Wrapper<IntegerContext>(b))); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<IntegerContext, expressions::Expression<IntegerContext, A, O>,
                                           expressions::PowerCOp<signed long>::impl> power(const expressions::Expression<IntegerContext, A, O> & a, signed long b)
            { return expressions::Expression<IntegerContext, expressions::Expression<IntegerContext, A, O>,
                                             expressions::PowerCOp<signed long>::impl>(a, b); }
        }
    }
}

#endif

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

#ifndef PLLL_INCLUDE_GUARD__ARITHMETIC_GMP_REAL_OPS_HPP
#define PLLL_INCLUDE_GUARD__ARITHMETIC_GMP_REAL_OPS_HPP

/**
   \file
   \brief Operator definitions for floating point numbers.
   
   This header contains templates and instantiations to implement all operations on floating point
   numbers.
*/
namespace plll
{
    namespace arithmetic
    {
        namespace expressions
        {
            // Adding real and (un)signed long, double or Integer
            template<class D>
            void do_assign(arithmetic::Real & x,
                           const AddOp<RealContext, std::pair<D, expressions::Expression<RealContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> > > &,
                           const std::pair<D, expressions::Expression<RealContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> > & data)
            {
                add_ui(x, data.first.evaluate(), data.second.data().evaluate());
            }
            
            template<class D>
            void do_assign(arithmetic::Real & x,
                           const AddOp<RealContext, std::pair<D, expressions::Expression<RealContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> > > &,
                           const std::pair<D, expressions::Expression<RealContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> > & data)
            {
                add_si(x, data.first.evaluate(), data.second.data().evaluate());
            }
            
            template<class D>
            void do_assign(arithmetic::Real & x,
                           const AddOp<RealContext, std::pair<D, expressions::Expression<RealContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> > > &,
                           const std::pair<D, expressions::Expression<RealContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> > & data)
            {
                add_d(x, data.first.evaluate(), data.second.data().evaluate());
            }
            
            template<class D1, class D2, template<typename, typename> class O2, class E>
            void do_assign(arithmetic::Real & x,
                           const AddOp<RealContext, std::pair<D1, Expression<RealContext, Expression<IntegerContext, D2, O2>, ConvertOp_Context> > > &,
                           const std::pair<D1, Expression<RealContext, Expression<IntegerContext, D2, O2>, ConvertOp_Context> > & data)
            {
                add_z(x, data.first.evaluate(), data.second.data().evaluate());
            }
            
            template<class D>
            void do_assign(arithmetic::Real & x,
                           const AddOp<RealContext, std::pair<expressions::Expression<RealContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context>, D> > &,
                           const std::pair<expressions::Expression<RealContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context>, D> & data)
            {
                add_ui(x, data.second.evaluate(), data.first.data().evaluate());
            }
            
            template<class D>
            void do_assign(arithmetic::Real & x,
                           const AddOp<RealContext, std::pair<expressions::Expression<RealContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context>, D> > &,
                           const std::pair<expressions::Expression<RealContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context>, D> & data)
            {
                add_si(x, data.second.evaluate(), data.first.data().evaluate());
            }
            
            template<class D>
            void do_assign(arithmetic::Real & x,
                           const AddOp<RealContext, std::pair<expressions::Expression<RealContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context>, D> > &,
                           const std::pair<expressions::Expression<RealContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context>, D> & data)
            {
                add_d(x, data.second.evaluate(), data.first.data().evaluate());
            }
            
            template<class D1, class D2, template<typename, typename> class O2, class E>
            void do_assign(arithmetic::Real & x,
                           const AddOp<RealContext, std::pair<Expression<RealContext, Expression<IntegerContext, D2, O2>, ConvertOp_Context>, D1> > &,
                           const std::pair<Expression<RealContext, Expression<IntegerContext, D2, O2>, ConvertOp_Context>, D1> & data)
            {
                add_z(x, data.second.evaluate(), data.first.data().evaluate());
            }
            
            // Subtracting real and (un)signed long, double or Integer
            template<class D>
            void do_assign(arithmetic::Real & x,
                           const SubOp<RealContext, std::pair<D, expressions::Expression<RealContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> > > &,
                           const std::pair<D, expressions::Expression<RealContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> > & data)
            {
                sub_ui(x, data.first.evaluate(), data.second.data().evaluate());
            }
            
            template<class D>
            void do_assign(arithmetic::Real & x,
                           const SubOp<RealContext, std::pair<D, expressions::Expression<RealContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> > > &,
                           const std::pair<D, expressions::Expression<RealContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> > & data)
            {
                sub_si(x, data.first.evaluate(), data.second.data().evaluate());
            }
            
            template<class D>
            void do_assign(arithmetic::Real & x,
                           const SubOp<RealContext, std::pair<D, expressions::Expression<RealContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> > > &,
                           const std::pair<D, expressions::Expression<RealContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> > & data)
            {
                sub_d(x, data.first.evaluate(), data.second.data().evaluate());
            }
            
            template<class D1, class D2, template<typename, typename> class O2, class E>
            void do_assign(arithmetic::Real & x,
                           const SubOp<RealContext, std::pair<D1, Expression<RealContext, Expression<IntegerContext, D2, O2>, ConvertOp_Context> > > &,
                           const std::pair<D1, Expression<RealContext, Expression<IntegerContext, D2, O2>, ConvertOp_Context> > & data)
            {
                sub_z(x, data.first.evaluate(), data.second.data().evaluate());
            }
            
            template<class D>
            void do_assign(arithmetic::Real & x,
                           const SubOp<RealContext, std::pair<expressions::Expression<RealContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context>, D> > &,
                           const std::pair<expressions::Expression<RealContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context>, D> & data)
            {
                ui_sub(x, data.first.data().evaluate(), data.second.evaluate());
            }
            
            template<class D>
            void do_assign(arithmetic::Real & x,
                           const SubOp<RealContext, std::pair<expressions::Expression<RealContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context>, D> > &,
                           const std::pair<expressions::Expression<RealContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context>, D> & data)
            {
                si_sub(x, data.first.data().evaluate(), data.second.evaluate());
            }
            
            template<class D>
            void do_assign(arithmetic::Real & x,
                           const SubOp<RealContext, std::pair<expressions::Expression<RealContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context>, D> > &,
                           const std::pair<expressions::Expression<RealContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context>, D> & data)
            {
                d_sub(x, data.first.data().evaluate(), data.second.evaluate());
            }
            
            template<class D1, class D2, template<typename, typename> class O2, class E>
            void do_assign(arithmetic::Real & x,
                           const SubOp<RealContext, std::pair<Expression<RealContext, Expression<IntegerContext, D2, O2>, ConvertOp_Context>, D1> > &,
                           const std::pair<Expression<RealContext, Expression<IntegerContext, D2, O2>, ConvertOp_Context>, D1> & data)
            {
                z_sub(x, data.first.data().evaluate(), data.second.evaluate());
            }
            
            // Multiplyign real and (un)signed long, double or Integer
            template<class D>
            void do_assign(arithmetic::Real & x,
                           const MulOp<RealContext, std::pair<D, expressions::Expression<RealContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> > > &,
                           const std::pair<D, expressions::Expression<RealContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> > & data)
            {
                mul_ui(x, data.first.evaluate(), data.second.data().evaluate());
            }
            
            template<class D>
            void do_assign(arithmetic::Real & x,
                           const MulOp<RealContext, std::pair<D, expressions::Expression<RealContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> > > &,
                           const std::pair<D, expressions::Expression<RealContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> > & data)
            {
                mul_si(x, data.first.evaluate(), data.second.data().evaluate());
            }
            
            template<class D>
            void do_assign(arithmetic::Real & x,
                           const MulOp<RealContext, std::pair<D, expressions::Expression<RealContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> > > &,
                           const std::pair<D, expressions::Expression<RealContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> > & data)
            {
                mul_d(x, data.first.evaluate(), data.second.data().evaluate());
            }
            
            template<class D1, class D2, template<typename, typename> class O2, class E>
            void do_assign(arithmetic::Real & x,
                           const MulOp<RealContext, std::pair<D1, Expression<RealContext, Expression<IntegerContext, D2, O2>, ConvertOp_Context> > > &,
                           const std::pair<D1, Expression<RealContext, Expression<IntegerContext, D2, O2>, ConvertOp_Context> > & data)
            {
                mul_z(x, data.first.evaluate(), data.second.data().evaluate());
            }
            
            template<class D>
            void do_assign(arithmetic::Real & x,
                           const MulOp<RealContext, std::pair<expressions::Expression<RealContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context>, D> > &,
                           const std::pair<expressions::Expression<RealContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context>, D> & data)
            {
                mul_ui(x, data.second.evaluate(), data.first.data().evaluate());
            }
            
            template<class D>
            void do_assign(arithmetic::Real & x,
                           const MulOp<RealContext, std::pair<expressions::Expression<RealContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context>, D> > &,
                           const std::pair<expressions::Expression<RealContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context>, D> & data)
            {
                mul_si(x, data.second.evaluate(), data.first.data().evaluate());
            }
            
            template<class D>
            void do_assign(arithmetic::Real & x,
                           const MulOp<RealContext, std::pair<expressions::Expression<RealContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context>, D> > &,
                           const std::pair<expressions::Expression<RealContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context>, D> & data)
            {
                mul_d(x, data.second.evaluate(), data.first.data().evaluate());
            }
            
            template<class D1, class D2, template<typename, typename> class O2, class E>
            void do_assign(arithmetic::Real & x,
                           const MulOp<RealContext, std::pair<Expression<RealContext, Expression<IntegerContext, D2, O2>, ConvertOp_Context>, D1> > &,
                           const std::pair<Expression<RealContext, Expression<IntegerContext, D2, O2>, ConvertOp_Context>, D1> & data)
            {
                mul_z(x, data.second.evaluate(), data.first.data().evaluate());
            }
            
            // Dividing real and (un)signed long, double or Integer
            template<class D>
            void do_assign(arithmetic::Real & x,
                           const DivOp<RealContext, std::pair<D, expressions::Expression<RealContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> > > &,
                           const std::pair<D, expressions::Expression<RealContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> > & data)
            {
                div_ui(x, data.first.evaluate(), data.second.data().evaluate());
            }
            
            template<class D>
            void do_assign(arithmetic::Real & x,
                           const DivOp<RealContext, std::pair<D, expressions::Expression<RealContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> > > &,
                           const std::pair<D, expressions::Expression<RealContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> > & data)
            {
                div_si(x, data.first.evaluate(), data.second.data().evaluate());
            }
            
            template<class D>
            void do_assign(arithmetic::Real & x,
                           const DivOp<RealContext, std::pair<D, expressions::Expression<RealContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> > > &,
                           const std::pair<D, expressions::Expression<RealContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> > & data)
            {
                div_d(x, data.first.evaluate(), data.second.data().evaluate());
            }
            
            template<class D1, class D2, template<typename, typename> class O2, class E>
            void do_assign(arithmetic::Real & x,
                           const DivOp<RealContext, std::pair<D1, Expression<RealContext, Expression<IntegerContext, D2, O2>, ConvertOp_Context> > > &,
                           const std::pair<D1, Expression<RealContext, Expression<IntegerContext, D2, O2>, ConvertOp_Context> > & data)
            {
                div_z(x, data.first.evaluate(), data.second.data().evaluate());
            }
            
            template<class D>
            void do_assign(arithmetic::Real & x,
                           const DivOp<RealContext, std::pair<expressions::Expression<RealContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context>, D> > &,
                           const std::pair<expressions::Expression<RealContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context>, D> & data)
            {
                ui_div(x, data.first.data().evaluate(), data.second.evaluate());
            }
            
            template<class D>
            void do_assign(arithmetic::Real & x,
                           const DivOp<RealContext, std::pair<expressions::Expression<RealContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context>, D> > &,
                           const std::pair<expressions::Expression<RealContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context>, D> & data)
            {
                si_div(x, data.first.data().evaluate(), data.second.evaluate());
            }
            
            template<class D>
            void do_assign(arithmetic::Real & x,
                           const DivOp<RealContext, std::pair<expressions::Expression<RealContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context>, D> > &,
                           const std::pair<expressions::Expression<RealContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context>, D> & data)
            {
                d_div(x, data.first.data().evaluate(), data.second.evaluate());
            }
            
            /*
              Disabled, since these interfere with other definitions from above, for example in the
              case: a = convert() + c * d
              
              // Special add/subtract and multiplication mixes
              template<class D1, class D2, template<typename, typename> class O2, class E>
              void do_assign(arithmetic::Real & x,
                             const AddOp<RealContext, std::pair<D1, Expression<RealContext, D2, MulOp> > > &,
                             const std::pair<D1, Expression<RealContext, D2, MulOp> > & data)
              {
                  addmul4(x, data.second.data().first.evaluate(), data.second.data().second.evaluate(), data.first.evaluate());
              }
              
              template<class D1, class D2, template<typename, typename> class O2, class E>
              void do_assign(arithmetic::Real & x,
                             const AddOp<RealContext, std::pair<Expression<RealContext, D2, MulOp>, D1> > &,
                             const std::pair<Expression<RealContext, D2, MulOp>, D1> & data)
              {
                  addmul4(x, data.first.data().first.evaluate(), data.first.data().second.evaluate(), data.second.evaluate());
              }
              
              template<class D1, class D2, template<typename, typename> class O2, class E>
              void do_assign(arithmetic::Real & x,
                             const SubOp<RealContext, std::pair<D1, Expression<RealContext, D2, MulOp> > > &,
                             const std::pair<D1, Expression<RealContext, D2, MulOp> > & data)
              {
                  submul4(x, data.second.data().first.evaluate(), data.second.data().second.evaluate(), data.first.evaluate());
                  neg(x, x);
              }
              
              template<class D1, class D2, template<typename, typename> class O2, class E>
              void do_assign(arithmetic::Real & x,
                             const SubOp<RealContext, std::pair<Expression<RealContext, D2, MulOp>, D1> > &,
                             const std::pair<Expression<RealContext, D2, MulOp>, D1> & data)
              {
                  submul4(x, data.first.data().first.evaluate(), data.first.data().second.evaluate(), data.second.evaluate());
              }
            */
        }
        
        /**@{
           \name Operators.
        */
        
        /**
           \brief Negates the floating point number.
           
           \param a The floating point number.
           \return The result.
        */
        inline expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                    expressions::NegOp> operator - (const Real & a)
        { return expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                      expressions::NegOp>(expressions::Wrapper<RealContext>(a)); }
        /**
           \brief Adds the two floating point numbers and returns the result.
           
           \param a The first operand.
           \param b The second operand.
           \return The result.
        */
        inline expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                              expressions::Wrapper<RealContext> >,
                                                    expressions::AddOp> operator + (const Real & a, const Real & b)
        { return expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                                expressions::Wrapper<RealContext> >,
                                                      expressions::AddOp>(std::make_pair(expressions::Wrapper<RealContext>(a),
                                                                                         expressions::Wrapper<RealContext>(b))); }
        /**
           \brief Subtracts the second from the first floating point number and returns the result.
           
           \param a The first operand.
           \param b The second operand.
           \return The result.
        */
        inline expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                              expressions::Wrapper<RealContext> >,
                                                    expressions::SubOp> operator - (const Real & a, const Real & b)
        { return expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                                expressions::Wrapper<RealContext> >,
                                                      expressions::SubOp>(std::make_pair(expressions::Wrapper<RealContext>(a),
                                                                                         expressions::Wrapper<RealContext>(b))); }
        /**
           \brief Multiplies the two floating point numbers and returns the result.
           
           \param a The first operand.
           \param b The second operand.
           \return The result.
        */
        inline expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                              expressions::Wrapper<RealContext> >,
                                                    expressions::MulOp> operator * (const Real & a, const Real & b)
        { return expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                                expressions::Wrapper<RealContext> >,
                                                      expressions::MulOp>(std::make_pair(expressions::Wrapper<RealContext>(a),
                                                                                         expressions::Wrapper<RealContext>(b))); }
        /**
           \brief Divides the first by the second floating point number and returns the result.
           
           \param a The first operand.
           \param b The second operand.
           \return The result.
        */
        inline expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                              expressions::Wrapper<RealContext> >,
                                                    expressions::DivOp> operator / (const Real & a, const Real & b)
        { return expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                                expressions::Wrapper<RealContext> >,
                                                      expressions::DivOp>(std::make_pair(expressions::Wrapper<RealContext>(a),
                                                                                         expressions::Wrapper<RealContext>(b))); }
        /**
           \brief Divides the first by the second floating point number and returns the remainder.
           
           \param a The first operand.
           \param b The second operand.
           \return The result.
        */
        inline expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                              expressions::Wrapper<RealContext> >,
                                                    expressions::ModOp> operator % (const Real & a, const Real & b)
        { return expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                                expressions::Wrapper<RealContext> >,
                                                      expressions::ModOp>(std::make_pair(expressions::Wrapper<RealContext>(a),
                                                                                         expressions::Wrapper<RealContext>(b))); }
        /**
           \brief Multiplies `a` with 2 to the power of `b` and returns the result.
           
           \param a The first operand.
           \param b The second operand.
           \return The result.
           
           \sa \ref arithctxts_sqrtfp
        */
        inline expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                              expressions::Wrapper<RealContext> >,
                                                    expressions::ShLOp> operator << (const Real & a, const Real & b)
        { return expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                                expressions::Wrapper<RealContext> >,
                                                      expressions::ShLOp>(std::make_pair(expressions::Wrapper<RealContext>(a),
                                                                                         expressions::Wrapper<RealContext>(b))); }
        /**
           \brief Divides `a` by 2 to the power of `b` and returns the
                  result.
           
           \param a The first operand.
           \param b The second operand.
           \return The result.
           
           \sa \ref arithctxts_sqrtfp
        */
        inline expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                              expressions::Wrapper<RealContext> >,
                                                    expressions::ShROp> operator >> (const Real & a, const Real & b)
        { return expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                                expressions::Wrapper<RealContext> >,
                                                      expressions::ShROp>(std::make_pair(expressions::Wrapper<RealContext>(a),
                                                                                         expressions::Wrapper<RealContext>(b))); }
        /**
           \brief Multiplies `a` with 2 to the power of `b` and returns the
                  result.
                  
           \param a The first operand.
           \param b The second operand.
           \return The result.
        */
        inline expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                    expressions::ShiftCOp> operator << (const Real & a, signed long b)
        { return expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                      expressions::ShiftCOp>(expressions::Wrapper<RealContext>(a), b); }
        /**
           \brief Divides `a` by 2 to the power of `b` and returns the
                  result.
           
           \param a The first operand.
           \param b The second operand.
           \return The result.
        */
        inline expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                    expressions::ShiftCOp> operator >> (const Real & a, signed long b)
        { return expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                      expressions::ShiftCOp>(expressions::Wrapper<RealContext>(a), -b); }
        ///@}
        
        // second operand is a expression
        template<class A2, template<typename, typename> class O2>
        inline expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                              expressions::Expression<RealContext, A2, O2> >,
                                                    expressions::AddOp> operator + (const Real & a, const expressions::Expression<RealContext, A2, O2> & b)
        { return expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                                expressions::Expression<RealContext, A2, O2> >,
                                                      expressions::AddOp>(std::make_pair(expressions::Wrapper<RealContext>(a), b)); }
        template<class A2, template<typename, typename> class O2>
        inline expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                              expressions::Expression<RealContext, A2, O2> >,
                                                    expressions::SubOp> operator - (const Real & a, const expressions::Expression<RealContext, A2, O2> & b)
        { return expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                                expressions::Expression<RealContext, A2, O2> >,
                                                      expressions::SubOp>(std::make_pair(expressions::Wrapper<RealContext>(a), b)); }
        template<class A2, template<typename, typename> class O2>
        inline expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                              expressions::Expression<RealContext, A2, O2> >,
                                                    expressions::MulOp> operator * (const Real & a, const expressions::Expression<RealContext, A2, O2> & b)
        { return expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                                expressions::Expression<RealContext, A2, O2> >,
                                                      expressions::MulOp>(std::make_pair(expressions::Wrapper<RealContext>(a), b)); }
        template<class A2, template<typename, typename> class O2>
        inline expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                              expressions::Expression<RealContext, A2, O2> >,
                                                    expressions::DivOp> operator / (const Real & a, const expressions::Expression<RealContext, A2, O2> & b)
        { return expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                                expressions::Expression<RealContext, A2, O2> >,
                                                      expressions::DivOp>(std::make_pair(expressions::Wrapper<RealContext>(a), b)); }
        template<class A2, template<typename, typename> class O2>
        inline expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                              expressions::Expression<RealContext, A2, O2> >,
                                                    expressions::ModOp> operator % (const Real & a, const expressions::Expression<RealContext, A2, O2> & b)
        { return expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                                expressions::Expression<RealContext, A2, O2> >,
                                                      expressions::ModOp>(std::make_pair(expressions::Wrapper<RealContext>(a), b)); }
        template<class A2, template<typename, typename> class O2>
        inline expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                              expressions::Expression<RealContext, A2, O2> >,
                                                    expressions::ShLOp> operator << (const Real & a, const expressions::Expression<RealContext, A2, O2> & b)
        { return expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                                expressions::Expression<RealContext, A2, O2> >,
                                                      expressions::ShLOp>(std::make_pair(expressions::Wrapper<RealContext>(a), b)); }
        template<class A2, template<typename, typename> class O2>
        inline expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                              expressions::Expression<RealContext, A2, O2> >,
                                                    expressions::ShROp> operator >> (const Real & a, const expressions::Expression<RealContext, A2, O2> & b)
        { return expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                                expressions::Expression<RealContext, A2, O2> >,
                                                      expressions::ShROp>(std::make_pair(expressions::Wrapper<RealContext>(a), b)); }
        
        namespace expressions
        {
            // first operand is a expression
            template<class A1, template<typename, typename> class O1>
            inline expressions::Expression<RealContext, expressions::Expression<RealContext, A1, O1>,
                                           expressions::NegOp> operator - (const expressions::Expression<RealContext, A1, O1> & a)
            { return expressions::Expression<RealContext, expressions::Expression<RealContext, A1, O1>,
                                             expressions::NegOp>(a); }
            template<class A1, template<typename, typename> class O1>
            inline expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A1, O1>,
                                                                  expressions::Wrapper<RealContext> >,
                                           expressions::AddOp> operator + (const expressions::Expression<RealContext, A1, O1> & a,
                                                                           const arithmetic::Real & b)
            { return expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A1, O1>,
                                                                    expressions::Wrapper<RealContext> >,
                                             expressions::AddOp>(std::make_pair(a, expressions::Wrapper<RealContext>(b))); }
            template<class A1, template<typename, typename> class O1>
            inline expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A1, O1>,
                                                                  expressions::Wrapper<RealContext> >,
                                           expressions::SubOp> operator - (const expressions::Expression<RealContext, A1, O1> & a,
                                                                           const arithmetic::Real & b)
            { return expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A1, O1>,
                                                                    expressions::Wrapper<RealContext> >,
                                             expressions::SubOp>(std::make_pair(a, expressions::Wrapper<RealContext>(b))); }
            template<class A1, template<typename, typename> class O1>
            inline expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A1, O1>,
                                                                  expressions::Wrapper<RealContext> >,
                                           expressions::MulOp> operator * (const expressions::Expression<RealContext, A1, O1> & a,
                                                                           const arithmetic::Real & b)
            { return expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A1, O1>,
                                                                    expressions::Wrapper<RealContext> >,
                                             expressions::MulOp>(std::make_pair(a, expressions::Wrapper<RealContext>(b))); }
            template<class A1, template<typename, typename> class O1>
            inline expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A1, O1>,
                                                                  expressions::Wrapper<RealContext> >,
                                           expressions::DivOp> operator / (const expressions::Expression<RealContext, A1, O1> & a,
                                                                           const arithmetic::Real & b)
            { return expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A1, O1>,
                                                                    expressions::Wrapper<RealContext> >,
                                             expressions::DivOp>(std::make_pair(a, expressions::Wrapper<RealContext>(b))); }
            template<class A1, template<typename, typename> class O1>
            inline expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A1, O1>,
                                                                  expressions::Wrapper<RealContext> >,
                                           expressions::ModOp> operator % (const expressions::Expression<RealContext, A1, O1> & a,
                                                                           const arithmetic::Real & b)
            { return expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A1, O1>,
                                                                    expressions::Wrapper<RealContext> >,
                                             expressions::ModOp>(std::make_pair(a, expressions::Wrapper<RealContext>(b))); }
            template<class A1, template<typename, typename> class O1>
            inline expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A1, O1>,
                                                                  expressions::Wrapper<RealContext> >,
                                           expressions::ShLOp> operator << (const expressions::Expression<RealContext, A1, O1> & a,
                                                                            const arithmetic::Real & b)
            { return expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A1, O1>,
                                                                    expressions::Wrapper<RealContext> >,
                                             expressions::ShLOp>(std::make_pair(a, expressions::Wrapper<RealContext>(b))); }
            template<class A1, template<typename, typename> class O1>
            inline expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A1, O1>,
                                                                  expressions::Wrapper<RealContext> >,
                                           expressions::ShROp> operator >> (const expressions::Expression<RealContext, A1, O1> & a,
                                                                            const arithmetic::Real & b)
            { return expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A1, O1>,
                                                                    expressions::Wrapper<RealContext> >,
                                             expressions::ShROp>(std::make_pair(a, expressions::Wrapper<RealContext>(b))); }
            template<class A1, template<typename, typename> class O1>
            inline expressions::Expression<RealContext, expressions::Expression<RealContext, A1, O1>,
                                           expressions::ShiftCOp> operator << (const expressions::Expression<RealContext, A1, O1> & a,
                                                                               signed long b)
            { return expressions::Expression<RealContext, expressions::Expression<RealContext, A1, O1>,
                                             expressions::ShiftCOp>(a, b); }
            template<class A1, template<typename, typename> class O1>
            inline expressions::Expression<RealContext, expressions::Expression<RealContext, A1, O1>,
                                           expressions::ShiftCOp> operator >> (const expressions::Expression<RealContext, A1, O1> & a,
                                                                               signed long b)
            { return expressions::Expression<RealContext, expressions::Expression<RealContext, A1, O1>,
                                             expressions::ShiftCOp>(a, -b); }
            
            // both operands are expressions
            template<class A1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
            inline expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A1, O1>,
                                                                  expressions::Expression<RealContext, A2, O2> >,
                                           expressions::AddOp> operator + (const expressions::Expression<RealContext, A1, O1> & a,
                                                                           const expressions::Expression<RealContext, A2, O2> & b)
            { return expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A1, O1>,
                                                                    expressions::Expression<RealContext, A2, O2> >,
                                             expressions::AddOp>(std::make_pair(a, b)); }
            template<class A1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
            inline expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A1, O1>,
                                                                  expressions::Expression<RealContext, A2, O2> >,
                                           expressions::SubOp> operator - (const expressions::Expression<RealContext, A1, O1> & a,
                                                                           const expressions::Expression<RealContext, A2, O2> & b)
            { return expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A1, O1>,
                                                                    expressions::Expression<RealContext, A2, O2> >,
                                             expressions::SubOp>(std::make_pair(a, b)); }
            template<class A1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
            inline expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A1, O1>,
                                                                  expressions::Expression<RealContext, A2, O2> >,
                                           expressions::MulOp> operator * (const expressions::Expression<RealContext, A1, O1> & a,
                                                                           const expressions::Expression<RealContext, A2, O2> & b)
            { return expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A1, O1>,
                                                                    expressions::Expression<RealContext, A2, O2> >,
                                             expressions::MulOp>(std::make_pair(a, b)); }
            template<class A1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
            inline expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A1, O1>,
                                                                  expressions::Expression<RealContext, A2, O2> >,
                                           expressions::DivOp> operator / (const expressions::Expression<RealContext, A1, O1> & a,
                                                                           const expressions::Expression<RealContext, A2, O2> & b)
            { return expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A1, O1>,
                                                                    expressions::Expression<RealContext, A2, O2> >,
                                             expressions::DivOp>(std::make_pair(a, b)); }
            template<class A1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
            inline expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A1, O1>,
                                                                  expressions::Expression<RealContext, A2, O2> >,
                                           expressions::ModOp> operator % (const expressions::Expression<RealContext, A1, O1> & a,
                                                                           const expressions::Expression<RealContext, A2, O2> & b)
            { return expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A1, O1>,
                                                                    expressions::Expression<RealContext, A2, O2> >,
                                             expressions::ModOp>(std::make_pair(a, b)); }
            template<class A1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
            inline expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A1, O1>,
                                                                  expressions::Expression<RealContext, A2, O2> >,
                                           expressions::ShLOp> operator << (const expressions::Expression<RealContext, A1, O1> & a,
                                                                            const expressions::Expression<RealContext, A2, O2> & b)
            { return expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A1, O1>,
                                                                    expressions::Expression<RealContext, A2, O2> >,
                                             expressions::ShLOp>(std::make_pair(a, b)); }
            template<class A1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
            inline expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A1, O1>,
                                                                  expressions::Expression<RealContext, A2, O2> >,
                                           expressions::ShROp> operator >> (const expressions::Expression<RealContext, A1, O1> & a,
                                                                            const expressions::Expression<RealContext, A2, O2> & b)
            { return expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A1, O1>,
                                                                    expressions::Expression<RealContext, A2, O2> >,
                                             expressions::ShROp>(std::make_pair(a, b)); }
            
            template<class A, template<typename, typename> class O>
            inline bool isZero(const expressions::Expression<RealContext, A, O> & r)
            {
                return isZero(r.evaluate());
            }
            
            template<class A, template<typename, typename> class O>
            inline bool isOne(const expressions::Expression<RealContext, A, O> & r)
            {
                return isOne(r.evaluate());
            }
            
            template<class A, template<typename, typename> class O>
            inline bool isPositive(const expressions::Expression<RealContext, A, O> & r)
            {
                return isPositive(r.evaluate());
            }
            
            template<class A, template<typename, typename> class O>
            inline bool isNonNegative(const expressions::Expression<RealContext, A, O> & r)
            {
                return isNonNegative(r.evaluate());
            }
            
            template<class A, template<typename, typename> class O>
            inline bool isNegative(const expressions::Expression<RealContext, A, O> & r)
            {
                return isNegative(r.evaluate());
            }
            
            template<class A, template<typename, typename> class O>
            inline bool isNonPositive(const expressions::Expression<RealContext, A, O> & r)
            {
                return isZero(r.evaluate());
            }
            
            template<class A, template<typename, typename> class O>
            inline std::ostream & operator << (std::ostream & s, const expressions::Expression<RealContext, A, O> & r)
            {
                return s << r.evaluate();
            }
        }
        
        /**@{
           \name Assignment-operation operators.
        */
        
        /**
           \brief Subtracts the floating point number `r` from `cur`.
           
           \param cur The floating point number to operate on.
           \param r The floating point number to subtract from `cur`.
           \return A reference to `cur` containing the result.
        */
        inline Real & operator -= (Real & cur, const Real & r)
        {
            sub(cur, cur, r);
            return cur;
        }
        
        /**
           \brief Subtracts the multiplication expression `E` from `cur`.
           
           Note that the faster `submul()` command is used for this special case where the product
           of two floating point numbers is subtracted from a third one.
           
           \param cur The floating point number to operate on.
           \param E The multiplication expression to subtract from `cur`.
           \return A reference to `cur` containing the result.
        */
        template<class D>
        inline Real & operator -= (Real & cur, const expressions::Expression<RealContext, D, expressions::MulOp> & E)
        {
            submul(cur, E.data().first.evaluate(), E.data().second.evaluate());
            return cur;
        }
        
        /**
           \brief Adds the floating point number `r` to `cur`.
           
           \param cur The floating point number to operate on.
           \param r The floating point number to add to `cur`.
           \return A reference to `cur` containing the result.
        */
        inline Real & operator += (Real & cur, const Real & r)
        {
            add(cur, cur, r);
            return cur;
        }
        
        /**
           \brief Adds the given multiplication expression to `cur`.
           
           Note that the faster `addmul()` command is used for this special case where the product
           of two floating point numbers is added to a third one.
           
           \param cur The floating point number to operate on.
           \param E The multiplication expression to add to `cur`.
           \return A reference to `cur` containing the result.
        */
        template<class D>
        inline Real & operator += (Real & cur, const expressions::Expression<RealContext, D, expressions::MulOp> & E)
        {
            addmul(cur, E.data().first.evaluate(), E.data().second.evaluate());
            return cur;
        }
        
        /**
           \brief Multiplies the given floating point number `r` with `cur`.
               
           \param cur The floating point number to operate on.
           \param r The floating point number to multiply to `cur`.
           \return A reference to `cur` containing the result.
        */
        inline Real & operator *= (Real & cur, const Real & r)
        {
            mul(cur, cur, r);
            return cur;
        }
        
        /**
           \brief Divides `cur` by the given floating point number `r`.
           
           \param cur The floating point number to operate on.
           \param r The floating point number to divide by.
           \return A reference to `cur` containing the result.
        */
        inline Real & operator /= (Real & cur, const Real & r)
        {
            div(cur, cur, r);
            return cur;
        }
        
        /**
           \brief Divides `cur` by the given floating point number `r` and
                  stores the remainder in `cur`.
           
           \param cur The floating point number to operate on.
           \param r The floating point number to divide by.
           \return A reference to `cur` containing the result.
        */
        inline Real & operator %= (Real & cur, const Real & r)
        {
            mod(cur, cur, r);
            return cur;
        }
        
        /**
           \brief Multiplies `cur` by 2 to the power of `r`.
           
           \param cur The floating point number to operate on.
           \param r The exponent.
           \return A reference to `cur` containing the result.
           
           \sa \ref arithctxts_sqrtfp
        */
        inline Real & operator <<= (Real & cur, const Real & r)
        {
            shl(cur, cur, r);
            return cur;
        }
        
        /**
           \brief Divides `cur` by 2 to the power of `r`.
           
           \param cur The floating point number to operate on.
           \param r The exponent.
           \return A reference to `cur` containing the result.
           
           \sa \ref arithctxts_sqrtfp
        */
        inline Real & operator >>= (Real & cur, const Real & r)
        {
            shr(cur, cur, r);
            return cur;
        }
        
        /**
           \brief Multiplies `cur` by 2 to the power of `r`.
           
           \param cur The floating point number to operate on.
           \param r The exponent.
           \return A reference to `cur` containing the result.
        */
        inline Real & operator <<= (Real & cur, long r)
        {
            shl(cur, cur, r);
            return cur;
        }
        
        /**
           \brief Divides `cur` by 2 to the power of `r`.
           
           \param cur The floating point number to operate on.
           \param r The exponent.
           \return A reference to `cur` containing the result.
        */
        inline Real & operator >>= (Real & cur, long r)
        {
            shr(cur, cur, r);
            return cur;
        }
        
        ///@}
        
        /**@{
           \name Operators.
        */
        
        /**
           \brief Increments `cur` by one and returns the previous value.
           
           \param cur The floating point number to operate on.
           \return The previous value before incrementing.
        */
        inline Real operator ++ (Real & cur, int)
        {
            Real r(cur, true);
            increment(cur, cur);
            return r;
        }
        
        /**
           \brief Decrements `cur` by one and returns the previous value.
           
           \param cur The floating point number to operate on.
           \return The previous value before decrementing.
        */
        inline Real operator -- (Real & cur, int)
        {
            Real r(cur, true);
            decrement(cur, cur);
            return r;
        }
        
        /**
           \brief Increments `cur` by one and returns the new value.
           
           \param cur The floating point number to operate on.
           \return The new value.
        */
        inline Real & operator ++ (Real & cur)
        {
            increment(cur, cur);
            return cur;
        }
        
        /**
           \brief Decrements `cur` by one and returns the new value.
           
           \param cur The floating point number to operate on.
           \return The new value.
        */
        inline Real & operator -- (Real & cur)
        {
            decrement(cur, cur);
            return cur;
        }
        
        ///@}
        
        /**@{
           \name Comparisons.
        */
        
        /**
           \brief Compares the two floating point numbers `a` and `b` for equality.
           
           \param a The first operand.
           \param b The second operand.
           \return `true` if and only if `a` equals `b`.
        */
        inline bool operator == (const Real & a, const Real & b)
        {
            return mpfr_equal_p(a.d_value, b.d_value);
        }
        
        /**
           \brief Compares the two floating point numbers `a` and `b` for inequality.
           
           \param a The first operand.
           \param b The second operand.
           \return `true` if and only if `a` does not equal `b`.
        */
        inline bool operator != (const Real & a, const Real & b)
        {
            return !mpfr_equal_p(a.d_value, b.d_value);
        }
        
        /**
           \brief Compares the two floating point numbers `a` and `b`.
               
           \param a The first operand.
           \param b The second operand.
           \return `true` if and only if `a` is less than or equal to `b`.
        */
        inline bool operator <= (const Real & a, const Real & b)
        {
            return mpfr_lessequal_p(a.d_value, b.d_value);
        }
        
        /**
           \brief Compares the two floating point numbers `a` and `b`.
               
           \param a The first operand.
           \param b The second operand.
           \return `true` if and only if `a` is greater than or equal to `b`.
        */
        inline bool operator >= (const Real & a, const Real & b)
        {
            return mpfr_greaterequal_p(a.d_value, b.d_value);
        }
        
        /**
           \brief Compares the two floating point numbers `a` and `b`.
               
           \param a The first operand.
           \param b The second operand.
           \return `true` if and only if `a` is less than `b`.
        */
        inline bool operator < (const Real & a, const Real & b)
        {
            return mpfr_less_p(a.d_value, b.d_value);
        }
        
        /**
           \brief Compares the two floating point numbers `a` and `b`.
               
           \param a The first operand.
           \param b The second operand.
           \return `true` if and only if `a` is greater than `b`.
        */
        inline bool operator > (const Real & a, const Real & b)
        {
            return mpfr_greater_p(a.d_value, b.d_value);
        }
        
        ///@}
        
        // Second operand is expression: also make the first one an expression
        template<class A2, template<typename, typename> class O2>
        inline bool operator == (const Real & a, const expressions::Expression<RealContext, A2, O2> & b)
        {
            return expressions::make_expression<RealContext>(a) == b;
        }
        
        template<class A2, template<typename, typename> class O2>
        inline bool operator != (const Real & a, const expressions::Expression<RealContext, A2, O2> & b)
        {
            return expressions::make_expression<RealContext>(a) != b;
        }
        
        template<class A2, template<typename, typename> class O2>
        inline bool operator <= (const Real & a, const expressions::Expression<RealContext, A2, O2> & b)
        {
            return expressions::make_expression<RealContext>(a) <= b;
        }
        
        template<class A2, template<typename, typename> class O2>
        inline bool operator >= (const Real & a, const expressions::Expression<RealContext, A2, O2> & b)
        {
            return expressions::make_expression<RealContext>(a) >= b;
        }
        
        template<class A2, template<typename, typename> class O2>
        inline bool operator < (const Real & a, const expressions::Expression<RealContext, A2, O2> & b)
        {
            return expressions::make_expression<RealContext>(a) < b;
        }
                
        template<class A2, template<typename, typename> class O2>
        inline bool operator > (const Real & a, const expressions::Expression<RealContext, A2, O2> & b)
        {
            return expressions::make_expression<RealContext>(a) > b;
        }
        
        namespace expressions
        {
            // First operand is expression: also make the second one an expression
            template<class A1, template<typename, typename> class O1>
            inline bool operator == (const expressions::Expression<RealContext, A1, O1> & a, const arithmetic::Real & b)
            {
                return a == expressions::make_expression<RealContext>(b);
            }
            
            template<class A1, template<typename, typename> class O1>
            inline bool operator != (const expressions::Expression<RealContext, A1, O1> & a, const arithmetic::Real & b)
            {
                return a != expressions::make_expression<RealContext>(b);
            }
            
            template<class A1, template<typename, typename> class O1>
            inline bool operator <= (const expressions::Expression<RealContext, A1, O1> & a, const arithmetic::Real & b)
            {
                return a <= expressions::make_expression<RealContext>(b);
            }
            
            template<class A1, template<typename, typename> class O1>
            inline bool operator >= (const expressions::Expression<RealContext, A1, O1> & a, const arithmetic::Real & b)
            {
                return a >= expressions::make_expression<RealContext>(b);
            }
            
            template<class A1, template<typename, typename> class O1>
            inline bool operator < (const expressions::Expression<RealContext, A1, O1> & a, const arithmetic::Real & b)
            {
                return a < expressions::make_expression<RealContext>(b);
            }
            
            template<class A1, template<typename, typename> class O1>
            inline bool operator > (const expressions::Expression<RealContext, A1, O1> & a, const arithmetic::Real & b)
            {
                return a > expressions::make_expression<RealContext>(b);
            }
            
            // Both operands are expressions
            template<class A1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
            inline bool operator == (const expressions::Expression<RealContext, A1, O1> & a, const expressions::Expression<RealContext, A2, O2> & b)
            {
                return a.evaluate() == b.evaluate();
            }
            
            template<class A1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
            inline bool operator != (const expressions::Expression<RealContext, A1, O1> & a, const expressions::Expression<RealContext, A2, O2> & b)
            {
                return a.evaluate() != b.evaluate();
            }
            
            template<class A1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
            inline bool operator <= (const expressions::Expression<RealContext, A1, O1> & a, const expressions::Expression<RealContext, A2, O2> & b)
            {
                return a.evaluate() <= b.evaluate();
            }
            
            template<class A1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
            inline bool operator >= (const expressions::Expression<RealContext, A1, O1> & a, const expressions::Expression<RealContext, A2, O2> & b)
            {
                return a.evaluate() >= b.evaluate();
            }
            
            template<class A1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
            inline bool operator < (const expressions::Expression<RealContext, A1, O1> & a, const expressions::Expression<RealContext, A2, O2> & b)
            {
                return a.evaluate() < b.evaluate();
            }
            
            template<class A1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
            inline bool operator > (const expressions::Expression<RealContext, A1, O1> & a, const expressions::Expression<RealContext, A2, O2> & b)
            {
                return a.evaluate() > b.evaluate();
            }
            
            // First operand is integer/double/... conversion expression, second operand arbitrary expression
            template<class A2, template<typename, typename> class O2>
            inline bool operator == (const expressions::Expression<RealContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> & a,
                                     const expressions::Expression<RealContext, A2, O2> & b)
            { return compare_ui(b.evaluate(), a.data().evaluate()) == 0; }
            template<class A2, template<typename, typename> class O2>
            inline bool operator != (const expressions::Expression<RealContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> & a,
                                     const expressions::Expression<RealContext, A2, O2> & b)
            { return compare_ui(b.evaluate(), a.data().evaluate()) != 0; }
            template<class A2, template<typename, typename> class O2>
            inline bool operator <= (const expressions::Expression<RealContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> & a,
                                     const expressions::Expression<RealContext, A2, O2> & b)
            { return compare_ui(b.evaluate(), a.data().evaluate()) >= 0; }
            template<class A2, template<typename, typename> class O2>
            inline bool operator >= (const expressions::Expression<RealContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> & a,
                                     const expressions::Expression<RealContext, A2, O2> & b)
            { return compare_ui(b.evaluate(), a.data().evaluate()) <= 0; }
            template<class A2, template<typename, typename> class O2>
            inline bool operator < (const expressions::Expression<RealContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> & a,
                                    const expressions::Expression<RealContext, A2, O2> & b)
            { return compare_ui(b.evaluate(), a.data().evaluate()) > 0; }
            template<class A2, template<typename, typename> class O2>
            inline bool operator > (const expressions::Expression<RealContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> & a,
                                    const expressions::Expression<RealContext, A2, O2> & b)
            { return compare_ui(b.evaluate(), a.data().evaluate()) < 0; }
            
            template<class A2, template<typename, typename> class O2>
            inline bool operator == (const expressions::Expression<RealContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> & a,
                                     const expressions::Expression<RealContext, A2, O2> & b)
            { return compare_si(b.evaluate(), a.data().evaluate()) == 0; }
            template<class A2, template<typename, typename> class O2>
            inline bool operator != (const expressions::Expression<RealContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> & a,
                                     const expressions::Expression<RealContext, A2, O2> & b)
            { return compare_si(b.evaluate(), a.data().evaluate()) != 0; }
            template<class A2, template<typename, typename> class O2>
            inline bool operator <= (const expressions::Expression<RealContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> & a,
                                     const expressions::Expression<RealContext, A2, O2> & b)
            { return compare_si(b.evaluate(), a.data().evaluate()) >= 0; }
            template<class A2, template<typename, typename> class O2>
            inline bool operator >= (const expressions::Expression<RealContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> & a,
                                     const expressions::Expression<RealContext, A2, O2> & b)
            { return compare_si(b.evaluate(), a.data().evaluate()) <= 0; }
            template<class A2, template<typename, typename> class O2>
            inline bool operator < (const expressions::Expression<RealContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> & a,
                                    const expressions::Expression<RealContext, A2, O2> & b)
            { return compare_si(b.evaluate(), a.data().evaluate()) > 0; }
            template<class A2, template<typename, typename> class O2>
            inline bool operator > (const expressions::Expression<RealContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> & a,
                                    const expressions::Expression<RealContext, A2, O2> & b)
            { return compare_si(b.evaluate(), a.data().evaluate()) < 0; }
            
            template<class A2, template<typename, typename> class O2>
            inline bool operator == (const expressions::Expression<RealContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> & a,
                                     const expressions::Expression<RealContext, A2, O2> & b)
            { return compare_d(b.evaluate(), a.data().evaluate()) == 0; }
            template<class A2, template<typename, typename> class O2>
            inline bool operator != (const expressions::Expression<RealContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> & a,
                                     const expressions::Expression<RealContext, A2, O2> & b)
            { return compare_d(b.evaluate(), a.data().evaluate()) != 0; }
            template<class A2, template<typename, typename> class O2>
            inline bool operator <= (const expressions::Expression<RealContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> & a,
                                     const expressions::Expression<RealContext, A2, O2> & b)
            { return compare_d(b.evaluate(), a.data().evaluate()) >= 0; }
            template<class A2, template<typename, typename> class O2>
            inline bool operator >= (const expressions::Expression<RealContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> & a,
                                     const expressions::Expression<RealContext, A2, O2> & b)
            { return compare_d(b.evaluate(), a.data().evaluate()) <= 0; }
            template<class A2, template<typename, typename> class O2>
            inline bool operator < (const expressions::Expression<RealContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> & a,
                                    const expressions::Expression<RealContext, A2, O2> & b)
            { return compare_d(b.evaluate(), a.data().evaluate()) > 0; }
            template<class A2, template<typename, typename> class O2>
            inline bool operator > (const expressions::Expression<RealContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> & a,
                                    const expressions::Expression<RealContext, A2, O2> & b)
            { return compare_d(b.evaluate(), a.data().evaluate()) < 0; }
            
            template<class A2, template<typename, typename> class O2>
            inline bool operator == (const expressions::Expression<RealContext, expressions::ConversionWrapper<long double>, expressions::ConvertOp_Context> & a,
                                     const expressions::Expression<RealContext, A2, O2> & b)
            { return compare_ld(b.evaluate(), a.data().evaluate()) == 0; }
            template<class A2, template<typename, typename> class O2>
            inline bool operator != (const expressions::Expression<RealContext, expressions::ConversionWrapper<long double>, expressions::ConvertOp_Context> & a,
                                     const expressions::Expression<RealContext, A2, O2> & b)
            { return compare_ld(b.evaluate(), a.data().evaluate()) != 0; }
            template<class A2, template<typename, typename> class O2>
            inline bool operator <= (const expressions::Expression<RealContext, expressions::ConversionWrapper<long double>, expressions::ConvertOp_Context> & a,
                                     const expressions::Expression<RealContext, A2, O2> & b)
            { return compare_ld(b.evaluate(), a.data().evaluate()) >= 0; }
            template<class A2, template<typename, typename> class O2>
            inline bool operator >= (const expressions::Expression<RealContext, expressions::ConversionWrapper<long double>, expressions::ConvertOp_Context> & a,
                                     const expressions::Expression<RealContext, A2, O2> & b)
            { return compare_ld(b.evaluate(), a.data().evaluate()) <= 0; }
            template<class A2, template<typename, typename> class O2>
            inline bool operator < (const expressions::Expression<RealContext, expressions::ConversionWrapper<long double>, expressions::ConvertOp_Context> & a,
                                    const expressions::Expression<RealContext, A2, O2> & b)
            { return compare_ld(b.evaluate(), a.data().evaluate()) > 0; }
            template<class A2, template<typename, typename> class O2>
            inline bool operator > (const expressions::Expression<RealContext, expressions::ConversionWrapper<long double>, expressions::ConvertOp_Context> & a,
                                    const expressions::Expression<RealContext, A2, O2> & b)
            { return compare_ld(b.evaluate(), a.data().evaluate()) < 0; }
            
            template<class D1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
            inline bool operator == (const Expression<RealContext, Expression<IntegerContext, D1, O1>, ConvertOp_Context> & a,
                                     const expressions::Expression<RealContext, A2, O2> & b)
            { return compare_z(b.evaluate(), a.data().evaluate()) == 0; }
            template<class D1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
            inline bool operator != (const Expression<RealContext, Expression<IntegerContext, D1, O1>, ConvertOp_Context> & a,
                                     const expressions::Expression<RealContext, A2, O2> & b)
            { return compare_z(b.evaluate(), a.data().evaluate()) != 0; }
            template<class D1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
            inline bool operator <= (const Expression<RealContext, Expression<IntegerContext, D1, O1>, ConvertOp_Context> & a,
                                     const expressions::Expression<RealContext, A2, O2> & b)
            { return compare_z(b.evaluate(), a.data().evaluate()) >= 0; }
            template<class D1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
            inline bool operator >= (const Expression<RealContext, Expression<IntegerContext, D1, O1>, ConvertOp_Context> & a,
                                     const expressions::Expression<RealContext, A2, O2> & b)
            { return compare_z(b.evaluate(), a.data().evaluate()) <= 0; }
            template<class D1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
            inline bool operator < (const Expression<RealContext, Expression<IntegerContext, D1, O1>, ConvertOp_Context> & a,
                                    const expressions::Expression<RealContext, A2, O2> & b)
            { return compare_z(b.evaluate(), a.data().evaluate()) > 0; }
            template<class D1, template<typename, typename> class O1, class A2, template<typename, typename> class O2>
            inline bool operator > (const Expression<RealContext, Expression<IntegerContext, D1, O1>, ConvertOp_Context> & a,
                                    const expressions::Expression<RealContext, A2, O2> & b)
            { return compare_z(b.evaluate(), a.data().evaluate()) < 0; }
            
            template<class A2, template<typename, typename> class O2>
            inline bool operator == (const expressions::Expression<RealContext, expressions::Expression<RealContext, expressions::ConversionWrapper<unsigned long>,
                                                                                                        expressions::ConvertOp_Context>, expressions::ShiftCOp> & a,
                                     const expressions::Expression<RealContext, A2, O2> & b)
            { return compare_ui2(b.evaluate(), a.data().data().evaluate(), a.op().shift()) == 0; }
            template<class A2, template<typename, typename> class O2>
            inline bool operator != (const expressions::Expression<RealContext, expressions::Expression<RealContext, expressions::ConversionWrapper<unsigned long>,
                                                                                                        expressions::ConvertOp_Context>, expressions::ShiftCOp> & a,
                                     const expressions::Expression<RealContext, A2, O2> & b)
            { return compare_ui2(b.evaluate(), a.data().data().evaluate(), a.op().shift()) != 0; }
            template<class A2, template<typename, typename> class O2>
            inline bool operator <= (const expressions::Expression<RealContext, expressions::Expression<RealContext, expressions::ConversionWrapper<unsigned long>,
                                                                                                        expressions::ConvertOp_Context>, expressions::ShiftCOp> & a,
                                     const expressions::Expression<RealContext, A2, O2> & b)
            { return compare_ui2(b.evaluate(), a.data().data().evaluate(), a.op().shift()) >= 0; }
            template<class A2, template<typename, typename> class O2>
            inline bool operator >= (const expressions::Expression<RealContext, expressions::Expression<RealContext, expressions::ConversionWrapper<unsigned long>,
                                                                                                        expressions::ConvertOp_Context>, expressions::ShiftCOp> & a,
                                     const expressions::Expression<RealContext, A2, O2> & b)
            { return compare_ui2(b.evaluate(), a.data().data().evaluate(), a.op().shift()) <= 0; }
            template<class A2, template<typename, typename> class O2>
            inline bool operator < (const expressions::Expression<RealContext, expressions::Expression<RealContext, expressions::ConversionWrapper<unsigned long>,
                                                                                                       expressions::ConvertOp_Context>, expressions::ShiftCOp> & a,
                                    const expressions::Expression<RealContext, A2, O2> & b)
            { return compare_ui2(b.evaluate(), a.data().data().evaluate(), a.op().shift()) > 0; }
            template<class A2, template<typename, typename> class O2>
            inline bool operator > (const expressions::Expression<RealContext, expressions::Expression<RealContext, expressions::ConversionWrapper<unsigned long>,
                                                                                                       expressions::ConvertOp_Context>, expressions::ShiftCOp> & a,
                                    const expressions::Expression<RealContext, A2, O2> & b)
            { return compare_ui2(b.evaluate(), a.data().data().evaluate(), a.op().shift()) < 0; }
            
            // First operand is arbitrary expression, Second operand is integer/double/... conversion expression
            template<class A1, template<typename, typename> class O1>
            inline bool operator == (const expressions::Expression<RealContext, A1, O1> & a,
                                     const expressions::Expression<RealContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> & b)
            { return compare_ui(a.evaluate(), b.data().evaluate()) == 0; }
            template<class A1, template<typename, typename> class O1>
            inline bool operator != (const expressions::Expression<RealContext, A1, O1> & a,
                                     const expressions::Expression<RealContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> & b)
            { return compare_ui(a.evaluate(), b.data().evaluate()) != 0; }
            template<class A1, template<typename, typename> class O1>
            inline bool operator <= (const expressions::Expression<RealContext, A1, O1> & a,
                                     const expressions::Expression<RealContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> & b)
            { return compare_ui(a.evaluate(), b.data().evaluate()) <= 0; }
            template<class A1, template<typename, typename> class O1>
            inline bool operator >= (const expressions::Expression<RealContext, A1, O1> & a,
                                     const expressions::Expression<RealContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> & b)
            { return compare_ui(a.evaluate(), b.data().evaluate()) >= 0; }
            template<class A1, template<typename, typename> class O1>
            inline bool operator < (const expressions::Expression<RealContext, A1, O1> & a,
                                    const expressions::Expression<RealContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> & b)
            { return compare_ui(a.evaluate(), b.data().evaluate()) < 0; }
            template<class A1, template<typename, typename> class O1>
            inline bool operator > (const expressions::Expression<RealContext, A1, O1> & a,
                                    const expressions::Expression<RealContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> & b)
            { return compare_ui(a.evaluate(), b.data().evaluate()) > 0; }
            
            template<class A1, template<typename, typename> class O1>
            inline bool operator == (const expressions::Expression<RealContext, A1, O1> & a,
                                     const expressions::Expression<RealContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> & b)
            { return compare_si(a.evaluate(), b.data().evaluate()) == 0; }
            template<class A1, template<typename, typename> class O1>
            inline bool operator != (const expressions::Expression<RealContext, A1, O1> & a,
                                     const expressions::Expression<RealContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> & b)
            { return compare_si(a.evaluate(), b.data().evaluate()) != 0; }
            template<class A1, template<typename, typename> class O1>
            inline bool operator <= (const expressions::Expression<RealContext, A1, O1> & a,
                                     const expressions::Expression<RealContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> & b)
            { return compare_si(a.evaluate(), b.data().evaluate()) <= 0; }
            template<class A1, template<typename, typename> class O1>
            inline bool operator >= (const expressions::Expression<RealContext, A1, O1> & a,
                                     const expressions::Expression<RealContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> & b)
            { return compare_si(a.evaluate(), b.data().evaluate()) >= 0; }
            template<class A1, template<typename, typename> class O1>
            inline bool operator < (const expressions::Expression<RealContext, A1, O1> & a,
                                    const expressions::Expression<RealContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> & b)
            { return compare_si(a.evaluate(), b.data().evaluate()) < 0; }
            template<class A1, template<typename, typename> class O1>
            inline bool operator > (const expressions::Expression<RealContext, A1, O1> & a,
                                    const expressions::Expression<RealContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> & b)
            { return compare_si(a.evaluate(), b.data().evaluate()) > 0; }
            
            template<class A1, template<typename, typename> class O1>
            inline bool operator == (const expressions::Expression<RealContext, A1, O1> & a,
                                     const expressions::Expression<RealContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> & b)
            { return compare_d(a.evaluate(), b.data().evaluate()) == 0; }
            template<class A1, template<typename, typename> class O1>
            inline bool operator != (const expressions::Expression<RealContext, A1, O1> & a,
                                     const expressions::Expression<RealContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> & b)
            { return compare_d(a.evaluate(), b.data().evaluate()) != 0; }
            template<class A1, template<typename, typename> class O1>
            inline bool operator <= (const expressions::Expression<RealContext, A1, O1> & a,
                                     const expressions::Expression<RealContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> & b)
            { return compare_d(a.evaluate(), b.data().evaluate()) <= 0; }
            template<class A1, template<typename, typename> class O1>
            inline bool operator >= (const expressions::Expression<RealContext, A1, O1> & a,
                                     const expressions::Expression<RealContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> & b)
            { return compare_d(a.evaluate(), b.data().evaluate()) >= 0; }
            template<class A1, template<typename, typename> class O1>
            inline bool operator < (const expressions::Expression<RealContext, A1, O1> & a,
                                    const expressions::Expression<RealContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> & b)
            { return compare_d(a.evaluate(), b.data().evaluate()) < 0; }
            template<class A1, template<typename, typename> class O1>
            inline bool operator > (const expressions::Expression<RealContext, A1, O1> & a,
                                    const expressions::Expression<RealContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> & b)
            { return compare_d(a.evaluate(), b.data().evaluate()) > 0; }
            
            template<class A1, template<typename, typename> class O1>
            inline bool operator == (const expressions::Expression<RealContext, A1, O1> & a,
                                     const expressions::Expression<RealContext, expressions::ConversionWrapper<long double>, expressions::ConvertOp_Context> & b)
            { return compare_ld(a.evaluate(), b.data().evaluate()) == 0; }
            template<class A1, template<typename, typename> class O1>
            inline bool operator != (const expressions::Expression<RealContext, A1, O1> & a,
                                     const expressions::Expression<RealContext, expressions::ConversionWrapper<long double>, expressions::ConvertOp_Context> & b)
            { return compare_ld(a.evaluate(), b.data().evaluate()) != 0; }
            template<class A1, template<typename, typename> class O1>
            inline bool operator <= (const expressions::Expression<RealContext, A1, O1> & a,
                                     const expressions::Expression<RealContext, expressions::ConversionWrapper<long double>, expressions::ConvertOp_Context> & b)
            { return compare_ld(a.evaluate(), b.data().evaluate()) <= 0; }
            template<class A1, template<typename, typename> class O1>
            inline bool operator >= (const expressions::Expression<RealContext, A1, O1> & a,
                                     const expressions::Expression<RealContext, expressions::ConversionWrapper<long double>, expressions::ConvertOp_Context> & b)
            { return compare_ld(a.evaluate(), b.data().evaluate()) >= 0; }
            template<class A1, template<typename, typename> class O1>
            inline bool operator < (const expressions::Expression<RealContext, A1, O1> & a,
                                    const expressions::Expression<RealContext, expressions::ConversionWrapper<long double>, expressions::ConvertOp_Context> & b)
            { return compare_ld(a.evaluate(), b.data().evaluate()) < 0; }
            template<class A1, template<typename, typename> class O1>
            inline bool operator > (const expressions::Expression<RealContext, A1, O1> & a,
                                    const expressions::Expression<RealContext, expressions::ConversionWrapper<long double>, expressions::ConvertOp_Context> & b)
            { return compare_ld(a.evaluate(), b.data().evaluate()) > 0; }
            
            template<class A1, template<typename, typename> class O1, class D2, template<typename, typename> class O2>
            inline bool operator == (const expressions::Expression<RealContext, A1, O1> & a,
                                     const Expression<RealContext, Expression<IntegerContext, D2, O2>, ConvertOp_Context> & b)
            { return compare_z(a.evaluate(), b.data().evaluate()) == 0; }
            template<class A1, template<typename, typename> class O1, class D2, template<typename, typename> class O2>
            inline bool operator != (const expressions::Expression<RealContext, A1, O1> & a,
                                     const Expression<RealContext, Expression<IntegerContext, D2, O2>, ConvertOp_Context> & b)
            { return compare_z(a.evaluate(), b.data().evaluate()) != 0; }
            template<class A1, template<typename, typename> class O1, class D2, template<typename, typename> class O2>
            inline bool operator <= (const expressions::Expression<RealContext, A1, O1> & a,
                                     const Expression<RealContext, Expression<IntegerContext, D2, O2>, ConvertOp_Context> & b)
            { return compare_z(a.evaluate(), b.data().evaluate()) <= 0; }
            template<class A1, template<typename, typename> class O1, class D2, template<typename, typename> class O2>
            inline bool operator >= (const expressions::Expression<RealContext, A1, O1> & a,
                                     const Expression<RealContext, Expression<IntegerContext, D2, O2>, ConvertOp_Context> & b)
            { return compare_z(a.evaluate(), b.data().evaluate()) >= 0; }
            template<class A1, template<typename, typename> class O1, class D2, template<typename, typename> class O2>
            inline bool operator < (const expressions::Expression<RealContext, A1, O1> & a,
                                    const Expression<RealContext, Expression<IntegerContext, D2, O2>, ConvertOp_Context> & b)
            { return compare_z(a.evaluate(), b.data().evaluate()) < 0; }
            template<class A1, template<typename, typename> class O1, class D2, template<typename, typename> class O2>
            inline bool operator > (const expressions::Expression<RealContext, A1, O1> & a,
                                    const Expression<RealContext, Expression<IntegerContext, D2, O2>, ConvertOp_Context> & b)
            { return compare_z(a.evaluate(), b.data().evaluate()) > 0; }
            
            template<class A2, template<typename, typename> class O2>
            inline bool operator == (const expressions::Expression<RealContext, A2, O2> & a,
                                     const expressions::Expression<RealContext, expressions::Expression<RealContext, expressions::ConversionWrapper<unsigned long>,
                                                                                                        expressions::ConvertOp_Context>, expressions::ShiftCOp> & b)
            { return compare_ui2(a.evaluate(), b.data().data().evaluate(), b.op().shift()) == 0; }
            template<class A2, template<typename, typename> class O2>
            inline bool operator != (const expressions::Expression<RealContext, A2, O2> & a,
                                     const expressions::Expression<RealContext, expressions::Expression<RealContext, expressions::ConversionWrapper<unsigned long>,
                                                                                                        expressions::ConvertOp_Context>, expressions::ShiftCOp> & b)
            { return compare_ui2(a.evaluate(), b.data().data().evaluate(), b.op().shift()) != 0; }
            template<class A2, template<typename, typename> class O2>
            inline bool operator <= (const expressions::Expression<RealContext, A2, O2> & a,
                                     const expressions::Expression<RealContext, expressions::Expression<RealContext, expressions::ConversionWrapper<unsigned long>,
                                                                                                        expressions::ConvertOp_Context>, expressions::ShiftCOp> & b)
            { return compare_ui2(a.evaluate(), b.data().data().evaluate(), b.op().shift()) <= 0; }
            template<class A2, template<typename, typename> class O2>
            inline bool operator >= (const expressions::Expression<RealContext, A2, O2> & a,
                                     const expressions::Expression<RealContext, expressions::Expression<RealContext, expressions::ConversionWrapper<unsigned long>,
                                                                                                        expressions::ConvertOp_Context>, expressions::ShiftCOp> & b)
            { return compare_ui2(a.evaluate(), b.data().data().evaluate(), b.op().shift()) >= 0; }
            template<class A2, template<typename, typename> class O2>
            inline bool operator < (const expressions::Expression<RealContext, A2, O2> & a,
                                    const expressions::Expression<RealContext, expressions::Expression<RealContext, expressions::ConversionWrapper<unsigned long>,
                                                                                                       expressions::ConvertOp_Context>, expressions::ShiftCOp> & b)
            { return compare_ui2(a.evaluate(), b.data().data().evaluate(), b.op().shift()) < 0; }
            template<class A2, template<typename, typename> class O2>
            inline bool operator > (const expressions::Expression<RealContext, A2, O2> & a,
                                    const expressions::Expression<RealContext, expressions::Expression<RealContext, expressions::ConversionWrapper<unsigned long>,
                                                                                                       expressions::ConvertOp_Context>, expressions::ShiftCOp> & b)
            { return compare_ui2(a.evaluate(), b.data().data().evaluate(), b.op().shift()) > 0; }
        }
            
        /**@{
           \name Operators.
        */
        /**
           \brief Returns the absolute value of the given floating point number.
           
           \param i The operand.
           \return The result.
        */
        inline expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                    expressions::AbsOp> abs(const Real & i)
        { return expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                      expressions::AbsOp>(expressions::Wrapper<RealContext>(i)); }
        /**
           \brief Returns the square of the given floating point number.
           
           \param i The operand.
           \return The result.
        */
        inline expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                    expressions::SquareOp> square(const Real & i)
        { return expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                      expressions::SquareOp>(expressions::Wrapper<RealContext>(i)); }
        ///@}
        /**@{
           \name Trigonometric functions.
        */
        /**
           \brief Returns the sine of the given floating point number.
           
           \param i The operand.
           \return The result.
           \sa \ref arithctxts_trig
        */
        inline expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                    expressions::SinOp> sin(const Real & i)
        { return expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                      expressions::SinOp>(expressions::Wrapper<RealContext>(i)); }
        /**
           \brief Returns the cosine of the given floating point number.
           
           \param i The operand.
           \return The result.
           \sa \ref arithctxts_trig
        */
        inline expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                    expressions::CosOp> cos(const Real & i)
        { return expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                      expressions::CosOp>(expressions::Wrapper<RealContext>(i)); }
        /**
           \brief Returns the tangent of the given floating point number.
           
           \param i The operand.
           \return The result.
           \sa \ref arithctxts_trig
        */
        inline expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                    expressions::TanOp> tan(const Real & i)
        { return expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                      expressions::TanOp>(expressions::Wrapper<RealContext>(i)); }
        /**
           \brief Returns the arcsine of the given floating point number.
           
           \param i The operand.
           \return The result.
           \sa \ref arithctxts_trig
        */
        inline expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                    expressions::ASinOp> asin(const Real & i)
        { return expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                      expressions::ASinOp>(expressions::Wrapper<RealContext>(i)); }
        /**
           \brief Returns the arccosine of the given floating point number.
           
           \param i The operand.
           \return The result.
           \sa \ref arithctxts_trig
        */
        inline expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                    expressions::ACosOp> acos(const Real & i)
        { return expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                      expressions::ACosOp>(expressions::Wrapper<RealContext>(i)); }
        /**
           \brief Returns the arctangent of the given floating point number.
           
           \param i The operand.
           \return The result.
           \sa \ref arithctxts_trig
        */
        inline expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                    expressions::ATanOp> atan(const Real & i)
        { return expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                      expressions::ATanOp>(expressions::Wrapper<RealContext>(i)); }
        /**
           \brief Returns the arctangent of \f$\tfrac{y}{x}\f$, using the signs of `x` and `y` to
                  determine the quadrant.
                  
           \param y The numerator.
           \param x The denominator.
           \return The result.
           \sa \ref arithctxts_trig
        */
        inline expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                              expressions::Wrapper<RealContext> >,
                                                    expressions::ATan2Op> atan2(const Real & y, const Real & x)
        { return expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                                expressions::Wrapper<RealContext> >,
                                                      expressions::ATan2Op>(std::make_pair(expressions::Wrapper<RealContext>(y),
                                                                                           expressions::Wrapper<RealContext>(x))); }
        ///@}
        /**@{
           \name Exponential and logarithmic functions.
        */
        /**
           \brief Returns the exponential function evaluated at the given floating point number.
           
           \param i The operand.
           \return The result.
           \sa \ref arithctxts_specfns
        */
        inline expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                    expressions::ExpOp> exp(const Real & i)
        { return expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                      expressions::ExpOp>(expressions::Wrapper<RealContext>(i)); }
        /**
           \brief Returns the natural logarithm evaluated at the given floating point number.
           
           \param i The operand.
           \return The result.
           \sa \ref arithctxts_specfns
        */
        inline expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                    expressions::LogOp> log(const Real & i)
        { return expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                      expressions::LogOp>(expressions::Wrapper<RealContext>(i)); }
        /**
           \brief Returns the base-2 logarithm evaluated at the given floating point number.
           
           \param i The operand.
           \return The result.
           \sa \ref arithctxts_specfns
        */
        inline expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                    expressions::Log2Op> log2(const Real & i)
        { return expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                      expressions::Log2Op>(expressions::Wrapper<RealContext>(i)); }
        /**
           \brief Returns the base-10 logarithm evaluated at the given floating point number.
           
           \param i The operand.
           \return The result.
           \sa \ref arithctxts_specfns
        */
        inline expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                    expressions::Log10Op> log10(const Real & i)
        { return expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                      expressions::Log10Op>(expressions::Wrapper<RealContext>(i)); }
        /**
           \brief Returns the square root of the given floating point number.
           
           \param i The operand.
           \return The result.
           \sa \ref arithctxts_sqrtfp
        */
        inline expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                    expressions::SqrtOp> sqrt(const Real & i)
        { return expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                      expressions::SqrtOp>(expressions::Wrapper<RealContext>(i)); }
        ///@}
        /**@{
           \name Special functions.
        */
        /**
           \brief Returns the Gamma function evaluated at the given floating point number.
           
           \param i The operand.
           \return The result.
           \sa \ref arithctxts_specfns
        */
        inline expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                    expressions::GammaOp> gamma(const Real & i)
        { return expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                      expressions::GammaOp>(expressions::Wrapper<RealContext>(i)); }
        /**
           \brief Returns the natural logarithm of the absolute value of the Gamma function
                  evaluated at the given floating point number.
           
           \param i The operand.
           \return The result.
           \sa \ref arithctxts_specfns
        */
        inline expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                    expressions::LGammaOp> lgamma(const Real & i)
        { return expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                      expressions::LGammaOp>(expressions::Wrapper<RealContext>(i)); }
        /**
           \brief Returns the natural logarithm of the absolute value of the Gamma function
                  evaluated at the given floating point number, as well as the sign of the Gamma
                  function.
                  
           \param sign The integer variable to store the sign of \f$\Gamma(i)\f$ in.
           \param i The operand.
           \return The result.
           \sa \ref arithctxts_specfns
        */
        inline expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                    expressions::LGamma2Op> lgamma(int & sign, const Real & i)
        { return expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                      expressions::LGamma2Op>(expressions::Wrapper<RealContext>(i), sign); }
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
        inline expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                    expressions::PowerCOp<signed long>::impl> power(const Real & a, signed long b)
        { return expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                      expressions::PowerCOp<signed long>::impl>(expressions::Wrapper<RealContext>(a), b); }
        /**
           \brief Returns `a` raised to the power of `b`.
           
           \param a The base.
           \param b The exponent.
           \return The result.
        */
        inline expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                    expressions::PowerCOp<unsigned long>::impl> power(const Real & a, unsigned long b)
        { return expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                      expressions::PowerCOp<unsigned long>::impl>(expressions::Wrapper<RealContext>(a), b); }
        /**
           \brief Returns `a` raised to the power of `b`.
           
           \param a The base.
           \param b The exponent.
           \return The result.
        */
        inline expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                              expressions::Wrapper<IntegerContext> >,
                                                    expressions::PowerOp> power(const Real & a, const Integer & b)
        { return expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                                expressions::Wrapper<IntegerContext> >,
                                                      expressions::PowerOp>(std::make_pair(expressions::Wrapper<RealContext>(a),
                                                                                           expressions::Wrapper<IntegerContext>(b))); }
        /**
           \brief Returns `a` raised to the power of `b`.
           
           \param a The base.
           \param b The exponent.
           \return The result.
           \sa \ref arithctxts_sqrtfp
        */
        inline expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                              expressions::Wrapper<RealContext> >,
                                                    expressions::PowerOp> power(const Real & a, const Real & b)
        { return expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                                expressions::Wrapper<RealContext> >,
                                                      expressions::PowerOp>(std::make_pair(expressions::Wrapper<RealContext>(a),
                                                                                           expressions::Wrapper<RealContext>(b))); }
        ///@}
        
        template<class A2, template<typename, typename> class O2>
        inline expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                              expressions::Expression<RealContext, A2, O2> >,
                                                    expressions::ATan2Op> atan2(const Real & a, const expressions::Expression<RealContext, A2, O2> & b)
        { return expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                                expressions::Expression<RealContext, A2, O2> >,
                                                      expressions::ATan2Op>(std::make_pair(expressions::Wrapper<RealContext>(a), b)); }
        template<class A2, template<typename, typename> class O2>
        
        inline expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                              expressions::Expression<RealContext, A2, O2> >,
                                                    expressions::PowerOp> power(const Real & a, const expressions::Expression<RealContext, A2, O2> & b)
        { return expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                                expressions::Expression<RealContext, A2, O2> >,
                                                      expressions::PowerOp>(std::make_pair(expressions::Wrapper<RealContext>(a), b)); }
        
        namespace expressions
        {
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                           expressions::AbsOp> abs(const expressions::Expression<RealContext, A, O> & b)
            { return expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                             expressions::AbsOp>(b); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                           expressions::SquareOp> square(const expressions::Expression<RealContext, A, O> & b)
            { return expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                             expressions::SquareOp>(b); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                           expressions::SinOp> sin(const expressions::Expression<RealContext, A, O> & b)
            { return expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                             expressions::SinOp>(b); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                           expressions::CosOp> cos(const expressions::Expression<RealContext, A, O> & b)
            { return expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                             expressions::CosOp>(b); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                           expressions::TanOp> tan(const expressions::Expression<RealContext, A, O> & b)
            { return expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                             expressions::TanOp>(b); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                           expressions::ASinOp> asin(const expressions::Expression<RealContext, A, O> & b)
            { return expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                             expressions::ASinOp>(b); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                           expressions::ACosOp> acos(const expressions::Expression<RealContext, A, O> & b)
            { return expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                             expressions::ACosOp>(b); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                           expressions::ATanOp> atan(const expressions::Expression<RealContext, A, O> & b)
            { return expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                             expressions::ATanOp>(b); }
            template<class A, template<typename, typename> class O, class A2, template<typename, typename> class O2>
            inline expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A, O>,
                                                                  expressions::Expression<RealContext, A2, O2> >,
                                           expressions::ATan2Op> atan2(const expressions::Expression<RealContext, A, O> & a,
                                                                       const expressions::Expression<RealContext, A2, O2> & b)
            { return expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A, O>,
                                                                    expressions::Expression<RealContext, A2, O2> >,
                                             expressions::ATan2Op>(std::make_pair(a, b)); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A, O>,
                                                                  expressions::Wrapper<RealContext> >,
                                           expressions::ATan2Op> atan2(const expressions::Expression<RealContext, A, O> & a, const arithmetic::Real & b)
            { return expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A, O>,
                                                                    expressions::Wrapper<RealContext> >,
                                             expressions::ATan2Op>(std::make_pair(a, expressions::Wrapper<RealContext>(b))); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                           expressions::ExpOp> exp(const expressions::Expression<RealContext, A, O> & b)
            { return expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                             expressions::ExpOp>(b); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                           expressions::LogOp> log(const expressions::Expression<RealContext, A, O> & b)
            { return expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                             expressions::LogOp>(b); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                           expressions::Log2Op> log2(const expressions::Expression<RealContext, A, O> & b)
            { return expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                             expressions::Log2Op>(b); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                           expressions::Log10Op> log10(const expressions::Expression<RealContext, A, O> & b)
            { return expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                             expressions::Log10Op>(b); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                           expressions::SqrtOp> sqrt(const expressions::Expression<RealContext, A, O> & b)
            { return expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                             expressions::SqrtOp>(b); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                           expressions::GammaOp> gamma(const expressions::Expression<RealContext, A, O> & b)
            { return expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                             expressions::GammaOp>(b); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                           expressions::LGammaOp> lgamma(const expressions::Expression<RealContext, A, O> & b)
            { return expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                             expressions::LGammaOp>(b); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                           expressions::LGamma2Op> lgamma(int & sign, const expressions::Expression<RealContext, A, O> & b)
            { return expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                             expressions::LGamma2Op>(b, sign); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                           expressions::PowerCOp<signed long>::impl> power(const expressions::Expression<RealContext, A, O> & a, signed long b)
            { return expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                             expressions::PowerCOp<signed long>::impl>(a, b); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                           expressions::PowerCOp<unsigned long>::impl> power(const expressions::Expression<RealContext, A, O> & a, unsigned long b)
            { return expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                             expressions::PowerCOp<unsigned long>::impl>(a, b); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A, O>,
                                                                  expressions::Wrapper<IntegerContext> >,
                                           expressions::PowerOp> power(const expressions::Expression<RealContext, A, O> & a, const Integer & b)
            { return expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A, O>,
                                                                    expressions::Wrapper<IntegerContext> >,
                                             expressions::PowerOp>(std::make_pair(a, expressions::Wrapper<IntegerContext>(b))); }
            template<class A, template<typename, typename> class O, class A2, template<typename, typename> class O2>
            inline expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A, O>,
                                                                  expressions::Expression<RealContext, A2, O2> >,
                                           expressions::PowerOp> power(const expressions::Expression<RealContext, A, O> & a,
                                                                       const expressions::Expression<RealContext, A2, O2> & b)
            { return expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A, O>,
                                                                    expressions::Expression<RealContext, A2, O2> >,
                                             expressions::PowerOp>(std::make_pair(a, b)); }
            template<class A, template<typename, typename> class O>
            inline expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A, O>,
                                                                  expressions::Wrapper<RealContext> >,
                                           expressions::PowerOp> power(const expressions::Expression<RealContext, A, O> & a, const arithmetic::Real & b)
            { return expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A, O>,
                                                                    expressions::Wrapper<RealContext> >,
                                             expressions::PowerOp>(std::make_pair(a, expressions::Wrapper<RealContext>(b))); }
        }
        
        /**@{
           \name Operators.
        */
        /**
           \brief Returns the absolute value of the given floating point number.
           
           \param i The operand.
           \param rc The context whose precision to use for the result.
           \return The result.
        */
        inline expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                    expressions::AbsOp_Context> abs(const Real & i, const RealContext & rc)
        { return expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                      expressions::AbsOp_Context>(expressions::Wrapper<RealContext>(i), rc); }
        /**
           \brief Returns the square of the given floating point number.
           
           \param i The operand.
           \param rc The context whose precision to use for the result.
           \return The result.
        */
        inline expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                    expressions::SquareOp_Context> square(const Real & i, const RealContext & rc)
        { return expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                      expressions::SquareOp_Context>(expressions::Wrapper<RealContext>(i), rc); }
        ///@}
        /**@{
           \name Trigonometric functions.
        */
        /**
           \brief Returns the sine of the given floating point number.
           
           \param i The operand.
           \param rc The context whose precision to use for the result.
           \return The result.
           \sa \ref arithctxts_trig
        */
        inline expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                    expressions::SinOp_Context> sin(const Real & i, const RealContext & rc)
        { return expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                      expressions::SinOp_Context>(expressions::Wrapper<RealContext>(i), rc); }
        /**
           \brief Returns the cosine of the given floating point number.
           
           \param i The operand.
           \param rc The context whose precision to use for the result.
           \return The result.
           \sa \ref arithctxts_trig
        */
        inline expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                    expressions::CosOp_Context> cos(const Real & i, const RealContext & rc)
        { return expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                      expressions::CosOp_Context>(expressions::Wrapper<RealContext>(i), rc); }
        /**
           \brief Returns the tangent of the given floating point number.
           
           \param i The operand.
           \param rc The context whose precision to use for the result.
           \return The result.
           \sa \ref arithctxts_trig
        */
        inline expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                    expressions::TanOp_Context> tan(const Real & i, const RealContext & rc)
        { return expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                      expressions::TanOp_Context>(expressions::Wrapper<RealContext>(i), rc); }
        /**
           \brief Returns the arcsine of the given floating point number.
           
           \param i The operand.
           \param rc The context whose precision to use for the result.
           \return The result.
           \sa \ref arithctxts_trig
        */
        inline expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                    expressions::ASinOp_Context> asin(const Real & i, const RealContext & rc)
        { return expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                      expressions::ASinOp_Context>(expressions::Wrapper<RealContext>(i), rc); }
        /**
           \brief Returns the arccosine of the given floating point number.
           
           \param i The operand.
           \param rc The context whose precision to use for the result.
           \return The result.
           \sa \ref arithctxts_trig
        */
        inline expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                    expressions::ACosOp_Context> acos(const Real & i, const RealContext & rc)
        { return expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                      expressions::ACosOp_Context>(expressions::Wrapper<RealContext>(i), rc); }
        /**
           \brief Returns the arctangent of the given floating point number.
           
           \param i The operand.
           \param rc The context whose precision to use for the result.
           \return The result.
           \sa \ref arithctxts_trig
        */
        inline expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                    expressions::ATanOp_Context> atan(const Real & i, const RealContext & rc)
        { return expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                      expressions::ATanOp_Context>(expressions::Wrapper<RealContext>(i), rc); }
        /**
           \brief Returns the arctangent of \f$\tfrac{y}{x}\f$, using the signs of `x` and `y`
           to determine the quadrant.
           
           \param y The numerator.
           \param x The denominator.
           \param rc The context whose precision to use for the result.
           \return The result.
           \sa \ref arithctxts_trig
        */
        inline expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                              expressions::Wrapper<RealContext> >,
                                                    expressions::ATan2Op_Context> atan2(const Real & y, const Real & x, const RealContext & rc)
        { return expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                                expressions::Wrapper<RealContext> >,
                                                      expressions::ATan2Op_Context>(std::make_pair(expressions::Wrapper<RealContext>(y),
                                                                                             expressions::Wrapper<RealContext>(x)), rc); }
        ///@}
        /**@{
           \name Exponential and logarithmic functions.
        */
        /**
           \brief Returns the exponential function evaluated at the given floating point number.
           
           \param i The operand.
           \param rc The context whose precision to use for the result.
           \param rc The context whose precision to use for the result.
           \return The result.
           \sa \ref arithctxts_specfns
        */
        inline expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                    expressions::ExpOp_Context> exp(const Real & i, const RealContext & rc)
        { return expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                      expressions::ExpOp_Context>(expressions::Wrapper<RealContext>(i), rc); }
        /**
           \brief Returns the natural logarithm evaluated at the given floating point number.
           
           \param i The operand.
           \param rc The context whose precision to use for the result.
           \return The result.
           \sa \ref arithctxts_specfns
        */
        inline expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                    expressions::LogOp_Context> log(const Real & i, const RealContext & rc)
        { return expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                      expressions::LogOp_Context>(expressions::Wrapper<RealContext>(i), rc); }
        /**
           \brief Returns the base-2 logarithm evaluated at the given floating point number.
           
           \param i The operand.
           \param rc The context whose precision to use for the result.
           \return The result.
           \sa \ref arithctxts_specfns
        */
        inline expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                    expressions::Log2Op_Context> log2(const Real & i, const RealContext & rc)
        { return expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                      expressions::Log2Op_Context>(expressions::Wrapper<RealContext>(i), rc); }
        /**
           \brief Returns the base-10 logarithm evaluated at the given floating point number.
           
           \param i The operand.
           \param rc The context whose precision to use for the result.
           \return The result.
           \sa \ref arithctxts_specfns
        */
        inline expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                    expressions::Log10Op_Context> log10(const Real & i, const RealContext & rc)
        { return expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                      expressions::Log10Op_Context>(expressions::Wrapper<RealContext>(i), rc); }
        /**
           \brief Returns the square root of the given floating point number.
           
           \param i The operand.
           \param rc The context whose precision to use for the result.
           \return The result.
           \sa \ref arithctxts_sqrtfp
        */
        inline expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                    expressions::SqrtOp_Context> sqrt(const Real & i, const RealContext & rc)
        { return expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                      expressions::SqrtOp_Context>(expressions::Wrapper<RealContext>(i), rc); }
        ///@}
        /**@{
           \name Special functions.
        */
        /**
           \brief Returns the Gamma function evaluated at the given floating point number.
           
           \param i The operand.
           \param rc The context whose precision to use for the result.
           \return The result.
           \sa \ref arithctxts_specfns
        */
        inline expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                    expressions::GammaOp_Context> gamma(const Real & i, const RealContext & rc)
        { return expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                      expressions::GammaOp_Context>(expressions::Wrapper<RealContext>(i), rc); }
        /**
           \brief Returns the natural logarithm of the absolute value of the Gamma function
                  evaluated at the given floating point number.
           
           \param i The operand.
           \param rc The context whose precision to use for the result.
           \return The result.
           \sa \ref arithctxts_specfns
        */
        inline expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                    expressions::LGammaOp_Context> lgamma(const Real & i, const RealContext & rc)
        { return expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                      expressions::LGammaOp_Context>(expressions::Wrapper<RealContext>(i), rc); }
        /**
           \brief Returns the natural logarithm of the absolute value of the Gamma function
                  evaluated at the given floating point number, as well as the sign of the Gamma
                  function.
           
           \param sign The integer variable to store the sign of \f$\Gamma(i)\f$ in.
           \param i The operand.
           \param rc The context whose precision to use for the result.
           \return The result.
           \sa \ref arithctxts_specfns
        */
        inline expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                    expressions::LGamma2Op_Context> lgamma(int & sign, const Real & i, const RealContext & rc)
        { return expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                      expressions::LGamma2Op_Context>(expressions::Wrapper<RealContext>(i),
                                                                                      expressions::LGamma2Op_Context<RealContext, expressions::Wrapper<RealContext> >(sign, rc)); }
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
        inline expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                    expressions::PowerCOp_Context<signed long>::impl> power(const Real & a, signed long b, const RealContext & rc)
        { return expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                      expressions::PowerCOp_Context<signed long>::impl>(expressions::Wrapper<RealContext>(a),
                                                                                                        expressions::PowerCOp_Context<signed long>::impl<RealContext, expressions::Wrapper<RealContext> >(b, rc)); }
        /**
           \brief Returns `a` raised to the power of `b`.
           
           \param a The base.
           \param b The exponent.
           \param rc The context whose precision to use for the result.
           \return The result.
        */
        inline expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                    expressions::PowerCOp_Context<unsigned long>::impl> power(const Real & a, unsigned long b, const RealContext & rc)
        { return expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                      expressions::PowerCOp_Context<unsigned long>::impl>(expressions::Wrapper<RealContext>(a),
                                                                                                          expressions::PowerCOp_Context<unsigned long>::impl<RealContext, expressions::Wrapper<RealContext> >(b, rc)); }
        /**
           \brief Returns `a` raised to the power of `b`.
           
           \param a The base.
           \param b The exponent.
           \param rc The context whose precision to use for the result.
           \return The result.
        */
        inline expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                              expressions::Wrapper<IntegerContext> >,
                                                    expressions::PowerOp_Context> power(const Real & a, const Integer & b, const RealContext & rc)
        { return expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                                expressions::Wrapper<IntegerContext> >,
                                                      expressions::PowerOp_Context>(std::make_pair(expressions::Wrapper<RealContext>(a),
                                                                                                   expressions::Wrapper<IntegerContext>(b)), rc); }
        /**
           \brief Returns `a` raised to the power of `b`.
           
           \param a The base.
           \param b The exponent.
           \param rc The context whose precision to use for the result.
           \return The result.
           \sa \ref arithctxts_sqrtfp
        */
        inline expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                              expressions::Wrapper<RealContext> >,
                                                    expressions::PowerOp_Context> power(const Real & a, const Real & b, const RealContext & rc)
        { return expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                                expressions::Wrapper<RealContext> >,
                                                      expressions::PowerOp_Context>(std::make_pair(expressions::Wrapper<RealContext>(a),
                                                                                                   expressions::Wrapper<RealContext>(b)), rc); }
        ///@}
        
        template<class A2, template<typename, typename> class O2>
        inline expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                              expressions::Expression<RealContext, A2, O2> >,
                                                    expressions::ATan2Op_Context> atan2(const Real & a, const expressions::Expression<RealContext, A2, O2> & b,
                                                                                  const RealContext & rc)
        { return expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                                expressions::Expression<RealContext, A2, O2> >,
                                                      expressions::ATan2Op_Context>(std::make_pair(expressions::Wrapper<RealContext>(a), b), rc); }
        template<class A2, template<typename, typename> class O2>
        inline expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                              expressions::Expression<RealContext, A2, O2> >,
                                                    expressions::PowerOp_Context> power(const Real & a, const expressions::Expression<RealContext, A2, O2> & b, const RealContext & rc)
        { return expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,
                                                                 expressions::Expression<RealContext, A2, O2> >,
                                                      expressions::PowerOp_Context>(std::make_pair(expressions::Wrapper<RealContext>(a), b), rc); }
        
        namespace implementation
        {
            namespace Real_impl
            {
                template<class A, template<typename, typename> class O>
                inline expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                                            expressions::AbsOp_Context> abs(const expressions::Expression<RealContext, A, O> & b,
                                                                                            const RealContext & rc)
                { return expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                                              expressions::AbsOp_Context>(b, rc); }
                template<class A, template<typename, typename> class O>
                inline expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                                            expressions::SquareOp_Context> square(const expressions::Expression<RealContext, A, O> & b,
                                                                                                  const RealContext & rc)
                { return expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                                              expressions::SquareOp_Context>(b, rc); }
                template<class A, template<typename, typename> class O>
                inline expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                                            expressions::SinOp_Context> sin(const expressions::Expression<RealContext, A, O> & b,
                                                                                            const RealContext & rc)
                { return expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                                              expressions::SinOp_Context>(b, rc); }
                template<class A, template<typename, typename> class O>
                inline expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                                            expressions::CosOp_Context> cos(const expressions::Expression<RealContext, A, O> & b,
                                                                                            const RealContext & rc)
                { return expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                                              expressions::CosOp_Context>(b, rc); }
                template<class A, template<typename, typename> class O>
                inline expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                                            expressions::TanOp_Context> tan(const expressions::Expression<RealContext, A, O> & b,
                                                                                            const RealContext & rc)
                { return expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                                              expressions::TanOp_Context>(b, rc); }
                template<class A, template<typename, typename> class O>
                inline expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                                            expressions::ASinOp_Context> asin(const expressions::Expression<RealContext, A, O> & b,
                                                                                              const RealContext & rc)
                { return expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                                              expressions::ASinOp_Context>(b, rc); }
                template<class A, template<typename, typename> class O>
                inline expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                                            expressions::ACosOp_Context> acos(const expressions::Expression<RealContext, A, O> & b,
                                                                                              const RealContext & rc)
                { return expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                                              expressions::ACosOp_Context>(b, rc); }
                template<class A, template<typename, typename> class O>
                inline expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                                            expressions::ATanOp_Context> atan(const expressions::Expression<RealContext, A, O> & b,
                                                                                              const RealContext & rc)
                { return expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                                              expressions::ATanOp_Context>(b, rc); }
                template<class A, template<typename, typename> class O, class A2, template<typename, typename> class O2>
                inline expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A, O>,
                                                                      expressions::Expression<RealContext, A2, O2> >,
                                                            expressions::ATan2Op_Context> atan2(const expressions::Expression<RealContext, A, O> & a,
                                                                                                const expressions::Expression<RealContext, A2, O2> & b,
                                                                                                const RealContext & rc)
                { return expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A, O>,
                                                                        expressions::Expression<RealContext, A2, O2> >,
                                                              expressions::ATan2Op_Context>(std::make_pair(a, b), rc); }
                template<class A, template<typename, typename> class O>
                inline expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A, O>,
                                                                      expressions::Wrapper<RealContext> >,
                                                            expressions::ATan2Op_Context> atan2(const expressions::Expression<RealContext, A, O> & a,
                                                                                                const arithmetic::Real & b, const RealContext & rc)
                { return expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A, O>,
                                                                        expressions::Wrapper<RealContext> >,
                                                              expressions::ATan2Op_Context>(std::make_pair(a, expressions::Wrapper<RealContext>(b)), rc); }
                template<class A, template<typename, typename> class O>
                inline expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                                            expressions::ExpOp_Context> exp(const expressions::Expression<RealContext, A, O> & b,
                                                                                            const RealContext & rc)
                { return expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                                              expressions::ExpOp_Context>(b, rc); }
                template<class A, template<typename, typename> class O>
                inline expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                                            expressions::LogOp_Context> log(const expressions::Expression<RealContext, A, O> & b,
                                                                                            const RealContext & rc)
                { return expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                                              expressions::LogOp_Context>(b, rc); }
                template<class A, template<typename, typename> class O>
                inline expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                                            expressions::Log2Op_Context> log2(const expressions::Expression<RealContext, A, O> & b,
                                                                                              const RealContext & rc)
                { return expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                                              expressions::Log2Op_Context>(b, rc); }
                template<class A, template<typename, typename> class O>
                inline expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                                            expressions::Log10Op_Context> log10(const expressions::Expression<RealContext, A, O> & b,
                                                                                                const RealContext & rc)
                { return expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                                              expressions::Log10Op_Context>(b, rc); }
                template<class A, template<typename, typename> class O>
                inline expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                                            expressions::SqrtOp_Context> sqrt(const expressions::Expression<RealContext, A, O> & b,
                                                                                              const RealContext & rc)
                { return expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                                              expressions::SqrtOp_Context>(b, rc); }
                template<class A, template<typename, typename> class O>
                inline expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                                            expressions::GammaOp_Context> gamma(const expressions::Expression<RealContext, A, O> & b,
                                                                                                const RealContext & rc)
                { return expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                                              expressions::GammaOp_Context>(b, rc); }
                template<class A, template<typename, typename> class O>
                inline expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                                            expressions::LGammaOp_Context> lgamma(const expressions::Expression<RealContext, A, O> & b,
                                                                                                  const RealContext & rc)
                { return expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                                              expressions::LGammaOp_Context>(b, rc); }
                template<class A, template<typename, typename> class O>
                inline expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                                            expressions::PowerCOp_Context<signed long>::impl> power(const expressions::Expression<RealContext, A, O> & a,
                                                                                                                    signed long b, const RealContext & rc)
                { return expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                                              expressions::PowerCOp_Context<signed long>::impl>(a, expressions::PowerCOp_Context<signed long>::impl<RealContext, expressions::Expression<RealContext, A, O> >(b, rc)); }
                template<class A, template<typename, typename> class O>
                inline expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                                            expressions::PowerCOp_Context<unsigned long>::impl> power(const expressions::Expression<RealContext, A, O> & a,
                                                                                                                      unsigned long b, const RealContext & rc)
                { return expressions::Expression<RealContext, expressions::Expression<RealContext, A, O>,
                                                              expressions::PowerCOp_Context<unsigned long>::impl>(a, expressions::PowerCOp_Context<unsigned long>::impl<RealContext, expressions::Expression<RealContext, A, O> >(b, rc)); }
                template<class A, template<typename, typename> class O>
                inline expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A, O>,
                                                                      expressions::Wrapper<IntegerContext> >,
                                                            expressions::PowerOp_Context> power(const expressions::Expression<RealContext, A, O> & a,
                                                                                                const Integer & b, const RealContext & rc)
                { return expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A, O>,
                                                                        expressions::Wrapper<IntegerContext> >,
                                                              expressions::PowerOp_Context>(std::make_pair(a, expressions::Wrapper<IntegerContext>(b)), rc); }
                template<class A, template<typename, typename> class O, class A2, template<typename, typename> class O2>
                inline expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A, O>,
                                                                      expressions::Expression<RealContext, A2, O2> >,
                                                            expressions::PowerOp_Context> power(const expressions::Expression<RealContext, A, O> & a,
                                                                                                const expressions::Expression<RealContext, A2, O2> & b, const RealContext & rc)
                { return expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A, O>,
                                                                        expressions::Expression<RealContext, A2, O2> >,
                                                              expressions::PowerOp_Context>(std::make_pair(a, b), rc); }
                template<class A, template<typename, typename> class O>
                inline expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A, O>,
                                                                      expressions::Wrapper<RealContext> >,
                                                            expressions::PowerOp_Context> power(const expressions::Expression<RealContext, A, O> & a,
                                                                                                const arithmetic::Real & b, const RealContext & rc)
                { return expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, A, O>,
                                                                        expressions::Wrapper<RealContext> >,
                                                              expressions::PowerOp_Context>(std::make_pair(a, expressions::Wrapper<RealContext>(b)), rc); }
            }
        }
    }
}

#endif

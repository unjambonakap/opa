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

#ifndef PLLL_INCLUDE_GUARD__RATIONAL_CONV_HPP
#define PLLL_INCLUDE_GUARD__RATIONAL_CONV_HPP

/**
   \file
   \brief Conversion definitions for rational numbers.
   
   This header contains abstract templates used to implement most conversions involving rational
   numbers.
*/
namespace plll
{
    namespace arithmetic
    {
        namespace implementation
        {
            // Implementation for RationalContext
            
            template<> class conversion_impl<float, RationalContext>
            {
            public:
                typedef expressions::Expression<RationalContext, expressions::ConversionWrapper<float>, expressions::ConvertOp_Context> RetVal;
                typedef Rational RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(Rational & d, float v, const RationalContext & c)
                {
                    d.initFromDouble(v);
                }
                
                static RetVal convert(float v, const RationalContext & c)
                {
                    return RetVal(v, c);
                }
                
                static void convert_frac(Rational & d, float v1, float v2, const RationalContext & c)
                {
                    convert(d, v1, c);
                    Rational tmp;
                    convert(tmp, v2, c);
                    d /= tmp;
                }
                
                static RetVal_Frac convert_frac(float v1, float v2, const RationalContext & c)
                {
                    Rational res, tmp;
                    convert(res, v1, c);
                    convert(tmp, v2, c);
                    return res /= tmp;
                }
            };
            
            template<> class conversion_impl<double, RationalContext>
            {
            public:
                typedef expressions::Expression<RationalContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> RetVal;
                typedef Rational RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(Rational & d, double v, const RationalContext & c)
                {
                    d.initFromDouble(v);
                }
                
                static RetVal convert(double v, const RationalContext & c)
                {
                    return RetVal(v, c);
                }
                
                static void convert_frac(Rational & d, double v1, double v2, const RationalContext & c)
                {
                    convert(d, v1, c);
                    Rational tmp;
                    convert(tmp, v2, c);
                    d /= tmp;
                }
                
                static RetVal_Frac convert_frac(double v1, double v2, const RationalContext & c)
                {
                    Rational res, tmp;
                    convert(res, v1, c);
                    convert(tmp, v2, c);
                    return res /= tmp;
                }
            };
            
            template<> class conversion_impl<long double, RationalContext>
            {
            public:
                typedef expressions::Expression<RationalContext, expressions::ConversionWrapper<long double>, expressions::ConvertOp_Context> RetVal;
                typedef Rational RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(Rational & d, long double v, const RationalContext & c)
                {
                    d.initFromDouble(v);
                }
                
                static RetVal convert(long double v, const RationalContext & c)
                {
                    return RetVal(v, c);
                }
                
                static void convert_frac(Rational & d, long double v1, long double v2, const RationalContext & c)
                {
                    convert(d, v1, c);
                    Rational tmp;
                    convert(tmp, v2, c);
                    d /= tmp;
                }
                
                static RetVal_Frac convert_frac(long double v1, long double v2, const RationalContext & c)
                {
                    Rational res, tmp;
                    convert(res, v1, c);
                    convert(tmp, v2, c);
                    return res /= tmp;
                }
            };
            
            template<> class conversion_impl<signed int, RationalContext>
            {
            private:
                typedef expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> Expr;
                
            public:
                typedef expressions::Expression<RationalContext, Expr, expressions::ConvertOp_Context> RetVal;
                typedef Rational RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(Rational & d, signed int v, const RationalContext & c)
                {
                    IntegerContext ic;
                    arithmetic::convert(d.d_num, v, ic);
                    setOne(d.d_denom);
                }
                
                static RetVal convert(signed int v, const RationalContext & c)
                {
                    return RetVal(Expr(v, g_intcontext), c);
                }
                
                static void convert_frac(Rational & d, signed int v1, signed int v2, const RationalContext & c)
                {
                    arithmetic::convert(d.d_num, v1, IntegerContext());
                    arithmetic::convert(d.d_denom, v2, IntegerContext());
                    d.normalize();
                    // Make sign correct
                    if (sign(d.d_denom) < 0)
                    {
                        neg(d.d_denom, d.d_denom);
                        neg(d.d_num, d.d_num);
                    }
                }
                
                static RetVal_Frac convert_frac(signed int v1, signed int v2, const RationalContext & c)
                {
                    Rational r;
                    convert_frac(r, v1, v2, c);
                    return r;
                }
            };
            
            template<> class conversion_impl<unsigned int, RationalContext>
            {
            private:
                typedef expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> Expr;
                
            public:
                typedef expressions::Expression<RationalContext, Expr, expressions::ConvertOp_Context> RetVal;
                typedef Rational RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(Rational & d, unsigned int v, const RationalContext & c)
                {
                    IntegerContext ic;
                    arithmetic::convert(d.d_num, v, ic);
                    setOne(d.d_denom);
                }
                
                static RetVal convert(unsigned int v, const RationalContext & c)
                {
                    return RetVal(Expr(v, g_intcontext), c);
                }
                
                static void convert_frac(Rational & d, unsigned int v1, unsigned int v2, const RationalContext & c)
                {
                    arithmetic::convert(d.d_num, v1, IntegerContext());
                    arithmetic::convert(d.d_denom, v2, IntegerContext());
                    d.normalize();
                    // Make sign correct
                    if (sign(d.d_denom) < 0)
                    {
                        neg(d.d_denom, d.d_denom);
                        neg(d.d_num, d.d_num);
                    }
                }
                
                static RetVal_Frac convert_frac(unsigned int v1, unsigned int v2, const RationalContext & c)
                {
                    Rational r;
                    convert_frac(r, v1, v2, c);
                    return r;
                }
            };
            
            template<> class conversion_impl<signed long, RationalContext>
            {
            private:
                typedef expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> Expr;
                
            public:
                typedef expressions::Expression<RationalContext, Expr, expressions::ConvertOp_Context> RetVal;
                typedef Rational RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(Rational & d, signed long v, const RationalContext & c)
                {
                    IntegerContext ic;
                    arithmetic::convert(d.d_num, v, ic);
                    setOne(d.d_denom);
                }
                
                static RetVal convert(signed long v, const RationalContext & c)
                {
                    return RetVal(Expr(v, g_intcontext), c);
                }
                
                static void convert_frac(Rational & d, signed long v1, signed long v2, const RationalContext & c)
                {
                    arithmetic::convert(d.d_num, v1, IntegerContext());
                    arithmetic::convert(d.d_denom, v2, IntegerContext());
                    d.normalize();
                    // Make sign correct
                    if (sign(d.d_denom) < 0)
                    {
                        neg(d.d_denom, d.d_denom);
                        neg(d.d_num, d.d_num);
                    }
                }
                
                static RetVal_Frac convert_frac(signed long v1, signed long v2, const RationalContext & c)
                {
                    Rational r;
                    convert_frac(r, v1, v2, c);
                    return r;
                }
            };
            
            template<> class conversion_impl<unsigned long, RationalContext>
            {
            private:
                typedef expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> Expr;
                
            public:
                typedef expressions::Expression<RationalContext, Expr, expressions::ConvertOp_Context> RetVal;
                typedef Rational RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(Rational & d, unsigned long v, const RationalContext & c)
                {
                    IntegerContext ic;
                    arithmetic::convert(d.d_num, v, ic);
                    setOne(d.d_denom);
                }
                
                static RetVal convert(unsigned long v, const RationalContext & c)
                {
                    return RetVal(Expr(v, g_intcontext), c);
                }
                
                static void convert_frac(Rational & d, unsigned long v1, unsigned long v2, const RationalContext & c)
                {
                    arithmetic::convert(d.d_num, v1, IntegerContext());
                    arithmetic::convert(d.d_denom, v2, IntegerContext());
                    d.normalize();
                    // Make sign correct
                    if (sign(d.d_denom) < 0)
                    {
                        neg(d.d_denom, d.d_denom);
                        neg(d.d_num, d.d_num);
                    }
                }
                
                static RetVal_Frac convert_frac(unsigned long v1, unsigned long v2, const RationalContext & c)
                {
                    Rational r;
                    convert_frac(r, v1, v2, c);
                    return r;
                }
            };
            
            template<> class conversion_impl<long long, RationalContext>
            {
            private:
                typedef expressions::Expression<IntegerContext, expressions::ConversionWrapper<long long>, expressions::ConvertOp_Context> Expr;
                
            public:
                typedef expressions::Expression<RationalContext, Expr, expressions::ConvertOp_Context> RetVal;
                typedef Rational RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(Rational & d, long long v, const RationalContext & c)
                {
                    IntegerContext ic;
                    arithmetic::convert(d.d_num, v, ic);
                    setOne(d.d_denom);
                }
                
                static RetVal convert(long long v, const RationalContext & c)
                {
                    return RetVal(Expr(v, g_intcontext), c);
                }
                
                static void convert_frac(Rational & d, long long v1, long long v2, const RationalContext & c)
                {
                    arithmetic::convert(d.d_num, v1, IntegerContext());
                    arithmetic::convert(d.d_denom, v2, IntegerContext());
                    d.normalize();
                    // Make sign correct
                    if (sign(d.d_denom) < 0)
                    {
                        neg(d.d_denom, d.d_denom);
                        neg(d.d_num, d.d_num);
                    }
                }
                
                static RetVal_Frac convert_frac(long long v1, long long v2, const RationalContext & c)
                {
                    Rational r;
                    convert_frac(r, v1, v2, c);
                    return r;
                }
            };
            
            template<> class conversion_impl<Rational, RealContext>
            {
            public:
                typedef Real RetVal;
                typedef Real RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(Real & res, const Rational & r, const RealContext & rc)
                {
                    arithmetic::convert(res, r.numerator(), rc);
                    Real tmp(r.denominator(), rc);
                    res /= tmp;
                }
                
                static Real convert(const Rational & r, const RealContext & rc)
                {
                    Real res(r.numerator(), rc), tmp(r.denominator(), rc);
                    res /= tmp;
                    return res;
                }
                
                static void convert_frac(Real & res, const Rational & r1, const Rational & r2, const RealContext & rc)
                {
                    arithmetic::convert(res, r1.numerator() * r2.denominator(), rc);
                    Real tmp(r1.denominator() * r2.numerator(), rc);
                    res /= tmp;
                }
                
                static Real convert_frac(const Rational & r1, const Rational & r2, const RealContext & rc)
                {
                    Real res(r1.numerator() * r2.denominator(), rc), tmp(r1.denominator() * r2.numerator(), rc);
                    return res /= tmp;
                }
            };
            
            template<> class conversion_impl<Real, RationalContext>
            {
            public:
                typedef Rational RetVal;
                typedef Rational RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(Rational & res, const Real & r, const RationalContext & rc)
                {
                    setOne(res.d_denom);
                    mpfr_exp_t exp = mpfr_get_z_2exp(res.d_num.d_value, r.d_value);
                    if (exp > 0)
                        res.d_num <<= exp;
                    else if (exp < 0)
                        res.d_denom <<= -exp;
                }
                
                static Rational convert(const Real & r, const RationalContext & rc)
                {
                    Rational res;
                    convert(res, r, rc);
                    return res;
                }
                
                static void convert_frac(Rational & res, const Real & r1, const Real & r2, const RationalContext & rc)
                {
                    setOne(res.d_denom);
                    mpfr_exp_t exp1 = mpfr_get_z_2exp(res.d_num.d_value, r1.d_value);
                    mpfr_exp_t exp2 = mpfr_get_z_2exp(res.d_denom.d_value, r2.d_value);
                    if (exp1 > exp2)
                        res.d_num <<= exp1 - exp2;
                    else if (exp1 < exp2)
                        res.d_denom <<= exp2 - exp1;
                    res.normalize();
                    // Make sign correct
                    if (sign(res.d_denom) < 0)
                    {
                        neg(res.d_denom, res.d_denom);
                        neg(res.d_num, res.d_num);
                    }
                }
                
                static Rational convert_frac(const Real & r1, const Real & r2, const RationalContext & rc)
                {
                    Rational res;
                    convert_frac(res, r1, r2, rc);
                    return res;
                }
            };
            
            template<> class conversion_impl<Rational, IntegerContext>
            {
#ifndef PLLL_INTERNAL_NO_TEMPLATE_FRIENDS
                friend class nativeconversion_impl<Rational>;
                
            private:
#else
            public:
#endif
                static inline void toint(Integer & r, const Rational & v)
                {
                    r = v.d_num / v.d_denom;
                }
                
                static inline void toint_floor(Integer & r, const Rational & v)
                {
                    floorDiv(r, v.d_num, v.d_denom);
                }
                
                static inline void toint_round(Integer & r, const Rational & v)
                {
                    Integer rr;
                    euclideanDivisionPos(r, rr, v.d_num, v.d_denom);
                    rr <<= 1;
                    if (rr >= v.d_denom)
                        ++r;
                }
                
                static inline void toint_round(Integer & r, const Rational & v, bool & up)
                {
                    Integer rr;
                    euclideanDivisionPos(r, rr, v.d_num, v.d_denom);
                    rr <<= 1;
                    if ((up = (rr >= v.d_denom)) == true)
                        ++r;
                }
                
                static inline void toint_ceil(Integer & r, const Rational & v)
                {
                    ceilDiv(r, v.d_num, v.d_denom);
                }
                
            public:
                typedef typename expressions::Expression<IntegerContext, expressions::Expression<RationalContext, expressions::Wrapper<RationalContext>,
                                                                                                 expressions::NoneOp>, expressions::ConvertFloorOp_Context> RetVal;
                typedef void RetVal_Frac;
                typedef typename expressions::Expression<IntegerContext, expressions::Expression<RationalContext, expressions::Wrapper<RationalContext>,
                                                                                                 expressions::NoneOp>, expressions::ConvertFloorOp_Context> RetVal_Floor;
                typedef typename expressions::Expression<IntegerContext, expressions::Expression<RationalContext, expressions::Wrapper<RationalContext>,
                                                                                                 expressions::NoneOp>, expressions::ConvertRoundOp_Context> RetVal_Round;
                typedef typename expressions::Expression<IntegerContext, expressions::Expression<RationalContext, expressions::Wrapper<RationalContext>,
                                                                                                 expressions::NoneOp>, expressions::ConvertRound2Op_Context> RetVal_Round2;
                typedef typename expressions::Expression<IntegerContext, expressions::Expression<RationalContext, expressions::Wrapper<RationalContext>,
                                                                                                 expressions::NoneOp>, expressions::ConvertCeilOp_Context> RetVal_Ceil;
                
                static void convert(Integer & res, const Rational & v, const IntegerContext &)
                {
                    toint(res, v);
                }
                
                static RetVal convert(const Rational & v, const IntegerContext & c)
                {
                    return RetVal(expressions::make_expression<RationalContext>(v), c);
                }
                
                static void floor(Integer & res, const Rational & v, const IntegerContext &)
                {
                    toint_floor(res, v);
                }
                
                static RetVal_Floor floor(const Rational & v, const IntegerContext & c)
                {
                    return RetVal_Floor(expressions::make_expression<RationalContext>(v), c);
                }
                
                static void round(Integer & res, const Rational & v, const IntegerContext &)
                {
                    toint_round(res, v);
                }
                
                static RetVal_Round round(const Rational & v, const IntegerContext & c)
                {
                    return RetVal_Round(expressions::make_expression<RationalContext>(v), c);
                }
                
                static void round(Integer & res, const Rational & v, bool & up, const IntegerContext &)
                {
                    toint_round(res, v, up);
                }
                
                static RetVal_Round2 round(const Rational & v, bool & up, const IntegerContext & c)
                {
                    return RetVal_Round2(expressions::make_expression<RationalContext>(v),
                                         expressions::ConvertRound2Op_Context<IntegerContext,
                                                                              expressions::Expression<RationalContext, expressions::Wrapper<RationalContext>,
                                                                                                      expressions::NoneOp> >(up, c));
                }
                
                static void ceil(Integer & res, const Rational & v, const IntegerContext &)
                {
                    toint_ceil(res, v);
                }
                
                static RetVal_Ceil ceil(const Rational & v, const IntegerContext & c)
                {
                    return RetVal_Ceil(expressions::make_expression<RationalContext>(v), c);
                }
            };
            
            template<> class conversion_impl<Integer, RationalContext>
            {
            public:
                typedef expressions::Expression<RationalContext, expressions::Expression<IntegerContext, expressions::Wrapper<IntegerContext>,
                                                                                         expressions::NoneOp>, expressions::ConvertOp_Context> RetVal;
                typedef Rational RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(Rational & res, const Integer & v, const RationalContext & rc)
                {
                    res.d_num = v;
                    setOne(res.d_denom);
                }
                
                static RetVal convert(const arithmetic::Integer & v, const RationalContext & c)
                {
                    return RetVal(expressions::make_expression<IntegerContext>(v), c);
                }
                
                static void convert_frac(Rational & res, const Integer & v1, const Integer & v2, const RationalContext & rc)
                {
                    if (&res.d_num == &v2)
                    {
                        if (&res.d_denom == &v1)
                        {
                            swap(res.d_num, res.d_denom);
                        }
                        else
                        {
                            res.d_denom = v2; // v2 == res.d_num
                            res.d_num = v1;
                            res.normalize();
                        }
                    }
                    else
                    {
                        res.d_num = v1;
                        res.d_denom = v2;
                        res.normalize();
                    }
                    // Make sign correct
                    if (sign(res.d_denom) < 0)
                    {
                        neg(res.d_denom, res.d_denom);
                        neg(res.d_num, res.d_num);
                    }
                }
                
                static RetVal_Frac convert_frac(const Integer & v1, const Integer & v2, const RationalContext & rc)
                {
                    return Rational(v1, v2);
                }
            };
            
            template<class Data, template<typename, typename> class Op>
            class conversion_impl<expressions::Expression<IntegerContext, Data, Op>, RationalContext>
            {
            private:
                typedef expressions::Expression<IntegerContext, Data, Op> Expr;
                
            public:
                typedef expressions::Expression<RationalContext, Expr, expressions::ConvertOp_Context> RetVal;
                typedef Rational RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(Rational & res, const Expr & v, const RationalContext & rc)
                {
                    res.d_num = v;
                    setOne(res.d_denom);
                }
                
                static RetVal convert(const Expr & v, const RationalContext & c)
                {
                    return RetVal(v);
                }
                
                static void convert_frac(Rational & res, const Integer & v1, const Integer & v2, const RationalContext & rc)
                {
                    conversion_impl<arithmetic::Integer, RationalContext>::convert_frac(res, v1, v2, rc);
                }
                
                static RetVal_Frac convert_frac(const Integer & v1, const Integer & v2, const RationalContext & rc)
                {
                    return conversion_impl<arithmetic::Integer, RationalContext>::convert_frac(v1, v2, rc);
                }
            };
            
            template<> class nativeconversion_impl<Rational>
            {
            public:
                static int toInt(const Rational & v)
                {
                    Integer r;
                    conversion_impl<Rational, IntegerContext>::toint(r, v);
                    return convert<int>(r);
                }
                
                static int toInt_Floor(const Rational & v)
                {
                    Integer r;
                    conversion_impl<Rational, IntegerContext>::toint_floor(r, v);
                    return convert<int>(r);
                }
                
                static int toInt_Round(const Rational & v)
                {
                    Integer r;
                    conversion_impl<Rational, IntegerContext>::toint_round(r, v);
                    return convert<int>(r);
                }
                
                static int toInt_Round(const Rational & v, bool & up)
                {
                    Integer r;
                    conversion_impl<Rational, IntegerContext>::toint_round(r, v, up);
                    return convert<int>(r);
                }
                
                static int toInt_Ceil(const Rational & v)
                {
                    Integer r;
                    conversion_impl<Rational, IntegerContext>::toint_ceil(r, v);
                    return convert<int>(r);
                }
                
                static unsigned int toUInt(const Rational & v)
                {
                    Integer r;
                    conversion_impl<Rational, IntegerContext>::toint(r, v);
                    return convert<unsigned int>(r);
                }
                
                static unsigned int toUInt_Floor(const Rational & v)
                {
                    Integer r;
                    conversion_impl<Rational, IntegerContext>::toint_floor(r, v);
                    return convert<unsigned int>(r);
                }
                
                static unsigned int toUInt_Round(const Rational & v)
                {
                    Integer r;
                    conversion_impl<Rational, IntegerContext>::toint_round(r, v);
                    return convert<unsigned int>(r);
                }
                
                static unsigned int toUInt_Round(const Rational & v, bool & up)
                {
                    Integer r;
                    conversion_impl<Rational, IntegerContext>::toint_round(r, v, up);
                    return convert<unsigned int>(r);
                }
                
                static unsigned int toUInt_Ceil(const Rational & v)
                {
                    Integer r;
                    conversion_impl<Rational, IntegerContext>::toint_ceil(r, v);
                    return convert<unsigned int>(r);
                }
                
                static long toLong(const Rational & v)
                {
                    Integer r;
                    conversion_impl<Rational, IntegerContext>::toint(r, v);
                    return convert<long>(r);
                }
                
                static long toLong_Floor(const Rational & v)
                {
                    Integer r;
                    conversion_impl<Rational, IntegerContext>::toint_floor(r, v);
                    return convert<long>(r);
                }
                
                static long toLong_Round(const Rational & v)
                {
                    Integer r;
                    conversion_impl<Rational, IntegerContext>::toint_round(r, v);
                    return convert<long>(r);
                }
                
                static long toLong_Round(const Rational & v, bool & up)
                {
                    Integer r;
                    conversion_impl<Rational, IntegerContext>::toint_round(r, v, up);
                    return convert<long>(r);
                }
                
                static long toLong_Ceil(const Rational & v)
                {
                    Integer r;
                    conversion_impl<Rational, IntegerContext>::toint_ceil(r, v);
                    return convert<long>(r);
                }
                
                static unsigned long toULong(const Rational & v)
                {
                    Integer r;
                    conversion_impl<Rational, IntegerContext>::toint(r, v);
                    return convert<unsigned long>(r);
                }
                
                static unsigned long toULong_Floor(const Rational & v)
                {
                    Integer r;
                    conversion_impl<Rational, IntegerContext>::toint_floor(r, v);
                    return convert<unsigned long>(r);
                }
                
                static unsigned long toULong_Round(const Rational & v)
                {
                    Integer r;
                    conversion_impl<Rational, IntegerContext>::toint_round(r, v);
                    return convert<unsigned long>(r);
                }
                
                static unsigned long toULong_Round(const Rational & v, bool & up)
                {
                    Integer r;
                    conversion_impl<Rational, IntegerContext>::toint_round(r, v, up);
                    return convert<unsigned long>(r);
                }
                
                static unsigned long toULong_Ceil(const Rational & v)
                {
                    Integer r;
                    conversion_impl<Rational, IntegerContext>::toint_ceil(r, v);
                    return convert<unsigned long>(r);
                }
                
                static long long toLongLong(const Rational & v)
                {
                    Integer r;
                    conversion_impl<Rational, IntegerContext>::toint(r, v);
                    return convert<long long>(r);
                }
                
                static long long toLongLong_Floor(const Rational & v)
                {
                    Integer r;
                    conversion_impl<Rational, IntegerContext>::toint_floor(r, v);
                    return convert<long long>(r);
                }
                
                static long long toLongLong_Round(const Rational & v)
                {
                    Integer r;
                    conversion_impl<Rational, IntegerContext>::toint_round(r, v);
                    return convert<long long>(r);
                }
                
                static long long toLongLong_Round(const Rational & v, bool & up)
                {
                    Integer r;
                    conversion_impl<Rational, IntegerContext>::toint_round(r, v, up);
                    return convert<long long>(r);
                }
                
                static long long toLongLong_Ceil(const Rational & v)
                {
                    Integer r;
                    conversion_impl<Rational, IntegerContext>::toint_ceil(r, v);
                    return convert<long long>(r);
                }
                
                static float toFloat(const Rational & v)
                {
                    return toDouble(v);
                }
                
                static double toDouble(const Rational & v)
                {
                    long n = approxLog2(v.d_num), d = approxLog2(v.d_denom);
                    if ((n < std::numeric_limits<double>::max_exponent) && (d < std::numeric_limits<double>::max_exponent))
                        return convert<double>(v.d_num) / convert<double>(v.d_denom);
                    else
                    {
                        long nr = n - std::numeric_limits<double>::digits - 4, dr = d - std::numeric_limits<double>::digits - 4;
                        if (nr < 0) nr = 0;
                        if (dr < 0) dr = 0;
                        Integer N = v.d_num >> nr, D = v.d_denom >> dr;
                        double r = convert<double>(N) / convert<double>(D);
                        return ldexp(r, nr - dr);
                    }
                }
                
                static long double toLongDouble(const Rational & v)
                {
                    long n = approxLog2(v.d_num), d = approxLog2(v.d_denom);
                    if ((n < std::numeric_limits<long double>::max_exponent) && (d < std::numeric_limits<long double>::max_exponent))
                        return convert<long double>(v.d_num) / convert<long double>(v.d_denom);
                    else
                    {
                        long nr = n - std::numeric_limits<long double>::digits - 4, dr = d - std::numeric_limits<long double>::digits - 4;
                        if (nr < 0) nr = 0;
                        if (dr < 0) dr = 0;
                        Integer N = v.d_num >> nr, D = v.d_denom >> dr;
                        long double r = convert<long double>(N) / convert<long double>(D);
                        return ldexp(r, nr - dr);
                    }
                }
            };
        }
        
        namespace traits
        {
            template<>
            struct type_traits<Rational>
            {
                enum
                {
                    is_number = true,
                    is_cputype = false,
                    is_realtype = true,
                    is_inttype = false,
                    is_string = false,
                    is_cpp_string = false,
                    is_c_string = false,
                    is_exact = true,
                    has_context = true,
                    is_native = false,
                    is_modulo = false,
                    has_infinity = false,
                    is_variable_precision = false,
                    has_squareroot = false,
                    has_full_power = false,
                    has_special_fns = false,
                    has_huge_exponent = true,
                    has_uniform_rng = false,
                    has_constants = false,
                    has_trigonometric = false
                };
                
                typedef RationalContext Context;
                typedef Rational PromoteType;
                typedef const Rational & ConstReferenceType;
            };
            
            template<class Data, template<typename, typename> class Op>
            struct type_traits<expressions::Expression<RationalContext, Data, Op> >
            {
                enum
                {
                    is_number = true,
                    is_cputype = false,
                    is_realtype = true,
                    is_inttype = false,
                    is_string = false,
                    is_cpp_string = false,
                    is_c_string = false,
                    is_exact = true,
                    has_context = true,
                    is_native = false,
                    is_modulo = false,
                    has_infinity = false,
                    is_variable_precision = false,
                    has_squareroot = false,
                    has_full_power = false,
                    has_special_fns = false,
                    has_huge_exponent = true,
                    has_uniform_rng = false,
                    has_constants = false,
                    has_trigonometric = false
                };
                
                typedef RationalContext Context;
                typedef Rational PromoteType;
                typedef expressions::Expression<RationalContext, Data, Op> ConstReferenceType;
            };
        }
        
        namespace implementation
        {
            #define MAKE_RATIONAL_DEFS(Op1, Op2)                                                                                                            \
                template<>                                                                                                                                  \
                struct binary_operation_impl<Rational, Rational, Op1>                                                                                       \
                {                                                                                                                                           \
                    enum { supported = true, intermediate_expression = true };                                                                              \
                    typedef expressions::Expression<RationalContext, std::pair<expressions::Wrapper<RationalContext>,                                       \
                                                                               expressions::Wrapper<RationalContext> >, Op2> IntermediateType;              \
                    typedef Rational ResultType;                                                                                                            \
                };                                                                                                                                          \
                                                                                                                                                            \
                template<template<typename, typename> class O1, typename D1>                                                                                \
                struct binary_operation_impl<expressions::Expression<RationalContext, D1, O1>, Rational, Op1>                                               \
                {                                                                                                                                           \
                    enum { supported = true, intermediate_expression = true };                                                                              \
                    typedef expressions::Expression<RationalContext, std::pair<expressions::Expression<RationalContext, D1, O1>,                            \
                                                                               expressions::Wrapper<RationalContext> >, Op2> IntermediateType;              \
                    typedef Rational ResultType;                                                                                                            \
                };                                                                                                                                          \
                                                                                                                                                            \
                template<template<typename, typename> class O2, typename D2>                                                                                \
                struct binary_operation_impl<Rational, expressions::Expression<RationalContext, D2, O2>, Op1>                                               \
                {                                                                                                                                           \
                    enum { supported = true, intermediate_expression = true };                                                                              \
                    typedef expressions::Expression<RationalContext, std::pair<expressions::Wrapper<RationalContext>,                                       \
                                                                               expressions::Expression<RationalContext, D2, O2> >, Op2> IntermediateType;   \
                    typedef Rational ResultType;                                                                                                            \
                };                                                                                                                                          \
                                                                                                                                                            \
                template<template<typename, typename> class O1, typename D1, template<typename, typename> class O2, typename D2>                            \
                struct binary_operation_impl<expressions::Expression<RationalContext, D1, O1>, expressions::Expression<RationalContext, D2, O2>, Op1>       \
                {                                                                                                                                           \
                    enum { supported = true, intermediate_expression = true };                                                                              \
                    typedef expressions::Expression<RationalContext, std::pair<expressions::Expression<RationalContext, D1, O1>,                            \
                                                                               expressions::Expression<RationalContext, D2, O2> >, Op2> IntermediateType;   \
                    typedef Rational ResultType;                                                                                                            \
                };                                                                                                                                          \
                
            MAKE_RATIONAL_DEFS(arithmetic::op::addition, expressions::AddOp)
            MAKE_RATIONAL_DEFS(arithmetic::op::subtraction, expressions::SubOp)
            MAKE_RATIONAL_DEFS(arithmetic::op::multiplication, expressions::MulOp)
            MAKE_RATIONAL_DEFS(arithmetic::op::division, expressions::DivOp)
            MAKE_RATIONAL_DEFS(arithmetic::op::modulo, expressions::ModOp)
            #undef MAKE_RATIONAL_DEFS
            
            template<>
            struct unary_operation_impl<Rational, arithmetic::op::negation>
            {
                enum { supported = true, intermediate_expression = true };
                typedef expressions::Expression<RationalContext, expressions::Wrapper<RationalContext>,
                                                                 expressions::NegOp> IntermediateType;
                typedef Rational ResultType;
            };
            
            template<template<typename, typename> class O, class D>
            struct unary_operation_impl<expressions::Expression<RationalContext, D, O>, arithmetic::op::negation>
            {
                enum { supported = true, intermediate_expression = true };
                typedef expressions::Expression<RationalContext, expressions::Expression<RationalContext, D, O>,
                                                                 expressions::NegOp> IntermediateType;
                typedef Rational ResultType;
            };
        }
    }
}

#endif

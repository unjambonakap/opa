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

#ifndef PLLL_INCLUDE_GUARD__ARITHMETIC_GMP_CONVERSIONS_HPP
#define PLLL_INCLUDE_GUARD__ARITHMETIC_GMP_CONVERSIONS_HPP

/**
   \file
   \brief Conversion definitions for integers and floating point numbers.
   
   This header contains abstract templates used to implement most conversions involving integers and
   floating point numbers.
*/
namespace plll
{
    namespace arithmetic
    {
        namespace implementation
        {
            // Implementation for IntegerContext
    
            template<class Data, template<typename, typename> class Op, class Dest>
            class conversion_impl<expressions::Expression<IntegerContext, Data, Op>, Dest>
            {
            public:
                typedef typename conversion_impl<Integer, Dest>::RetVal RetVal;
                typedef void RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(typename Dest::Type & d, const expressions::Expression<IntegerContext, Data, Op> & v, const Dest & c)
                {
                    conversion_impl<Integer, Dest>::convert(d, static_cast<Integer>(v), c);
                }
        
                static RetVal convert(const expressions::Expression<IntegerContext, Data, Op> & v, const Dest & c)
                {
                    return conversion_impl<Integer, Dest>::convert(static_cast<Integer>(v), c);
                }
            };
    
            template<> class conversion_impl<signed long, IntegerContext>
            {
            public:
                typedef expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> RetVal;
                typedef void RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(Integer & d, const signed long & v, const IntegerContext & c)
                {
                    mpz_set_si(d.d_value, v);
                }
                
                static RetVal convert(const signed long & v, const IntegerContext & c)
                {
                    return RetVal(v, c);
                }
            };
    
            template<> class conversion_impl<unsigned long, IntegerContext>
            {
            public:
                typedef expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> RetVal;
                typedef void RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(Integer & d, const unsigned long & v, const IntegerContext & c)
                {
                    mpz_set_ui(d.d_value, v);
                }
                
                static RetVal convert(const unsigned long & v, const IntegerContext & c)
                {
                    return RetVal(v, c);
                }
            };
            
            template<> class conversion_impl<signed int, IntegerContext>
            {
            public:
                typedef expressions::Expression<IntegerContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> RetVal;
                typedef void RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(Integer & d, const signed int & v, const IntegerContext & c)
                {
                    mpz_set_si(d.d_value, v);
                }
                
                static RetVal convert(const signed int & v, const IntegerContext & c)
                {
                    return RetVal(v, c);
                }
            };
            
            template<> class conversion_impl<unsigned int, IntegerContext>
            {
            public:
                typedef expressions::Expression<IntegerContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> RetVal;
                typedef void RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(Integer & d, const unsigned int & v, const IntegerContext & c)
                {
                    mpz_set_ui(d.d_value, v);
                }
                
                static RetVal convert(const unsigned int & v, const IntegerContext & c)
                {
                    return RetVal(v, c);
                }
            };
            
            template<> class conversion_impl<long long, IntegerContext>
            {
            public:
                typedef expressions::Expression<IntegerContext, expressions::ConversionWrapper<long long>, expressions::ConvertOp_Context> RetVal;
                typedef void RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(Integer & d, const long long & v, const IntegerContext & c)
                {
                    Integer::mpz_set_ll(d.d_value, v);
                }
                
                static RetVal convert(const long long & v, const IntegerContext & c)
                {
                    return RetVal(v, c);
                }
            };
            
            template<> class conversion_impl<float, IntegerContext>
            {
            public:
                typedef expressions::Expression<IntegerContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> RetVal;
                typedef void RetVal_Frac;
                typedef expressions::Expression<IntegerContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> RetVal_Floor;
                typedef expressions::Expression<IntegerContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> RetVal_Round;
                typedef expressions::Expression<IntegerContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> RetVal_Round2;
                typedef expressions::Expression<IntegerContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> RetVal_Ceil;
            
                static void convert(Integer & d, const float & v, const IntegerContext & c)
                {
                    mpz_set_d(d.d_value, v);
                }
            
                static RetVal convert(const float & v, const IntegerContext & c)
                {
                    return RetVal(v, c);
                }
            
                static void floor(Integer & d, const float & v, const IntegerContext & c)
                {
                    mpz_set_d(d.d_value, std::floor(v));
                }
            
                static RetVal floor(const float & v, const IntegerContext & c)
                {
                    return RetVal(std::floor(v), c);
                }
            
                static void round(Integer & d, const float & v, const IntegerContext & c)
                {
                    mpz_set_d(d.d_value, rintf(v));
                }
            
                static RetVal round(const float & v, const IntegerContext & c)
                {
                    return RetVal(rintf(v), c);
                }
            
                static void round(Integer & d, const float & v, bool & roundUp, const IntegerContext & c)
                {
                    float rv = rintf(v);
                    roundUp = rv > v;
                    mpz_set_d(d.d_value, rv);
                }
            
                static RetVal round(const float & v, bool & roundUp, const IntegerContext & c)
                {
                    float rv = rintf(v);
                    roundUp = rv > v;
                    return RetVal(rv, c);
                }
            
                static void ceil(Integer & d, const float & v, const IntegerContext & c)
                {
                    mpz_set_d(d.d_value, std::ceil(v));
                }
            
                static RetVal ceil(const float & v, const IntegerContext & c)
                {
                    return RetVal(std::ceil(v), c);
                }
            };
        
            template<> class conversion_impl<double, IntegerContext>
            {
            public:
                typedef expressions::Expression<IntegerContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> RetVal;
                typedef void RetVal_Frac;
                typedef expressions::Expression<IntegerContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> RetVal_Floor;
                typedef expressions::Expression<IntegerContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> RetVal_Round;
                typedef expressions::Expression<IntegerContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> RetVal_Round2;
                typedef expressions::Expression<IntegerContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> RetVal_Ceil;
            
                static void convert(Integer & d, const double & v, const IntegerContext & c)
                {
                    mpz_set_d(d.d_value, v);
                }
            
                static RetVal convert(const double & v, const IntegerContext & c)
                {
                    return RetVal(v, c);
                }
            
                static void floor(Integer & d, const double & v, const IntegerContext & c)
                {
                    mpz_set_d(d.d_value, std::floor(v));
                }
            
                static RetVal floor(const double & v, const IntegerContext & c)
                {
                    return RetVal(std::floor(v), c);
                }
            
                static void round(Integer & d, const double & v, const IntegerContext & c)
                {
                    mpz_set_d(d.d_value, rint(v));
                }
                
                static RetVal round(const double & v, const IntegerContext & c)
                {
                    return RetVal(rint(v), c);
                }
                
                static void round(Integer & d, const double & v, bool & roundUp, const IntegerContext & c)
                {
                    double rv = rint(v);
                    roundUp = rv > v;
                    mpz_set_d(d.d_value, rv);
                }
                
                static RetVal round(const double & v, bool & roundUp, const IntegerContext & c)
                {
                    double rv = rint(v);
                    roundUp = rv > v;
                    return RetVal(rv, c);
                }
                
                static void ceil(Integer & d, const double & v, const IntegerContext & c)
                {
                    mpz_set_d(d.d_value, std::ceil(v));
                }
                
                static RetVal ceil(const double & v, const IntegerContext & c)
                {
                    return RetVal(v, c);
                }
            };
            
            template<> class conversion_impl<long double, IntegerContext>
            {
            public:
                typedef expressions::Expression<IntegerContext, expressions::ConversionWrapper<long double>, expressions::ConvertOp_Context> RetVal;
                typedef void RetVal_Frac;
                typedef expressions::Expression<IntegerContext, expressions::ConversionWrapper<long double>, expressions::ConvertOp_Context> RetVal_Floor;
                typedef expressions::Expression<IntegerContext, expressions::ConversionWrapper<long double>, expressions::ConvertOp_Context> RetVal_Round;
                typedef expressions::Expression<IntegerContext, expressions::ConversionWrapper<long double>, expressions::ConvertOp_Context> RetVal_Round2;
                typedef expressions::Expression<IntegerContext, expressions::ConversionWrapper<long double>, expressions::ConvertOp_Context> RetVal_Ceil;
                
                static void convert(Integer & d, const long double & v, const IntegerContext & c)
                {
                    Integer::mpz_set_ld(d.d_value, v);
                }
                
                static RetVal convert(const long double & v, const IntegerContext & c)
                {
                    return RetVal(v, c);
                }
                
                static void floor(Integer & d, const long double & v, const IntegerContext & c)
                {
                    Integer::mpz_set_ld(d.d_value, std::floor(v));
                }
                
                static RetVal floor(const long double & v, const IntegerContext & c)
                {
                    return RetVal(v, c);
                }
                
                static void round(Integer & d, const long double & v, const IntegerContext & c)
                {
                    Integer::mpz_set_ld(d.d_value, rintl(v));
                }
                
                static RetVal round(const long double & v, const IntegerContext & c)
                {
                    return RetVal(rintl(v), c);
                }
                
                static void round(Integer & d, const long double & v, bool & roundUp, const IntegerContext & c)
                {
                    long double rv = rintl(v);
                    roundUp = rv > v;
                    Integer::mpz_set_ld(d.d_value, rv);
                }
                
                static RetVal round(const long double & v, bool & roundUp, const IntegerContext & c)
                {
                    long double rv = rintl(v);
                    roundUp = rv > v;
                    return RetVal(rv, c);
                }
                
                static void ceil(Integer & d, const long double & v, const IntegerContext & c)
                {
                    Integer::mpz_set_ld(d.d_value, std::ceil(v));
                }
                
                static RetVal ceil(const long double & v, const IntegerContext & c)
                {
                    return RetVal(std::ceil(v), c);
                }
            };
            
            template<> class nativeconversion_impl<Integer>
            {
            public:
                static int toInt(const Integer & v)
                {
                    return mpz_get_si(v.d_value);
                }
        
                static unsigned int toUInt(const Integer & v)
                {
                    return mpz_get_ui(v.d_value);
                }
        
                static long toLong(const Integer & v)
                {
                    return mpz_get_si(v.d_value);
                }
        
                static unsigned long toULong(const Integer & v)
                {
                    return mpz_get_ui(v.d_value);
                }
        
                static long long toLongLong(const Integer & v)
                {
                    return arithmetic::Integer::mpz_get_ll(v.d_value);
                }
        
                static float toFloat(const Integer & v)
                {
                    return mpz_get_d(v.d_value);
                }
        
                static double toDouble(const Integer & v)
                {
                    return mpz_get_d(v.d_value);
                }
        
                static long double toLongDouble(const Integer & v)
                {
                    return arithmetic::Integer::mpz_get_ld(v.d_value);
                }
            };
    
            template<class Data, template<typename, typename> class Op> class nativeconversion_impl<expressions::Expression<IntegerContext, Data, Op> >
            {
            public:
                static int toInt(const expressions::Expression<IntegerContext, Data, Op> & v)
                {
                    return nativeconversion_impl<Integer>::toInt(static_cast<Integer>(v));
                }
        
                static unsigned int toUInt(const expressions::Expression<IntegerContext, Data, Op> & v)
                {
                    return nativeconversion_impl<Integer>::toUInt(static_cast<Integer>(v));
                }
        
                static long toLong(const expressions::Expression<IntegerContext, Data, Op> & v)
                {
                    return nativeconversion_impl<Integer>::toLong(static_cast<Integer>(v));
                }
        
                static unsigned long toULong(const expressions::Expression<IntegerContext, Data, Op> & v)
                {
                    return nativeconversion_impl<Integer>::toULong(static_cast<Integer>(v));
                }
        
                static long long toLongLong(const expressions::Expression<IntegerContext, Data, Op> & v)
                {
                    return nativeconversion_impl<Integer>::toLongLong(static_cast<Integer>(v));
                }
        
                static float toFloat(const expressions::Expression<IntegerContext, Data, Op> & v)
                {
                    return nativeconversion_impl<Integer>::toFloat(static_cast<Integer>(v));
                }
        
                static double toDouble(const expressions::Expression<IntegerContext, Data, Op> & v)
                {
                    return nativeconversion_impl<Integer>::toDouble(static_cast<Integer>(v));
                }
        
                static long double toLongDouble(const expressions::Expression<IntegerContext, Data, Op> & v)
                {
                    return nativeconversion_impl<Integer>::toLongDouble(static_cast<Integer>(v));
                }
            };
            
            // Implementation for RealContext
            
            template<> class conversion_impl<signed int, RealContext>
            {
            public:
                typedef expressions::Expression<RealContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> RetVal;
                typedef Real RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(Real & d, signed int v, const RealContext & c)
                {
                    mpfr_set_si(d.d_value, v, MPFR_RNDN);
                }
                
                static RetVal convert(signed int v, const RealContext & c)
                {
                    return RetVal(v, c);
                }
                
                static void convert_frac(Real & d, signed int v1, signed int v2, const RealContext & c)
                {
                    mpfr_set_si(d.d_value, v1, MPFR_RNDN);
                    Real t((signed long)v2, c);
                    d /= t;
                }
                
                static RetVal_Frac convert_frac(signed int v1, signed int v2, const RealContext & c)
                {
                    Real t1((signed long)v1, c), t2((signed long)v2, c);
                    t1 /= t2;
                    return t1;
                }
            };
    
            template<> class conversion_impl<unsigned int, RealContext>
            {
            public:
                typedef expressions::Expression<RealContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> RetVal;
                typedef Real RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(Real & d, unsigned int v, const RealContext & c)
                {
                    mpfr_set_ui(d.d_value, v, MPFR_RNDN);
                }
                
                static RetVal convert(unsigned int v, const RealContext & c)
                {
                    return RetVal(v, c);
                }
                
                static void convert_frac(Real & d, unsigned int v1, unsigned int v2, const RealContext & c)
                {
                    mpfr_set_ui(d.d_value, v1, MPFR_RNDN);
                    Real t((unsigned long)v2, c);
                    d /= t;
                }
                
                static RetVal_Frac convert_frac(unsigned int v1, unsigned int v2, const RealContext & c)
                {
                    Real t1((unsigned long)v1, c), t2((unsigned long)v2, c);
                    t1 /= t2;
                    return t1;
                }
            };
            
            template<> class conversion_impl<signed long, RealContext>
            {
            public:
                typedef expressions::Expression<RealContext, expressions::ConversionWrapper<signed long>, expressions::ConvertOp_Context> RetVal;
                typedef Real RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(Real & d, signed long v, const RealContext & c)
                {
                    mpfr_set_si(d.d_value, v, MPFR_RNDN);
                }
                
                static RetVal convert(signed long v, const RealContext & c)
                {
                    return RetVal(v, c);
                }
                
                static void convert_frac(Real & d, signed long v1, signed long v2, const RealContext & c)
                {
                    mpfr_set_si(d.d_value, v1, MPFR_RNDN);
                    Real t(v2, c);
                    d /= t;
                }
                
                static RetVal_Frac convert_frac(signed long v1, signed long v2, const RealContext & c)
                {
                    Real t1(v1, c), t2(v2, c);
                    t1 /= t2;
                    return t1;
                }
            };
            
            template<> class conversion_impl<unsigned long, RealContext>
            {
            public:
                typedef expressions::Expression<RealContext, expressions::ConversionWrapper<unsigned long>, expressions::ConvertOp_Context> RetVal;
                typedef Real RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
        
                static void convert(Real & d, unsigned long v, const RealContext & c)
                {
                    mpfr_set_ui(d.d_value, v, MPFR_RNDN);
                }
        
                static RetVal convert(unsigned long v, const RealContext & c)
                {
                    return RetVal(v, c);
                }
        
                static void convert_frac(Real & d, unsigned long v1, unsigned long v2, const RealContext & c)
                {
                    mpfr_set_ui(d.d_value, v1, MPFR_RNDN);
                    Real t(v2, c);
                    d /= t;
                }
        
                static RetVal_Frac convert_frac(unsigned long v1, unsigned long v2, const RealContext & c)
                {
                    Real t1(v1, c), t2(v2, c);
                    t1 /= t2;
                    return t1;
                }
            };
            
            template<> class conversion_impl<long long, RealContext>
            {
            public:
                typedef expressions::Expression<RealContext, expressions::ConversionWrapper<long long>, expressions::ConvertOp_Context> RetVal;
                typedef Real RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(Real & d, long long v, const RealContext & c)
                {
                    Real::mpfr_set_ll(d.d_value, v, MPFR_RNDN);
                }
                
                static RetVal convert(long long v, const RealContext & c)
                {
                    return RetVal(v, c);
                }
                
                static void convert_frac(Real & d, long long v1, long long v2, const RealContext & c)
                {
                    Real::mpfr_set_ll(d.d_value, v1, MPFR_RNDN);
                    Real t(v2, c);
                    d /= t;
                }
                
                static RetVal_Frac convert_frac(long long v1, long long v2, const RealContext & c)
                {
                    Real t1(v1, c), t2(v2, c);
                    t1 /= t2;
                    return t1;
                }
            };
            
            template<> class conversion_impl<float, RealContext>
            {
            public:
                typedef expressions::Expression<RealContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> RetVal;
                typedef Real RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(Real & d, float v, const RealContext & c)
                {
                    mpfr_set_d(d.d_value, v, MPFR_RNDN);
                }
                
                static RetVal convert(float v, const RealContext & c)
                {
                    return RetVal(v, c);
                }
                
                static void convert_frac(Real & d, float v1, float v2, const RealContext & c)
                {
                    mpfr_set_d(d.d_value, v1, MPFR_RNDN);
                    Real t(v2, c);
                    d /= t;
                }
                
                static RetVal_Frac convert_frac(float v1, float v2, const RealContext & c)
                {
                    Real t1(v1, c), t2(v2, c);
                    t1 /= t2;
                    return t1;
                }
            };
            
            template<> class conversion_impl<double, RealContext>
            {
            public:
                typedef expressions::Expression<RealContext, expressions::ConversionWrapper<double>, expressions::ConvertOp_Context> RetVal;
                typedef Real RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(Real & d, double v, const RealContext & c)
                {
                    mpfr_set_d(d.d_value, v, MPFR_RNDN);
                }
                
                static RetVal convert(double v, const RealContext & c)
                {
                    return RetVal(v, c);
                }
                
                static void convert_frac(Real & d, double v1, double v2, const RealContext & c)
                {
                    mpfr_set_d(d.d_value, v1, MPFR_RNDN);
                    Real t(v2, c);
                    d /= t;
                }
                
                static RetVal_Frac convert_frac(double v1, double v2, const RealContext & c)
                {
                    Real t1(v1, c), t2(v2, c);
                    t1 /= t2;
                    return t1;
                }
            };
            
            template<> class conversion_impl<long double, RealContext>
            {
            public:
                typedef expressions::Expression<RealContext, expressions::ConversionWrapper<long double>, expressions::ConvertOp_Context> RetVal;
                typedef Real RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(Real & d, long double v, const RealContext & c)
                {
                    mpfr_set_ld(d.d_value, v, MPFR_RNDN);
                }
                
                static RetVal convert(long double v, const RealContext & c)
                {
                    return RetVal(v, c);
                }
                
                static void convert_frac(Real & d, long double v1, long double v2, const RealContext & c)
                {
                    mpfr_set_ld(d.d_value, v1, MPFR_RNDN);
                    Real t(v2, c);
                    d /= t;
                }
                
                static RetVal_Frac convert_frac(long double v1, long double v2, const RealContext & c)
                {
                    Real t1(v1, c), t2(v2, c);
                    t1 /= t2;
                    return t1;
                }
            };
    
            template<> class conversion_impl<Integer, RealContext>
            {
            public:
                typedef expressions::Expression<RealContext, expressions::Expression<IntegerContext, expressions::Wrapper<IntegerContext>,
                                                                                     expressions::NoneOp>, expressions::ConvertOp_Context> RetVal;
                typedef Real RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(Real & d, const Integer & v, const RealContext & c)
                {
                    mpfr_set_z(d.d_value, v.d_value, MPFR_RNDN);
                }
                
                static RetVal convert(const Integer & v, const RealContext & c)
                {
                    return RetVal(expressions::make_expression<IntegerContext>(v), c);
                }
                
                static void convert_frac(Real & d, const Integer & v1, const Integer & v2, const RealContext & c)
                {
                    mpfr_set_z(d.d_value, v1.d_value, MPFR_RNDN);
                    Real t(v2, c);
                    d /= t;
                }
                
                static RetVal_Frac convert_frac(const Integer & v1, const Integer & v2, const RealContext & c)
                {
                    Real t1(v1, c), t2(v2, c);
                    t1 /= t2;
                    return t1;
                }
            };
            
            template<> class conversion_impl<Real, IntegerContext>
            {
            public:
                typedef typename expressions::Expression<IntegerContext, expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                                                                 expressions::NoneOp>, expressions::ConvertFloorOp_Context> RetVal;
                typedef void RetVal_Frac;
                typedef typename expressions::Expression<IntegerContext, expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                                                                 expressions::NoneOp>, expressions::ConvertFloorOp_Context> RetVal_Floor;
                typedef typename expressions::Expression<IntegerContext, expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                                                                 expressions::NoneOp>, expressions::ConvertRoundOp_Context> RetVal_Round;
                typedef typename expressions::Expression<IntegerContext, expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                                                                 expressions::NoneOp>, expressions::ConvertRound2Op_Context> RetVal_Round2;
                typedef typename expressions::Expression<IntegerContext, expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                                                                 expressions::NoneOp>, expressions::ConvertCeilOp_Context> RetVal_Ceil;
                
                static void convert(Integer & d, const Real & v, const IntegerContext & c)
                {
                    mpfr_get_z(d.d_value, v.d_value, MPFR_RNDD);
                }
                
                static RetVal convert(const Real & v, const IntegerContext & c)
                {
                    return RetVal(expressions::make_expression<RealContext>(v), implementation::g_intcontext);
                }
                
                static void floor(Integer & d, const Real & v, const IntegerContext & c)
                {
                    mpfr_get_z(d.d_value, v.d_value, MPFR_RNDD);
                }
                
                static RetVal_Floor floor(const Real & v, const IntegerContext & c)
                {
                    return RetVal_Floor(expressions::make_expression<RealContext>(v), implementation::g_intcontext);
                }
                
                static void round(Integer & d, const Real & v, const IntegerContext & c)
                {
                    mpfr_get_z(d.d_value, v.d_value, MPFR_RNDN);
                }
                
                static RetVal_Round round(const Real & v, const IntegerContext & c)
                {
                    return RetVal_Round(expressions::make_expression<RealContext>(v), implementation::g_intcontext);
                }
                
                static void round(Integer & d, const Real & v, bool & up, const IntegerContext & c)
                {
                    mpfr_get_z(d.d_value, v.d_value, MPFR_RNDN);
                    mpz_t tmp;
                    mpz_init(tmp);
                    mpfr_get_z(tmp, v.d_value, MPFR_RNDD);
                    up = mpz_cmp(tmp, d.d_value) != 0;
                    mpz_clear(tmp);
                }
                
                static RetVal_Round2 round(const Real & v, bool & up, const IntegerContext & c)
                {
                    return RetVal_Round2(expressions::make_expression<RealContext>(v),
                                         expressions::ConvertRound2Op_Context<IntegerContext,
                                                                              expressions::Expression<RealContext, expressions::Wrapper<RealContext>,
                                                                                                      expressions::NoneOp> >(up, implementation::g_intcontext));
                }
                
                static void ceil(Integer & d, const Real & v, const IntegerContext & c)
                {
                    mpfr_get_z(d.d_value, v.d_value, MPFR_RNDU);
                }
                
                static RetVal_Ceil ceil(const Real & v, const IntegerContext & c)
                {
                    return RetVal_Ceil(expressions::make_expression<RealContext>(v), implementation::g_intcontext);
                }
            };
            
            template<> class conversion_impl<Real, RealContext>
            /**
               \brief Provides an identity conversion from a type to itself for the `Real`
                      type. Note that the context is _not_ ignored.
             */
            {
            public:
                typedef Real RetVal;
                typedef Real RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(Real & d, const Real & v, const RealContext & c)
                // Nothing special: no precision change
                {
                    d = v;
                }
                
                static RetVal convert(const Real & v, const RealContext & c)
                {
                    return Real(v, c);
                }
                
                static void convert_frac(Real & d, const Real & v1, const Real & v2, const RealContext & c)
                // Nothing special: no precision change
                {
                    d = v1 / v2;
                }
                
                static RetVal_Frac convert_frac(const Real & v1, const Real & v2, const RealContext & c)
                // Adjust precision of result!
                {
                    Real r(c);
                    r = v1 / v2;
                    return r;
                }
            };
            
            template<> class nativeconversion_impl<Real>
            {
            public:
                static int toInt(const Real & v)
                {
                    return mpfr_get_si(v.d_value, MPFR_RNDD);
                }
                
                static int toInt_Floor(const Real & v)
                {
                    return mpfr_get_si(v.d_value, MPFR_RNDD);
                }
                
                static int toInt_Round(const Real & v)
                {
                    return mpfr_get_si(v.d_value, MPFR_RNDN);
                }
                
                static int toInt_Round(const Real & v, bool & up)
                {
                    int r = mpfr_get_si(v.d_value, MPFR_RNDN);
                    up = mpfr_cmp_si(v.d_value, r) < 0;
                    return r;
                }
                
                static int toInt_Ceil(const Real & v)
                {
                    return mpfr_get_si(v.d_value, MPFR_RNDU);
                }
                
                static unsigned int toUInt(const Real & v)
                {
                    return mpfr_get_ui(v.d_value, MPFR_RNDD);
                }
                
                static unsigned int toUInt_Floor(const Real & v)
                {
                    return mpfr_get_ui(v.d_value, MPFR_RNDD);
                }
                
                static unsigned int toUInt_Round(const Real & v)
                {
                    return mpfr_get_ui(v.d_value, MPFR_RNDN);
                }
                
                static unsigned int toUInt_Round(const Real & v, bool & up)
                {
                    unsigned int r = mpfr_get_ui(v.d_value, MPFR_RNDN);
                    up = mpfr_cmp_ui(v.d_value, r) < 0;
                    return r;
                }
                
                static unsigned int toUInt_Ceil(const Real & v)
                {
                    return mpfr_get_ui(v.d_value, MPFR_RNDU);
                }
                
                static long toLong(const Real & v)
                {
                    return mpfr_get_si(v.d_value, MPFR_RNDD);
                }
                
                static long toLong_Floor(const Real & v)
                {
                    return mpfr_get_si(v.d_value, MPFR_RNDD);
                }
                
                static long toLong_Round(const Real & v)
                {
                    return mpfr_get_si(v.d_value, MPFR_RNDN);
                }
                
                static long toLong_Round(const Real & v, bool & up)
                {
                    long r = mpfr_get_si(v.d_value, MPFR_RNDN);
                    up = mpfr_cmp_si(v.d_value, r) < 0;
                    return r;
                }
                
                static long toLong_Ceil(const Real & v)
                {
                    return mpfr_get_si(v.d_value, MPFR_RNDU);
                }
                
                static unsigned long toULong(const Real & v)
                {
                    return mpfr_get_ui(v.d_value, MPFR_RNDD);
                }
                
                static unsigned long toULong_Floor(const Real & v)
                {
                    return mpfr_get_ui(v.d_value, MPFR_RNDD);
                }
                
                static unsigned long toULong_Round(const Real & v)
                {
                    return mpfr_get_ui(v.d_value, MPFR_RNDN);
                }
                
                static unsigned long toULong_Round(const Real & v, bool & up)
                {
                    unsigned long r = mpfr_get_ui(v.d_value, MPFR_RNDN);
                    up = mpfr_cmp_ui(v.d_value, r) < 0;
                    return r;
                }
                
                static unsigned long toULong_Ceil(const Real & v)
                {
                    return mpfr_get_ui(v.d_value, MPFR_RNDU);
                }
                
                static long long toLongLong(const Real & v)
                {
                    return Real::mpfr_get_ll(v.d_value, MPFR_RNDD);
                }
                
                static long long toLongLong_Floor(const Real & v)
                {
                    return Real::mpfr_get_ll(v.d_value, MPFR_RNDD);
                }
                
                static long long toLongLong_Round(const Real & v)
                {
                    return Real::mpfr_get_ll(v.d_value, MPFR_RNDN);
                }
                
                static long long toLongLong_Round(const Real & v, bool & up)
                {
                    return Real::mpfr_get_ll(v.d_value, MPFR_RNDN, up);
                }
                
                static long long toLongLong_Ceil(const Real & v)
                {
                    return Real::mpfr_get_ll(v.d_value, MPFR_RNDU);
                }
                
                static float toFloat(const Real & v)
                {
                    return mpfr_get_d(v.d_value, MPFR_RNDN);
                }
        
                static double toDouble(const Real & v)
                {
                    return mpfr_get_d(v.d_value, MPFR_RNDN);
                }
        
                static long double toLongDouble(const Real & v)
                {
                    return mpfr_get_ld(v.d_value, MPFR_RNDN);
                }
            };
        }
        
        namespace traits
        {
            template<class Data, template<typename, typename> class Op>
            struct type_traits<expressions::Expression<IntegerContext, Data, Op> >
            {
                enum
                {
                    is_number = true,
                    is_cputype = false,
                    is_realtype = false,
                    is_inttype = true,
                    is_string = false,
                    is_cpp_string = false,
                    is_c_string = false,
                    is_exact = true,
                    has_uniform_rng = true,
                    has_context = true,
                    is_native = false,
                    is_modulo = false,
                    has_infinity = false,
                    is_variable_precision = false,
                    has_squareroot = false,
                    has_full_power = false,
                    has_special_fns = false,
                    has_huge_exponent = false,
                    has_constants = false,
                    has_trigonometric = false
                };
                
                typedef IntegerContext Context;
                typedef Integer PromoteType;
                typedef expressions::Expression<IntegerContext, Data, Op> ConstReferenceType;
            };
            
            template<class Data, template<typename, typename> class Op>
            struct type_traits<expressions::Expression<RealContext, Data, Op> >
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
                    is_exact = false,
                    has_uniform_rng = true,
                    has_context = true,
                    is_native = false,
                    is_modulo = false,
                    has_infinity = true,
                    is_variable_precision = true,
                    has_squareroot = true,
                    has_full_power = true,
                    has_special_fns = true,
                    has_huge_exponent = true,
                    has_constants = true,
                    has_trigonometric = true
                };
                
                typedef RealContext Context;
                typedef Real PromoteType;
                typedef expressions::Expression<RealContext, Data, Op> ConstReferenceType;
            };
        }
        
        namespace implementation
        {
            template<>
            class from_string_conversion<IntegerContext>
            {
            public:
                typedef Integer RetVal1;
                typedef Integer RetVal2;
            
                static bool convert(Integer & res, const std::string & s, const IntegerContext & c);
                static bool convert(Integer & res, const char * s, const IntegerContext & c);
            
                static RetVal1 convert(const std::string & s, const IntegerContext & c)
                {
                    Integer res;
                    convert(res, s, c);
                    return res;
                }
            
                static RetVal2 convert(const char * s, const IntegerContext & c)
                {
                    Integer res;
                    convert(res, s, c);
                    return res;
                }
            };
        
            template<>
            class to_string_conversion<Integer>
            {
            public:
                typedef std::string RetVal;
            
                static RetVal convert(const Integer &);
                static RetVal convert(const Integer &, unsigned);
                static void convert(std::string &, const Integer &);
                static void convert(std::string &, const Integer &, unsigned);
            };
        
            template<>
            class from_string_conversion<RealContext>
            {
            public:
                typedef Real RetVal1;
                typedef Real RetVal2;
            
                static bool convert(Real & res, const std::string & s, const RealContext & c);
                static bool convert(Real & res, const char * s, const RealContext & c);
            
                static RetVal1 convert(const std::string & s, const RealContext & c)
                {
                    Real res(c);
                    convert(res, s, c);
                    return res;
                }
            
                static RetVal2 convert(const char * s, const RealContext & c)
                {
                    Real res(c);
                    convert(res, s, c);
                    return res;
                }
            };
        
            template<>
            class to_string_conversion<Real>
            {
            public:
                typedef std::string RetVal;
            
                static RetVal convert(const Real &);
                static RetVal convert(const Real &, unsigned);
                static void convert(std::string &, const Real &);
                static void convert(std::string &, const Real &, unsigned);
            };
            
            #define MAKE_INTEGER_DEFS(Op1, Op2)                                                                                                           \
                template<>                                                                                                                                \
                struct binary_operation_impl<Integer, Integer, Op1>                                                                                       \
                {                                                                                                                                         \
                    enum { supported = true, intermediate_expression = true };                                                                            \
                    typedef expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,                                       \
                                                                              expressions::Wrapper<IntegerContext> >, Op2> IntermediateType;              \
                    typedef Integer ResultType;                                                                                                           \
                };                                                                                                                                        \
                                                                                                                                                          \
                template<template<typename, typename> class O1, typename D1>                                                                              \
                struct binary_operation_impl<expressions::Expression<IntegerContext, D1, O1>, Integer, Op1>                                               \
                {                                                                                                                                         \
                    enum { supported = true, intermediate_expression = true };                                                                            \
                    typedef expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, D1, O1>,                            \
                                                                              expressions::Wrapper<IntegerContext> >, Op2> IntermediateType;              \
                    typedef Integer ResultType;                                                                                                           \
                };                                                                                                                                        \
                                                                                                                                                          \
                template<template<typename, typename> class O2, typename D2>                                                                              \
                struct binary_operation_impl<Integer, expressions::Expression<IntegerContext, D2, O2>, Op1>                                               \
                {                                                                                                                                         \
                    enum { supported = true, intermediate_expression = true };                                                                            \
                    typedef expressions::Expression<IntegerContext, std::pair<expressions::Wrapper<IntegerContext>,                                       \
                                                                              expressions::Expression<IntegerContext, D2, O2> >, Op2> IntermediateType;   \
                    typedef Integer ResultType;                                                                                                           \
                };                                                                                                                                        \
                                                                                                                                                          \
                template<template<typename, typename> class O1, typename D1, template<typename, typename> class O2, typename D2>                          \
                struct binary_operation_impl<expressions::Expression<IntegerContext, D1, O1>, expressions::Expression<IntegerContext, D2, O2>, Op1>       \
                {                                                                                                                                         \
                    enum { supported = true, intermediate_expression = true };                                                                            \
                    typedef expressions::Expression<IntegerContext, std::pair<expressions::Expression<IntegerContext, D1, O1>,                            \
                                                                              expressions::Expression<IntegerContext, D2, O2> >, Op2> IntermediateType;   \
                    typedef Integer ResultType;                                                                                                           \
                };                                                                                                                                        \
            
            MAKE_INTEGER_DEFS(arithmetic::op::addition, expressions::AddOp)
            MAKE_INTEGER_DEFS(arithmetic::op::subtraction, expressions::SubOp)
            MAKE_INTEGER_DEFS(arithmetic::op::multiplication, expressions::MulOp)
            MAKE_INTEGER_DEFS(arithmetic::op::division, expressions::DivOp)
            MAKE_INTEGER_DEFS(arithmetic::op::modulo, expressions::ModOp)
            #undef MAKE_INTEGER_DEFS
            
            template<>
            struct unary_operation_impl<Integer, arithmetic::op::negation>
            {
                enum { supported = true, intermediate_expression = true };
                typedef expressions::Expression<IntegerContext, expressions::Wrapper<IntegerContext>,
                                                                expressions::NegOp> IntermediateType;
                typedef Integer ResultType;
            };
            
            template<template<typename, typename> class O, class D>
            struct unary_operation_impl<expressions::Expression<IntegerContext, D, O>, arithmetic::op::negation>
            {
                enum { supported = true, intermediate_expression = true };
                typedef expressions::Expression<IntegerContext, expressions::Expression<IntegerContext, D, O>,
                                                                expressions::NegOp> IntermediateType;
                typedef Integer ResultType;
            };
            
            #define MAKE_REAL_DEFS(Op1, Op2)                                                                                                        \
                template<>                                                                                                                          \
                struct binary_operation_impl<Real, Real, Op1>                                                                                       \
                {                                                                                                                                   \
                    enum { supported = true, intermediate_expression = true };                                                                      \
                    typedef expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,                                       \
                                                                           expressions::Wrapper<RealContext> >, Op2> IntermediateType;              \
                    typedef Real ResultType;                                                                                                        \
                };                                                                                                                                  \
                                                                                                                                                    \
                template<template<typename, typename> class O1, typename D1>                                                                        \
                struct binary_operation_impl<expressions::Expression<RealContext, D1, O1>, Real, Op1>                                               \
                {                                                                                                                                   \
                    enum { supported = true, intermediate_expression = true };                                                                      \
                    typedef expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, D1, O1>,                            \
                                                                           expressions::Wrapper<RealContext> >, Op2> IntermediateType;              \
                    typedef Real ResultType;                                                                                                        \
                };                                                                                                                                  \
                                                                                                                                                    \
                template<template<typename, typename> class O2, typename D2>                                                                        \
                struct binary_operation_impl<Real, expressions::Expression<RealContext, D2, O2>, Op1>                                               \
                {                                                                                                                                   \
                    enum { supported = true, intermediate_expression = true };                                                                      \
                    typedef expressions::Expression<RealContext, std::pair<expressions::Wrapper<RealContext>,                                       \
                                                                           expressions::Expression<RealContext, D2, O2> >, Op2> IntermediateType;   \
                    typedef Real ResultType;                                                                                                        \
                };                                                                                                                                  \
                                                                                                                                                    \
                template<template<typename, typename> class O1, typename D1, template<typename, typename> class O2, typename D2>                    \
                struct binary_operation_impl<expressions::Expression<RealContext, D1, O1>, expressions::Expression<RealContext, D2, O2>, Op1>       \
                {                                                                                                                                   \
                    enum { supported = true, intermediate_expression = true };                                                                      \
                    typedef expressions::Expression<RealContext, std::pair<expressions::Expression<RealContext, D1, O1>,                            \
                                                                           expressions::Expression<RealContext, D2, O2> >, Op2> IntermediateType;   \
                    typedef Real ResultType;                                                                                                        \
                };                                                                                                                                  \
            
            MAKE_REAL_DEFS(arithmetic::op::addition, expressions::AddOp)
            MAKE_REAL_DEFS(arithmetic::op::subtraction, expressions::SubOp)
            MAKE_REAL_DEFS(arithmetic::op::multiplication, expressions::MulOp)
            MAKE_REAL_DEFS(arithmetic::op::division, expressions::DivOp)
            MAKE_REAL_DEFS(arithmetic::op::modulo, expressions::ModOp)
            #undef MAKE_REAL_DEFS
            
            template<>
            struct unary_operation_impl<Real, arithmetic::op::negation>
            {
                enum { supported = true, intermediate_expression = true };
                typedef expressions::Expression<RealContext, expressions::Wrapper<RealContext>, expressions::NegOp> IntermediateType;
                typedef Real ResultType;
            };
            
            template<template<typename, typename> class O, class D>
            struct unary_operation_impl<expressions::Expression<RealContext, D, O>, arithmetic::op::negation>
            {
                enum { supported = true, intermediate_expression = true };
                typedef expressions::Expression<RealContext, expressions::Expression<RealContext, D, O>, expressions::NegOp> IntermediateType;
                typedef Real ResultType;
            };
        }
    }
}

#endif

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

#ifndef PLLL_INCLUDE_GUARD__NINT_WRAPPER_CONV
#define PLLL_INCLUDE_GUARD__NINT_WRAPPER_CONV

/**
   \file
   \brief Conversion definitions for native integers.
   
   This header contains abstract templates used to implement most conversions involving native
   integers.
*/

#include <plll/rational.hpp>

namespace plll
{
    namespace arithmetic
    {
        namespace implementation
        {
            // Implementation for NIntContext<Type>
            
            // To avoid that conversion_impl<NInt<Type>, NIntContext<Type> > is defined twice, we
            // instanciate these conversions explicitly for distinct types.
            template<> class conversion_impl<NInt<long>, NIntContext<int> >
            {
            public:
                typedef NInt<int> RetVal;
                typedef void RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
        
                static void convert(NInt<int> & d, const NInt<long> & v, const NIntContext<int> & c)
                {
                    d.d_value = v.d_value;
                }
        
                static RetVal convert(const NInt<long> & v, const NIntContext<int> & c)
                {
                    return NInt<int>(false, v.d_value);
                }
            };
        
            template<> class conversion_impl<NInt<long long>, NIntContext<int> >
            {
            public:
                typedef NInt<int> RetVal;
                typedef void RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
        
                static void convert(NInt<int> & d, const NInt<long long> & v, const NIntContext<int> & c)
                {
                    d.d_value = v.d_value;
                }
        
                static RetVal convert(const NInt<long long> & v, const NIntContext<int> & c)
                {
                    return NInt<int>(false, v.d_value);
                }
            };
        
            template<> class conversion_impl<NInt<int>, NIntContext<long> >
            {
            public:
                typedef NInt<long> RetVal;
                typedef void RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
        
                static void convert(NInt<long> & d, const NInt<int> & v, const NIntContext<long> & c)
                {
                    d.d_value = v.d_value;
                }
        
                static RetVal convert(const NInt<int> & v, const NIntContext<long> & c)
                {
                    return NInt<long>(false, v.d_value);
                }
            };
        
            template<> class conversion_impl<NInt<long long>, NIntContext<long> >
            {
            public:
                typedef NInt<long> RetVal;
                typedef void RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
        
                static void convert(NInt<long> & d, const NInt<long long> & v, const NIntContext<long> & c)
                {
                    d.d_value = v.d_value;
                }
        
                static RetVal convert(const NInt<long long> & v, const NIntContext<long> & c)
                {
                    return NInt<long>(false, v.d_value);
                }
            };
        
            template<> class conversion_impl<NInt<int>, NIntContext<long long> >
            {
            public:
                typedef NInt<long long> RetVal;
                typedef void RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
        
                static void convert(NInt<long long> & d, const NInt<int> & v, const NIntContext<long long> & c)
                {
                    d.d_value = v.d_value;
                }
        
                static RetVal convert(const NInt<int> & v, const NIntContext<long long> & c)
                {
                    return NInt<long long>(false, v.d_value);
                }
            };
        
            template<> class conversion_impl<NInt<long>, NIntContext<long long> >
            {
            public:
                typedef NInt<long long> RetVal;
                typedef void RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
        
                static void convert(NInt<long long> & d, const NInt<long> & v, const NIntContext<long long> & c)
                {
                    d.d_value = v.d_value;
                }
        
                static RetVal convert(const NInt<long> & v, const NIntContext<long long> & c)
                {
                    return NInt<long long>(false, v.d_value);
                }
            };
    
            template<typename Type> class conversion_impl<signed int, NIntContext<Type> >
            {
            public:
                typedef NInt<Type> RetVal;
                typedef void RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
        
                static void convert(NInt<Type> & d, signed int v, const NIntContext<Type> & c)
                {
                    d.d_value = v;
                }
        
                static RetVal convert(signed int v, const NIntContext<Type> & c)
                {
                    return NInt<Type>(false, v);
                }
            };
    
            template<typename Type> class conversion_impl<unsigned int, NIntContext<Type> >
            {
            public:
                typedef NInt<Type> RetVal;
                typedef void RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
        
                static void convert(NInt<Type> & d, unsigned int v, const NIntContext<Type> & c)
                {
                    d.d_value = v;
                }
        
                static RetVal convert(unsigned int v, const NIntContext<Type> & c)
                {
                    return NInt<Type>(false, v);
                }
            };
    
            template<typename Type> class conversion_impl<signed long, NIntContext<Type> >
            {
            public:
                typedef NInt<Type> RetVal;
                typedef void RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
        
                static void convert(NInt<Type> & d, signed long v, const NIntContext<Type> & c)
                {
                    d.d_value = v;
                }
        
                static RetVal convert(signed long v, const NIntContext<Type> & c)
                {
                    return NInt<Type>(false, v);
                }
            };
    
            template<typename Type> class conversion_impl<unsigned long, NIntContext<Type> >
            {
            public:
                typedef NInt<Type> RetVal;
                typedef void RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
        
                static void convert(NInt<Type> & d, unsigned long v, const NIntContext<Type> & c)
                {
                    d.d_value = v;
                }
        
                static RetVal convert(unsigned long v, const NIntContext<Type> & c)
                {
                    return NInt<Type>(false, v);
                }
            };
    
            template<typename Type> class conversion_impl<long long, NIntContext<Type> >
            {
            public:
                typedef NInt<Type> RetVal;
                typedef void RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
        
                static void convert(NInt<Type> & d, long long v, const NIntContext<Type> & c)
                {
                    d.d_value = v;
                }
        
                static RetVal convert(long long v, const NIntContext<Type> & c)
                {
                    return NInt<Type>(false, v);
                }
            };
            
            template<typename Type> // since there is no std::round for float, double and long double, we have
            class RoundingNFP;   // to implement our own
            
            template<>
            class RoundingNFP<float>
            {
            public:
                static float round(float f) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    return roundf(f);
                }
                
                static float round(float f, bool & up) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    float v = roundf(f);
                    up = v > f;
                    return v;
                }
            };
            
            template<>
            class RoundingNFP<double>
            {
            public:
                static double round(double d) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    return ::round(d);
                }
                
                static double round(double d, bool & up) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    double v = ::round(d);
                    up = v > d;
                    return v;
                }
            };
            
            template<>
            class RoundingNFP<long double>
            {
            public:
                static long double round(long double ld) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    return roundl(ld);
                }
                
                static long double round(long double ld, bool & up) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    long double v = roundl(ld);
                    up = v > ld;
                    return v;
                }
            };
            
            template<typename Type> class conversion_impl<float, NIntContext<Type> >
            {
            public:
                typedef NInt<Type> RetVal;
                typedef void RetVal_Frac;
                typedef NInt<Type> RetVal_Floor;
                typedef NInt<Type> RetVal_Round;
                typedef NInt<Type> RetVal_Round2;
                typedef NInt<Type> RetVal_Ceil;
        
                static void convert(NInt<Type> & d, float v, const NIntContext<Type> & c)
                {
                    d.d_value = v;
                }
        
                static RetVal convert(float v, const NIntContext<Type> & c)
                {
                    return NInt<Type>(false, v);
                }
        
                static void floor(NInt<Type> & d, float v, const NIntContext<Type> & c)
                {
                    d.d_value = std::floor(v);
                }
        
                static RetVal_Floor floor(float v, const NIntContext<Type> & c)
                {
                    return NInt<Type>(false, std::floor(v));
                }
        
                static void round(NInt<Type> & d, float v, const NIntContext<Type> & c)
                {
                    d.d_value = RoundingNFP<float>::round(v);
                }
        
                static RetVal_Round round(float v, const NIntContext<Type> & c)
                {
                    return NInt<Type>(false, RoundingNFP<float>::round(v));
                }
        
                static void round(NInt<Type> & d, float v, bool & up, const NIntContext<Type> & c)
                {
                    d.d_value = RoundingNFP<float>::round(v, up);
                }
        
                static RetVal_Round2 round(float v, bool & up, const NIntContext<Type> & c)
                {
                    return NInt<Type>(false, RoundingNFP<float>::round(v, up));
                }
        
                static void ceil(NInt<Type> & d, float v, const NIntContext<Type> & c)
                {
                    d.d_value = std::ceil(v);
                }
        
                static RetVal_Ceil ceil(float v, const NIntContext<Type> & c)
                {
                    return NInt<Type>(false, std::ceil(v));
                }
            };
    
            template<typename Type> class conversion_impl<double, NIntContext<Type> >
            {
            public:
                typedef NInt<Type> RetVal;
                typedef void RetVal_Frac;
                typedef NInt<Type> RetVal_Floor;
                typedef NInt<Type> RetVal_Round;
                typedef NInt<Type> RetVal_Round2;
                typedef NInt<Type> RetVal_Ceil;
        
                static void convert(NInt<Type> & d, double v, const NIntContext<Type> & c)
                {
                    d.d_value = v;
                }
        
                static RetVal convert(double v, const NIntContext<Type> & c)
                {
                    return NInt<Type>(false, v);
                }
        
                static void floor(NInt<Type> & d, double v, const NIntContext<Type> & c)
                {
                    d.d_value = std::floor(v);
                }
        
                static RetVal_Floor floor(double v, const NIntContext<Type> & c)
                {
                    return NInt<Type>(false, std::floor(v));
                }
        
                static void round(NInt<Type> & d, double v, const NIntContext<Type> & c)
                {
                    d.d_value = RoundingNFP<double>::round(v);
                }
        
                static RetVal_Round round(double v, const NIntContext<Type> & c)
                {
                    return NInt<Type>(false, RoundingNFP<double>::round(v));
                }
        
                static void round(NInt<Type> & d, double v, bool & up, const NIntContext<Type> & c)
                {
                    d.d_value = RoundingNFP<double>::round(v, up);
                }
        
                static RetVal_Round2 round(double v, bool & up, const NIntContext<Type> & c)
                {
                    return NInt<Type>(false, RoundingNFP<double>::round(v, up));
                }
        
                static void ceil(NInt<Type> & d, double v, const NIntContext<Type> & c)
                {
                    d.d_value = std::ceil(v);
                }
        
                static RetVal_Ceil ceil(double v, const NIntContext<Type> & c)
                {
                    return NInt<Type>(false, std::ceil(v));
                }
            };
    
            template<typename Type> class conversion_impl<long double, NIntContext<Type> >
            {
            public:
                typedef NInt<Type> RetVal;
                typedef void RetVal_Frac;
                typedef NInt<Type> RetVal_Floor;
                typedef NInt<Type> RetVal_Round;
                typedef NInt<Type> RetVal_Round2;
                typedef NInt<Type> RetVal_Ceil;
        
                static void convert(NInt<Type> & d, long double v, const NIntContext<Type> & c)
                {
                    d.d_value = v;
                }
        
                static RetVal convert(long double v, const NIntContext<Type> & c)
                {
                    return NInt<Type>(false, v);
                }
        
                static void floor(NInt<Type> & d, long double v, const NIntContext<Type> & c)
                {
                    d.d_value = std::floor(v);
                }
        
                static RetVal_Floor floor(long double v, const NIntContext<Type> & c)
                {
                    return NInt<Type>(false, std::floor(v));
                }
        
                static void round(NInt<Type> & d, long double v, const NIntContext<Type> & c)
                {
                    d.d_value = RoundingNFP<long double>::round(v);
                }
        
                static RetVal_Round round(long double v, const NIntContext<Type> & c)
                {
                    return NInt<Type>(false, RoundingNFP<long double>::round(v));
                }
        
                static void round(NInt<Type> & d, long double v, bool & up, const NIntContext<Type> & c)
                {
                    d.d_value = RoundingNFP<long double>::round(v, up);
                }
        
                static RetVal_Round2 round(long double v, bool & up, const NIntContext<Type> & c)
                {
                    return NInt<Type>(false, RoundingNFP<long double>::round(v, up));
                }
        
                static void ceil(NInt<Type> & d, long double v, const NIntContext<Type> & c)
                {
                    d.d_value = std::ceil(v);
                }
        
                static RetVal_Ceil ceil(long double v, const NIntContext<Type> & c)
                {
                    return NInt<Type>(false, std::ceil(v));
                }
            };
            
            template<typename Type> class conversion_impl<Integer, NIntContext<Type> >
            {
            public:
                typedef NInt<Type> RetVal;
                typedef void RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
        
                static void convert(NInt<Type> & d, const Integer & v, const NIntContext<Type> & c)
                {
                    d.d_value = arithmetic::convert<Type>(v);
                }
        
                static RetVal convert(const Integer & v, const NIntContext<Type> & c)
                {
                    return NInt<Type>(false, arithmetic::convert<Type>(v));
                }
            };
            
            template<typename Type> class conversion_impl<NInt<Type>, IntegerContext>
            {
            public:
                typedef typename conversion_impl<Type, IntegerContext>::RetVal RetVal;
                typedef void RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
        
                static void convert(Integer & d, const NInt<Type> & v, const IntegerContext & c)
                {
                    conversion_impl<Type, IntegerContext>::convert(d, v.d_value, c);
                }
        
                static RetVal convert(const NInt<Type> & v, const IntegerContext & c)
                {
                    return conversion_impl<Type, IntegerContext>::convert(v.d_value, c);
                }
            };
            
            template<typename Type> class conversion_impl<Real, NIntContext<Type> >
            {
            public:
                typedef NInt<Type> RetVal;
                typedef void RetVal_Frac;
                typedef NInt<Type> RetVal_Floor;
                typedef NInt<Type> RetVal_Round;
                typedef NInt<Type> RetVal_Round2;
                typedef NInt<Type> RetVal_Ceil;
        
                static void convert(NInt<Type> & d, const Real & v, const NIntContext<Type> & c)
                {
                    d.d_value = arithmetic::convert<Type>(v);
                }
        
                static RetVal convert(const Real & v, const NIntContext<Type> & c)
                {
                    return NInt<Type>(false, arithmetic::convert<Type>(v));
                }
        
                static void floor(NInt<Type> & d, const Real & v, const NIntContext<Type> & c)
                {
                    d.d_value = convert_floor<Type>(v);
                }
        
                static RetVal_Floor floor(const Real & v, const NIntContext<Type> & c)
                {
                    return NInt<Type>(false, convert_floor<Type>(v));
                }
        
                static void round(NInt<Type> & d, const Real & v, const NIntContext<Type> & c)
                {
                    d.d_value = convert_round<Type>(v);
                }
        
                static RetVal_Round round(const Real & v, const NIntContext<Type> & c)
                {
                    return NInt<Type>(false, convert_round<Type>(v));
                }
        
                static void round(NInt<Type> & d, const Real & v, bool & up, const NIntContext<Type> & c)
                {
                    d.d_value = convert_round<Type>(v, up);
                }
        
                static RetVal_Round2 round(const Real & v, bool & up, const NIntContext<Type> & c)
                {
                    return NInt<Type>(false, convert_round<Type>(v, up));
                }
        
                static void ceil(NInt<Type> & d, const Real & v, const NIntContext<Type> & c)
                {
                    d.d_value = convert_ceil<Type>(v);
                }
        
                static RetVal_Ceil ceil(const Real & v, const NIntContext<Type> & c)
                {
                    return NInt<Type>(false, convert_ceil<Type>(v));
                }
            };
            
            template<typename Type> class conversion_impl<NInt<Type>, RealContext>
            {
            public:
                typedef typename conversion_impl<Type, RealContext>::RetVal RetVal;
                typedef typename conversion_impl<Type, RealContext>::RetVal_Frac RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
        
                static void convert(Real & d, const NInt<Type> & v, const RealContext & c)
                {
                    conversion_impl<Type, RealContext>::convert(d, v.d_value, c);
                }
        
                static RetVal convert(const NInt<Type> & v, const RealContext & c)
                {
                    return conversion_impl<Type, RealContext>::convert(v.d_value, c);
                }
        
                static void convert_frac(Real & d, const NInt<Type> & v1, const NInt<Type> & v2, const RealContext & c)
                {
                    conversion_impl<Type, RealContext>::convert_frac(d, v1.d_value, v2.d_value, c);
                }
        
                static RetVal_Frac convert_frac(const NInt<Type> & v1, const NInt<Type> & v2, const RealContext & c)
                {
                    return conversion_impl<Type, RealContext>::convert_frac(v1.d_value, v2.d_value, c);
                }
            };
            
            template<typename Type> class conversion_impl<Rational, NIntContext<Type> >
            {
            public:
                typedef NInt<Type> RetVal;
                typedef void RetVal_Frac;
                typedef NInt<Type> RetVal_Floor;
                typedef NInt<Type> RetVal_Round;
                typedef NInt<Type> RetVal_Round2;
                typedef NInt<Type> RetVal_Ceil;
        
                static void convert(NInt<Type> & d, const Rational & v, const NIntContext<Type> & c)
                {
                    d.d_value = arithmetic::convert<Type>(v);
                }
        
                static RetVal convert(const Rational & v, const NIntContext<Type> & c)
                {
                    return NInt<Type>(false, arithmetic::convert<Type>(v));
                }
        
                static void floor(NInt<Type> & d, const Rational & v, const NIntContext<Type> & c)
                {
                    d.d_value = convert_floor<Type>(v);
                }
        
                static RetVal_Floor floor(const Rational & v, const NIntContext<Type> & c)
                {
                    return NInt<Type>(false, convert_floor<Type>(v));
                }
        
                static void round(NInt<Type> & d, const Rational & v, const NIntContext<Type> & c)
                {
                    d.d_value = convert_round<Type>(v);
                }
        
                static RetVal_Round round(const Rational & v, const NIntContext<Type> & c)
                {
                    return NInt<Type>(false, convert_round<Type>(v));
                }
        
                static void round(NInt<Type> & d, const Rational & v, bool & up, const NIntContext<Type> & c)
                {
                    d.d_value = convert_round<Type>(v, up);
                }
        
                static RetVal_Round2 round(const Rational & v, bool & up, const NIntContext<Type> & c)
                {
                    return NInt<Type>(false, convert_round<Type>(v, up));
                }
        
                static void ceil(NInt<Type> & d, const Rational & v, const NIntContext<Type> & c)
                {
                    d.d_value = convert_ceil<Type>(v);
                }
        
                static RetVal_Ceil ceil(const Rational & v, const NIntContext<Type> & c)
                {
                    return NInt<Type>(false, convert_ceil<Type>(v));
                }
            };
            
            template<typename Type> class conversion_impl<NInt<Type>, RationalContext>
            {
            public:
                typedef typename conversion_impl<Type, RationalContext>::RetVal RetVal;
                typedef typename conversion_impl<Type, RationalContext>::RetVal_Frac RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
        
                static void convert(Rational & d, const NInt<Type> & v, const RationalContext & c)
                {
                    conversion_impl<Type, RationalContext>::convert(d, v.d_value, c);
                }
        
                static RetVal convert(const NInt<Type> & v, const RationalContext & c)
                {
                    return conversion_impl<Type, RationalContext>::convert(v.d_value, c);
                }
                
                static void convert_frac(Rational & d, const NInt<Type> & v1, const NInt<Type> & v2, const RationalContext & c)
                {
                    conversion_impl<Type, RationalContext>::convert_frac(d, v1.d_value, v2.d_value, c);
                }
        
                static RetVal_Frac convert_frac(const NInt<Type> & v1, const NInt<Type> & v2, const RationalContext & c)
                {
                    return conversion_impl<Type, RationalContext>::convert_frac(v1.d_value, v2.d_value, c);
                }
            };
            
            template<typename Type> class nativeconversion_impl<NInt<Type> >
            {
            public:
                static int toInt(const NInt<Type> & v)
                {
                    return v.d_value;
                }
        
                static unsigned int toUInt(const NInt<Type> & v)
                {
                    return v.d_value;
                }
        
                static long toLong(const NInt<Type> & v)
                {
                    return v.d_value;
                }
        
                static unsigned long toULong(const NInt<Type> & v)
                {
                    return v.d_value;
                }
        
                static long long toLongLong(const NInt<Type> & v)
                {
                    return v.d_value;
                }
        
                static float toFloat(const NInt<Type> & v)
                {
                    return v.d_value;
                }
        
                static double toDouble(const NInt<Type> & v)
                {
                    return v.d_value;
                }
        
                static long double toLongDouble(const NInt<Type> & v)
                {
                    return v.d_value;
                }
            };
        
            template<typename Type, typename Op>
            struct binary_operation_impl<NInt<Type>, NInt<Type>, Op>
            {
                enum { supported = true, intermediate_expression = false };
                typedef NInt<Type> IntermediateType;
                typedef NInt<Type> ResultType;
            };
        
            template<typename Type>
            struct unary_operation_impl<NInt<Type>, arithmetic::op::negation>
            {
                enum { supported = true, intermediate_expression = false };
                typedef NInt<Type> IntermediateType;
                typedef NInt<Type> ResultType;
            };
        }
    }
}

#endif

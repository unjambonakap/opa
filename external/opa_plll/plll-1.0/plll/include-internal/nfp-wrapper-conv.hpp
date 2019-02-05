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

#ifndef PLLL_INCLUDE_GUARD__NFP_WRAPPER_CONVERSION
#define PLLL_INCLUDE_GUARD__NFP_WRAPPER_CONVERSION

namespace plll
{
    namespace arithmetic
    {
        namespace implementation
        {
            // Implementation for NFPContext<Type>
            
            template<class Type, class IType> class conversion_impl<NFP<Type>, NIntContext<IType> >
            {
            public:
                typedef NInt<IType> RetVal;
                typedef void RetVal_Frac;
                typedef NInt<IType> RetVal_Floor;
                typedef NInt<IType> RetVal_Round;
                typedef NInt<IType> RetVal_Round2;
                typedef NInt<IType> RetVal_Ceil;
        
                static void convert(NInt<IType> & d, const NFP<Type> & v, const NIntContext<IType> & c)
                {
                    d.d_value = arithmetic::convert<IType>(v);
                }
        
                static RetVal convert(const NFP<Type> & v, const NIntContext<IType> & c)
                {
                    return NInt<IType>(false, arithmetic::convert<IType>(v));
                }
        
                static void floor(NInt<IType> & d, const NFP<Type> & v, const NIntContext<IType> & c)
                {
                    d.d_value = arithmetic::convert_floor<IType>(v);
                }
        
                static RetVal_Floor floor(const NFP<Type> & v, const NIntContext<IType> & c)
                {
                    return NInt<IType>(false, arithmetic::convert_floor<IType>(v));
                }
        
                static void round(NInt<IType> & d, const NFP<Type> & v, const NIntContext<IType> & c)
                {
                    d.d_value = arithmetic::convert_round<IType>(v);
                }
        
                static RetVal_Round round(const NFP<Type> & v, const NIntContext<IType> & c)
                {
                    return NInt<IType>(false, arithmetic::convert_round<IType>(v));
                }
        
                static void round(NInt<IType> & d, const NFP<Type> & v, bool & up, const NIntContext<IType> & c)
                {
                    d.d_value = arithmetic::convert_round<IType>(v, up);
                }
        
                static RetVal_Round2 round(const NFP<Type> & v, bool & up, const NIntContext<IType> & c)
                {
                    return NInt<IType>(false, arithmetic::convert_round<IType>(v, up));
                }
        
                static void ceil(NInt<IType> & d, const NFP<Type> & v, const NIntContext<IType> & c)
                {
                    d.d_value = arithmetic::convert_ceil<IType>(v);
                }
        
                static RetVal_Ceil ceil(const NFP<Type> & v, const NIntContext<IType> & c)
                {
                    return NInt<IType>(false, arithmetic::convert_ceil<IType>(v));
                }
            };
            
            // To avoid that conversion_impl<NInt<Type>, NIntContext<Type> > is defined twice, we
            // instanciate these conversions explicitly for distinct types.
            template<> class conversion_impl<NFP<double>, NFPContext<float> >
            {
            public:
                typedef NFP<float> RetVal;
                typedef NFP<float> RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(NFP<float> & d, const NFP<double> & v, const NFPContext<float> & c)
                {
                    d.d_value = v.d_value;
                }
                
                static RetVal convert(const NFP<double> & v, const NFPContext<float> & c)
                {
                    return NFP<float>(false, v.d_value);
                }
                
                static void convert_frac(NFP<float> & d, const NFP<double> & v1, const NFP<double> & v2, const NFPContext<float> & c)
                {
                    d.d_value = v1.d_value / v2.d_value;
                }
                
                static RetVal_Frac convert_frac(const NFP<double> & v1, const NFP<double> & v2, const NFPContext<float> & c)
                {
                    return NFP<float>(false, v1.d_value / v2.d_value);
                }
            };
        
            template<> class conversion_impl<NFP<long double>, NFPContext<float> >
            {
            public:
                typedef NFP<float> RetVal;
                typedef NFP<float> RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
        
                static void convert(NFP<float> & d, const NFP<long double> & v, const NFPContext<float> & c)
                {
                    d.d_value = v.d_value;
                }
        
                static RetVal convert(const NFP<long double> & v, const NFPContext<float> & c)
                {
                    return NFP<float>(false, v.d_value);
                }
                
                static void convert_frac(NFP<float> & d, const NFP<long double> & v1, const NFP<long double> & v2, const NFPContext<float> & c)
                {
                    d.d_value = v1.d_value / v2.d_value;
                }
                
                static RetVal_Frac convert_frac(const NFP<long double> & v1, const NFP<long double> & v2, const NFPContext<float> & c)
                {
                    return NFP<float>(false, v1.d_value / v2.d_value);
                }
            };
        
            template<> class conversion_impl<NFP<float>, NFPContext<double> >
            {
            public:
                typedef NFP<double> RetVal;
                typedef NFP<double> RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
        
                static void convert(NFP<double> & d, const NFP<float> & v, const NFPContext<double> & c)
                {
                    d.d_value = v.d_value;
                }
        
                static RetVal convert(const NFP<float> & v, const NFPContext<double> & c)
                {
                    return NFP<double>(false, v.d_value);
                }
                
                static void convert_frac(NFP<double> & d, const NFP<float> & v1, const NFP<float> & v2, const NFPContext<double> & c)
                {
                    d.d_value = static_cast<double>(v1.d_value) / static_cast<double>(v2.d_value);
                }
                
                static RetVal_Frac convert_frac(const NFP<float> & v1, const NFP<float> & v2, const NFPContext<double> & c)
                {
                    return NFP<double>(false, static_cast<double>(v1.d_value) / static_cast<double>(v2.d_value));
                }
            };
        
            template<> class conversion_impl<NFP<long double>, NFPContext<double> >
            {
            public:
                typedef NFP<double> RetVal;
                typedef NFP<double> RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
        
                static void convert(NFP<double> & d, const NFP<long double> & v, const NFPContext<double> & c)
                {
                    d.d_value = v.d_value;
                }
        
                static RetVal convert(const NFP<long double> & v, const NFPContext<double> & c)
                {
                    return NFP<double>(false, v.d_value);
                }
                
                static void convert_frac(NFP<double> & d, const NFP<long double> & v1, const NFP<long double> & v2, const NFPContext<double> & c)
                {
                    d.d_value = v1.d_value / v2.d_value;
                }
                
                static RetVal_Frac convert_frac(const NFP<long double> & v1, const NFP<long double> & v2, const NFPContext<double> & c)
                {
                    return NFP<double>(false, v1.d_value / v2.d_value);
                }
            };
        
            template<> class conversion_impl<NFP<float>, NFPContext<long double> >
            {
            public:
                typedef NFP<long double> RetVal;
                typedef NFP<long double> RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
        
                static void convert(NFP<long double> & d, const NFP<float> & v, const NFPContext<long double> & c)
                {
                    d.d_value = v.d_value;
                }
        
                static RetVal convert(const NFP<float> & v, const NFPContext<long double> & c)
                {
                    return NFP<long double>(false, v.d_value);
                }
                
                static void convert_frac(NFP<long double> & d, const NFP<float> & v1, const NFP<float> & v2, const NFPContext<long double> & c)
                {
                    d.d_value = static_cast<long double>(v1.d_value) / static_cast<long double>(v2.d_value);
                }
                
                static RetVal_Frac convert_frac(const NFP<float> & v1, const NFP<float> & v2, const NFPContext<long double> & c)
                {
                    return NFP<long double>(false, static_cast<long double>(v1.d_value) / static_cast<long double>(v2.d_value));
                }
            };
        
            template<> class conversion_impl<NFP<double>, NFPContext<long double> >
            {
            public:
                typedef NFP<long double> RetVal;
                typedef NFP<long double> RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
        
                static void convert(NFP<long double> & d, const NFP<double> & v, const NFPContext<long double> & c)
                {
                    d.d_value = v.d_value;
                }
        
                static RetVal convert(const NFP<double> & v, const NFPContext<long double> & c)
                {
                    return NFP<long double>(false, v.d_value);
                }
                
                static void convert_frac(NFP<long double> & d, const NFP<double> & v1, const NFP<double> & v2, const NFPContext<long double> & c)
                {
                    d.d_value = static_cast<long double>(v1.d_value) / static_cast<long double>(v2.d_value);
                }
                
                static RetVal_Frac convert_frac(const NFP<double> & v1, const NFP<double> & v2, const NFPContext<long double> & c)
                {
                    return NFP<long double>(false, static_cast<long double>(v1.d_value) / static_cast<long double>(v2.d_value));
                }
            };
            
            template<class Type, class IType> class conversion_impl<NInt<IType>, NFPContext<Type> >
            {
            public:
                typedef NFP<Type> RetVal;
                typedef NFP<Type> RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
        
                static void convert(NFP<Type> & d, const NInt<IType> & v, const NFPContext<Type> & c)
                {
                    d.d_value = arithmetic::convert<Type>(v);
                }
        
                static RetVal convert(const NInt<IType> & v, const NFPContext<Type> & c)
                {
                    return NFP<Type>(false, arithmetic::convert<Type>(v));
                }
                
                static void convert_frac(NFP<Type> & d, const NInt<IType> & v1, const NInt<IType> & v2, const NFPContext<Type> & c)
                {
                    d.d_value = static_cast<Type>(v1.d_value) / static_cast<Type>(v2.d_value);
                }
                
                static RetVal_Frac convert_frac(const NInt<IType> & v1, const NInt<IType> & v2, const NFPContext<Type> & c)
                {
                    return NFP<Type>(false, static_cast<Type>(v1.d_value) / static_cast<Type>(v2.d_value));
                }
            };
            
            template<class Type> class conversion_impl<float, NFPContext<Type> >
            {
            public:
                typedef NFP<Type> RetVal;
                typedef NFP<Type> RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
        
                static void convert(NFP<Type> & d, float v, const NFPContext<Type> & c)
                {
                    d.d_value = v;
                }
        
                static RetVal convert(float v, const NFPContext<Type> & c)
                {
                    return NFP<Type>(false, v);
                }
                
                static void convert_frac(NFP<Type> & d, float v1, float v2, const NFPContext<Type> & c)
                {
                    d.d_value = static_cast<Type>(v1) / static_cast<Type>(v2);
                }
                
                static RetVal_Frac convert_frac(float v1, float v2, const NFPContext<Type> & c)
                {
                    return NFP<Type>(false, static_cast<Type>(v1) / static_cast<Type>(v2));
                }
            };

            template<class Type> class conversion_impl<double, NFPContext<Type> >
            {
            public:
                typedef NFP<Type> RetVal;
                typedef NFP<Type> RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
        
                static void convert(NFP<Type> & d, double v, const NFPContext<Type> & c)
                {
                    d.d_value = v;
                }
        
                static RetVal convert(double v, const NFPContext<Type> & c)
                {
                    return NFP<Type>(false, v);
                }
                
                static void convert_frac(NFP<Type> & d, double v1, double v2, const NFPContext<Type> & c)
                {
                    d.d_value = static_cast<Type>(v1) / static_cast<Type>(v2);
                }
                
                static RetVal_Frac convert_frac(double v1, double v2, const NFPContext<Type> & c)
                {
                    return NFP<Type>(false, static_cast<Type>(v1) / static_cast<Type>(v2));
                }
            };

            template<class Type> class conversion_impl<long double, NFPContext<Type> >
            {
            public:
                typedef NFP<Type> RetVal;
                typedef NFP<Type> RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
        
                static void convert(NFP<Type> & d, long double v, const NFPContext<Type> & c)
                {
                    d.d_value = v;
                }
        
                static RetVal convert(long double v, const NFPContext<Type> & c)
                {
                    return NFP<Type>(false, v);
                }
                
                static void convert_frac(NFP<Type> & d, long double v1, long double v2, const NFPContext<Type> & c)
                {
                    d.d_value = static_cast<Type>(v1) / static_cast<Type>(v2);
                }
                
                static RetVal_Frac convert_frac(long double v1, long double v2, const NFPContext<Type> & c)
                {
                    return NFP<Type>(false, static_cast<Type>(v1) / static_cast<Type>(v2));
                }
            };

            template<class Type> class conversion_impl<signed int, NFPContext<Type> >
            {
            public:
                typedef NFP<Type> RetVal;
                typedef NFP<Type> RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(NFP<Type> & d, signed int v, const NFPContext<Type> & c)
                {
                    d.d_value = v;
                }
                
                static RetVal convert(signed int v, const NFPContext<Type> & c)
                {
                    return NFP<Type>(false, v);
                }
                
                static void convert_frac(NFP<Type> & d, signed int v1, signed int v2, const NFPContext<Type> & c)
                {
                    d.d_value = static_cast<Type>(v1) / static_cast<Type>(v2);
                }
                
                static RetVal_Frac convert_frac(signed int v1, signed int v2, const NFPContext<Type> & c)
                {
                    return NFP<Type>(false, static_cast<Type>(v1) / static_cast<Type>(v2));
                }
            };

            template<class Type> class conversion_impl<unsigned int, NFPContext<Type> >
            {
            public:
                typedef NFP<Type> RetVal;
                typedef NFP<Type> RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
        
                static void convert(NFP<Type> & d, unsigned int v, const NFPContext<Type> & c)
                {
                    d.d_value = v;
                }
        
                static RetVal convert(unsigned int v, const NFPContext<Type> & c)
                {
                    return NFP<Type>(false, v);
                }
                
                static void convert_frac(NFP<Type> & d, unsigned int v1, unsigned int v2, const NFPContext<Type> & c)
                {
                    d.d_value = static_cast<Type>(v1) / static_cast<Type>(v2);
                }
                
                static RetVal_Frac convert_frac(unsigned int v1, unsigned int v2, const NFPContext<Type> & c)
                {
                    return NFP<Type>(false, static_cast<Type>(v1) / static_cast<Type>(v2));
                }
            };

            template<class Type> class conversion_impl<signed long, NFPContext<Type> >
            {
            public:
                typedef NFP<Type> RetVal;
                typedef NFP<Type> RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
        
                static void convert(NFP<Type> & d, signed long v, const NFPContext<Type> & c)
                {
                    d.d_value = v;
                }
        
                static RetVal convert(signed long v, const NFPContext<Type> & c)
                {
                    return NFP<Type>(false, v);
                }
                
                static void convert_frac(NFP<Type> & d, signed long v1, signed long v2, const NFPContext<Type> & c)
                {
                    d.d_value = static_cast<Type>(v1) / static_cast<Type>(v2);
                }
                
                static RetVal_Frac convert_frac(signed long v1, signed long v2, const NFPContext<Type> & c)
                {
                    return NFP<Type>(false, static_cast<Type>(v1) / static_cast<Type>(v2));
                }
            };

            template<class Type> class conversion_impl<unsigned long, NFPContext<Type> >
            {
            public:
                typedef NFP<Type> RetVal;
                typedef NFP<Type> RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
        
                static void convert(NFP<Type> & d, unsigned long v, const NFPContext<Type> & c)
                {
                    d.d_value = v;
                }
        
                static RetVal convert(unsigned long v, const NFPContext<Type> & c)
                {
                    return NFP<Type>(false, v);
                }
                
                static void convert_frac(NFP<Type> & d, unsigned long v1, unsigned long v2, const NFPContext<Type> & c)
                {
                    d.d_value = static_cast<Type>(v1) / static_cast<Type>(v2);
                }
                
                static RetVal_Frac convert_frac(unsigned long v1, unsigned long v2, const NFPContext<Type> & c)
                {
                    return NFP<Type>(false, static_cast<Type>(v1) / static_cast<Type>(v2));
                }
            };

            template<class Type> class conversion_impl<long long, NFPContext<Type> >
            {
            public:
                typedef NFP<Type> RetVal;
                typedef NFP<Type> RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
        
                static void convert(NFP<Type> & d, long long v, const NFPContext<Type> & c)
                {
                    d.d_value = v;
                }
        
                static RetVal convert(long long v, const NFPContext<Type> & c)
                {
                    return NFP<Type>(false, v);
                }
                
                static void convert_frac(NFP<Type> & d, long long v1, long long v2, const NFPContext<Type> & c)
                {
                    d.d_value = static_cast<Type>(v1) / static_cast<Type>(v2);
                }
                
                static RetVal_Frac convert_frac(long long v1, long long v2, const NFPContext<Type> & c)
                {
                    return NFP<Type>(false, static_cast<Type>(v1) / static_cast<Type>(v2));
                }
            };
            
            template<class Type> class conversion_impl<Rational, NFPContext<Type> >
            {
            public:
                typedef NFP<Type> RetVal;
                typedef NFP<Type> RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
        
                static void convert(NFP<Type> & res, const Rational & r, const NFPContext<Type> & rc)
                {
                    arithmetic::convert_fraction(res, r.numerator(), r.denominator(), rc);
                }
                
                inline static NFP<Type> convert(const Rational & r, const NFPContext<Type> & rc)
                {
                    NFP<Type> res(rc);
                    convert(res, r, rc);
                    return res;
                }
                
                inline static void convert_frac(NFP<Type> & res, const Rational & r1, const Rational & r2, const NFPContext<Type> & rc)
                {
                    long nl1 = ceilOfLog2(r1.numerator()  ) - NFPContext<Type>::d_precision - 8,
                         dl1 = ceilOfLog2(r1.denominator()) - NFPContext<Type>::d_precision - 8;
                    long nl2 = ceilOfLog2(r2.numerator()  ) - NFPContext<Type>::d_precision - 8,
                         dl2 = ceilOfLog2(r2.denominator()) - NFPContext<Type>::d_precision - 8;
                    if (nl1 < 0) nl1 = 0;
                    if (nl2 < 0) nl2 = 0;
                    if (dl1 < 0) dl1 = 0;
                    if (dl2 < 0) dl2 = 0;
                    res.d_value = ldexp(IntRealConversion<Type>::intToType(r1.numerator() >> nl1) *
                                        IntRealConversion<Type>::intToType(r1.numerator() >> dl2) /
                                        IntRealConversion<Type>::intToType(r1.denominator() >> dl1) /
                                        IntRealConversion<Type>::intToType(r1.denominator() >> nl2), nl1 + dl2 - dl1 - nl2);
                }
                
                inline static NFP<Type> convert_frac(const Rational & r1, const Rational & r2, const NFPContext<Type> & rc)
                {
                    NFP<Type> res(rc);
                    convert_frac(res, r1, r2, rc);
                    return res;
                }
            };
            
            template<> class conversion_impl<Rational, NFPContext<long double> >
            {
            public:
                typedef NFP<long double> RetVal;
                typedef NFP<long double> RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                inline static void convert(NFP<long double> & res, const Rational & r, const NFPContext<long double> & rc)
                {
                    arithmetic::convert_fraction(res, r.numerator(), r.denominator(), rc);
                }
                
                inline static NFP<long double> convert(const Rational & r, const NFPContext<long double> & rc)
                {
                    NFP<long double> res(rc);
                    convert(res, r, rc);
                    return res;
                }
                
                inline static void convert_frac(NFP<long double> & res, const Rational & r1, const Rational & r2, const NFPContext<long double> & rc)
                {
                    long nl1 = ceilOfLog2(r1.numerator()  ) - NFPContext<long double>::d_precision - 8,
                         dl1 = ceilOfLog2(r1.denominator()) - NFPContext<long double>::d_precision - 8;
                    long nl2 = ceilOfLog2(r2.numerator()  ) - NFPContext<long double>::d_precision - 8,
                         dl2 = ceilOfLog2(r2.denominator()) - NFPContext<long double>::d_precision - 8;
                    if (nl1 < 0) nl1 = 0;
                    if (nl2 < 0) nl2 = 0;
                    if (dl1 < 0) dl1 = 0;
                    if (dl2 < 0) dl2 = 0;
                    res.d_value = ldexpl(IntRealConversion<long double>::intToType(r1.numerator() >> nl1) *
                                         IntRealConversion<long double>::intToType(r1.numerator() >> dl2) /
                                         IntRealConversion<long double>::intToType(r1.denominator() >> dl1) /
                                         IntRealConversion<long double>::intToType(r1.denominator() >> nl2), nl1 + dl2 - dl1 - nl2);
                }
                
                inline static NFP<long double> convert_frac(const Rational & r1, const Rational & r2, const NFPContext<long double> & rc)
                {
                    NFP<long double> res(rc);
                    convert_frac(res, r1, r2, rc);
                    return res;
                }
            };
            
            template<class Type> class conversion_impl<Real, NFPContext<Type> >
            {
            public:
                typedef NFP<Type> RetVal;
                typedef NFP<Type> RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(NFP<Type> & res, const Real & r, const NFPContext<Type> & rc)
                {
                    res.d_value = IntRealConversion<Type>::realToType(r);
                }
                
                inline static NFP<Type> convert(const Real & r, const NFPContext<Type> & rc)
                {
                    return NFP<Type>(false, IntRealConversion<Type>::realToType(r));
                }
                
                static void convert_frac(NFP<Type> & res, const Real & r1, const Real & r2, const NFPContext<Type> & rc)
                {
                    res.d_value = IntRealConversion<Type>::realToType(r1) / IntRealConversion<Type>::realToType(r2);
                }
                
                inline static NFP<Type> convert_frac(const Real & r1, const Real & r2, const NFPContext<Type> & rc)
                {
                    return NFP<Type>(false, IntRealConversion<Type>::realToType(r1) / IntRealConversion<Type>::realToType(r2));
                }
            };
            
            template<class Type> class conversion_impl<Integer, NFPContext<Type> >
            {
            public:
                typedef NFP<Type> RetVal;
                typedef NFP<Type> RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(NFP<Type> & res, const Integer & r, const NFPContext<Type> & rc)
                {
                    res.d_value = IntRealConversion<Type>::intToType(r);
                }
                
                inline static NFP<Type> convert(const Integer & r, const NFPContext<Type> & rc)
                {
                    return NFP<Type>(false, IntRealConversion<Type>::intToType(r));
                }
                
                inline static void convert_frac(NFP<Type> & res, const Integer & n, const Integer & d, const NFPContext<Type> & rc)
                {
                    long nl = ceilOfLog2(n) - NFPContext<Type>::d_precision - 8, dl = ceilOfLog2(d) - NFPContext<Type>::d_precision - 8;
                    if (nl < 0) nl = 0;
                    if (dl < 0) dl = 0;
                    res.d_value = ldexp(IntRealConversion<Type>::intToType(n >> nl) / IntRealConversion<Type>::intToType(d >> dl), nl - dl);
                }
                
                inline static NFP<Type> convert_frac(const Integer & n, const Integer & d, const NFPContext<Type> & rc)
                {
                    NFP<Type> res(rc);
                    convert_frac(res, n, d, rc);
                    return res;
                }
            };
            
            template<> class conversion_impl<Integer, NFPContext<long double> >
            // In C++98, frexp and ldexp do not have to be overloaded for long doubles.
            {
            public:
                typedef NFP<long double> RetVal;
                typedef NFP<long double> RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(NFP<long double> & res, const Integer & r, const NFPContext<long double> & rc)
                {
                    res.d_value = IntRealConversion<long double>::intToType(r);
                }
                
                inline static NFP<long double> convert(const Integer & r, const NFPContext<long double> & rc)
                {
                    return NFP<long double>(false, IntRealConversion<long double>::intToType(r));
                }
                
                inline static void convert_frac(NFP<long double> & res, const Integer & n, const Integer & d, const NFPContext<long double> & rc)
                {
                    long nl = ceilOfLog2(n) - NFPContext<long double>::d_precision - 8, dl = ceilOfLog2(d) - NFPContext<long double>::d_precision - 8;
                    if (nl < 0) nl = 0;
                    if (dl < 0) dl = 0;
                    res.d_value = ldexpl(IntRealConversion<long double>::intToType(n >> nl) / IntRealConversion<long double>::intToType(d >> dl), nl - dl);
                }
                
                inline static NFP<long double> convert_frac(const Integer & n, const Integer & d, const NFPContext<long double> & rc)
                {
                    NFP<long double> res(rc);
                    convert_frac(res, n, d, rc);
                    return res;
                }
            };
            
            template<class Type> class conversion_impl<NFP<Type>, RealContext>
            {
            public:
                typedef Real RetVal;
                typedef Real RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(Real & res, const NFP<Type> & r, const RealContext &)
                {
                    IntRealConversion<Type>::typeToReal(res, r.d_value);
                }
                
                inline static Real convert(const NFP<Type> & r, const RealContext & rc)
                {
                    Real res(rc);
                    IntRealConversion<Type>::typeToReal(res, r.d_value);
                    return res;
                }
                
                static void convert_frac(Real & res, const NFP<Type> & r1, const NFP<Type> & r2, const RealContext & rc)
                {
                    IntRealConversion<Type>::typeToReal(res, r1.d_value);
                    Real tmp(rc);
                    IntRealConversion<Type>::typeToReal(tmp, r2.d_value);
                    res /= tmp;
                }
                
                inline static Real convert_frac(const NFP<Type> & r1, const NFP<Type> & r2, const RealContext & rc)
                {
                    Real res(rc), tmp(rc);
                    IntRealConversion<Type>::typeToReal(res, r1.d_value);
                    IntRealConversion<Type>::typeToReal(tmp, r2.d_value);
                    return res /= tmp;
                }
            };
            
            template<class Type> class conversion_impl<NFP<Type>, IntegerContext>
            {
            public:
                typedef Integer RetVal;
                typedef void RetVal_Frac;
                typedef Integer RetVal_Floor;
                typedef Integer RetVal_Round;
                typedef Integer RetVal_Round2;
                typedef Integer RetVal_Ceil;
        
                static void convert(Integer & res, const NFP<Type> & r, const IntegerContext &)
                {
                    IntRealConversion<Type>::typeToInt(res, r.d_value);
                }
        
                inline static RetVal convert(const NFP<Type> & r, const IntegerContext & ic)
                {
                    Integer res(ic);
                    IntRealConversion<Type>::typeToInt(res, r.d_value);
                    return res;
                }
        
                static void floor(Integer & res, const NFP<Type> & r, const IntegerContext &)
                {
                    IntRealConversion<Type>::typeToInt(res, std::floor(r.d_value));
                }
        
                inline static RetVal_Floor floor(const NFP<Type> & r, const IntegerContext & ic)
                {
                    Integer res(ic);
                    IntRealConversion<Type>::typeToInt(res, std::floor(r.d_value));
                    return res;
                }
        
                static void round(Integer & res, const NFP<Type> & r, const IntegerContext &)
                {
                    IntRealConversion<Type>::typeToInt(res, implementation::RoundingNFP<Type>::round(r.d_value));
                }
        
                inline static RetVal_Round round(const NFP<Type> & r, const IntegerContext & ic)
                {
                    Integer res(ic);
                    IntRealConversion<Type>::typeToInt(res, implementation::RoundingNFP<Type>::round(r.d_value));
                    return res;
                }
        
                static void round(Integer & res, const NFP<Type> & r, bool & up, const IntegerContext &)
                {
                    IntRealConversion<Type>::typeToInt(res, implementation::RoundingNFP<Type>::round(r.d_value, up));
                }
        
                inline static RetVal_Round2 round(const NFP<Type> & r, bool & up, const IntegerContext & ic)
                {
                    Integer res(ic);
                    IntRealConversion<Type>::typeToInt(res, implementation::RoundingNFP<Type>::round(r.d_value, up));
                    return res;
                }
        
                static void ceil(Integer & res, const NFP<Type> & r, const IntegerContext &)
                {
                    IntRealConversion<Type>::typeToInt(res, std::ceil(r.d_value));
                }
        
                inline static RetVal_Ceil ceil(const NFP<Type> & r, const IntegerContext & ic)
                {
                    Integer res(ic);
                    IntRealConversion<Type>::typeToInt(res, std::ceil(r.d_value));
                    return res;
                }
            };
            
            template<class Type> class nativeconversion_impl<NFP<Type> >
            {
            public:
                static int toInt(const NFP<Type> & v)
                {
                    return v.d_value;
                }
        
                static int toInt_Floor(const NFP<Type> & v)
                {
                    return std::floor(v.d_value);
                }
        
                static int toInt_Round(const NFP<Type> & v)
                {
                    return implementation::RoundingNFP<Type>::round(v.d_value);
                }
        
                static int toInt_Round(const NFP<Type> & v, bool & up)
                {
                    int r = implementation::RoundingNFP<Type>::round(v.d_value);
                    up = v.d_value < r;
                    return r;
                }
        
                static int toInt_Ceil(const NFP<Type> & v)
                {
                    return std::ceil(v.d_value);
                }
            
                static unsigned int toUInt(const NFP<Type> & v)
                {
                    return v.d_value;
                }
        
                static unsigned int toUInt_Floor(const NFP<Type> & v)
                {
                    return std::floor(v.d_value);
                }
        
                static unsigned int toUInt_Round(const NFP<Type> & v)
                {
                    return implementation::RoundingNFP<Type>::round(v.d_value);
                }
        
                static unsigned int toUInt_Round(const NFP<Type> & v, bool & up)
                {
                    unsigned int r = implementation::RoundingNFP<Type>::round(v.d_value);
                    up = v.d_value < r;
                    return r;
                }
        
                static unsigned int toUInt_Ceil(const NFP<Type> & v)
                {
                    return std::ceil(v.d_value);
                }
        
                static long toLong(const NFP<Type> & v)
                {
                    return v.d_value;
                }
        
                static long toLong_Floor(const NFP<Type> & v)
                {
                    return std::floor(v.d_value);
                }
        
                static long toLong_Round(const NFP<Type> & v)
                {
                    return implementation::RoundingNFP<Type>::round(v.d_value);
                }
        
                static long toLong_Round(const NFP<Type> & v, bool & up)
                {
                    long r = implementation::RoundingNFP<Type>::round(v.d_value);
                    up = v.d_value < r;
                    return r;
                }
        
                static long toLong_Ceil(const NFP<Type> & v)
                {
                    return std::ceil(v.d_value);
                }
        
                static unsigned long toULong(const NFP<Type> & v)
                {
                    return v.d_value;
                }
        
                static unsigned long toULong_Floor(const NFP<Type> & v)
                {
                    return std::floor(v.d_value);
                }
        
                static unsigned long toULong_Round(const NFP<Type> & v)
                {
                    return implementation::RoundingNFP<Type>::round(v.d_value);
                }
        
                static unsigned long toULong_Round(const NFP<Type> & v, bool & up)
                {
                    unsigned long r = implementation::RoundingNFP<Type>::round(v.d_value);
                    up = v.d_value < r;
                    return r;
                }
        
                static unsigned long toULong_Ceil(const NFP<Type> & v)
                {
                    return std::ceil(v.d_value);
                }
        
                static long long toLongLong(const NFP<Type> & v)
                {
                    return v.d_value;
                }
        
                static long long toLongLong_Floor(const NFP<Type> & v)
                {
                    return std::floor(v.d_value);
                }
        
                static long long toLongLong_Round(const NFP<Type> & v)
                {
                    return implementation::RoundingNFP<Type>::round(v.d_value);
                }
        
                static long long toLongLong_Round(const NFP<Type> & v, bool & up)
                {
                    long long r = implementation::RoundingNFP<Type>::round(v.d_value);
                    up = v.d_value < r;
                    return r;
                }
        
                static long long toLongLong_Ceil(const NFP<Type> & v)
                {
                    return std::ceil(v.d_value);
                }
        
                static float toFloat(const NFP<Type> & v)
                {
                    return v.d_value;
                }
        
                static double toDouble(const NFP<Type> & v)
                {
                    return v.d_value;
                }
        
                static long double toLongDouble(const NFP<Type> & v)
                {
                    return v.d_value;
                }
            };
            
            template<typename Type, typename Op>
            struct binary_operation_impl<NFP<Type>, NFP<Type>, Op>
            {
                enum { supported = true, intermediate_expression = false };
                typedef NFP<Type> IntermediateType;
                typedef NFP<Type> ResultType;
            };
            
            template<typename Type>
            struct unary_operation_impl<NFP<Type>, arithmetic::op::negation>
            {
                enum { supported = true, intermediate_expression = false };
                typedef NFP<Type> IntermediateType;
                typedef NFP<Type> ResultType;
            };
        }
    }
}

#endif

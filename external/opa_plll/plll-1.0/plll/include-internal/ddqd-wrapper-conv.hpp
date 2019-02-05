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

#ifndef PLLL_INCLUDE_GUARD__DDQD_WRAPPER_CONVERSION
#define PLLL_INCLUDE_GUARD__DDQD_WRAPPER_CONVERSION

namespace plll
{
    namespace arithmetic
    {
        namespace implementation
        {
            // Implementation for DDQDContext<Type>
            
            template<class Type, class IType> class conversion_impl<NInt<IType>, DDQDContext<Type> >
            {
            public:
                typedef DDQD<Type> RetVal;
                typedef DDQD<Type> RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(DDQD<Type> & d, const NInt<IType> & v, const DDQDContext<Type> & c)
                {
                    arithmetic::convert(d, arithmetic::convert<IType>(v), c);
                }
                
                static RetVal convert(const NInt<IType> & v, const DDQDContext<Type> & c)
                {
                    return arithmetic::convert(arithmetic::convert<IType>(v), c);
                }
                
                static void convert_frac(DDQD<Type> & d, const NInt<IType> & v1, const NInt<IType> & v2, const DDQDContext<Type> & c)
                {
                    arithmetic::convert_fraction(d, arithmetic::convert<IType>(v1), arithmetic::convert<IType>(v2), c);
                }
                
                static RetVal_Frac convert_frac(const NInt<IType> & v1, const NInt<IType> & v2, const DDQDContext<Type> & c)
                {
                    return arithmetic::convert_fraction(arithmetic::convert<IType>(v1), arithmetic::convert<IType>(v2), c);
                }
            };
            
            template<class Type, class IType> class conversion_impl<DDQD<Type>, NIntContext<IType> >
            {
            public:
                typedef NInt<IType> RetVal;
                typedef void RetVal_Frac;
                typedef NInt<IType> RetVal_Floor;
                typedef NInt<IType> RetVal_Round;
                typedef NInt<IType> RetVal_Round2;
                typedef NInt<IType> RetVal_Ceil;
                
                static void convert(NInt<IType> & d, const DDQD<Type> & v, const NIntContext<IType> & c)
                {
                    arithmetic::convert(d, arithmetic::convert<IType>(v), c);
                }
                
                static RetVal convert(const DDQD<Type> & v, const NIntContext<IType> & c)
                {
                    return arithmetic::convert(arithmetic::convert<IType>(v), c);
                }
                
                static void floor(NInt<IType> & d, const DDQD<Type> & v, const NIntContext<IType> & c)
                {
                    arithmetic::convert(d, arithmetic::convert_floor<IType>(v), c);
                }
                
                static RetVal_Floor floor(const DDQD<Type> & v, const NIntContext<IType> & c)
                {
                    return arithmetic::convert(arithmetic::convert_floor<IType>(v), c);
                }
                
                static void round(NInt<IType> & d, const DDQD<Type> & v, const NIntContext<IType> & c)
                {
                    arithmetic::convert(d, arithmetic::convert_round<IType>(v), c);
                }
                
                static RetVal_Round round(const DDQD<Type> & v, const NIntContext<IType> & c)
                {
                    return arithmetic::convert(arithmetic::convert_round<IType>(v), c);
                }
                
                static void round(NInt<IType> & d, const DDQD<Type> & v, bool & up, const NIntContext<IType> & c)
                {
                    arithmetic::convert(d, arithmetic::convert_round<IType>(v, up), c);
                }
                
                static RetVal_Round2 round(const DDQD<Type> & v, bool & up, const NIntContext<IType> & c)
                {
                    return arithmetic::convert(arithmetic::convert_round<IType>(v, up), c);
                }
                
                static void ceil(NInt<IType> & d, const DDQD<Type> & v, const NIntContext<IType> & c)
                {
                    arithmetic::convert(d, arithmetic::convert_ceil<IType>(v), c);
                }
                
                static RetVal_Ceil ceil(const DDQD<Type> & v, const NIntContext<IType> & c)
                {
                    return arithmetic::convert(arithmetic::convert_ceil<IType>(v), c);
                }
            };
            
            template<class Type> class conversion_impl<float, DDQDContext<Type> >
            {
            public:
                typedef DDQD<Type> RetVal;
                typedef DDQD<Type> RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(DDQD<Type> & d, float v, const DDQDContext<Type> & c)
                {
                    d.d_value = IntRealConversion<Type>::doubleToType(v);
                }
                
                static RetVal convert(float v, const DDQDContext<Type> & c)
                {
                    return DDQD<Type>(false, IntRealConversion<Type>::doubleToType(v));
                }
                
                static void convert_frac(DDQD<Type> & d, float v1, float v2, const DDQDContext<Type> & c)
                {
                    d.d_value = IntRealConversion<Type>::doubleToType(v1) / IntRealConversion<Type>::doubleToType(v2);
                }
                
                static RetVal convert_frac(float v1, float v2, const DDQDContext<Type> & c)
                {
                    return DDQD<Type>(false, IntRealConversion<Type>::doubleToType(v1) / IntRealConversion<Type>::doubleToType(v2));
                }
            };
            
            template<class Type> class conversion_impl<double, DDQDContext<Type> >
            {
            public:
                typedef DDQD<Type> RetVal;
                typedef DDQD<Type> RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(DDQD<Type> & d, double v, const DDQDContext<Type> & c)
                {
                    d.d_value = IntRealConversion<Type>::doubleToType(v);
                }
                
                static RetVal convert(double v, const DDQDContext<Type> & c)
                {
                    return DDQD<Type>(false, IntRealConversion<Type>::doubleToType(v));
                }
                
                static void convert_frac(DDQD<Type> & d, double v1, double v2, const DDQDContext<Type> & c)
                {
                    d.d_value = IntRealConversion<Type>::doubleToType(v1) / IntRealConversion<Type>::doubleToType(v2);
                }
                
                static RetVal convert_frac(double v1, double v2, const DDQDContext<Type> & c)
                {
                    return DDQD<Type>(false, IntRealConversion<Type>::doubleToType(v1) / IntRealConversion<Type>::doubleToType(v2));
                }
            };
            
            template<class Type> class conversion_impl<long double, DDQDContext<Type> >
            {
            public:
                typedef DDQD<Type> RetVal;
                typedef DDQD<Type> RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(DDQD<Type> & d, long double v, const DDQDContext<Type> & c)
                {
                    d.d_value = IntRealConversion<Type>::lDoubleToType(v);
                }
                
                static RetVal convert(long double v, const DDQDContext<Type> & c)
                {
                    return DDQD<Type>(false, IntRealConversion<Type>::lDoubleToType(v));
                }
                
                static void convert_frac(DDQD<Type> & d, long double v1, long double v2, const DDQDContext<Type> & c)
                {
                    d.d_value = IntRealConversion<Type>::lDoubleToType(v1) / IntRealConversion<Type>::lDoubleToType(v2);
                }
                
                static RetVal convert_frac(long double v1, long double v2, const DDQDContext<Type> & c)
                {
                    return DDQD<Type>(false, IntRealConversion<Type>::lDoubleToType(v1) / IntRealConversion<Type>::lDoubleToType(v2));
                }
            };
            
            template<class Type> class conversion_impl<signed int, DDQDContext<Type> >
            {
            public:
                typedef DDQD<Type> RetVal;
                typedef DDQD<Type> RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(DDQD<Type> & d, signed int v, const DDQDContext<Type> & c)
                {
                    d.d_value = IntRealConversion<Type>::longintToType(v);
                }
                
                static RetVal convert(signed int v, const DDQDContext<Type> & c)
                {
                    return DDQD<Type>(false, IntRealConversion<Type>::longintToType(v));
                }
                
                static void convert_frac(DDQD<Type> & d, signed int v1, signed int v2, const DDQDContext<Type> & c)
                {
                    d.d_value = IntRealConversion<Type>::longintToType(v1) / IntRealConversion<Type>::longintToType(v2);
                }
                
                static RetVal convert_frac(signed int v1, signed int v2, const DDQDContext<Type> & c)
                {
                    return DDQD<Type>(false, IntRealConversion<Type>::longintToType(v1) / IntRealConversion<Type>::longintToType(v2));
                }
            };
            
            template<class Type> class conversion_impl<unsigned int, DDQDContext<Type> >
            {
            public:
                typedef DDQD<Type> RetVal;
                typedef DDQD<Type> RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(DDQD<Type> & d, unsigned int v, const DDQDContext<Type> & c)
                {
                    d.d_value = IntRealConversion<Type>::ulongintToType(v);
                }
                
                static RetVal convert(unsigned int v, const DDQDContext<Type> & c)
                {
                    return DDQD<Type>(false, IntRealConversion<Type>::ulongintToType(v));
                }
                
                static void convert_frac(DDQD<Type> & d, unsigned int v1, unsigned int v2, const DDQDContext<Type> & c)
                {
                    d.d_value = IntRealConversion<Type>::ulongintToType(v1) / IntRealConversion<Type>::ulongintToType(v2);
                }
                
                static RetVal convert_frac(unsigned int v1, unsigned int v2, const DDQDContext<Type> & c)
                {
                    return DDQD<Type>(false, IntRealConversion<Type>::ulongintToType(v1) / IntRealConversion<Type>::ulongintToType(v2));
                }
            };
            
            template<class Type> class conversion_impl<signed long, DDQDContext<Type> >
            {
            public:
                typedef DDQD<Type> RetVal;
                typedef DDQD<Type> RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(DDQD<Type> & d, signed long v, const DDQDContext<Type> & c)
                {
                    d.d_value = IntRealConversion<Type>::longintToType(v);
                }
                
                static RetVal convert(signed long v, const DDQDContext<Type> & c)
                {
                    return DDQD<Type>(false, IntRealConversion<Type>::longintToType(v));
                }
                
                static void convert_frac(DDQD<Type> & d, signed long v1, signed long v2, const DDQDContext<Type> & c)
                {
                    d.d_value = IntRealConversion<Type>::longintToType(v1) / IntRealConversion<Type>::longintToType(v2);
                }
                
                static RetVal convert_frac(signed long v1, signed long v2, const DDQDContext<Type> & c)
                {
                    return DDQD<Type>(false, IntRealConversion<Type>::longintToType(v1) / IntRealConversion<Type>::longintToType(v2));
                }
            };
            
            template<class Type> class conversion_impl<unsigned long, DDQDContext<Type> >
            {
            public:
                typedef DDQD<Type> RetVal;
                typedef DDQD<Type> RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
        
                static void convert(DDQD<Type> & d, unsigned long v, const DDQDContext<Type> & c)
                {
                    d.d_value = IntRealConversion<Type>::ulongintToType(v);
                }
        
                static RetVal convert(unsigned long v, const DDQDContext<Type> & c)
                {
                    return DDQD<Type>(false, IntRealConversion<Type>::ulongintToType(v));
                }
                
                static void convert_frac(DDQD<Type> & d, unsigned long v1, unsigned long v2, const DDQDContext<Type> & c)
                {
                    d.d_value = IntRealConversion<Type>::ulongintToType(v1) / IntRealConversion<Type>::ulongintToType(v2);
                }
                
                static RetVal convert_frac(unsigned long v1, unsigned long v2, const DDQDContext<Type> & c)
                {
                    return DDQD<Type>(false, IntRealConversion<Type>::ulongintToType(v1) / IntRealConversion<Type>::ulongintToType(v2));
                }
            };
            
            template<class Type> class conversion_impl<long long, DDQDContext<Type> >
            {
            public:
                typedef DDQD<Type> RetVal;
                typedef DDQD<Type> RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
        
                static void convert(DDQD<Type> & d, long long v, const DDQDContext<Type> & c)
                {
                    d.d_value = IntRealConversion<Type>::longlongToType(v);
                }
        
                static RetVal convert(long long v, const DDQDContext<Type> & c)
                {
                    return DDQD<Type>(false, IntRealConversion<Type>::longlongToType(v));
                }
                
                static void convert_frac(DDQD<Type> & d, long long v1, long long v2, const DDQDContext<Type> & c)
                {
                    d.d_value = IntRealConversion<Type>::longlongToType(v1) / IntRealConversion<Type>::longlongToType(v2);
                }
                
                static RetVal convert_frac(long long v1, long long v2, const DDQDContext<Type> & c)
                {
                    return DDQD<Type>(false, IntRealConversion<Type>::longlongToType(v1) / IntRealConversion<Type>::longlongToType(v2));
                }
            };
            
            template<class Type> class conversion_impl<Rational, DDQDContext<Type> >
            {
            public:
                typedef DDQD<Type> RetVal;
                typedef DDQD<Type> RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                inline static void convert(DDQD<Type> & res, const Rational & r, const DDQDContext<Type> & rc)
                {
                    arithmetic::convert_fraction(res, r.numerator(), r.denominator(), rc);
                }
                
                inline static DDQD<Type> convert(const Rational & r, const DDQDContext<Type> & rc)
                {
                    DDQD<Type> res(rc);
                    convert(res, r, rc);
                    return res;
                }
                
                inline static void convert_frac(DDQD<Type> & res, const Rational & r1, const Rational & r2, const DDQDContext<Type> & rc)
                {
                    long nl1 = ceilOfLog2(r1.numerator()  ) - DDQDContext<Type>::d_precision - 8,
                         dl1 = ceilOfLog2(r1.denominator()) - DDQDContext<Type>::d_precision - 8;
                    long nl2 = ceilOfLog2(r2.numerator()  ) - DDQDContext<Type>::d_precision - 8,
                         dl2 = ceilOfLog2(r2.denominator()) - DDQDContext<Type>::d_precision - 8;
                    if (nl1 < 0) nl1 = 0;
                    if (nl2 < 0) nl2 = 0;
                    if (dl1 < 0) dl1 = 0;
                    if (dl2 < 0) dl2 = 0;
                    res.d_value = ldexp(IntRealConversion<Type>::intToType(r1.numerator() >> nl1) *
                                        IntRealConversion<Type>::intToType(r1.numerator() >> dl2) /
                                        IntRealConversion<Type>::intToType(r1.denominator() >> dl1) /
                                        IntRealConversion<Type>::intToType(r1.denominator() >> nl2), nl1 + dl2 - dl1 - nl2);
                }
                
                inline static DDQD<Type> convert_frac(const Rational & r1, const Rational & r2, const DDQDContext<Type> & rc)
                {
                    DDQD<Type> res(rc);
                    convert_frac(res, r1, r2, rc);
                    return res;
                }
            };
            
            template<class Type> class conversion_impl<Real, DDQDContext<Type> >
            {
            public:
                typedef DDQD<Type> RetVal;
                typedef DDQD<Type> RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(DDQD<Type> & res, const Real & r, const DDQDContext<Type> & rc)
                {
                    res.d_value = IntRealConversion<Type>::realToType(r);
                }
                
                inline static DDQD<Type> convert(const Real & r, const DDQDContext<Type> & rc)
                {
                    return DDQD<Type>(false, IntRealConversion<Type>::realToType(r));
                }
                
                static void convert_frac(DDQD<Type> & res, const Real & r1, const Real & r2, const DDQDContext<Type> & rc)
                {
                    res.d_value = IntRealConversion<Type>::realToType(r1) / IntRealConversion<Type>::realToType(r2);
                }
                
                inline static DDQD<Type> convert_frac(const Real & r1, const Real & r2, const DDQDContext<Type> & rc)
                {
                    return DDQD<Type>(false, IntRealConversion<Type>::realToType(r1) / IntRealConversion<Type>::realToType(r2));
                }
            };
            
            template<class Type> class conversion_impl<Integer, DDQDContext<Type> >
            {
            public:
                typedef DDQD<Type> RetVal;
                typedef DDQD<Type> RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(DDQD<Type> & res, const Integer & r, const DDQDContext<Type> & rc)
                {
                    res.d_value = IntRealConversion<Type>::intToType(r);
                }
                
                inline static DDQD<Type> convert(const Integer & r, const DDQDContext<Type> & rc)
                {
                    return DDQD<Type>(false, IntRealConversion<Type>::intToType(r));
                }
                
                inline static void convert_frac(DDQD<Type> & res, const Integer & n, const Integer & d, const DDQDContext<Type> & rc)
                {
                    long nl = ceilOfLog2(n) - DDQDContext<Type>::d_precision - 8,
                         dl = ceilOfLog2(d) - DDQDContext<Type>::d_precision - 8;
                    if (nl < 0) nl = 0;
                    if (dl < 0) dl = 0;
                    res.d_value = ldexp(IntRealConversion<Type>::intToType(n >> nl) /
                                        IntRealConversion<Type>::intToType(d >> dl), nl - dl);
                }
                
                inline static DDQD<Type> convert_frac(const Integer & n, const Integer & d, const DDQDContext<Type> & rc)
                {
                    DDQD<Type> res(rc);
                    convert_frac(res, n, d, rc);
                    return res;
                }
            };
            
            template<class Type> class conversion_impl<DDQD<Type>, RealContext>
            {
            public:
                typedef Real RetVal;
                typedef Real RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(Real & res, const DDQD<Type> & r, const RealContext &)
                {
                    IntRealConversion<Type>::typeToReal(res, r.d_value);
                }
                
                inline static Real convert(const DDQD<Type> & r, const RealContext & rc)
                {
                    Real res(rc);
                    IntRealConversion<Type>::typeToReal(res, r.d_value);
                    return res;
                }
                
                static void convert_frac(Real & res, const DDQD<Type> & r1, const DDQD<Type> & r2, const RealContext & rc)
                {
                    IntRealConversion<Type>::typeToReal(res, r1.d_value);
                    Real tmp(rc);
                    IntRealConversion<Type>::typeToReal(tmp, r2.d_value);
                    res /= tmp;
                }
                
                inline static Real convert_frac(const DDQD<Type> & r1, const DDQD<Type> & r2, const RealContext & rc)
                {
                    Real res(rc), tmp(rc);
                    IntRealConversion<Type>::typeToReal(res, r1.d_value);
                    IntRealConversion<Type>::typeToReal(res, r2.d_value);
                    return res /= tmp;
                }
            };
            
            template<class Type> class conversion_impl<DDQD<Type>, IntegerContext>
            {
            public:
                typedef Integer RetVal;
                typedef void RetVal_Frac;
                typedef Integer RetVal_Floor;
                typedef Integer RetVal_Round;
                typedef Integer RetVal_Round2;
                typedef Integer RetVal_Ceil;
        
                static void convert(Integer & res, const DDQD<Type> & r, const IntegerContext &)
                {
                    IntRealConversion<Type>::typeToInt(res, r.d_value);
                }
        
                inline static RetVal convert(const DDQD<Type> & r, const IntegerContext & ic)
                {
                    Integer res(ic);
                    IntRealConversion<Type>::typeToInt(res, r.d_value);
                    return res;
                }
        
                static void floor(Integer & res, const DDQD<Type> & r, const IntegerContext &)
                {
                    IntRealConversion<Type>::typeToInt(res, ::floor(r.d_value));
                }
        
                inline static RetVal_Floor floor(const DDQD<Type> & r, const IntegerContext & ic)
                {
                    Integer res(ic);
                    IntRealConversion<Type>::typeToInt(res, ::floor(r.d_value));
                    return res;
                }
        
                static void round(Integer & res, const DDQD<Type> & r, const IntegerContext &)
                {
                    IntRealConversion<Type>::typeToInt(res, nint(r.d_value));
                }
        
                inline static RetVal_Round round(const DDQD<Type> & r, const IntegerContext & ic)
                {
                    Integer res(ic);
                    IntRealConversion<Type>::typeToInt(res, nint(r.d_value));
                    return res;
                }
        
                static void round(Integer & res, const DDQD<Type> & r, bool & up, const IntegerContext &)
                {
                    Type v = nint(r.d_value);
                    up = r.d_value < v;
                    IntRealConversion<Type>::typeToInt(res, v);
                }
        
                inline static RetVal_Round2 round(const DDQD<Type> & r, bool & up, const IntegerContext & ic)
                {
                    Integer res(ic);
                    Type v = nint(r.d_value);
                    up = r.d_value < v;
                    IntRealConversion<Type>::typeToInt(res, v);
                    return res;
                }
        
                static void ceil(Integer & res, const DDQD<Type> & r, const IntegerContext & rc)
                {
                    IntRealConversion<Type>::typeToInt(res, ::ceil(r.d_value));
                }
        
                inline static RetVal_Ceil ceil(const DDQD<Type> & r, const IntegerContext & ic)
                {
                    Integer res(ic);
                    IntRealConversion<Type>::typeToInt(res, ::ceil(r.d_value));
                    return res;
                }
            };
            
            template<class Type> class nativeconversion_impl<DDQD<Type> >
            {
            public:
                static int toInt(const DDQD<Type> & v)
                {
                    return IntRealConversion<Type>::typeToLongint(v.d_value);
                }
        
                static int toInt_Floor(const DDQD<Type> & v)
                {
                    return IntRealConversion<Type>::typeToLongint(::floor(v.d_value));
                }
        
                static int toInt_Round(const DDQD<Type> & v)
                {
                    return IntRealConversion<Type>::typeToLongint(nint(v.d_value));
                }
        
                static int toInt_Round(const DDQD<Type> & v, bool & up)
                {
                    Type r = nint(v.d_value);
                    up = v.d_value < r;
                    return IntRealConversion<Type>::typeToLongint(r);
                }
        
                static int toInt_Ceil(const DDQD<Type> & v)
                {
                    return IntRealConversion<Type>::typeToLongint(::ceil(v.d_value));
                }
            
                static unsigned int toUInt(const DDQD<Type> & v)
                {
                    return IntRealConversion<Type>::typeToUlongint(v.d_value);
                }
        
                static unsigned int toUInt_Floor(const DDQD<Type> & v)
                {
                    return IntRealConversion<Type>::typeToUlongint(::floor(v.d_value));
                }
        
                static unsigned int toUInt_Round(const DDQD<Type> & v)
                {
                    return IntRealConversion<Type>::typeToUlongint(nint(v.d_value));
                }
        
                static unsigned int toUInt_Round(const DDQD<Type> & v, bool & up)
                {
                    Type r = nint(v.d_value);
                    up = v.d_value < r;
                    return IntRealConversion<Type>::typeToUlongint(r);
                }
        
                static unsigned int toUInt_Ceil(const DDQD<Type> & v)
                {
                    return IntRealConversion<Type>::typeToUlongint(::ceil(v.d_value));
                }
        
                static long toLong(const DDQD<Type> & v)
                {
                    return IntRealConversion<Type>::typeToLongint(nint(v.d_value));
                }
        
                static long toLong_Floor(const DDQD<Type> & v)
                {
                    return IntRealConversion<Type>::typeToLongint(::floor(v.d_value));
                }
        
                static long toLong_Round(const DDQD<Type> & v)
                {
                    return IntRealConversion<Type>::typeToLongint(nint(v.d_value));
                }
        
                static long toLong_Round(const DDQD<Type> & v, bool & up)
                {
                    Type r = nint(v.d_value);
                    up = v.d_value < r;
                    return IntRealConversion<Type>::typeToLongint(r);
                }
        
                static long toLong_Ceil(const DDQD<Type> & v)
                {
                    return IntRealConversion<Type>::typeToLongint(::ceil(v.d_value));
                }
        
                static unsigned long toULong(const DDQD<Type> & v)
                {
                    return IntRealConversion<Type>::typeToUlongint(nint(v.d_value));
                }
        
                static unsigned long toULong_Floor(const DDQD<Type> & v)
                {
                    return IntRealConversion<Type>::typeToUlongint(::floor(v.d_value));
                }
        
                static unsigned long toULong_Round(const DDQD<Type> & v)
                {
                    return IntRealConversion<Type>::typeToUlongint(nint(v.d_value));
                }
        
                static unsigned long toULong_Round(const DDQD<Type> & v, bool & up)
                {
                    Type r = nint(v.d_value);
                    up = v.d_value < r;
                    return IntRealConversion<Type>::typeToUlongint(r);
                }
        
                static unsigned long toULong_Ceil(const DDQD<Type> & v)
                {
                    return IntRealConversion<Type>::typeToUlongint(::ceil(v.d_value));
                }
        
                static long long toLongLong(const DDQD<Type> & v)
                {
                    return IntRealConversion<Type>::typeToLonglong(nint(v.d_value));
                }
        
                static long long toLongLong_Floor(const DDQD<Type> & v)
                {
                    return IntRealConversion<Type>::typeToLonglong(::floor(v.d_value));
                }
        
                static long long toLongLong_Round(const DDQD<Type> & v)
                {
                    return IntRealConversion<Type>::typeToLonglong(nint(v.d_value));
                }
        
                static long long toLongLong_Round(const DDQD<Type> & v, bool & up)
                {
                    Type r = nint(v.d_value);
                    up = v.d_value < r;
                    return IntRealConversion<Type>::typeToLonglong(r);
                }
        
                static long long toLongLong_Ceil(const DDQD<Type> & v)
                {
                    return IntRealConversion<Type>::typeToLonglong(::ceil(v.d_value));
                }
        
                static float toFloat(const DDQD<Type> & v)
                {
                    return IntRealConversion<Type>::typeToDouble(v.d_value);
                }
        
                static double toDouble(const DDQD<Type> & v)
                {
                    return IntRealConversion<Type>::typeToDouble(v.d_value);
                }
        
                static long double toLongDouble(const DDQD<Type> & v)
                {
                    return IntRealConversion<Type>::typeToLDouble(v.d_value);
                }
            };
            
            template<typename Type, typename Op>
            struct binary_operation_impl<DDQD<Type>, DDQD<Type>, Op>
            {
                enum { supported = true, intermediate_expression = false };
                typedef DDQD<Type> IntermediateType;
                typedef DDQD<Type> ResultType;
            };
            
            template<typename Type>
            struct unary_operation_impl<DDQD<Type>, arithmetic::op::negation>
            {
                enum { supported = true, intermediate_expression = false };
                typedef DDQD<Type> IntermediateType;
                typedef DDQD<Type> ResultType;
            };
        }
    }
}
        
#endif

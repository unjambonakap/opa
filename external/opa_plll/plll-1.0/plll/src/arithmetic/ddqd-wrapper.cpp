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

#include <plll/config.hpp>

#if !defined(PLLL_CONFIG_NO_ARITHMETIC_QUADDOUBLE) || !defined(PLLL_CONFIG_NO_ARITHMETIC_DOUBLEDOUBLE)

#include <plll/helper.hpp>
#include "ddqd-wrapper.hpp"

namespace plll
{
    namespace arithmetic
    {
        namespace implementation
        {
            void IntRealConversion<dd_real>::typeToInt(Integer & i, dd_real t)
            {
                i = convert<Integer>(t._hi());
                i += convert<Integer>(t._lo());
            }
            
            dd_real IntRealConversion<dd_real>::intToType(const Integer & i)
            {
                // TODO: CREATE PROPER CONVERSION !!! ??? ...
                return lDoubleToType(convert<long double>(i));
            }
            
            void IntRealConversion<dd_real>::typeToReal(Real & r, dd_real t)
            {
                RealContext rc((long)r.precision());
                convert(r, t._hi(), rc);
                Real tmp(rc);
                convert(tmp, t._lo(), rc);
                r += tmp;
            }
            
            dd_real IntRealConversion<dd_real>::realToType(const Real & r)
            {
                // TODO: CREATE PROPER CONVERSION !!! ??? ...
                return lDoubleToType(convert<long double>(r));
            }
            
            long double IntRealConversion<dd_real>::typeToLDouble(dd_real t)
            {
                return (long double)t._hi() + (long double)t._lo();
            }
            
            dd_real IntRealConversion<dd_real>::lDoubleToType(long double d)
            {
                // TODO: CREATE PROPER CONVERSION !!! ??? ...
                dd_real r = (double)d;
                d -= (double)d; // should yield the "rest" of d
                r += (double)d;
                return r;
            }
            
            long int IntRealConversion<dd_real>::typeToLongint(dd_real t)
            {
                return (long int)t._hi() + (long int)t._lo();
            }
            
            dd_real IntRealConversion<dd_real>::longintToType(long int i)
            {
                // TODO: CREATE PROPER CONVERSION !!! ??? ...
                if (i < 0)
                    return -IntRealConversion<dd_real>::ulongintToType(-i);
                else
                    return IntRealConversion<dd_real>::ulongintToType(i);
            }
            
            unsigned long int IntRealConversion<dd_real>::typeToUlongint(dd_real t)
            {
                return (unsigned long int)t._hi() + (unsigned long int)t._lo();
            }
            
            dd_real IntRealConversion<dd_real>::ulongintToType(unsigned long int i)
            {
                // TODO: CREATE PROPER CONVERSION !!! ??? ...
                if (sizeof(int) == sizeof(long int))
                    return (int)i;
                else
                    // ASSUME: sizeof(int) = 4, sizeof(long int) = 8 !!! ??? ...
                    return (dd_real)(int)(i >> 62) * pow((dd_real)2, (dd_real)62)
                        + (dd_real)(int)(i >> 31) * pow((dd_real)2, (dd_real)31)
                        + (dd_real)(int)(i & ((1l << 31) - 1));
            }
            
            long long IntRealConversion<dd_real>::typeToLonglong(dd_real t)
            {
                return (long long)t.x[0] + (long long)t.x[1];
            }
            
            dd_real IntRealConversion<dd_real>::longlongToType(long long i)
            {
                PLLL_INTERNAL_STATIC_CHECK(sizeof(long long) == 8, Implementation_requires_longlong_to_be_64bits);
                
                // TODO: CREATE PROPER CONVERSION !!! ??? ...
                bool sign = i < 0;
                if (sign)
                    i = -i;
                dd_real r = ulongintToType((unsigned long long)(i & 0xFFFFFFFFull)) + ulongintToType((unsigned long long)(i >> 32));
                if (sign)
                    r = -r;
                return r;
            }
            
            void IntRealConversion<qd_real>::typeToInt(Integer & i, qd_real t)
            {
                i = convert<Integer>(t[0]);
                i += convert<Integer>(t[1]);
                i += convert<Integer>(t[2]);
                i += convert<Integer>(t[3]);
            }
            
            qd_real IntRealConversion<qd_real>::intToType(const Integer & i)
            {
                // TODO: CREATE PROPER CONVERSION !!! ??? ...
                return lDoubleToType(convert<long double>(i));
            }
            
            void IntRealConversion<qd_real>::typeToReal(Real & r, qd_real t)
            {
                RealContext rc((long)r.precision());
                convert(r, t[0], rc);
                Real tmp(rc);
                convert(tmp, t[1], rc);
                r += tmp;
                convert(tmp, t[2], rc);
                r += tmp;
                convert(tmp, t[3], rc);
                r += tmp;
            }
            
            qd_real IntRealConversion<qd_real>::realToType(const Real & r)
            {
                // TODO: CREATE PROPER CONVERSION !!! ??? ...
                return lDoubleToType(convert<long double>(r));
            }
            
            long double IntRealConversion<qd_real>::typeToLDouble(qd_real t)
            {
                return (long double)t[0] + (long double)t[1] + (long double)t[2] + (long double)t[3];
            }
            
            qd_real IntRealConversion<qd_real>::lDoubleToType(long double d)
            {
                // TODO: CREATE PROPER CONVERSION !!! ??? ...
                qd_real r = (double)d;
                d -= (double)d; // should yield the "rest" of d
                r += (double)d;
                return r;
            }
            
            long int IntRealConversion<qd_real>::typeToLongint(qd_real t)
            {
                return (long int)t[0] + (long int)t[1] + (long int)t[2] + (long int)t[3];
            }
            
            qd_real IntRealConversion<qd_real>::longintToType(long int i)
            {
                // TODO: CREATE PROPER CONVERSION !!! ??? ...
                if (i < 0)
                    return -IntRealConversion<qd_real>::ulongintToType(-i);
                else
                    return IntRealConversion<qd_real>::ulongintToType(i);
            }
            
            unsigned long int IntRealConversion<qd_real>::typeToUlongint(qd_real t)
            {
                return (unsigned long int)t[0] + (unsigned long int)t[1] + (unsigned long int)t[2] + (unsigned long int)t[3];
            }
            
            qd_real IntRealConversion<qd_real>::ulongintToType(unsigned long int i)
            {
                // TODO: CREATE PROPER CONVERSION !!! ??? ...
                if (sizeof(int) == sizeof(long int))
                    return (int)i;
                else
                    // ASSUME: sizeof(int) = 4, sizeof(long int) = 8 !!! ??? ...
                    return (qd_real)(int)(i >> 62) * pow((qd_real)2, (qd_real)62)
                        + (qd_real)(int)(i >> 31) * pow((qd_real)2, (qd_real)31)
                        + (qd_real)(int)(i & ((1l << 31) - 1));
            }
            
            long long IntRealConversion<qd_real>::typeToLonglong(qd_real t)
            {
                return (long long)t.x[0] + (long long)t.x[1] + (long long)t.x[2] + (long long)t.x[3];
            }
            
            qd_real IntRealConversion<qd_real>::longlongToType(long long i)
            {
                PLLL_INTERNAL_STATIC_CHECK(sizeof(long long) == 8, Implementation_requires_longlong_to_be_64bits);
                
                // TODO: CREATE PROPER CONVERSION !!! ??? ...
                bool sign = i < 0;
                if (sign)
                    i = -i;
                qd_real r = ulongintToType((unsigned long long)(i & 0xFFFFFFFFull)) + ulongintToType((unsigned long long)(i >> 32));
                if (sign)
                    r = -r;
                return r;
            }
        }
        
        long ExponentArithmetic<dd_real>::GetExponent(const dd_real & v)
        {
            int e;
            frexp(v._hi(), &e);
            return e;
        }
    
        void ExponentArithmetic<dd_real>::SetExponent(dd_real & v, long e)
        {
            v = ldexp(v, e - GetExponent(v));
        }
    
        long ExponentArithmetic<dd_real>::Normalize(dd_real & v)
        {
            long e = GetExponent(v);
            SetExponent(v, 0);
            return e;
        }
    
        long ExponentArithmetic<qd_real>::GetExponent(const qd_real & v)
        {
            int e;
            frexp(v[0], &e);
            return e;
        }
    
        void ExponentArithmetic<qd_real>::SetExponent(qd_real & v, long e)
        {
            v = ldexp(v, e - GetExponent(v));
        }
    
        long ExponentArithmetic<qd_real>::Normalize(qd_real & v)
        {
            long e = GetExponent(v);
            SetExponent(v, 0);
            return e;
        }
    
        void randomUniform(double & r, RandomNumberGenerator & rng); // defined in bigint-gmp.cpp
    
        void randomUniform(dd_real & r, RandomNumberGenerator & rng)
        {
            double d;
            randomUniform(d, rng);
            r.x[0] = d;
            randomUniform(d, rng);
            r.x[1] = ldexp(d, -std::numeric_limits<dd_real>::digits / 2);
        }
    
        void randomUniform(qd_real & r, RandomNumberGenerator & rng)
        {
            double d;
            randomUniform(d, rng);
            r.x[0] = d;
            randomUniform(d, rng);
            r.x[1] = ldexp(d, -std::numeric_limits<qd_real>::digits / 4);
            randomUniform(d, rng);
            r.x[2] = ldexp(d, -2 * std::numeric_limits<qd_real>::digits / 4);
            randomUniform(d, rng);
            r.x[3] = ldexp(d, -3 * std::numeric_limits<qd_real>::digits / 4);
            r.renorm();
        }
    }
}

#endif

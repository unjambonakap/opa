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

#ifndef PLLL_INCLUDE_GUARD__DDQD_WRAPPER
#define PLLL_INCLUDE_GUARD__DDQD_WRAPPER

#include <sstream>
#include <limits>
#include <cmath>
#include <plll/arithmetic.hpp>
#include <plll/arithmetic-nint.hpp>
#include <plll/rational.hpp>

#include <qd/dd_real.h>
#include <qd/qd_real.h>

namespace plll
{
    namespace arithmetic
    {
        template<class Type>
        // Type should be dd_real or qd_real.
        class DDQD;
        // DDQD stands for "Native Floating Point" (i.e. something native to the CPU/FPU)
    }
}

namespace plll
{
    namespace arithmetic
    {
        template<class Type>
        class NInt;
        
        namespace implementation
        {
            template<class T> class IntRealConversion;
            // Casting between type T and types Integer/Real. Will be specialized for T = dd_real, qd_real
        }
        
        template<class Type_>
        class DDQDContext
        // Emulates RealContext
        {
#ifdef PLLL_INTERNAL_NO_TEMPLATE_FRIENDS
        public:
#else
        private:
            template<class X, class Y>
            friend class implementation::conversion_impl;
            
            friend class implementation::nativeconversion_impl<DDQD<Type_> >;
#endif
            
            enum { d_precision = std::numeric_limits<Type_>::digits };
            
        public:
            typedef arithmetic::DDQD<Type_> Real;
            typedef arithmetic::DDQD<Type_> Type;
            
            enum { is_cputype = false, is_realtype = true, is_inttype = false, is_exact = false, is_variable_precision = false,
                   has_squareroot = true, has_full_power = true, has_special_fns = true, has_huge_exponent = false,
                   has_infinity = true, has_uniform_rng = true, has_constants = true, has_trigonometric = true };
            
            inline DDQDContext() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
            }
            
            inline explicit DDQDContext(long) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
            }
            
            inline ~DDQDContext() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
            }
            
            static inline void setRealPrecision(unsigned long) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
            }
            
            static inline unsigned long getRealPrecision() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return d_precision;
            }
            
            static inline unsigned long getMinRealPrecision() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE // returns minimal value for precision
            {
                return d_precision;
            }
            
            static inline unsigned long getMaxRealPrecision() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE // returns maximal value for precision
            {
                return d_precision;
            }
            
            static inline DDQD<Type_> getEpsilon() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return DDQD<Type_>(false, std::numeric_limits<Type_>::epsilon());
            }
            
            static inline void getEpsilon(DDQD<Type_> & x) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                x.assignType(std::numeric_limits<Type_>::epsilon());
            }
            
            static inline DDQD<Type_> getPi() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return DDQD<Type_>(false, implementation::IntRealConversion<Type_>::lDoubleToType(3.1415926535897932384626433832795028841971693993751l));
            }
            
            static inline void getPi(DDQD<Type_> & x) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                x = getPi();
            }
            
            static inline DDQD<Type_> getEuler() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return DDQD<Type_>(false, implementation::IntRealConversion<Type_>::lDoubleToType(2.7182818284590452353602874713526624977572470937000l));
            }
            
            static inline void getEuler(DDQD<Type_> & x) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                x = getEuler();
            }
            
            static inline DDQD<Type_> getLog2() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return DDQD<Type_>(false, implementation::IntRealConversion<Type_>::lDoubleToType(0.69314718055994530941723212145817656807550013436026l));
            }
            
            static inline void getLog2(DDQD<Type_> & x) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                x = getLog2();
            }
            
            class UniformRNG
            {
            private:
                RandomNumberGenerator & d_rng;
                
            public:
                UniformRNG(RandomNumberGenerator & rng) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    : d_rng(rng)
                {
                }
                
                inline DDQD<Type_> randomUniform(const DDQDContext<Type_> & rc)
                {
                    DDQD<Type_> r;
                    randomUniform(r);
                    return r;
                }
                
                inline void randomUniform(DDQD<Type_> & r);
            };
        };
        
        namespace traits
        {
            template<class Type>
            struct type_traits<DDQD<Type> >
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
                    is_variable_precision = false,
                    has_squareroot = true,
                    has_full_power = true,
                    has_special_fns = true,
                    has_huge_exponent = false,
                    has_constants = true,
                    has_trigonometric = true
                };
                
                typedef DDQDContext<Type> Context;
                typedef DDQD<Type> PromoteType;
                typedef DDQD<Type> ConstReferenceType;
            };
        }
        
        template<class Type>
        inline DDQD<Type> operator << (const DDQD<Type> & a, long e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline DDQD<Type> operator >> (const DDQD<Type> & a, long e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline DDQD<Type> operator << (const DDQD<Type> & a, const DDQD<Type> & e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline DDQD<Type> operator >> (const DDQD<Type> & a, const DDQD<Type> & e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline bool operator == (const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline bool operator != (const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline bool operator <= (const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline bool operator >= (const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline bool operator < (const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline bool operator > (const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        
        template<class Type>
        inline bool isZero(const DDQD<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        // Tests whether the given value is = 0.
        
        template<class Type>
        inline bool isOne(const DDQD<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        // Tests whether the given value is = 1.
        
        template<class Type>
        inline bool isPositive(const DDQD<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        // Tests whether the given value is > 0.
        
        template<class Type>
        inline bool isNonNegative(const DDQD<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        // Tests whether the given value is >= 0.
        
        template<class Type>
        inline bool isNegative(const DDQD<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        // Tests whether the given value is < 0.
        
        template<class Type>
        inline bool isNonPositive(const DDQD<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        // Tests whether the given value is <= 0.
        
        template<class Type>
        inline void add(DDQD<Type> & r, const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void sub(DDQD<Type> & r, const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void addmul(DDQD<Type> & r, const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void submul(DDQD<Type> & r, const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void mul(DDQD<Type> & r, const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void div(DDQD<Type> & r, const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void mod(DDQD<Type> & r, const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void divmod(DDQD<Type> & q, DDQD<Type> & r, const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void shl(DDQD<Type> & r, const DDQD<Type> & a, long b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void shr(DDQD<Type> & r, const DDQD<Type> & a, long b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void shl(DDQD<Type> & r, const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void shr(DDQD<Type> & r, const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void increment(DDQD<Type> & r, const DDQD<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void decrement(DDQD<Type> & r, const DDQD<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void neg(DDQD<Type> & r, const DDQD<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void abs(DDQD<Type> & r, const DDQD<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        
        template<class Type>
        inline DDQD<Type> abs(const DDQD<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline DDQD<Type> abs(const DDQD<Type> &, const DDQDContext<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void makeAbs(DDQD<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        
        template<class Type>
        std::ostream & operator << (std::ostream &, const DDQD<Type> &);
        template<class Type>
        std::istream & operator >> (std::istream &, DDQD<Type> &);
        // Output to stream or input from stream
        
        template<class Type>
        inline void setZero(DDQD<Type> &, bool positive = true) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void setOne(DDQD<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        // Set to zero/one
        
        template<class Type>
        inline int compare(const DDQD<Type> &, const DDQD<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        // Tests whether the first integer is < (-1), = (0) or > (1) than the second.
        template<class Type>
        inline int compareAbsValues(const DDQD<Type> &, const DDQD<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        // Tests whether the absolute value of the first integer is < (-1), = (0) or > (1) than the
        // absolute value of the second integer.
        template<class Type>
        inline int sign(const DDQD<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        // Returns the sign (i.e. -1, 0, 1).
        
        template<class Type>
        inline DDQD<Type> square(const DDQD<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline DDQD<Type> square(const DDQD<Type> &, const DDQDContext<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void square(DDQD<Type> &, const DDQD<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        // Some basic arithmetic operations
        
        template<class Type>
        inline DDQD<Type> sin(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline DDQD<Type> cos(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline DDQD<Type> tan(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline DDQD<Type> asin(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline DDQD<Type> acos(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline DDQD<Type> atan(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline DDQD<Type> atan2(const DDQD<Type> & r, const DDQD<Type> & r2) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline DDQD<Type> exp(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline DDQD<Type> log(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline DDQD<Type> log2(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline DDQD<Type> log10(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline DDQD<Type> sqrt(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline DDQD<Type> sin(const DDQD<Type> & r, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline DDQD<Type> cos(const DDQD<Type> & r, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline DDQD<Type> tan(const DDQD<Type> & r, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline DDQD<Type> asin(const DDQD<Type> & r, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline DDQD<Type> acos(const DDQD<Type> & r, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline DDQD<Type> atan(const DDQD<Type> & r, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline DDQD<Type> atan2(const DDQD<Type> & r, const DDQD<Type> & r2, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline DDQD<Type> exp(const DDQD<Type> & r, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline DDQD<Type> log(const DDQD<Type> & r, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline DDQD<Type> log2(const DDQD<Type> & r, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline DDQD<Type> log10(const DDQD<Type> & r, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline DDQD<Type> sqrt(const DDQD<Type> & r, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void sin(DDQD<Type> & res, const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void cos(DDQD<Type> & res, const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void tan(DDQD<Type> & res, const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void asin(DDQD<Type> & res, const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void acos(DDQD<Type> & res, const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void atan(DDQD<Type> & res, const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void atan2(DDQD<Type> & res, const DDQD<Type> & r, const DDQD<Type> & r2) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void exp(DDQD<Type> & res, const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void log(DDQD<Type> & res, const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void log2(DDQD<Type> & res, const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void log10(DDQD<Type> & res, const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void sqrt(DDQD<Type> & res, const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline DDQD<Type> power(const DDQD<Type> &, signed long) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline DDQD<Type> power(const DDQD<Type> &, unsigned long) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline DDQD<Type> power(const DDQD<Type> &, const Integer &);
        template<class Type>
        inline DDQD<Type> power(const DDQD<Type> &, const DDQD<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline DDQD<Type> power(const DDQD<Type> &, signed long, const DDQDContext<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline DDQD<Type> power(const DDQD<Type> &, unsigned long, const DDQDContext<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline DDQD<Type> power(const DDQD<Type> &, const Integer &, const DDQDContext<Type> &);
        template<class Type>
        inline DDQD<Type> power(const DDQD<Type> &, const DDQD<Type> &, const DDQDContext<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void power(DDQD<Type> &, const DDQD<Type> &, signed long) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void power(DDQD<Type> &, const DDQD<Type> &, unsigned long) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void power(DDQD<Type> &, const DDQD<Type> &, const Integer &);
        template<class Type>
        inline void power(DDQD<Type> &, const DDQD<Type> &, const DDQD<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        // Some extended functions for floats
        template<class Type>
        inline void swap(DDQD<Type> &, DDQD<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        
        namespace implementation
        {
            template<class T> class IntRealConversion;
            // Casting between type T and types Integer/Real. Will be specialized for T = dd_real, qd_real
            
            template<> class IntRealConversion<dd_real>
            {
            public:
                static void typeToInt(Integer & i, dd_real t);
                static dd_real intToType(const Integer & i);
                static void typeToReal(Real & r, dd_real t);
                static dd_real realToType(const Real & r);
                inline static double typeToDouble(dd_real t) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE { return t._hi(); }
                inline static dd_real doubleToType(double d) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE { return d; }
                static long double typeToLDouble(dd_real t);
                static dd_real lDoubleToType(long double d);
                static long int typeToLongint(dd_real t);
                static dd_real longintToType(long int i);
                static unsigned long int typeToUlongint(dd_real t);
                static dd_real ulongintToType(unsigned long int i);
                static long long typeToLonglong(dd_real t);
                static dd_real longlongToType(long long i);
            };
            
            template<> class IntRealConversion<qd_real>
            {
            public:
                static void typeToInt(Integer & i, qd_real t);
                static qd_real intToType(const Integer & i);
                static void typeToReal(Real & r, qd_real t);
                static qd_real realToType(const Real & r);
                inline static double typeToDouble(qd_real t) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE { return t[0]; }
                inline static qd_real doubleToType(double d) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE { return d; }
                static long double typeToLDouble(qd_real t);
                static qd_real lDoubleToType(long double d);
                static long int typeToLongint(qd_real t);
                static qd_real longintToType(long int i);
                static unsigned long int typeToUlongint(qd_real t);
                static qd_real ulongintToType(unsigned long int i);
                static long long typeToLonglong(qd_real t);
                static qd_real longlongToType(long long i);
            };
        }
        
        template<class T> class ExponentArithmetic;
        
        template<> class ExponentArithmetic<dd_real>
        {
        public:
            static long GetExponent(const dd_real &);
            static void SetExponent(dd_real &, long e);
            static long Normalize(dd_real & f);
        };
        
        template<> class ExponentArithmetic<qd_real>
        {
        public:
            static long GetExponent(const qd_real &);
            static void SetExponent(qd_real &, long e);
            static long Normalize(qd_real & f);
        };
        
        template<class Type>
        class DDQD
        {
#ifdef PLLL_INTERNAL_NO_TEMPLATE_FRIENDS
        public:
#else
        private:
            friend class DDQDContext<Type>;
            
            template<class X, class Y>
            friend class implementation::conversion_impl;
            
            friend class implementation::nativeconversion_impl<DDQD<Type> >;
            
#endif
            Type d_value;
            
            inline DDQD(bool, Type v) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_value(v)
            {
            }
            
            inline void assignType(Type v) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                d_value = v;
            }
            
            inline static void typeToInt(Integer & res, Type r)
            // Promise: r represents an integer
            {
                implementation::IntRealConversion<Type>::typeToInt(res, r);
            }
            
            inline static Type intToType(const Integer & i)
            {
                return implementation::IntRealConversion<Type>::intToType(i);
            }
            
        public:
            typedef DDQDContext<Type> Context;
            
            inline DDQD() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
            }
            
            inline DDQD(const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
            }
            
            inline DDQD(const DDQD & r, bool = false) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_value(r.d_value)
            {
            }
            
            inline DDQD(const DDQD & r, const DDQDContext<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_value(r.d_value)
            {
            }
            
            static inline unsigned long precision() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return DDQDContext<Type>::d_precision;
            }
            
            static inline void setContext(const DDQDContext<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
            }
            
            inline explicit DDQD(double d) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_value(implementation::IntRealConversion<Type>::doubleToType(d))
            {
            }
            
            inline explicit DDQD(double d, const DDQDContext<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_value(implementation::IntRealConversion<Type>::doubleToType(d))
            {
            }
            
            inline explicit DDQD(long double d) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_value(implementation::IntRealConversion<Type>::lDoubleToType(d))
            {
            }
            
            inline explicit DDQD(long double d, const DDQDContext<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_value(implementation::IntRealConversion<Type>::lDoubleToType(d))
            {
            }
            
            inline explicit DDQD(long i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_value(implementation::IntRealConversion<Type>::longintToType(i))
            {
            }
            
            inline explicit DDQD(long i, const DDQDContext<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_value(implementation::IntRealConversion<Type>::longintToType(i))
            {
            }
            
            inline explicit DDQD(unsigned long i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_value(implementation::IntRealConversion<Type>::ulongintToType(i))
            {
            }
            
            inline explicit DDQD(unsigned long i, const DDQDContext<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_value(implementation::IntRealConversion<Type>::ulongintToType(i))
            {
            }
            
            inline explicit DDQD(long long i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_value(implementation::IntRealConversion<Type>::longlongToType(i))
            {
            }
            
            inline explicit DDQD(long long i, const DDQDContext<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_value(implementation::IntRealConversion<Type>::longlongToType(i))
            {
            }
            
            inline explicit DDQD(const Integer & i)
                : d_value(intToType(i))
            {
            }
            
            inline explicit DDQD(const Integer & i, const DDQDContext<Type> &)
                : d_value(intToType(i))
            {
            }
            
            inline DDQD operator - () const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return DDQD(true, -d_value);
            }
        
            inline DDQD operator + (const DDQD<Type> & b) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return DDQD(true, d_value + b.d_value);
            }
        
            inline DDQD operator - (const DDQD<Type> & b) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return DDQD(true, d_value - b.d_value);
            }
            
            inline DDQD operator * (const DDQD<Type> & b) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return DDQD(true, d_value * b.d_value);
            }
            
            inline DDQD operator / (const DDQD<Type> & b) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return DDQD(true, d_value / b.d_value);
            }
            
            inline DDQD operator % (const DDQD<Type> & b) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return DDQD(true, fmod(d_value, b.d_value));
            }
            
#ifndef PLLL_INTERNAL_NO_TEMPLATE_FRIENDS
            friend DDQD<Type> operator << <>(const DDQD<Type> & a, long e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend DDQD<Type> operator >> <>(const DDQD<Type> & a, long e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend DDQD<Type> operator << <>(const DDQD<Type> & a, const DDQD<Type> & e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend DDQD<Type> operator >> <>(const DDQD<Type> & a, const DDQD<Type> & e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
#endif
            
            inline DDQD & operator = (const DDQD & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                d_value = r.d_value;
                return *this;
            }
            
            inline DDQD & operator -= (const DDQD & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                d_value -= r.d_value;
                return *this;
            }
            
            inline DDQD & operator += (const DDQD & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                d_value += r.d_value;
                return *this;
            }
            
            inline DDQD & operator *= (const DDQD & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                d_value *= r.d_value;
                return *this;
            }
            
            inline DDQD & operator /= (const DDQD & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                d_value /= r.d_value;
                return *this;
            }
            
            inline DDQD & operator %= (const DDQD & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                d_value = fmod(d_value, r.d_value);
                return *this;
            }
            
            inline DDQD & operator <<= (long r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                shl(*this, *this, r);
                return *this;
            }
            
            inline DDQD & operator >>= (long r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                shr(*this, *this, r);
                return *this;
            }
            
            inline DDQD & operator <<= (const DDQD & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                shl(*this, *this, r);
                return *this;
            }
            
            inline DDQD & operator >>= (const DDQD & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                shr(*this, *this, r);
                return *this;
            }
            
            inline DDQD operator ++ (int) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return DDQD(true, d_value + 1);
            }
            
            inline DDQD operator -- (int) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return DDQD(true, d_value - 1);
            }
            
            inline DDQD & operator ++ () PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                d_value += 1;
                return *this;
            }
            
            inline DDQD & operator -- () PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                d_value -= 1;
                return *this;
            }
            
#ifndef PLLL_INTERNAL_NO_TEMPLATE_FRIENDS
            friend bool operator == <>(const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend bool operator != <>(const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend bool operator <= <>(const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend bool operator >= <>(const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend bool operator < <>(const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend bool operator > <>(const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend bool isZero<>(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            // Tests whether the given value is = 0.
            friend bool isOne<>(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            // Tests whether the given value is = 1.
            friend bool isPositive<>(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            // Tests whether the given value is > 0.
            friend bool isNonNegative<>(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            // Tests whether the given value is >= 0.
            friend bool isNegative<>(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            // Tests whether the given value is < 0.
            friend bool isNonPositive<>(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            // Tests whether the given value is <= 0.
            friend void add<>(DDQD<Type> & r, const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void sub<>(DDQD<Type> & r, const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void addmul<>(DDQD<Type> & r, const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void submul<>(DDQD<Type> & r, const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void mul<>(DDQD<Type> & r, const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void div<>(DDQD<Type> & r, const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void mod<>(DDQD<Type> & r, const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void shl<>(DDQD<Type> & r, const DDQD<Type> & a, long b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void shr<>(DDQD<Type> & r, const DDQD<Type> & a, long b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void shl<>(DDQD<Type> & r, const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void shr<>(DDQD<Type> & r, const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void increment<>(DDQD<Type> & r, const DDQD<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void decrement<>(DDQD<Type> & r, const DDQD<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void neg<>(DDQD<Type> & r, const DDQD<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void abs<>(DDQD<Type> & r, const DDQD<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend DDQD<Type> abs<>(const DDQD<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend DDQD<Type> abs<>(const DDQD<Type> & a, const DDQDContext<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void makeAbs<>(DDQD<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend std::ostream & operator << <>(std::ostream & s, const DDQD<Type> & r);
            // Output to stream
            friend std::istream & operator >> <>(std::istream & s, DDQD<Type> & r);
            // Input from stream
            friend void setZero<>(DDQD<Type> & r, bool) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            // Set to zero
            friend void setOne<>(DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            // Set to one
            friend int compare<>(const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            // Tests whether the first real is < (-1), = (0) or > (1) than the second.
            friend int compareAbsValues<>(const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            // Tests whether the absolute value of the first DDQD is < (-1), = (0) or > (1) than the
            // absolute value of the second DDQD.
            friend int sign<>(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            // Returns the sign (i.e. -1, 0, 1).
            friend DDQD<Type> square<>(const DDQD<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend DDQD<Type> square<>(const DDQD<Type> & a, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void square<>(DDQD<Type> & r, const DDQD<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend DDQD<Type> sin<>(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend DDQD<Type> cos<>(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend DDQD<Type> tan<>(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend DDQD<Type> asin<>(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend DDQD<Type> acos<>(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend DDQD<Type> atan<>(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend DDQD<Type> atan2<>(const DDQD<Type> & r, const DDQD<Type> & r2) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend DDQD<Type> exp<>(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend DDQD<Type> log<>(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend DDQD<Type> log2<>(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend DDQD<Type> log10<>(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend DDQD<Type> sqrt<>(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend DDQD<Type> sin<>(const DDQD<Type> & r, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend DDQD<Type> cos<>(const DDQD<Type> & r, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend DDQD<Type> tan<>(const DDQD<Type> & r, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend DDQD<Type> asin<>(const DDQD<Type> & r, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend DDQD<Type> acos<>(const DDQD<Type> & r, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend DDQD<Type> atan<>(const DDQD<Type> & r, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend DDQD<Type> atan2<>(const DDQD<Type> & r, const DDQD<Type> & r2, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend DDQD<Type> exp<>(const DDQD<Type> & r, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend DDQD<Type> log<>(const DDQD<Type> & r, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend DDQD<Type> log2<>(const DDQD<Type> & r, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend DDQD<Type> log10<>(const DDQD<Type> & r, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend DDQD<Type> sqrt<>(const DDQD<Type> & r, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void sin<>(DDQD<Type> & res, const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void cos<>(DDQD<Type> & res, const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void tan<>(DDQD<Type> & res, const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void asin<>(DDQD<Type> & res, const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void acos<>(DDQD<Type> & res, const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void atan<>(DDQD<Type> & res, const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void atan2<>(DDQD<Type> & res, const DDQD<Type> & r, const DDQD<Type> & r2) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void exp<>(DDQD<Type> & res, const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void log<>(DDQD<Type> & res, const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void log2<>(DDQD<Type> & res, const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void log10<>(DDQD<Type> & res, const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void sqrt<>(DDQD<Type> & res, const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend DDQD<Type> power<>(const DDQD<Type> & a, signed long e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend DDQD<Type> power<>(const DDQD<Type> & a, unsigned long e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend DDQD<Type> power<>(const DDQD<Type> & a, const DDQD<Type> & e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend DDQD<Type> power<>(const DDQD<Type> & a, const Integer & e);
            friend DDQD<Type> power<>(const DDQD<Type> & a, signed long e, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend DDQD<Type> power<>(const DDQD<Type> & a, unsigned long e, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend DDQD<Type> power<>(const DDQD<Type> & a, const DDQD<Type> & e, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend DDQD<Type> power<>(const DDQD<Type> & a, const Integer & e, const DDQDContext<Type> & rc);
            friend void power<>(DDQD<Type> & r, const DDQD<Type> & a, signed long e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void power<>(DDQD<Type> & r, const DDQD<Type> & a, unsigned long e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void power<>(DDQD<Type> & r, const DDQD<Type> & a, const DDQD<Type> & e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void power<>(DDQD<Type> & r, const DDQD<Type> & a, const Integer & e);
            friend void arithmetic::swap<>(DDQD<Type> & a, DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
#endif
            
            long getExponent() const // number is +-f * 2^e, where e is the exponent and 0.5 <= f = 1
            {
                return ExponentArithmetic<Type>::GetExponent(d_value);
            }
            
            void setExponent(long e)
            {
                ExponentArithmetic<Type>::SetExponent(d_value, e);
            }
            
            long exponentNormalize() // moves to range [0.5, 1), returns exponent
            {
                return ExponentArithmetic<Type>::Normalize(d_value);
            }
            
            inline long getApproxExponent() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            // Returns approximate e such that |x| is approximately 2^e. The value of e should be OK up
            // to +-2.
            {
                return d_value.is_zero() ? std::numeric_limits<long>::min() : getExponent();
            }
        };
        
        template<class Type>
        inline DDQD<Type> operator << (const DDQD<Type> & a, long e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            DDQD<Type> r;
            shl(r, a, e);
            return r;
        }
        
        template<class Type>
        inline DDQD<Type> operator >> (const DDQD<Type> & a, long e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            DDQD<Type> r;
            shr(r, a, e);
            return r;
        }
        
        template<class Type>
        inline DDQD<Type> operator << (const DDQD<Type> & a, const DDQD<Type> & e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            DDQD<Type> r;
            shl(r, a, e);
            return r;
        }
        
        template<class Type>
        inline DDQD<Type> operator >> (const DDQD<Type> & a, const DDQD<Type> & e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            DDQD<Type> r;
            shr(r, a, e);
            return r;
        }
    
        template<class Type>
        inline bool operator == (const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return a.d_value == b.d_value;
        }
        
        template<class Type>
        inline bool operator != (const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return a.d_value != b.d_value;
        }
        
        template<class Type>
        inline bool operator <= (const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return a.d_value <= b.d_value;
        }
        
        template<class Type>
        inline bool operator >= (const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return a.d_value >= b.d_value;
        }
        
        template<class Type>
        inline bool operator < (const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return a.d_value < b.d_value;
        }
        
        template<class Type>
        inline bool operator > (const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return a.d_value > b.d_value;
        }
        
        template<class Type>
        inline bool isZero(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Tests whether the given value is = 0.
        {
            return r.d_value.is_zero();
        }
        
        template<class Type>
        inline bool isOne(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Tests whether the given value is = 1.
        {
            return r.d_value.is_one();
        }
        
        template<class Type>
        inline bool isPositive(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Tests whether the given value is > 0.
        {
            return r.d_value.is_positive();
        }
        
        template<class Type>
        inline bool isNonNegative(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Tests whether the given value is >= 0.
        {
            return r.d_value.is_positive() || r.d_value.is_zero();
        }
        
        template<class Type>
        inline bool isNegative(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Tests whether the given value is < 0.
        {
            return r.d_value.is_negative();
        }
        
        template<class Type>
        inline bool isNonPositive(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Tests whether the given value is <= 0.
        {
            return r.d_value.is_negative() || r.d_value.is_zero();
        }
        
        template<class Type>
        inline void add(DDQD<Type> & r, const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = a.d_value + b.d_value;
        }
        
        template<class Type>
        inline void sub(DDQD<Type> & r, const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = a.d_value - b.d_value;
        }
        
        template<class Type>
        inline void addmul(DDQD<Type> & r, const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value += a.d_value * b.d_value;
        }
    
        template<class Type>
        inline void submul(DDQD<Type> & r, const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value -= a.d_value * b.d_value;
        }
        
        template<class Type>
        inline void mul(DDQD<Type> & r, const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = a.d_value * b.d_value;
        }
        
        template<class Type>
        inline void div(DDQD<Type> & r, const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = a.d_value / b.d_value;
        }
        
        template<class Type>
        inline void mod(DDQD<Type> & r, const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = fmod(a.d_value, b.d_value);
        }
        
        template<class Type>
        inline void divmod(DDQD<Type> & q, DDQD<Type> & r, const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            Type qq = a.d_value / b.d_value, rr = a.d_value % b.d_value;
            q.d_value = qq;
            r.d_value = rr;
        }
        
        template<class Type>
        inline void shl(DDQD<Type> & r, const DDQD<Type> & a, long b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = a.d_value * pow((Type)2, implementation::IntRealConversion<Type>::longintToType(b));
        }
        
        template<class Type>
        inline void shr(DDQD<Type> & r, const DDQD<Type> & a, long b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = a.d_value * pow((Type)2, implementation::IntRealConversion<Type>::longintToType(-b));
        }
        
        template<class Type>
        inline void shl(DDQD<Type> & r, const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = a.d_value * pow((Type)2, b.d_value);
        }
        
        template<class Type>
        inline void shr(DDQD<Type> & r, const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = a.d_value * pow((Type)2, -b.d_value);
        }
        
        template<class Type>
        inline void increment(DDQD<Type> & r, const DDQD<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = a.d_value + (Type)1;
        }
        
        template<class Type>
        inline void decrement(DDQD<Type> & r, const DDQD<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = a.d_value - (Type)1;
        }
        
        template<class Type>
        inline void neg(DDQD<Type> & r, const DDQD<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = -a.d_value;
        }
        
        template<class Type>
        inline void abs(DDQD<Type> & r, const DDQD<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = fabs(a.d_value);
        }
        
        template<class Type>
        inline DDQD<Type> abs(const DDQD<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return DDQD<Type>(true, fabs(a.d_value));
        }
        
        template<class Type>
        inline DDQD<Type> abs(const DDQD<Type> & a, const DDQDContext<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return DDQD<Type>(true, fabs(a.d_value));
        }
        
        template<class Type>
        inline void makeAbs(DDQD<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            a.d_value = fabs(a.d_value);
        }
        
        template<class Type>
        std::ostream & operator << (std::ostream & s, const DDQD<Type> & r)
        // Output to stream
        {
            return s << r.d_value;
        }
        
        template<class Type>
        std::istream & operator >> (std::istream & s, DDQD<Type> & r)
        // Input from stream
        {
            return s >> r.d_value;
        }
        
        template<class Type>
        inline void setZero(DDQD<Type> & r, bool positive) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Set to zero
        {
            r.d_value = positive ? (Type)0 : -(Type)0;
        }
        
        template<class Type>
        inline void setOne(DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Set to one
        {
            r.d_value = (Type)1;
        }
        
        template<class Type>
        inline int compare(const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Tests whether the first real is < (-1), = (0) or > (1) than the second.
        {
            return a.d_value < b.d_value ? -1 : (a.d_value > b.d_value ? 1 : 0);
        }
        
        template<class Type>
        inline int compareAbsValues(const DDQD<Type> & a, const DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Tests whether the absolute value of the first DDQD is < (-1), = (0) or > (1) than the
        // absolute value of the second DDQD.
        {
            return fabs(a.d_value) < fabs(b.d_value) ? -1 : (fabs(a.d_value) > fabs(b.d_value) ? 1 : 0);
        }
        
        template<class Type>
        inline int sign(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Returns the sign (i.e. -1, 0, 1).
        {
            return r.d_value < 0 ? -1 : (r.d_value > 0 ? 1 : 0);
        }
        
        template<class Type>
        inline DDQD<Type> square(const DDQD<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return DDQD<Type>(true, a.d_value * a.d_value);
        }
        
        template<class Type>
        inline DDQD<Type> square(const DDQD<Type> & a, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return DDQD<Type>(true, a.d_value * a.d_value);
        }
        
        template<class Type>
        inline void square(DDQD<Type> & r, const DDQD<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = a.d_value * a.d_value;
        }
        
        template<class Type>
        inline DDQD<Type> sin(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return DDQD<Type>(true, sin(r.d_value));
        }
        
        template<class Type>
        inline DDQD<Type> cos(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return DDQD<Type>(true, cos(r.d_value));
        }
        
        template<class Type>
        inline DDQD<Type> tan(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return DDQD<Type>(true, tan(r.d_value));
        }
        
        template<class Type>
        inline DDQD<Type> asin(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return DDQD<Type>(true, asin(r.d_value));
        }
        
        template<class Type>
        inline DDQD<Type> acos(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return DDQD<Type>(true, acos(r.d_value));
        }
        
        template<class Type>
        inline DDQD<Type> atan(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return DDQD<Type>(true, atan(r.d_value));
        }
        
        template<class Type>
        inline DDQD<Type> atan2(const DDQD<Type> & r, const DDQD<Type> & r2) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return DDQD<Type>(true, atan2(r.d_value, r2.d_value));
        }
        
        template<class Type>
        inline DDQD<Type> exp(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return DDQD<Type>(true, exp(r.d_value));
        }
        
        template<class Type>
        inline DDQD<Type> log(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return DDQD<Type>(true, log(r.d_value));
        }
        
        template<class Type>
        inline DDQD<Type> log2(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return DDQD<Type>(true, log(r.d_value) / log((Type)2));
        }
        
        template<class Type>
        inline DDQD<Type> log10(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return DDQD<Type>(true, log10(r.d_value));
        }
        
        template<class Type>
        inline DDQD<Type> sqrt(const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return DDQD<Type>(true, sqrt(r.d_value));
        }
        
        template<class Type>
        inline DDQD<Type> sin(const DDQD<Type> & r, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return DDQD<Type>(true, sin(r.d_value));
        }
        
        template<class Type>
        inline DDQD<Type> cos(const DDQD<Type> & r, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return DDQD<Type>(true, cos(r.d_value));
        }
        
        template<class Type>
        inline DDQD<Type> tan(const DDQD<Type> & r, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return DDQD<Type>(true, tan(r.d_value));
        }
        
        template<class Type>
        inline DDQD<Type> asin(const DDQD<Type> & r, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return DDQD<Type>(true, asin(r.d_value));
        }
        
        template<class Type>
        inline DDQD<Type> acos(const DDQD<Type> & r, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return DDQD<Type>(true, acos(r.d_value));
        }
        
        template<class Type>
        inline DDQD<Type> atan(const DDQD<Type> & r, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return DDQD<Type>(true, atan(r.d_value));
        }
        
        template<class Type>
        inline DDQD<Type> atan2(const DDQD<Type> & r, const DDQD<Type> & r2, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return DDQD<Type>(true, atan2(r.d_value, r2.d_value));
        }
        
        template<class Type>
        inline DDQD<Type> exp(const DDQD<Type> & r, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return DDQD<Type>(true, exp(r.d_value));
        }
        
        template<class Type>
        inline DDQD<Type> log(const DDQD<Type> & r, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return DDQD<Type>(true, log(r.d_value));
        }
        
        template<class Type>
        inline DDQD<Type> log2(const DDQD<Type> & r, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return DDQD<Type>(true, log(r.d_value) / log((Type)2));
        }
        
        template<class Type>
        inline DDQD<Type> log10(const DDQD<Type> & r, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return DDQD<Type>(true, log10(r.d_value));
        }
        
        template<class Type>
        inline DDQD<Type> sqrt(const DDQD<Type> & r, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return DDQD<Type>(true, sqrt(r.d_value));
        }
        
        template<class Type>
        inline void sin(DDQD<Type> & res, const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            res.d_value = sin(r.d_value);
        }
        
        template<class Type>
        inline void cos(DDQD<Type> & res, const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            res.d_value = cos(r.d_value);
        }
        
        template<class Type>
        inline void tan(DDQD<Type> & res, const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            res.d_value = tan(r.d_value);
        }
        
        template<class Type>
        inline void asin(DDQD<Type> & res, const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            res.d_value = asin(r.d_value);
        }
        
        template<class Type>
        inline void acos(DDQD<Type> & res, const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            res.d_value = acos(r.d_value);
        }
        
        template<class Type>
        inline void atan(DDQD<Type> & res, const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            res.d_value = atan(r.d_value);
        }
        
        template<class Type>
        inline void atan2(DDQD<Type> & res, const DDQD<Type> & r, const DDQD<Type> & r2) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            res.d_value = atan2(r.d_value, r2.d_value);
        }
        
        template<class Type>
        inline void exp(DDQD<Type> & res, const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            res.d_value = exp(r.d_value);
        }
        
        template<class Type>
        inline void log(DDQD<Type> & res, const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            res.d_value = log(r.d_value);
        }
        
        template<class Type>
        inline void log2(DDQD<Type> & res, const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            res.d_value = log(r.d_value) / log((Type)2);
        }
        
        template<class Type>
        inline void log10(DDQD<Type> & res, const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            res.d_value = log10(r.d_value);
        }
        
        template<class Type>
        inline void sqrt(DDQD<Type> & res, const DDQD<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            res.d_value = sqrt(r.d_value);
        }
        
        template<class Type>
        inline DDQD<Type> power(const DDQD<Type> & a, signed long e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return DDQD<Type>(true, pow(a.d_value, (Type)e));
        }
        
        template<class Type>
        inline DDQD<Type> power(const DDQD<Type> & a, unsigned long e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return DDQD<Type>(true, pow(a.d_value, (Type)(int)e));
        }
        
        template<class Type>
        inline DDQD<Type> power(const DDQD<Type> & a, const DDQD<Type> & e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return DDQD<Type>(true, pow(a.d_value, e.d_value));
        }
        
        template<class Type>
        inline DDQD<Type> power(const DDQD<Type> & a, const Integer & e)
        {
            return DDQD<Type>(true, pow(a.d_value, convert<long>(e)));
        }
        
        template<class Type>
        inline DDQD<Type> power(const DDQD<Type> & a, signed long e, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return DDQD<Type>(true, pow(a.d_value, (Type)e));
        }
        
        template<class Type>
        inline DDQD<Type> power(const DDQD<Type> & a, unsigned long e, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return DDQD<Type>(true, pow(a.d_value, (Type)e));
        }
        
        template<class Type>
        inline DDQD<Type> power(const DDQD<Type> & a, const DDQD<Type> & e, const DDQDContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return DDQD<Type>(true, pow(a.d_value, e.d_value));
        }
        
        template<class Type>
        inline DDQD<Type> power(const DDQD<Type> & a, const Integer & e, const DDQDContext<Type> & rc)
        {
            return DDQD<Type>(true, pow(a.d_value, convert<long>(e)));
        }
        
        template<class Type>
        inline void power(DDQD<Type> & r, const DDQD<Type> & a, signed long e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = pow(a.d_value, (Type)(int)e);
        }
        
        template<class Type>
        inline void power(DDQD<Type> & r, const DDQD<Type> & a, unsigned long e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = pow(a.d_value, (Type)(int)e);
        }
        
        template<class Type>
        inline void power(DDQD<Type> & r, const DDQD<Type> & a, const DDQD<Type> & e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = pow(a.d_value, e.d_value);
        }
        
        template<class Type>
        inline void power(DDQD<Type> & r, const DDQD<Type> & a, const Integer & e)
        {
            r.d_value = pow(a.d_value, convert<long>(e));
        }
        
        void randomUniform(dd_real & r, RandomNumberGenerator & rng);
        void randomUniform(qd_real & r, RandomNumberGenerator & rng);
    
        template<class Type_>
        inline void DDQDContext<Type_>::UniformRNG::randomUniform(DDQD<Type_> & r)
        {
            arithmetic::randomUniform(r.d_value, d_rng);
        }
        
        template<class Type>
        inline void swap(DDQD<Type> & a, DDQD<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            std::swap(a.d_value, b.d_value);
        }
    }
}

#include "ddqd-wrapper-conv.hpp"

#endif

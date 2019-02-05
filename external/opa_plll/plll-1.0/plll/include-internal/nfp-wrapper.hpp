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

#ifndef PLLL_INCLUDE_GUARD__NFP_WRAPPER
#define PLLL_INCLUDE_GUARD__NFP_WRAPPER

#include <sstream>
#include <limits>
#include <cmath>
#include <plll/arithmetic.hpp>
#include <plll/arithmetic-nint.hpp>
#include <plll/rational.hpp>

namespace plll
{
    namespace arithmetic
    {
        template<class Type>
        // Type should be long double, double or float.
        class NFP;
        // NFP stands for "Native Floating Point" (i.e. something native to the CPU/FPU)
    }
}

namespace plll
{
    namespace arithmetic
    {
        template<class IType>
        class NInt;
    
        template<class Type_>
        class NFPContext
        // Emulates RealContext
        {
#ifdef PLLL_INTERNAL_NO_TEMPLATE_FRIENDS
        public:
#else
        private:
            template<class X, class Y>
            friend class implementation::conversion_impl;
        
            friend class implementation::nativeconversion_impl<NFP<Type_> >;
#endif
        
            enum { d_precision = std::numeric_limits<Type_>::digits }; // we assume that the radix is 2...
        
        public:
            typedef arithmetic::NFP<Type_> Real;
            typedef arithmetic::NFP<Type_> Type;
        
            enum { is_cputype = true, is_realtype = true, is_inttype = false, is_exact = false, is_variable_precision = false,
                   has_squareroot = true, has_full_power = true, has_special_fns = true, has_huge_exponent = false,
                   has_infinity = true, has_uniform_rng = true, has_constants = true, has_trigonometric = true };
        
            inline NFPContext() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
            }
            
            inline explicit NFPContext(long) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
            }
            
            inline ~NFPContext() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
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
            
            static inline NFP<Type_> getEpsilon() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return NFP<Type_>(false, std::numeric_limits<Type_>::epsilon());
            }
            
            static inline void getEpsilon(NFP<Type_> & x) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                x.assignType(std::numeric_limits<Type_>::epsilon());
            }
            
            static inline NFP<Type_> getPi() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return NFP<Type_>(false, 3.1415926535897932384626433832795028841971693993751l);
            }
            
            static inline void getPi(NFP<Type_> & x) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                x = getPi();
            }
            
            static inline NFP<Type_> getEuler() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return NFP<Type_>(false, 2.7182818284590452353602874713526624977572470937000l);
            }
            
            static inline void getEuler(NFP<Type_> & x) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                x = getEuler();
            }
            
            static inline NFP<Type_> getLog2() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return NFP<Type_>(false, 0.69314718055994530941723212145817656807550013436026l);
            }
            
            static inline void getLog2(NFP<Type_> & x) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
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
                
                inline NFP<Type_> randomUniform(const NFPContext<Type_> & rc)
                {
                    NFP<Type_> r;
                    randomUniform(r);
                    return r;
                }
                
                inline void randomUniform(NFP<Type_> & r);
            };
        };
        
        namespace traits
        {
            template<class Type>
            struct type_traits<NFP<Type> >
            {
                enum
                {
                    is_number = true,
                    is_cputype = true,
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
                
                typedef NFPContext<Type> Context;
                typedef NFP<typename type_traits<Type>::PromoteType> PromoteType;
                typedef NFP<Type> ConstReferenceType;
            };
        }
        
        template<class Type>
        inline NFP<Type> operator << (const NFP<Type> & a, long e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline NFP<Type> operator >> (const NFP<Type> & a, long e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline NFP<Type> operator << (const NFP<Type> & a, const NFP<Type> & e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline NFP<Type> operator >> (const NFP<Type> & a, const NFP<Type> & e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline bool operator == (const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline bool operator != (const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline bool operator <= (const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline bool operator >= (const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline bool operator < (const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline bool operator > (const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        
        template<class Type>
        inline bool isZero(const NFP<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        // Tests whether the given value is = 0.
        
        template<class Type>
        inline bool isOne(const NFP<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        // Tests whether the given value is = 1.
        
        template<class Type>
        inline bool isPositive(const NFP<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        // Tests whether the given value is > 0.
        
        template<class Type>
        inline bool isNonNegative(const NFP<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        // Tests whether the given value is >= 0.
        
        template<class Type>
        inline bool isNegative(const NFP<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        // Tests whether the given value is < 0.
        
        template<class Type>
        inline bool isNonPositive(const NFP<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        // Tests whether the given value is <= 0.
        
        template<class Type>
        inline void add(NFP<Type> & r, const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void sub(NFP<Type> & r, const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void addmul(NFP<Type> & r, const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void submul(NFP<Type> & r, const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void mul(NFP<Type> & r, const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void div(NFP<Type> & r, const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void mod(NFP<Type> & r, const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void divmod(NFP<Type> & q, NFP<Type> & r, const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void shl(NFP<Type> & r, const NFP<Type> & a, long b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void shr(NFP<Type> & r, const NFP<Type> & a, long b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void shl(NFP<Type> & r, const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void shr(NFP<Type> & r, const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void increment(NFP<Type> & r, const NFP<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void decrement(NFP<Type> & r, const NFP<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void neg(NFP<Type> & r, const NFP<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void abs(NFP<Type> & r, const NFP<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        
        template<class Type>
        inline NFP<Type> abs(const NFP<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline NFP<Type> abs(const NFP<Type> &, const NFPContext<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void makeAbs(NFP<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        
        template<class Type>
        std::ostream & operator << (std::ostream &, const NFP<Type> &);
        template<class Type>
        std::istream & operator >> (std::istream &, NFP<Type> &);
        // Output to stream or input from stream
        
        template<class Type>
        inline void setZero(NFP<Type> &, bool positive = true) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void setOne(NFP<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        // Set to zero/one
        
        template<class Type>
        inline int compare(const NFP<Type> &, const NFP<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        // Tests whether the first integer is < (-1), = (0) or > (1) than the second.
        template<class Type>
        inline int compareAbsValues(const NFP<Type> &, const NFP<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        // Tests whether the absolute value of the first integer is < (-1), = (0) or > (1) than the
        // absolute value of the second integer.
        template<class Type>
        inline int sign(const NFP<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        // Returns the sign (i.e. -1, 0, 1).
        
        template<class Type>
        inline NFP<Type> square(const NFP<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline NFP<Type> square(const NFP<Type> &, const NFPContext<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void square(NFP<Type> &, const NFP<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        // Some basic arithmetic operations
        
        template<class Type>
        inline NFP<Type> exp(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline NFP<Type> log(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline NFP<Type> log2(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline NFP<Type> log10(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline NFP<Type> sqrt(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline NFP<Type> gamma(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline NFP<Type> lgamma(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline NFP<Type> sin(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline NFP<Type> cos(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline NFP<Type> tan(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline NFP<Type> asin(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline NFP<Type> acos(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline NFP<Type> atan(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline NFP<Type> atan2(const NFP<Type> & r, const NFP<Type> & r2) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline NFP<Type> exp(const NFP<Type> & r, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline NFP<Type> log(const NFP<Type> & r, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline NFP<Type> log2(const NFP<Type> & r, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline NFP<Type> log10(const NFP<Type> & r, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline NFP<Type> sqrt(const NFP<Type> & r, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline NFP<Type> gamma(const NFP<Type> & r, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline NFP<Type> lgamma(const NFP<Type> & r, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline NFP<Type> sin(const NFP<Type> & r, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline NFP<Type> cos(const NFP<Type> & r, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline NFP<Type> tan(const NFP<Type> & r, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline NFP<Type> asin(const NFP<Type> & r, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline NFP<Type> acos(const NFP<Type> & r, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline NFP<Type> atan(const NFP<Type> & r, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline NFP<Type> atan2(const NFP<Type> & r, const NFP<Type> & r2, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void exp(NFP<Type> & res, const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void log(NFP<Type> & res, const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void log2(NFP<Type> & res, const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void log10(NFP<Type> & res, const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void sqrt(NFP<Type> & res, const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void gamma(NFP<Type> & res, const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void lgamma(NFP<Type> & res, const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void lgamma(NFP<Type> & res, int & sign, const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void sin(NFP<Type> & res, const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void cos(NFP<Type> & res, const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void tan(NFP<Type> & res, const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void asin(NFP<Type> & res, const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void acos(NFP<Type> & res, const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void atan(NFP<Type> & res, const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void atan2(NFP<Type> & res, const NFP<Type> & r, const NFP<Type> & r2) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline NFP<Type> power(const NFP<Type> &, signed long) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline NFP<Type> power(const NFP<Type> &, unsigned long) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline NFP<Type> power(const NFP<Type> &, const Integer &);
        template<class Type>
        inline NFP<Type> power(const NFP<Type> &, const NFP<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline NFP<Type> power(const NFP<Type> &, signed long, const NFPContext<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline NFP<Type> power(const NFP<Type> &, unsigned long, const NFPContext<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline NFP<Type> power(const NFP<Type> &, const Integer &, const NFPContext<Type> &);
        template<class Type>
        inline NFP<Type> power(const NFP<Type> &, const NFP<Type> &, const NFPContext<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void power(NFP<Type> &, const NFP<Type> &, signed long) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void power(NFP<Type> &, const NFP<Type> &, unsigned long) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        template<class Type>
        inline void power(NFP<Type> &, const NFP<Type> &, const Integer &);
        template<class Type>
        inline void power(NFP<Type> &, const NFP<Type> &, const NFP<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        // Some extended functions for floats
        template<class Type>
        inline void swap(NFP<Type> &, NFP<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        
        namespace implementation
        {
            template<class T> class IntRealConversion;
            // Casting between type T and types Integer/Real. Will be specialized for T = float, double, long double
            
            template<> class IntRealConversion<float>
            {
            public:
                inline static void typeToInt(Integer & i, float t) { i = convert<Integer>(t); }
                inline static float intToType(const Integer & i) { return convert<float>(i); }
                inline static void typeToReal(Real & r, float t) { mpfr_set_d(r.getInternal(), t, MPFR_RNDN); }
                inline static float realToType(const Real & r) { return convert<float>(r); }
            };
            
            template<> class IntRealConversion<double>
            {
            public:
                inline static void typeToInt(Integer & i, double t) { i = convert<Integer>(t); }
                inline static double intToType(const Integer & i) { return convert<double>(i); }
                inline static void typeToReal(Real & r, double t) { mpfr_set_d(r.getInternal(), t, MPFR_RNDN); }
                inline static double realToType(const Real & r) { return convert<double>(r); }
            };
            
            template<> class IntRealConversion<long double>
            {
            public:
                inline static void typeToInt(Integer & i, long double t) { i = convert<Integer>(t); }
                inline static long double intToType(const Integer & i) { return convert<long double>(i); }
                inline static void typeToReal(Real & r, long double t) { mpfr_set_ld(r.getInternal(), t, MPFR_RNDN); }
                inline static long double realToType(const Real & r) { return convert<long double>(r); }
            };
        }
        
        template<class T> class ExponentArithmetic
        {
        public:
            static long GetExponent(const T & f) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                int e;
                frexp(f, &e);
                return e;
            }
            
            static void SetExponent(T & f, long e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                int ee;
                frexp(f, &ee);
                f = ldexp(f, e - ee);
            }
            
            static long Normalize(T & f) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                int e;
                f = frexp(f, &e);
                return e;
            }
        };
        
        template<>
        class ExponentArithmetic<long double>
        // In C++98, frexp and ldexp do not have to be overloaded for long doubles.
        {
        public:
            static long GetExponent(const long double & f) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                int e;
                frexpl(f, &e);
                return e;
            }
            
            static void SetExponent(long double & f, long e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                int ee;
                frexpl(f, &ee);
                f = ldexpl(f, e - ee);
            }
            
            static long Normalize(long double & f) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                int e;
                f = frexpl(f, &e);
                return e;
            }
        };
        
        template<class Type>
        class NFP
        {
#ifdef PLLL_INTERNAL_NO_TEMPLATE_FRIENDS
        public:
#else
        private:
            friend class NFPContext<Type>;
        
            template<class X, class Y>
            friend class implementation::conversion_impl;
        
            friend class implementation::nativeconversion_impl<NFP<Type> >;
        
#endif
            Type d_value;
            
            inline NFP(bool, Type v) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
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
            typedef NFPContext<Type> Context;
            
            inline NFP() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
            }
            
            inline NFP(const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
            }
            
            inline NFP(const NFP & r, bool = false) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_value(r.d_value)
            {
            }
            
            template<typename T>
            inline NFP(const NFP<T> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_value(r.d_value)
            {
            }
            
            inline NFP(const NFP & r, const NFPContext<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_value(r.d_value)
            {
            }
        
            static inline unsigned long precision() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return NFPContext<Type>::d_precision;
            }
            
            static inline void setContext(const NFPContext<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
            }
            
            inline explicit NFP(double d) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_value(d)
            {
            }
            
            inline explicit NFP(double d, const NFPContext<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_value(d)
            {
            }
            
            inline explicit NFP(long double d) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_value(d)
            {
            }
            
            inline explicit NFP(long double d, const NFPContext<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_value(d)
            {
            }
            
            inline explicit NFP(long i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_value(i)
            {
            }
            
            inline explicit NFP(long i, const NFPContext<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_value(i)
            {
            }
            
            inline explicit NFP(unsigned long i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_value(i)
            {
            }
            
            inline explicit NFP(unsigned long i, const NFPContext<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_value(i)
            {
            }
            
            inline explicit NFP(const Integer & i)
                : d_value(intToType(i))
            {
            }
            
            inline explicit NFP(const Integer & i, const NFPContext<Type> &)
                : d_value(intToType(i))
            {
            }
            
            inline NFP operator - () const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return NFP(true, -d_value);
            }
            
            inline NFP operator + (const NFP<Type> & b) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return NFP(true, d_value + b.d_value);
            }
            
            inline NFP operator - (const NFP<Type> & b) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return NFP(true, d_value - b.d_value);
            }
            
            inline NFP operator * (const NFP<Type> & b) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return NFP(true, d_value * b.d_value);
            }
            
            inline NFP operator / (const NFP<Type> & b) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return NFP(true, d_value / b.d_value);
            }
            
            inline NFP operator % (const NFP<Type> & b) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return NFP(true, std::fmod(d_value, b.d_value));
            }
            
#ifndef PLLL_INTERNAL_NO_TEMPLATE_FRIENDS
            friend NFP<Type> operator << <>(const NFP<Type> & a, long e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NFP<Type> operator >> <>(const NFP<Type> & a, long e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NFP<Type> operator << <>(const NFP<Type> & a, const NFP<Type> & e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NFP<Type> operator >> <>(const NFP<Type> & a, const NFP<Type> & e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
#endif
            
            inline NFP & operator = (const NFP & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                d_value = r.d_value;
                return *this;
            }
            
            inline NFP & operator -= (const NFP & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                d_value -= r.d_value;
                return *this;
            }
            
            inline NFP & operator += (const NFP & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                d_value += r.d_value;
                return *this;
            }
            
            inline NFP & operator *= (const NFP & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                d_value *= r.d_value;
                return *this;
            }
            
            inline NFP & operator /= (const NFP & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                d_value /= r.d_value;
                return *this;
            }
            
            inline NFP & operator %= (const NFP & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                d_value = std::fmod(d_value, r.d_value);
                return *this;
            }
            
            inline NFP & operator <<= (long r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                shl(*this, *this, r);
                return *this;
            }
            
            inline NFP & operator >>= (long r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                shr(*this, *this, r);
                return *this;
            }
            
            inline NFP & operator <<= (const NFP & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                shl(*this, *this, r);
                return *this;
            }
            
            inline NFP & operator >>= (const NFP & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                shr(*this, *this, r);
                return *this;
            }
            
            inline NFP operator ++ (int) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return NFP(true, d_value++);
            }
            
            inline NFP operator -- (int) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return NFP(true, d_value--);
            }
            
            inline NFP & operator ++ () PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                ++d_value;
                return *this;
            }
            
            inline NFP & operator -- () PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                --d_value;
                return *this;
            }
            
#ifndef PLLL_INTERNAL_NO_TEMPLATE_FRIENDS
            friend bool operator == <>(const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend bool operator != <>(const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend bool operator <= <>(const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend bool operator >= <>(const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend bool operator < <>(const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend bool operator > <>(const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend bool isZero<>(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            // Tests whether the given value is = 0.
            friend bool isOne<>(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            // Tests whether the given value is = 1.
            friend bool isPositive<>(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            // Tests whether the given value is > 0.
            friend bool isNonNegative<>(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            // Tests whether the given value is >= 0.
            friend bool isNegative<>(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            // Tests whether the given value is < 0.
            friend bool isNonPositive<>(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            // Tests whether the given value is <= 0.
            friend void add<>(NFP<Type> & r, const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void sub<>(NFP<Type> & r, const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void addmul<>(NFP<Type> & r, const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void submul<>(NFP<Type> & r, const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void mul<>(NFP<Type> & r, const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void div<>(NFP<Type> & r, const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void mod<>(NFP<Type> & r, const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void shl<>(NFP<Type> & r, const NFP<Type> & a, long b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void shr<>(NFP<Type> & r, const NFP<Type> & a, long b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void shl<>(NFP<Type> & r, const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void shr<>(NFP<Type> & r, const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void increment<>(NFP<Type> & r, const NFP<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void decrement<>(NFP<Type> & r, const NFP<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void neg<>(NFP<Type> & r, const NFP<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void abs<>(NFP<Type> & r, const NFP<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NFP<Type> abs<>(const NFP<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NFP<Type> abs<>(const NFP<Type> & a, const NFPContext<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void makeAbs<>(NFP<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend std::ostream & operator << <>(std::ostream & s, const NFP<Type> & r);
            // Output to stream
            friend std::istream & operator >> <>(std::istream & s, NFP<Type> & r);
            // Input from stream
            friend void setZero<>(NFP<Type> & r, bool) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            // Set to zero
            friend void setOne<>(NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            // Set to one
            friend int compare<>(const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            // Tests whether the first real is < (-1), = (0) or > (1) than the second.
            friend int compareAbsValues<>(const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            // Tests whether the absolute value of the first NFP is < (-1), = (0) or > (1) than the
            // absolute value of the second NFP.
            friend int sign<>(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            // Returns the sign (i.e. -1, 0, 1).
            friend NFP<Type> square<>(const NFP<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NFP<Type> square<>(const NFP<Type> & a, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void square<>(NFP<Type> & r, const NFP<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NFP<Type> exp<>(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NFP<Type> log<>(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NFP<Type> log2<>(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NFP<Type> log10<>(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NFP<Type> sqrt<>(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NFP<Type> gamma<>(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NFP<Type> lgamma<>(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NFP<Type> sin<>(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NFP<Type> cos<>(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NFP<Type> tan<>(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NFP<Type> asin<>(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NFP<Type> acos<>(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NFP<Type> atan<>(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NFP<Type> atan2<>(const NFP<Type> & r, const NFP<Type> & r2) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NFP<Type> exp<>(const NFP<Type> & r, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NFP<Type> log<>(const NFP<Type> & r, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NFP<Type> log2<>(const NFP<Type> & r, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NFP<Type> log10<>(const NFP<Type> & r, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NFP<Type> sqrt<>(const NFP<Type> & r, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NFP<Type> gamma<>(const NFP<Type> & r, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NFP<Type> lgamma<>(const NFP<Type> & r, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NFP<Type> sin<>(const NFP<Type> & r, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NFP<Type> cos<>(const NFP<Type> & r, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NFP<Type> tan<>(const NFP<Type> & r, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NFP<Type> asin<>(const NFP<Type> & r, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NFP<Type> acos<>(const NFP<Type> & r, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NFP<Type> atan<>(const NFP<Type> & r, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NFP<Type> atan2<>(const NFP<Type> & r, const NFP<Type> & r2, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void exp<>(NFP<Type> & res, const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void log<>(NFP<Type> & res, const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void log2<>(NFP<Type> & res, const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void log10<>(NFP<Type> & res, const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void sqrt<>(NFP<Type> & res, const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void gamma<>(NFP<Type> & res, const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void lgamma<>(NFP<Type> & res, const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void lgamma<>(NFP<Type> & res, int & sign, const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void sin<>(NFP<Type> & res, const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void cos<>(NFP<Type> & res, const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void tan<>(NFP<Type> & res, const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void asin<>(NFP<Type> & res, const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void acos<>(NFP<Type> & res, const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void atan<>(NFP<Type> & res, const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void atan2<>(NFP<Type> & res, const NFP<Type> & r, const NFP<Type> & r2) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NFP<Type> power<>(const NFP<Type> & a, signed long e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NFP<Type> power<>(const NFP<Type> & a, unsigned long e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NFP<Type> power<>(const NFP<Type> & a, const NFP<Type> & e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NFP<Type> power<>(const NFP<Type> & a, const Integer & e);
            friend NFP<Type> power<>(const NFP<Type> & a, signed long e, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NFP<Type> power<>(const NFP<Type> & a, unsigned long e, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NFP<Type> power<>(const NFP<Type> & a, const NFP<Type> & e, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend NFP<Type> power<>(const NFP<Type> & a, const Integer & e, const NFPContext<Type> & rc);
            friend void power<>(NFP<Type> & r, const NFP<Type> & a, signed long e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void power<>(NFP<Type> & r, const NFP<Type> & a, unsigned long e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void power<>(NFP<Type> & r, const NFP<Type> & a, const NFP<Type> & e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            friend void power<>(NFP<Type> & r, const NFP<Type> & a, const Integer & e);
            friend void swap<>(NFP<Type> & a, NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
#endif
            
            long getExponent() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE // number is +-f * 2^e, where e is the exponent and 0.5 <= f < 1
            {
                return ExponentArithmetic<Type>::GetExponent(d_value);
            }
            
            void setExponent(long e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                ExponentArithmetic<Type>::SetExponent(d_value, e);
            }
            
            long exponentNormalize() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE // moves to range [0.5, 1), returns exponent
            {
                return ExponentArithmetic<Type>::Normalize(d_value);
            }
            
            inline long getApproxExponent() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            // Returns approximate e such that |x| is approximately 2^e. The value of e should be OK up
            // to +-2.
            {
                return (d_value == 0) ? std::numeric_limits<long>::min() : getExponent();
            }
        };
        
        template<class Type>
        inline NFP<Type> operator << (const NFP<Type> & a, long e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            NFP<Type> r;
            shl(r, a, e);
            return r;
        }
        
        template<class Type>
        inline NFP<Type> operator >> (const NFP<Type> & a, long e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            NFP<Type> r;
            shr(r, a, e);
            return r;
        }
        
        template<class Type>
        inline NFP<Type> operator << (const NFP<Type> & a, const NFP<Type> & e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            NFP<Type> r;
            shl(r, a, e);
            return r;
        }
        
        template<class Type>
        inline NFP<Type> operator >> (const NFP<Type> & a, const NFP<Type> & e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            NFP<Type> r;
            shr(r, a, e);
            return r;
        }
    
        template<class Type>
        inline bool operator == (const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return a.d_value == b.d_value;
        }
        
        template<class Type>
        inline bool operator != (const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return a.d_value != b.d_value;
        }
        
        template<class Type>
        inline bool operator <= (const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return a.d_value <= b.d_value;
        }
        
        template<class Type>
        inline bool operator >= (const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return a.d_value >= b.d_value;
        }
        
        template<class Type>
        inline bool operator < (const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return a.d_value < b.d_value;
        }
        
        template<class Type>
        inline bool operator > (const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return a.d_value > b.d_value;
        }
        
        template<class Type>
        inline bool isZero(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Tests whether the given value is = 0.
        {
            return r.d_value == (Type)0;
        }
        
        template<class Type>
        inline bool isOne(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Tests whether the given value is = 1.
        {
            return r.d_value == (Type)1;
        }
        
        template<class Type>
        inline bool isPositive(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Tests whether the given value is > 0.
        {
            return r.d_value > (Type)0;
        }
        
        template<class Type>
        inline bool isNonNegative(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Tests whether the given value is >= 0.
        {
            return r.d_value >= (Type)0;
        }
        
        template<class Type>
        inline bool isNegative(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Tests whether the given value is < 0.
        {
            return r.d_value < (Type)0;
        }
        
        template<class Type>
        inline bool isNonPositive(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Tests whether the given value is <= 0.
        {
            return r.d_value <= (Type)0;
        }
        
        template<class Type>
        inline void add(NFP<Type> & r, const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = a.d_value + b.d_value;
        }
    
        template<class Type>
        inline void sub(NFP<Type> & r, const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = a.d_value - b.d_value;
        }
        
        template<class Type>
        inline void addmul(NFP<Type> & r, const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value += a.d_value * b.d_value;
        }
    
        template<class Type>
        inline void submul(NFP<Type> & r, const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value -= a.d_value * b.d_value;
        }
        
        template<class Type>
        inline void mul(NFP<Type> & r, const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = a.d_value * b.d_value;
        }
        
        template<class Type>
        inline void div(NFP<Type> & r, const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = a.d_value / b.d_value;
        }
        
        template<class Type>
        inline void mod(NFP<Type> & r, const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = std::fmod(a.d_value, b.d_value);
        }
        
        template<class Type>
        inline void divmod(NFP<Type> & q, NFP<Type> & r, const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            Type qq = a.d_value / b.d_value, rr = a.d_value % b.d_value;
            q.d_value = qq;
            r.d_value = rr;
        }
        
        template<class Type>
        inline void shl(NFP<Type> & r, const NFP<Type> & a, long b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = a.d_value * std::pow((Type)2, (Type)b);
        }
        
        template<class Type>
        inline void shr(NFP<Type> & r, const NFP<Type> & a, long b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = a.d_value * std::pow((Type)2, (Type)-b);
        }
        
        template<class Type>
        inline void shl(NFP<Type> & r, const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = a.d_value * std::pow((Type)2, b.d_value);
        }
        
        template<class Type>
        inline void shr(NFP<Type> & r, const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = a.d_value * std::pow((Type)2, -b.d_value);
        }
        
        template<class Type>
        inline void increment(NFP<Type> & r, const NFP<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = a.d_value + (Type)1;
        }
        
        template<class Type>
        inline void decrement(NFP<Type> & r, const NFP<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = a.d_value - (Type)1;
        }
        
        template<class Type>
        inline void neg(NFP<Type> & r, const NFP<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = -a.d_value;
        }
        
        template<class Type>
        inline void abs(NFP<Type> & r, const NFP<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = std::fabs(a.d_value);
        }
        
        template<class Type>
        inline NFP<Type> abs(const NFP<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NFP<Type>(true, std::fabs(a.d_value));
        }
        
        template<class Type>
        inline NFP<Type> abs(const NFP<Type> & a, const NFPContext<Type> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NFP<Type>(true, std::fabs(a.d_value));
        }
        
        template<class Type>
        inline void makeAbs(NFP<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            a.d_value = std::fabs(a.d_value);
        }

        template<class Type>
        std::ostream & operator << (std::ostream & s, const NFP<Type> & r)
        // Output to stream
        {
            return s << r.d_value;
        }
        
        template<class Type>
        std::istream & operator >> (std::istream & s, NFP<Type> & r)
        // Input from stream
        {
            return s >> r.d_value;
        }
        
        template<class Type>
        inline void setZero(NFP<Type> & r, bool positive) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Set to zero
        {
            r.d_value = positive ? (Type)0 : -(Type)0;
        }
        
        template<class Type>
        inline void setOne(NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Set to one
        {
            r.d_value = (Type)1;
        }
        
        template<class Type>
        inline int compare(const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Tests whether the first real is < (-1), = (0) or > (1) than the second.
        {
            return a.d_value < b.d_value ? -1 : (a.d_value > b.d_value ? 1 : 0);
        }
        
        template<class Type>
        inline int compareAbsValues(const NFP<Type> & a, const NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Tests whether the absolute value of the first NFP is < (-1), = (0) or > (1) than the
        // absolute value of the second NFP.
        {
            return std::fabs(a.d_value) < std::fabs(b.d_value) ? -1 : (std::fabs(a.d_value) > std::fabs(b.d_value) ? 1 : 0);
        }
        
        template<class Type>
        inline int sign(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        // Returns the sign (i.e. -1, 0, 1).
        {
            return r.d_value < 0 ? -1 : (r.d_value > 0 ? 1 : 0);
        }
        
        template<class Type>
        inline NFP<Type> square(const NFP<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NFP<Type>(true, a.d_value * a.d_value);
        }
        
        template<class Type>
        inline NFP<Type> square(const NFP<Type> & a, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NFP<Type>(true, a.d_value * a.d_value);
        }
        
        template<class Type>
        inline void square(NFP<Type> & r, const NFP<Type> & a) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = a.d_value * a.d_value;
        }
        
        template<class Type>
        inline NFP<Type> sin(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NFP<Type>(true, std::sin(r.d_value));
        }
        
        template<class Type>
        inline NFP<Type> cos(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NFP<Type>(true, std::cos(r.d_value));
        }
        
        template<class Type>
        inline NFP<Type> tan(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NFP<Type>(true, std::tan(r.d_value));
        }
        
        template<class Type>
        inline NFP<Type> asin(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NFP<Type>(true, std::asin(r.d_value));
        }
        
        template<class Type>
        inline NFP<Type> acos(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NFP<Type>(true, std::acos(r.d_value));
        }
        
        template<class Type>
        inline NFP<Type> atan(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NFP<Type>(true, std::atan(r.d_value));
        }
        
        template<class Type>
        inline NFP<Type> atan2(const NFP<Type> & r, const NFP<Type> & r2) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NFP<Type>(true, std::atan2(r.d_value, r2.d_value));
        }
        
        template<class Type>
        inline NFP<Type> exp(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NFP<Type>(true, std::exp(r.d_value));
        }
        
        template<class Type>
        inline NFP<Type> log(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NFP<Type>(true, std::log(r.d_value));
        }
        
        template<class Type>
        inline NFP<Type> log2(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NFP<Type>(true, std::log(r.d_value) / std::log((Type)2));
        }
        
        template<class Type>
        inline NFP<Type> log10(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NFP<Type>(true, std::log10(r.d_value));
        }
        
        template<class Type>
        inline NFP<Type> sqrt(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NFP<Type>(true, std::sqrt(r.d_value));
        }
        
        template<class Type>
        class ComputeSpecialFunctionsNFP;
        
        template<>
        class ComputeSpecialFunctionsNFP<float>
        {
        public:
            static float gamma(float f) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return tgammaf(f);
            }
            
            static float lgamma(float f) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return lgammaf(f);
            }
            
            static float lgamma(float f, int & sign) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return lgammaf_r(f, &sign);
            }
        };
        
        template<>
        class ComputeSpecialFunctionsNFP<double>
        {
        public:
            static double gamma(double d) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return tgamma(d);
            }
            
            static double lgamma(double d) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return lgamma(d);
            }
            
            static double lgamma(double d, int & sign) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return lgamma_r(d, &sign);
            }
        };
        
        template<>
        class ComputeSpecialFunctionsNFP<long double>
        {
        public:
            static long double gamma(long double ld) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return tgammal(ld);
            }
            
            static long double lgamma(long double ld) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return lgammal(ld);
            }
            
            static long double lgamma(long double ld, int & sign) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return lgammal_r(ld, &sign);
            }
        };
        
        template<class Type>
        inline NFP<Type> gamma(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NFP<Type>(true, ComputeSpecialFunctionsNFP<Type>::gamma(r.d_value));
        }
        
        template<class Type>
        inline NFP<Type> lgamma(const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NFP<Type>(true, ComputeSpecialFunctionsNFP<Type>::lgamma(r.d_value));
        }
        
        template<class Type>
        inline NFP<Type> sin(const NFP<Type> & r, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NFP<Type>(true, std::sin(r.d_value));
        }
        
        template<class Type>
        inline NFP<Type> cos(const NFP<Type> & r, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NFP<Type>(true, std::cos(r.d_value));
        }
        
        template<class Type>
        inline NFP<Type> tan(const NFP<Type> & r, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NFP<Type>(true, std::tan(r.d_value));
        }
        
        template<class Type>
        inline NFP<Type> asin(const NFP<Type> & r, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NFP<Type>(true, std::asin(r.d_value));
        }
        
        template<class Type>
        inline NFP<Type> acos(const NFP<Type> & r, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NFP<Type>(true, std::acos(r.d_value));
        }
        
        template<class Type>
        inline NFP<Type> atan(const NFP<Type> & r, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NFP<Type>(true, std::atan(r.d_value));
        }
        
        template<class Type>
        inline NFP<Type> atan2(const NFP<Type> & r, const NFP<Type> & r2, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NFP<Type>(true, std::atan2(r.d_value, r2.d_value));
        }
        
        template<class Type>
        inline NFP<Type> exp(const NFP<Type> & r, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NFP<Type>(true, std::exp(r.d_value));
        }
        
        template<class Type>
        inline NFP<Type> log(const NFP<Type> & r, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NFP<Type>(true, std::log(r.d_value));
        }
        
        template<class Type>
        inline NFP<Type> log2(const NFP<Type> & r, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NFP<Type>(true, std::log(r.d_value) / std::log((Type)2));
        }
        
        template<class Type>
        inline NFP<Type> log10(const NFP<Type> & r, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NFP<Type>(true, std::log10(r.d_value));
        }
        
        template<class Type>
        inline NFP<Type> sqrt(const NFP<Type> & r, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NFP<Type>(true, std::sqrt(r.d_value));
        }
        
        template<class Type>
        inline NFP<Type> gamma(const NFP<Type> & r, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NFP<Type>(true, ComputeSpecialFunctionsNFP<Type>::gamma(r.d_value));
        }
        
        template<class Type>
        inline NFP<Type> lgamma(const NFP<Type> & r, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NFP<Type>(true, ComputeSpecialFunctionsNFP<Type>::lgamma(r.d_value));
        }
        
        template<class Type>
        inline void sin(NFP<Type> & res, const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            res.d_value = std::sin(r.d_value);
        }
        
        template<class Type>
        inline void cos(NFP<Type> & res, const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            res.d_value = std::cos(r.d_value);
        }
        
        template<class Type>
        inline void tan(NFP<Type> & res, const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            res.d_value = std::tan(r.d_value);
        }
        
        template<class Type>
        inline void asin(NFP<Type> & res, const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            res.d_value = std::asin(r.d_value);
        }
        
        template<class Type>
        inline void acos(NFP<Type> & res, const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            res.d_value = std::acos(r.d_value);
        }
        
        template<class Type>
        inline void atan(NFP<Type> & res, const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            res.d_value = std::atan(r.d_value);
        }
        
        template<class Type>
        inline void atan2(NFP<Type> & res, const NFP<Type> & r, const NFP<Type> & r2) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            res.d_value = std::atan2(r.d_value, r2.d_value);
        }
        
        template<class Type>
        inline void exp(NFP<Type> & res, const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            res.d_value = std::exp(r.d_value);
        }
        
        template<class Type>
        inline void log(NFP<Type> & res, const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            res.d_value = std::log(r.d_value);
        }
        
        template<class Type>
        inline void log2(NFP<Type> & res, const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            res.d_value = std::log(r.d_value) / std::log((Type)2);
        }
        
        template<class Type>
        inline void log10(NFP<Type> & res, const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            res.d_value = std::log10(r.d_value);
        }
        
        template<class Type>
        inline void sqrt(NFP<Type> & res, const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            res.d_value = std::sqrt(r.d_value);
        }
        
        template<class Type>
        inline void gamma(NFP<Type> & res, const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            res.d_value = ComputeSpecialFunctionsNFP<Type>::gamma(r.d_value);
        }
        
        template<class Type>
        inline void lgamma(NFP<Type> & res, const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            res.d_value = ComputeSpecialFunctionsNFP<Type>::lgamma(r.d_value);
        }
        
        template<class Type>
        inline void lgamma(NFP<Type> & res, int & sign, const NFP<Type> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            res.d_value = ComputeSpecialFunctionsNFP<Type>::lgamma(r.d_value, sign);
        }
        
        template<class Type>
        inline NFP<Type> power(const NFP<Type> & a, signed long e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NFP<Type>(true, std::pow(a.d_value, (Type)e));
        }
        
        template<class Type>
        inline NFP<Type> power(const NFP<Type> & a, unsigned long e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NFP<Type>(true, std::pow(a.d_value, (Type)e));
        }
        
        template<class Type>
        inline NFP<Type> power(const NFP<Type> & a, const NFP<Type> & e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NFP<Type>(true, std::pow(a.d_value, e.d_value));
        }
        
        template<class Type>
        inline NFP<Type> power(const NFP<Type> & a, const Integer & e)
        {
            return NFP<Type>(true, std::pow(a.d_value, convert<long>(e)));
        }
        
        template<class Type>
        inline NFP<Type> power(const NFP<Type> & a, signed long e, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NFP<Type>(true, std::pow(a.d_value, (Type)e));
        }
        
        template<class Type>
        inline NFP<Type> power(const NFP<Type> & a, unsigned long e, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NFP<Type>(true, std::pow(a.d_value, (Type)e));
        }
        
        template<class Type>
        inline NFP<Type> power(const NFP<Type> & a, const NFP<Type> & e, const NFPContext<Type> & rc) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return NFP<Type>(true, std::pow(a.d_value, e.d_value));
        }
        
        template<class Type>
        inline NFP<Type> power(const NFP<Type> & a, const Integer & e, const NFPContext<Type> & rc)
        {
            return NFP<Type>(true, std::pow(a.d_value, convert<long>(e)));
        }
        
        template<class Type>
        inline void power(NFP<Type> & r, const NFP<Type> & a, signed long e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = std::pow(a.d_value, (Type)e);
        }
        
        template<class Type>
        inline void power(NFP<Type> & r, const NFP<Type> & a, unsigned long e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = std::pow(a.d_value, (Type)e);
        }
        
        template<class Type>
        inline void power(NFP<Type> & r, const NFP<Type> & a, const NFP<Type> & e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            r.d_value = std::pow(a.d_value, e.d_value);
        }
        
        template<class Type>
        inline void power(NFP<Type> & r, const NFP<Type> & a, const Integer & e)
        {
            r.d_value = std::pow(a.d_value, convert<long>(e));
        }
        
        void randomUniform(float & r, RandomNumberGenerator & rng);
        void randomUniform(double & r, RandomNumberGenerator & rng);
        void randomUniform(long double & r, RandomNumberGenerator & rng);
        
        template<class Type_>
        inline void NFPContext<Type_>::UniformRNG::randomUniform(NFP<Type_> & r)
        {
            arithmetic::randomUniform(r.d_value, d_rng);
        }
        
        template<class Type>
        inline void swap(NFP<Type> & a, NFP<Type> & b) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            std::swap(a.d_value, b.d_value);
        }
    }
}

#include "nfp-wrapper-conv.hpp"

#endif

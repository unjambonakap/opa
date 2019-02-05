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

#include <plll/arithmetic.hpp>
#include <plll/rational.hpp>
#include <plll/arithmetic-nint.hpp>
#include <nfp-wrapper.hpp>
#if not(defined(PLLL_CONFIG_NO_ARITHMETIC_DOUBLEDOUBLE) and defined(PLLL_CONFIG_NO_ARITHMETIC_QUADDOUBLE))
  #include <ddqd-wrapper.hpp>
#endif
#include "profiling.hpp"

template<typename T>
void touch(T &)
{
}

template<class Context, class Native, bool context_float, bool native_float>
class test_native_conversions_impl;

template<class Context, class Native>
class test_native_conversions_impl<Context, Native, true, true>
{
public:
    void operator() (Context & context, typename Context::Type & var) const
    {
        // Conversion to native
        Native n = plll::arithmetic::convert<Native>(var);
        touch(n);
        // Conversion from native
        plll::arithmetic::convert(var, n, context);
        touch(var);
        var = plll::arithmetic::convert(n, context);
        touch(var);
        plll::arithmetic::convert_fraction(var, n, n, context);
        touch(var);
        var = plll::arithmetic::convert_fraction(n, n, context);
        touch(var);
    }
};

template<class Context, class Native>
class test_native_conversions_impl<Context, Native, true, false>
{
public:
    void operator() (Context & context, typename Context::Type & var) const
    {
        bool ru;
        // Conversion to native
        Native n = plll::arithmetic::convert<Native>(var);
        touch(n);
        n = plll::arithmetic::convert_floor<Native>(var);
        touch(n);
        n = plll::arithmetic::convert_round<Native>(var);
        touch(n);
        n = plll::arithmetic::convert_round<Native>(var, ru);
        touch(n);
        n = plll::arithmetic::convert_ceil<Native>(var);
        touch(n);
        // Conversion from native
        plll::arithmetic::convert(var, n, context);
        touch(var);
        var = plll::arithmetic::convert(n, context);
        touch(var);
        plll::arithmetic::convert_fraction(var, n, n, context);
        touch(var);
        var = plll::arithmetic::convert_fraction(n, n, context);
        touch(var);
    }
};

template<class Context, class Native>
class test_native_conversions_impl<Context, Native, false, true>
{
public:
    void operator() (Context & context, typename Context::Type & var) const
    {
        bool ru;
        // Conversion to native
        Native n = plll::arithmetic::convert<Native>(var);
        touch(n);
        // Conversion from native
        plll::arithmetic::convert(var, n, context);
        touch(var);
        var = plll::arithmetic::convert(n, context);
        touch(var);
        plll::arithmetic::convert_floor(var, n, context);
        touch(var);
        var = plll::arithmetic::convert_floor(n, context);
        touch(var);
        plll::arithmetic::convert_round(var, n, context);
        touch(var);
        var = plll::arithmetic::convert_round(n, context);
        touch(var);
        plll::arithmetic::convert_round(var, n, ru, context);
        touch(var);
        var = plll::arithmetic::convert_round(n, ru, context);
        touch(var);
        plll::arithmetic::convert_ceil(var, n, context);
        touch(var);
        var = plll::arithmetic::convert_ceil(n, context);
        touch(var);
    }
};

template<class Context, class Native>
class test_native_conversions_impl<Context, Native, false, false>
{
public:
    void operator() (Context & context, typename Context::Type & var) const
    {
        // Conversion to native
        Native n = plll::arithmetic::convert<Native>(var);
        touch(n);
        // Conversion from native
        plll::arithmetic::convert(var, n, context);
        touch(var);
        var = plll::arithmetic::convert(n, context);
        touch(var);
    }
};

template<class Context, class Native, bool native_float>
void test_native_conversions(Context & context, typename Context::Type & var)
{
    test_native_conversions_impl<Context, Native, Context::is_realtype, native_float>()(context, var);
}

template<class Context>
void test_native_conversions()
{
    Context context;
    typename Context::Type var(context);
    setZero(var);
    test_native_conversions<Context, signed int, false>(context, var);
    test_native_conversions<Context, unsigned int, false>(context, var);
    test_native_conversions<Context, signed long, false>(context, var);
    test_native_conversions<Context, unsigned long, false>(context, var);
    test_native_conversions<Context, long long, false>(context, var);
    test_native_conversions<Context, float, true>(context, var);
    test_native_conversions<Context, double, true>(context, var);
    test_native_conversions<Context, long double, true>(context, var);
    // Test string conversions
    std::string s1("0");
    const char * s2 = "0";
    char * s3 = (char*)new char[2];
    s3[0] = '0';
    s3[1] = 0;
    plll::arithmetic::convert(var, s1, context);
    var = plll::arithmetic::convert(s1, context);
    plll::arithmetic::convert(var, s2, context);
    var = plll::arithmetic::convert(s2, context);
    plll::arithmetic::convert(var, s3, context);
    var = plll::arithmetic::convert(s3, context);
    s1 = plll::arithmetic::convert<std::string>(var);
    s1 = plll::arithmetic::convert(var, plll::arithmetic::StringContext());
    plll::arithmetic::convert(s1, var, plll::arithmetic::StringContext());
    s1 = plll::arithmetic::convert(var, plll::arithmetic::StringContext(16));
    plll::arithmetic::convert(s1, var, plll::arithmetic::StringContext(16));
    s1 = plll::arithmetic::convert(var, plll::arithmetic::HexStringContext());
    plll::arithmetic::convert(s1, var, plll::arithmetic::HexStringContext());
    s1 = plll::arithmetic::convert(var, plll::arithmetic::OctalStringContext());
    plll::arithmetic::convert(s1, var, plll::arithmetic::OctalStringContext());
    delete[] s3;
}

template<class Native>
void test_native_conversions2_impl()
{
    Native var = 0;
    std::string s1 = "0";
    const char * s2 = "0";
    char * s3 = (char*)new char[2];
    s3[0] = '0';
    s3[1] = 0;
    var = plll::arithmetic::convert<Native>(s1);
    var = plll::arithmetic::convert<Native>(s2);
    var = plll::arithmetic::convert<Native>(s3);
    s1 = plll::arithmetic::convert<std::string>(var);
    s1 = plll::arithmetic::convert(var, plll::arithmetic::StringContext());
    plll::arithmetic::convert(s1, var, plll::arithmetic::StringContext());
    s1 = plll::arithmetic::convert(var, plll::arithmetic::StringContext(16));
    plll::arithmetic::convert(s1, var, plll::arithmetic::StringContext(16));
    s1 = plll::arithmetic::convert(var, plll::arithmetic::HexStringContext());
    plll::arithmetic::convert(s1, var, plll::arithmetic::HexStringContext());
    s1 = plll::arithmetic::convert(var, plll::arithmetic::OctalStringContext());
    plll::arithmetic::convert(s1, var, plll::arithmetic::OctalStringContext());
}

void test_native_conversions2()
{
    test_native_conversions2_impl<signed int>();
    test_native_conversions2_impl<unsigned int>();
    test_native_conversions2_impl<signed long>();
    test_native_conversions2_impl<unsigned long>();
    test_native_conversions2_impl<long long>();
    test_native_conversions2_impl<float>();
    test_native_conversions2_impl<double>();
    test_native_conversions2_impl<long double>();
}

template<class SourceContext, bool source_real, class DestinationContext, bool dest_real>
class test_conversions_impl;

template<class SourceContext, class DestinationContext>
class test_conversions_impl<SourceContext, false, DestinationContext, false>
{
public:
    void operator() (SourceContext & sc, typename SourceContext::Type & s, DestinationContext & dc, typename DestinationContext::Type & d)
    {
        convert(d, s, dc);
        touch(d);
        d = convert(s, dc);
        touch(d);
    }
};

template<class SourceContext, class DestinationContext>
class test_conversions_impl<SourceContext, false, DestinationContext, true>
{
public:
    void operator() (SourceContext & sc, typename SourceContext::Type & s, DestinationContext & dc, typename DestinationContext::Type & d)
    {
        convert(d, s, dc);
        touch(d);
        d = convert(s, dc);
        touch(d);
        convert_fraction(d, s, s, dc);
        touch(d);
        d = convert_fraction(s, s, dc);
        touch(d);
    }
};

template<class SourceContext, class DestinationContext>
class test_conversions_impl<SourceContext, true, DestinationContext, false>
{
public:
    void operator() (SourceContext & sc, typename SourceContext::Type & s, DestinationContext & dc, typename DestinationContext::Type & d)
    {
        bool ru;
        convert(d, s, dc);
        touch(d);
        d = convert(s, dc);
        touch(d);
        convert_floor(d, s, dc);
        touch(d);
        d = convert_floor(s, dc);
        touch(d);
        convert_round(d, s, dc);
        touch(d);
        d = convert_round(s, dc);
        touch(d);
        convert_round(d, s, ru, dc);
        touch(d);
        d = convert_round(s, ru, dc);
        touch(d);
        convert_ceil(d, s, dc);
        touch(d);
        d = convert_ceil(s, dc);
        touch(d);
    }
};

template<class SourceContext, class DestinationContext>
class test_conversions_impl<SourceContext, true, DestinationContext, true>
{
public:
    void operator() (SourceContext & sc, typename SourceContext::Type & s, DestinationContext & dc, typename DestinationContext::Type & d)
    {
        convert(d, s, dc);
        touch(d);
        d = convert(s, dc);
        touch(d);
        convert_fraction(d, s, s, dc);
        touch(d);
        d = convert_fraction(s, s, dc);
        touch(d);
    }
};

template<class SourceContext, class DestinationContext>
void test_conversions()
{
    SourceContext sc;
    typename SourceContext::Type s(sc);
    setZero(s);
    DestinationContext dc;
    typename DestinationContext::Type d(dc);
    test_conversions_impl<SourceContext, SourceContext::is_realtype, DestinationContext, DestinationContext::is_realtype>()(sc, s, dc, d);
}

template<class SourceContext>
void test_conversions_2()
// Conversions to the "elementary" types (integer types and Real)
{
    test_conversions<SourceContext, plll::arithmetic::IntegerContext>();
    test_conversions<SourceContext, plll::arithmetic::RealContext>();
    test_conversions<SourceContext, plll::arithmetic::NIntContext<int> >();
    test_conversions<SourceContext, plll::arithmetic::NIntContext<long> >();
    test_conversions<SourceContext, plll::arithmetic::NIntContext<long long> >();
}

template<class SourceContext>
void test_conversions_1()
// Conversions to all types
{
    test_conversions_2<SourceContext>();
    test_conversions<SourceContext, plll::arithmetic::RationalContext>();
    test_conversions<SourceContext, plll::arithmetic::NFPContext<float> >();
    test_conversions<SourceContext, plll::arithmetic::NFPContext<double> >();
    test_conversions<SourceContext, plll::arithmetic::NFPContext<long double> >();
#ifndef PLLL_CONFIG_NO_ARITHMETIC_DOUBLEDOUBLE
    test_conversions<SourceContext, plll::arithmetic::DDQDContext<dd_real> >();
#endif
#ifndef PLLL_CONFIG_NO_ARITHMETIC_QUADDOUBLE
    test_conversions<SourceContext, plll::arithmetic::DDQDContext<qd_real> >();
#endif
}

int main()
{
    plll::arithmetic::initArithmeticThreadAllocators();
    
    // Test native conversions
    test_native_conversions<plll::arithmetic::IntegerContext>();
    test_native_conversions<plll::arithmetic::RealContext>();
    test_native_conversions<plll::arithmetic::RationalContext>();
    test_native_conversions<plll::arithmetic::NFPContext<float> >();
    test_native_conversions<plll::arithmetic::NFPContext<double> >();
    test_native_conversions<plll::arithmetic::NFPContext<long double> >();
#ifndef PLLL_CONFIG_NO_ARITHMETIC_DOUBLEDOUBLE
    test_native_conversions<plll::arithmetic::DDQDContext<dd_real> >();
#endif
#ifndef PLLL_CONFIG_NO_ARITHMETIC_QUADDOUBLE
    test_native_conversions<plll::arithmetic::DDQDContext<qd_real> >();
#endif
    test_native_conversions2();
    
    test_conversions_1<plll::arithmetic::IntegerContext>();
    test_conversions_1<plll::arithmetic::RealContext>();
    test_conversions_1<plll::arithmetic::RationalContext>();
    test_conversions_1<plll::arithmetic::NIntContext<int> >();
    test_conversions_1<plll::arithmetic::NIntContext<long> >();
    test_conversions_1<plll::arithmetic::NIntContext<long long> >();
    test_conversions_2<plll::arithmetic::NFPContext<float> >();
    test_conversions_2<plll::arithmetic::NFPContext<double> >();
    test_conversions_2<plll::arithmetic::NFPContext<long double> >();
#ifndef PLLL_CONFIG_NO_ARITHMETIC_DOUBLEDOUBLE
    test_conversions_2<plll::arithmetic::DDQDContext<dd_real> >();
#endif
#ifndef PLLL_CONFIG_NO_ARITHMETIC_QUADDOUBLE
    test_conversions_2<plll::arithmetic::DDQDContext<qd_real> >();
#endif
}

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

#include <iostream> // for debugging!

#include <nfp-wrapper.hpp>
#include <plll/arithmetic-nint.hpp>
#include <plll/rational.hpp>
#include "profiling.hpp"

using namespace plll;

template<class IType>
void testInt()
{
    arithmetic::NIntContext<IType> ctx;
    arithmetic::IntegerContext ic;
    arithmetic::Integer i = arithmetic::convert<arithmetic::Integer>(1234);
    arithmetic::RealContext rc;
    arithmetic::Real r = arithmetic::convert(4321, rc);
    
    arithmetic::NInt<IType> a(ctx), b(ctx), c, d, e, f;
    arithmetic::convert(c, i, ctx);
    arithmetic::convert(i, c, ic);
    assert(arithmetic::convert<long>(i) == 1234);
    arithmetic::convert(c, r, ctx);
    arithmetic::convert(i, c, ic);
    assert(arithmetic::convert<long>(i) == 4321);

    a = arithmetic::convert(123, ctx);
    b = arithmetic::convert(234, ctx);
    c = arithmetic::convert(345, ctx);
    d = arithmetic::convert(11, ctx);
    std::cout << (a + b) * c << " " << (a - b) * (-c) << "\n";
    a = arithmetic::convert(3, ctx);
    std::cout << a << "^11 = " << power(a, 11) << "\n";
    std::cout << a << "^" << d << " = " << power(a, d) << "\n";
    a = arithmetic::convert(123, ctx);
    sqrtFloor(b, a);
    sqrtCeil(c, a);
    std::cout << a << ": " << b << " (" << b * b << ") " << c << " (" << c * c << ")\n";
    assert(b * b <= a);
    assert(c * c >= a);
    std::cout << a << ": " << floorOfLog2(a) << " " << ceilOfLog2(a) << " " << bitLength(a) << "\n";
    std::cout << a << " " << b << " " << c << " " << d << "\n";
    b = arithmetic::convert(23, ctx);
    floorDiv(c, a, b);
    ceilDiv(d, a, b);
    std::cout << c << " <= " << a << "/" << b << " <= " << d << "\n";
    assert(c * b <= a);
    assert(d * b >= a);
    euclideanDivision(e, f, a, b);
    std::cout << a << " = " << e << " * " << b << " + " << f << "\n";
    assert(isNonNegative(f));
    assert(compareAbsValues(f, b) < 0);
    assert(a == e * b + f);
    XGCD(c, d, e, a, b);
    GCD(f, a, b);
    std::cout << "GCD(" << a << ", " << b << ") = " << f << " = " << d << " * " << a << " + " << e << " * " << b << "\n";
    assert(f == c);
    assert(f == d * a + e * b);
    GCD(c, d, e);
    assert(isOne(c));
    LCM(c, a, b);
    std::cout << "LCM(" << a << ", " << b << ") = " << c << "\n";
    assert(c * f == a * b);
}

int main()
{
    arithmetic::initArithmeticThreadAllocators();

    // ...
    
    testInt<int>();
    testInt<long>();
    testInt<long long>();
}

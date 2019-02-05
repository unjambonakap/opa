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

#include <limits>
#include <iostream>
#include <iomanip>
#include <cmath>

#include <plll/arithmetic.hpp>
#include <ensurehugeexponent.hpp>
#include <nfp-wrapper.hpp>

using namespace plll;

typedef arithmetic::HugeExponent<arithmetic::NFPContext<long double> >::Type Huge;
typedef arithmetic::NFPContext<long double> LongDoubleContext;
typedef arithmetic::NFP<long double> LongDouble;

int main()
{
    LongDoubleContext rc;
    Huge a, b, c, d, e, f;
    long double
        A = 3.141592653589793238462643383279502884197169399375105820l,
        B = 2.718281828459045235360287471352662497757247093699959574l,
        C = 17*0.841470984807896506652502321630298999622563060798371065l;
    std::cout << std::setprecision(30);
    a = Huge(LongDouble(A));
    b = Huge(LongDouble(B));
    c = Huge(LongDouble(C));
    std::cout << "\nPI\n==\n";
    std::cout << A << "\n";
    std::cout << a << "\n";
    std::cout << "\nEXP(1)\n======\n";
    std::cout << B << "\n";
    std::cout << b << "\n";
    std::cout << "\n17 SIN(1)\n=========\n";
    std::cout << C << "\n";
    std::cout << c << "\n";
    std::cout << "\n";
    d = a * b;
    e = a / c;
    f = b * c;
    std::cout << d << " " << A * B << "\n";
    std::cout << e << " " << A / C << "\n";
    std::cout << f << " " << B * C << "\n";
    std::cout << (a < a) << " " << (a < b) << " " << (a < c) << " " << (a < d) << " " << (a < e) << " " << (a < f) << "\n";
    std::cout << (b < a) << " " << (b < b) << " " << (b < c) << " " << (b < d) << " " << (b < e) << " " << (b < f) << "\n";
    std::cout << (c < a) << " " << (c < b) << " " << (c < c) << " " << (c < d) << " " << (c < e) << " " << (c < f) << "\n";
    std::cout << (d < a) << " " << (d < b) << " " << (d < c) << " " << (d < d) << " " << (d < e) << " " << (d < f) << "\n";
    std::cout << (e < a) << " " << (e < b) << " " << (e < c) << " " << (e < d) << " " << (e < e) << " " << (e < f) << "\n";
    std::cout << (f < a) << " " << (f < b) << " " << (f < c) << " " << (f < d) << " " << (f < e) << " " << (f < f) << "\n";
    std::cout << a << " " << (a < b) << " " << b << "\n";
    std::cout << a << " " << (a > b) << " " << b << "\n";
    std::cout << a << " " << (a < d) << " " << d << "\n";
    std::cout << a << " " << (a > d) << " " << d << "\n";
    a *= square(Huge(LongDouble(123.0l)));
    a /= square(Huge(LongDouble(432.0l)));
    std::cout << a << "\n";
}

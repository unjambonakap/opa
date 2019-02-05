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
//#include <plll/matrix.hpp>
#include <plll/rational.hpp>
#include "profiling.hpp"

using namespace plll;

std::string binary(const arithmetic::Integer & I)
{
    std::ostringstream ss;
    if (sign(I) < 0)
        ss << '-';
    if (isZero(I))
        ss << '0';
    else
    {
        for (int b = bitLength(I) - 1; b >= 0; --b)
            ss << (bit(I, b) ? "1" : "0");
    }
    return ss.str();
}

int main()
{
    arithmetic::initArithmeticThreadAllocators();
    
    arithmetic::Integer x, y;
    setOne(x);
    x = x + y;
    std::cout << x << " " << y << "\n";
    std::cout << "\nINTEGER TO LONG DOUBLE:\n";
    std::cout << "\n";
    for (unsigned i = 0; i < 10; ++i)
    {
        std::cout << x << " " << arithmetic::convert<long double>(x) << "\n";
        ++x;
    }
    for (unsigned i = 0; i < 40; ++i)
    {
        x += x << 2;
        std::cout << x << " " << arithmetic::convert<long double>(x) << "\n";
        x = -x;
    }
    x <<= 100;
    std::cout << x << " " << arithmetic::convert<long double>(x) << "\n";
    y = x + arithmetic::convert<arithmetic::Integer>(123);
    x -= y;
    std::cout << x << " " << arithmetic::convert<long double>(x) << "\n";
    
    std::cout << "\nDOUBLE TO INTEGER:\n";
    {
        double z = (double)1/(double)16, z0 = z;
        for (int i = 0; i <= 99; ++i)
        {
            arithmetic::Integer I(z);
            std::cout << "2^" << i - 4 << "*x: " << z << " " << I << " " << binary(I) << "\n";
            z = 2*z + z0;
            z0 /= 2;
        }
        z = 0.0;
        arithmetic::Integer I(z);
        std::cout << z << " " << I << "\n";
        z = -0.0;
        I = arithmetic::convert<arithmetic::Integer>(z);
        std::cout << z << " " << I << "\n";
        z = -1.5;
        I = arithmetic::convert<arithmetic::Integer>(z);
        std::cout << z << " " << I << "\n";
    }
    std::cout << "\nLONG DOUBLE TO INTEGER:\n";
    {
        long double z = (long double)1/(long double)16, z0 = z;
        for (int i = 0; i < 130; ++i)
        {
            arithmetic::Integer I(z);
            std::cout << "2^" << i - 4 << "*x: " << z << " " << I << " " << binary(I) << "\n";
            z = 2*z + z0;
            z0 /= 2;
        }
        for (unsigned i = 130; i < 16387; ++i)
        {
            z = 2*z + z0;
            z0 /= 2;
        }
        for (unsigned i = 16387; i < 16389; ++i)
        {
            arithmetic::Integer I(z);
            std::cout << "2^" << i - 4 << "*x: " << z << " " << I << " " << binary(I) << "\n";
            z = 2*z + z0;
            z0 /= 2;
        }
        z = 0.0l;
        arithmetic::Integer I(z);
        std::cout << z << " " << I << "\n";
        z = -0.0l;
        I = arithmetic::convert<arithmetic::Integer>(z);
        std::cout << z << " " << I << "\n";
        z = -1.5l;
        I = arithmetic::convert<arithmetic::Integer>(z);
        std::cout << z << " " << I << "\n";
    }
//    linalg::math_rowvector<arithmetic::Integer> v(10, arithmetic::convert<arithmetic::Integer>(3)), w(10, arithmetic::convert<arithmetic::Integer>(2));
//    x -= x * y;
//    addmul(v, w, x);
//    x -= x * y;
//    std::cout << v.normSq() << "\n";
    arithmetic::Rational a, b, c;
    x -= x * y;
    a = b + c;
    x -= x * y;
    a = a + b;
    x -= x * y;
    a = b + a;
    x -= x * y;
    a = a + a;
    x -= x * y;
}

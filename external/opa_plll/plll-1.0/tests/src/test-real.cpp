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
#include "profiling.hpp"

using namespace plll;

int main()
{
    arithmetic::initArithmeticThreadAllocators();
    
    arithmetic::getThreadRealContext().setRealPrecision(1000);
    
    std::cout << std::setprecision(30);
    
//    std::cout << std::showpos;
//    std::cout << std::internal;
//    std::cout << std::left;
    std::cout << std::right;
    
//    std::cout << std::showpoint;
    
//    std::cout << std::showbase;
    std::cout << std::dec;
//    std::cout << std::hex;
//    std::cout << std::oct;
    
//    std::cout << std::fixed;
//    std::cout << std::scientific;

    /*
    arithmetic::Real p;
    arithmetic::convert(p, 3.141592654);
    p >>= 20;
    for (unsigned i = 0; i < 40; ++i)
    {
        std::cout << std::setw(30) << p << "    " << std::setw(30) << arithmetic::convert<long double>(p) << "\n";
        p <<= 1;
    }
    return 0;
    */
    
    arithmetic::RealContext rc;
    arithmetic::Real x, y, z;
    arithmetic::Integer xx, yy;
    double xxx, yyy, zzz;
    setOne(x);
    setZero(y);
    x = x + y;
    setOne(xx);
    setZero(yy);
    xx = xx + yy;
    xxx = 1;
    yyy = 0;
    xxx = xxx + yyy;
    std::cout << std::setw(30) << x << " " << std::setw(30) << y << " -- " << std::setw(30) << xx << " " << std::setw(30) << yy << " -- " << std::setw(30) << xxx << " " << std::setw(30) << yyy << "\n";
    for (unsigned i = 0; i < 60; ++i)
    {
        if (i & 1)
        {
            x = (x + y);
            x <<= 1;
        }
        else
        {
            y = (x + y);
            y <<= 1;
        }
        if (i & 1)
        {
            xx = (xx + yy);
            xx <<= 1;
        }
        else
        {
            yy = (xx + yy);
            yy <<= 1;
        }
        if (i & 1)
        {
            xxx = (xxx + yyy);
            xxx *= 2;
        }
        else
        {
            yyy = (xxx + yyy);
            yyy *= 2;
        }
        std::cout << std::setw(30) << x << " " << std::setw(30) << y << " -- " << std::setw(30) << xx << " " << std::setw(30) << yy << " -- " << std::setw(30) << xxx << " " << std::setw(30) << yyy << "\n";
    }
    setOne(x);
    setZero(y);
    setOne(z);
    xxx = 1;
    yyy = 0;
    zzz = 1;
    for (unsigned i = 0; i < 30; ++i)
    {
        if (i & 1)
        {
            x = (x + y);
            x >>= 1;
        }
        else
        {
            y = (x + y);
            y >>= 1;
        }
        z >>= 1;
        if (i & 1)
        {
            xxx = (xxx + yyy);
            xxx /= 2;
        }
        else
        {
            yyy = (xxx + yyy);
            yyy /= 2;
        }
        zzz /= 2;
        std::cout << std::setw(30) << x << " " << std::setw(30) << y << " " << std::setw(30) << z << " -- " << std::setw(30) << xxx << " " << std::setw(30) << yyy << " " << std::setw(30) << zzz << "\n";
    }
    std::cout << std::setw(30) << arithmetic::convert(123.4567, rc) << " " << std::setw(30) << arithmetic::convert(0.0123, rc) << "\n";
    std::cout << std::setw(30) << arithmetic::convert(0.123, rc) << " " << std::setw(30) << arithmetic::convert(-1.23, rc) << "\n";
    
    arithmetic::getThreadRealContext().getEpsilon(y);
    std::cout << std::setw(30) << x << " " << std::setw(30) << y << " " << std::setw(30) << (x + y) - x << " " << std::setw(30) << (x + (y >> 1)) - x << "\n";
}

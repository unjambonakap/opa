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
#include <vector>
#include <sstream>
#include <boost/math/distributions/chi_squared.hpp>

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

const char * RED = "\033[1;31m";
const char * DARKRED = "\033[0;31m";
const char * YELLOW = "\033[1;33m";
const char * DARKYELLOW = "\033[0;33m";
const char * GREEN = "\033[1;32m";
const char * DARKGREEN = "\033[0;32m";
const char * NORMAL = "\033[0m";

void test1part(arithmetic::RandomNumberGenerator & rng, unsigned long n)
{
    unsigned NN = n <= 10 ? 2500 + n * 250 : n <= 100 ? 500 : 100;
    signed long max = -1;
    std::vector<unsigned> pvalcnt(10);
    std::vector<unsigned long> gcnt(n);
    for (unsigned ii = 0; ii < NN; ++ii)
    {
        unsigned N = n <= 10 ? 5000 + 500 * n : n <= 100 ? 500 + 250 * n : 50 * n;
        std::vector<unsigned long> cnt(n);
        for (unsigned i = 0; i < N; ++i)
        {
            unsigned long x = rng.random((unsigned long)n);
            assert(x < n);
            ++cnt[x]; ++gcnt[x];
            if ((signed long)x > max)
                max = x;
        }
        if (n > 1)
        {
            double chis = 0.0;
            for (unsigned i = 0; i < n; ++i)
            {
                double s = (double)cnt[i] - (double)N / (double)n;
                chis += s * s / ((double)N / (double)n);
            }
            boost::math::chi_squared_distribution<double> dist = boost::math::chi_squared_distribution<double>(n - 1);
            double pval = chis > 0 ? 1.0 - boost::math::cdf(dist, chis) : 1.0;
//                std::cout << "chi^2 = " << chis << ", p-value = " << pval << "\n";
            
            unsigned pv = (unsigned)(pval * 10);
            if (pv == 10)
                --pv;
            ++pvalcnt[pv];
        }
    }
    if ((max < (signed long)n - 1) || (max < 0))
        std::cout << max << " " << n << "\n";
    double chis = 0.0;
    for (unsigned i = 0; i < 10; ++i)
    {
        double s = (double)pvalcnt[i] - (double)NN / (double)10;
        chis += s * s / ((double)NN / (double)10);
    }
    boost::math::chi_squared_distribution<double> dist = boost::math::chi_squared_distribution<double>(9);
    double pval = chis > 0 ? 1.0 - boost::math::cdf(dist, chis) : 1.0/+0.0;
    std::cout << n << ": chi^2 = " << chis << ", p-value = ";
    if (pval < 0.001)
        std::cout << RED;
    else if (pval < 0.01)
        std::cout << YELLOW;
    else if (pval < 0.1)
        std::cout << DARKYELLOW;
    else
        std::cout << GREEN;
    std::cout << pval << NORMAL << ":\t";
    for (unsigned i = 0; i < 10; ++i)
        std::cout << " " << pvalcnt[i];
    std::cout << "\n";
//    for (unsigned i = 0; i < n; ++i)
//        std::cout << " " << gcnt[i];
//    std::cout << "\n";
}

void test1(arithmetic::RandomNumberGenerator & rng)
{
    std::cout << "TEST 1 BEGIN\n";
    for (unsigned long n = 2; n < 50; ++n)
        test1part(rng, n);
    std::cout << "TEST 1 DONE\n";
}

void test2(arithmetic::RandomNumberGenerator & rng)
{
    std::cout << "TEST 2 BEGIN\n";
    for (unsigned N = 2; N < 190; ++N)
    {
        for (unsigned i = 0; i < 10; ++i)
        {
            arithmetic::Integer x = rng.randomBits(N);
            if ((bitLength(x) > N) || (sign(x) < 0))
                std::cout << "  " << N << " " << bitLength(x)  << " " << binary(x) << "\n";
            assert(bitLength(x) <= N);
        }
    }
    std::cout << "TEST 2 DONE\n";
}

void test3(arithmetic::RandomNumberGenerator & rng)
{
    std::cout << "TEST 3 BEGIN\n";
    for (unsigned N = 2; N < 190; ++N)
    {
        for (unsigned i = 0; i < 10; ++i)
        {
            arithmetic::Integer x = rng.randomLen(N);
            if (bitLength(x) != N)
                std::cout << "  " << N << " " << bitLength(x)  << " " << binary(x) << "\n";
            assert(bitLength(x) == N);
        }
    }
    std::cout << "TEST 3 DONE\n";
}

void test4part(arithmetic::RandomNumberGenerator & rng, unsigned long n)
{
    unsigned NN = n <= 10 ? 300 : n <= 100 ? 200 : 70;
    signed long max = -1;
    std::vector<unsigned> pvalcnt(10);
    for (unsigned ii = 0; ii < NN; ++ii)
    {
        unsigned N = n <= 10 ? 1000 + 100 * n : n <= 100 ? 1000 + 50 * n : 100 * n;
        std::vector<unsigned long> cnt(n);
        for (unsigned i = 0; i < N; ++i)
        {
            unsigned long x = arithmetic::convert<long>(rng.random(arithmetic::convert<arithmetic::Integer>(n)));
            assert(x < n);
            ++cnt[x];
            if ((signed long)x > max)
                max = x;
        }
        if (n > 1)
        {
            double chis = 0.0;
            for (unsigned i = 0; i < n; ++i)
            {
                double s = (double)cnt[i] - (double)N / (double)n;
                chis += s * s / ((double)N / (double)n);
            }
            boost::math::chi_squared_distribution<double> dist = boost::math::chi_squared_distribution<double>(n - 1);
            double pval = chis > 0 ? 1.0 - boost::math::cdf(dist, chis) : 1.0;
//                std::cout << "chi^2 = " << chis << ", p-value = " << pval << "\n";
            
            unsigned pv = (unsigned)(pval * 10);
            if (pv == 10)
                --pv;
            ++pvalcnt[pv];
        }
    }
    if ((max < (signed long)n - 1) || (max < 0))
        std::cout << max << " " << n << "\n";
    double chis = 0.0;
    for (unsigned i = 0; i < 10; ++i)
    {
        double s = (double)pvalcnt[i] - (double)NN / (double)10;
        chis += s * s / ((double)NN / (double)10);
    }
    boost::math::chi_squared_distribution<double> dist = boost::math::chi_squared_distribution<double>(9);
    double pval = chis > 0 ? 1.0 - boost::math::cdf(dist, chis) : 1.0/+0.0;
    std::cout << n << ": chi^2 = " << chis << ", p-value = ";
    if (pval < 0.001)
        std::cout << RED;
    else if (pval < 0.01)
        std::cout << YELLOW;
    else if (pval < 0.1)
        std::cout << DARKYELLOW;
    else
        std::cout << GREEN;
    std::cout << pval << NORMAL << ":\t";
    for (unsigned i = 0; i < 10; ++i)
        std::cout << " " << pvalcnt[i];
    std::cout << "\n";
}

void test4(arithmetic::RandomNumberGenerator & rng)
{
    std::cout << "TEST 4 BEGIN\n";
    for (unsigned long n = 2; n < 50; ++n)
        test4part(rng, n);
    std::cout << "TEST 4 DONE\n";
}

void test5(arithmetic::RandomNumberGenerator & rng)
{
}

void test6(arithmetic::RandomNumberGenerator & rng)
{
}

void test7(arithmetic::RandomNumberGenerator & rng)
{
}

void test8(arithmetic::RandomNumberGenerator & rng)
{
}

int main(int argc, char **argv)
{
    arithmetic::initArithmeticThreadAllocators();

    arithmetic::RandomNumberGenerator rng;
//    rng.randomizeTime();
    rng.randomizeDevRandom();

    if (argc > 1)
    {
        int n = atoi(argv[1]);
        if (n < 2)
            n = 2;
        for (int i = 0; i < 10; ++i)
            test1part(rng, n);
        return 0;
    }
    
    test1(rng);
    std::cout << "-----------------------------------------------------------\n";
//    test2(rng);
    std::cout << "-----------------------------------------------------------\n";
//    test3(rng);
    std::cout << "-----------------------------------------------------------\n";
    test4(rng);
    std::cout << "-----------------------------------------------------------\n";
    test5(rng);
    std::cout << "-----------------------------------------------------------\n";
    test6(rng);
    std::cout << "-----------------------------------------------------------\n";
    test7(rng);
    std::cout << "-----------------------------------------------------------\n";
    test8(rng);
}

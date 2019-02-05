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
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <plll/config.hpp>
#include "primes.hpp"

// TODO: Implement http://en.wikipedia.org/wiki/Sieve_of_Atkin or something similarly fast !!! ??? ...

namespace plll
{
    inline unsigned long bit(unsigned long n, unsigned long b)
    {
        return (n >> b) & 1;
    }
        
    inline long floorOfLog2(unsigned long n)
    {
        long r = 0;
        while (n)
        {
            n >>= 1;
            ++r;
        }
        return r;
    }
        
    inline void setOne(unsigned long & n)
    {
        n = 1;
    }
        
    inline bool isZero(unsigned long n)
    {
        return n == 0;
    }
        
    inline bool isOne(unsigned long n)
    {
        return n == 1;
    }
    
    template<class T> T fastexpmod(const T & x, const T & e, const T & n)
    {
        T res;
        if (isZero(e))
        {
            setOne(res);
            return res;
        }
        if (isOne(e))
            return res = x;
        long b = floorOfLog2(e);
        setOne(res);
        while (b >= 0)
        {
            res = (res * res) % n;
            if (bit(e, b))
                res = (res * x) % n;
            --b;
        }
        return res;
    }
    
    bool isPrimeMR(unsigned long n, unsigned long a, int r, unsigned long d)
    // Tests whether n is prime using NoTests Miller-Rabin strong pseudoprime tests. If false is
    // returned, the number is composite. If true is returned, the number is prime with probability
    // 4^{-NoTests}.
    {
        unsigned long x = fastexpmod(a, d, n), prevx = x;
        for (int e = 0; (e < r) && (x != 1); ++e)
        {
            prevx = x;
            x = (x * x) % n;
        }
        if (((prevx != 1) && (prevx != n - 1)) || (x != 1))
            return false;
        return true;
    }

    bool isPrimeMR(unsigned long n, unsigned long NoTests)
    // Tests whether n is prime using NoTests Miller-Rabin strong pseudoprime tests. If false is
    // returned, the number is composite. If true is returned, the number is prime with probability
    // 4^{-NoTests}.
    {
        if (n < 2)
            return false;
        if (n == 2)
            return true;
    
        if (n >= std::sqrt((double)std::numeric_limits<unsigned long>::max()))
            return isPrimeMR(arithmetic::convert<arithmetic::Integer>(n), NoTests);
    
        // Write n - 1 = 2^r * d
        int r = 0;
        unsigned long d = n - 1;
        while ((d & 1) == 0)
        {
            ++r;
            d >>= 1;
        }
    
        // Determine log_2(n)
        unsigned long nn = n, l = 0;
        while (nn)
        {
            nn >>= 1;
            ++l;
        }
        --l;
    
        // Makes the test fully deterministic! (see http://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test#Deterministic_variants_of_the_test)
        if (l <= 32)
        {
            if (l <= 23)
                return isPrimeMR(n, 31, r, d) && isPrimeMR(n, 73, r, d);
            else
                return isPrimeMR(n, 7, r, d) && isPrimeMR(n, 61, r, d);
        }
        else if (l <= 48)
        {
            if (!isPrimeMR(n, 3, r, d) || !isPrimeMR(n, 5, r, d) || !isPrimeMR(n, 7, r, d) || !isPrimeMR(n, 11, r, d))
                return false;
            if (l <= 40)
                return true;
            if (!isPrimeMR(n, 13, r, d))
                return false;
            if (l <= 41)
                return true;
            else
                return isPrimeMR(n, 17, r, d);
        }
        else
        {
            for (unsigned long test = 0; test < NoTests; ++test)
            {
                unsigned long a = std::rand() % n;
                if (a < 2)
                    a = 2;
                if (!isPrimeMR(n, a, r, d))
                    return false;
            }
            return true;
        }
    }

    bool isPrimeMR(const arithmetic::Integer & n, unsigned long NoTests)
    // Tests whether n is prime using NoTests Miller-Rabin strong pseudoprime tests. If false is
    // returned, the number is compoiste. If true is returned, the number is prime with probability
    // 4^{-NoTests}.
    {
        if (n < arithmetic::convert<arithmetic::Integer>(2))
            return false;
        if (n == arithmetic::convert<arithmetic::Integer>(2))
            return true;
        if (n < arithmetic::convert<arithmetic::Integer>(65535))
            return isPrimeMR(arithmetic::convert<long>(n), NoTests);
    
        // Check with some precomputed primes
        static arithmetic::Integer prod = arithmetic::convert<arithmetic::Integer>("4711930799906184953162487834760260422020574773409675520188634839616415335845034221205289256705544681972439104097777157991804380284218315038719444943990492579030720635990538452312528339864352999310398481791730017201031090"); // product(ithprime(i), i=1..100)
        if (!isOne(GCD(n, prod))) // the case that n is one of the first 100 primes cannot appear, since
            // for small values of n, we use the <unsigned long> version (see
            // above).
            return false;
    
        // Write n - 1 = 2^r * d
        int r = 0;
        arithmetic::Integer d = n - arithmetic::convert<arithmetic::Integer>(1);
        while (bit(d, 0) == 0)
        {
            ++r;
            d >>= 1;
        }
    
        for (unsigned long test = 0; test < NoTests; ++test)
        {
            arithmetic::Integer a = arithmetic::convert<arithmetic::Integer>(std::rand()) % n;
            if (a < arithmetic::convert<arithmetic::Integer>(2))
                a = arithmetic::convert<arithmetic::Integer>(2);
            arithmetic::Integer x = fastexpmod(a, d, n), prevx = x;
            for (int e = 0; (e < r) && !isOne(x); ++e)
            {
                prevx = x;
                x = (x * x) % n;
            }
            if ((!isOne(x) && !isOne(n - prevx)) || !isOne(x))
                return false;
        }
        return true;
    }

    void PrimesGenerator::sieve(unsigned long bound)
    {
        d_primes.clear();
        d_primes.push_back(2);
        std::vector<bool> sieve(bound);
        // Index Number
        // ------------
        //   0      3
        //   1      5
        //   2      7
        //   n    2n+3
        for (unsigned long i = 0; i < bound; ++i)
            if (!sieve[i])
            {
                // Found a prime!
                d_primes.push_back(2 * i + 3);
                // Mark the rest
                unsigned long m = 3 * i + 3; // index of 3*(2i+3)
                while (m < bound)
                {
                    sieve[m] = true;
                    m += 2 * i + 3; // go to the next odd multiple of 2i+3
                }
            }
    }

    void PrimesGenerator::generateAnotherPrime(unsigned long NoTests)
    {
        unsigned long n = d_primes.back() + 2;
        while (!isPrimeMR(n, NoTests))
            n += 2;
        d_primes.push_back(n);
    }

    void PrimesGenerator::generateUpTo(unsigned long upperBound)
    /* Generate all primes <= upperBound. Might generate more than that, but all entries in the d_primes
     * variable are primes. */
    {
        if (upperBound < 3)
            upperBound = 3;
        sieve(upperBound/2);
    }

    void PrimesGenerator::generateNOP(unsigned long numberOfPrimes)
    /* Generates at least numberOfPrimes primes. */
    {
        // Create bound for largest prime using PNT
        double l = std::log((double)numberOfPrimes), ll = std::log(l);
        double m = 1;
        // For these bounds, see P. Dusart (1999). "The kth prime is greater than k(ln k + ln ln k-1)
        // for k>=2". Mathematics of Computation 68: 411--415.
        // http://www.ams.org/mcom/1999-68-225/S0025-5718-99-01037-6/S0025-5718-99-01037-6.pdf.
        if (numberOfPrimes >= 27076)
            m = (ll - 1.8) / l;
        if ((numberOfPrimes >= 39017) && (m > 0.0516))
            m = 0.0516;
        unsigned long bound = (unsigned long)(numberOfPrimes * (l + ll - 1 + m));
        if (bound < 3)
            bound = 3;
//    std::cout << "bound = " << bound << "\n";
    
        // Generate all primes until this bound
        generateUpTo(bound);
//    std::cout << "found " << d_primes.size() << " primes\n";
//    std::cout << "the largest is " << d_primes.back() << "\n";

        // Generate missing primes
        if (d_primes.size() < numberOfPrimes)
        {
            std::cerr << "generating " << numberOfPrimes - d_primes.size() << " missing primes...\n";
            while (d_primes.size() < numberOfPrimes)
                generateAnotherPrime();
        }
    }

    arithmetic::Integer nextPrime(const arithmetic::Integer & x, unsigned long NoTests)
    {
        if (x < arithmetic::convert<arithmetic::Integer>(2))
            return arithmetic::convert<arithmetic::Integer>(2);
        arithmetic::Integer p = x;
        if (bit(p, 0) == 0)
            // Make odd
            --p;
        do
        {
            // Increase by 2
            p += arithmetic::convert<arithmetic::Integer>(2);
        } while (!isPrimeMR(p, NoTests));
        return p;
    }

    unsigned long nextPrime(unsigned long x, unsigned long NoTests)
    {
        if (x < 2)
            return 2;
        if ((x & 1) == 0)
            // Make odd
            --x;
        do
        {
            // Increase by 2
            x += 2;
        } while (!isPrimeMR(x, NoTests));
        return x;
    }
}

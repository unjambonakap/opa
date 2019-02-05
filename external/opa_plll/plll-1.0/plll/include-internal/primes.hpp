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

#ifndef PLLL_INCLUDE_GUARD__PRIMES_GENERATOR_HPP
#define PLLL_INCLUDE_GUARD__PRIMES_GENERATOR_HPP

#include <vector>
#include <plll/arithmetic.hpp>

namespace plll
{
    class PrimesGenerator
    {
    private:
        std::vector<unsigned long> d_primes;
    
        void sieve(unsigned long bound); // finds all primes until 2*bound+1
        void generateAnotherPrime(unsigned long NoTests = 100); // extends the list by one prime (with probability 1 - 4^{-NoTests})
    
    public:
        void generateUpTo(unsigned long upperBound);
        void generateNOP(unsigned long numberOfPrimes);
    
        void generateOneMore()
        {
            generateAnotherPrime();
        }
    
        unsigned long count() const
        {
            return d_primes.size();
        }
    
        unsigned long operator() (unsigned long i) const
        {
            return d_primes[i];
        }
    
        unsigned long last() const
        {
            return d_primes.back();
        }
    };
    
    bool isPrimeMR(unsigned long n, unsigned long NoTests = 100);
    // Tests whether n is prime using NoTests Miller-Rabin strong pseudoprime tests. If false is
    // returned, the number is compoiste. If true is returned, the number is prime with probability
    // 4^{-NoTests}.
    
    bool isPrimeMR(const arithmetic::Integer &, unsigned long NoTests = 100);
    // Tests whether n is prime using NoTests Miller-Rabin strong pseudoprime tests. If false is
    // returned, the number is compoiste. If true is returned, the number is prime with probability
    // 4^{-NoTests}.
    
    arithmetic::Integer nextPrime(const arithmetic::Integer &, unsigned long NoTests = 100);
    unsigned long nextPrime(unsigned long x, unsigned long NoTests = 100);
}

#endif

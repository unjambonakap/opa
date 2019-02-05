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

#ifndef PLLL_INCLUDE_GUARD__FACTOR_HPP
#define PLLL_INCLUDE_GUARD__FACTOR_HPP

#include <list>
#include <iostream>
#include <plll/arithmetic.hpp>

namespace plll
{
    class FactorizedNumber;

    void factorize(FactorizedNumber & result, const arithmetic::Integer & N);

    class FactorizedNumber
    // Stores a non-zero integer with its factorization
    {
        friend void factorize(FactorizedNumber & result, const arithmetic::Integer & N);
    
    public:
        struct PrimeExp
        {
            arithmetic::Integer p;
            long e;
        
            PrimeExp()
                : e(0)
            {
            }
            
            PrimeExp(const arithmetic::Integer & P)
                : p(P), e(1)
            {
            }
            
            PrimeExp(const arithmetic::Integer & P, long E)
                : p(P), e(E)
            {
            }
            
            bool operator < (const PrimeExp & b) const
            {
                return p < b.p;
            }
        };
    
    private:
        arithmetic::Integer d_N;
        std::list<PrimeExp> d_factorization;

        enum Mode { modMul, modGCD, modLCM };
    
        static void merge(std::list<PrimeExp> & res, const std::list<PrimeExp> & A, const std::list<PrimeExp> & B, Mode mode = modMul);
        static void optimize(std::list<PrimeExp> & res);
    
        FactorizedNumber(const arithmetic::Integer & N, const std::list<PrimeExp> & F)
            : d_N(N), d_factorization(F)
        {
        }
    
    public:
        FactorizedNumber()
        {
            setOne(d_N);
        }
    
        FactorizedNumber(const arithmetic::Integer & N)
        {
            factorize(*this, N);
        }
    
        FactorizedNumber(long N)
        {
            factorize(*this, arithmetic::convert<arithmetic::Integer>(N));
        }

        FactorizedNumber(const arithmetic::Integer & P, long e)
        {
            assert(e >= 0);
            // assert(isPrimeMR(P));
            d_factorization.push_back(PrimeExp(P, e));
            d_N = power(P, e);
        }
    
        FactorizedNumber(long P, long e)
        {
            assert(e >= 0);
            // assert(isPrimeMR(P));
            d_factorization.push_back(PrimeExp(arithmetic::convert<arithmetic::Integer>(P), e));
            d_N = power(arithmetic::convert<arithmetic::Integer>(P), e);
        }

        const arithmetic::Integer & number() const
        {
            return d_N;
        }

        const std::list<PrimeExp> & factorization() const
        {
            return d_factorization;
        }
    
        bool operator == (const FactorizedNumber & FN) const
        {
            return d_N == FN.d_N;
        }
    
        bool operator != (const FactorizedNumber & FN) const
        {
            return d_N != FN.d_N;
        }
    
        bool operator < (const FactorizedNumber & FN) const
        {
            return d_N < FN.d_N;
        }
    
        bool operator <= (const FactorizedNumber & FN) const
        {
            return d_N <= FN.d_N;
        }
    
        bool operator > (const FactorizedNumber & FN) const
        {
            return d_N > FN.d_N;
        }
    
        bool operator >= (const FactorizedNumber & FN) const
        {
            return d_N >= FN.d_N;
        }
    
        friend int sign(const FactorizedNumber & FN)
        {
            return sign(FN.d_N);
        }
    
        friend FactorizedNumber abs(const FactorizedNumber & FN)
        {
            return FactorizedNumber(abs(FN.d_N), FN.d_factorization);
        }
    
        friend void makeAbs(FactorizedNumber & FN)
        {
            makeAbs(FN.d_N);
        }
    
        FactorizedNumber operator * (const FactorizedNumber & FN) const
        {
            FactorizedNumber res;
            merge(res.d_factorization, d_factorization, FN.d_factorization, modMul);
            res.d_N = d_N * FN.d_N;
            return res;
        }
    
        FactorizedNumber & operator *= (const FactorizedNumber & FN)
        {
            std::list<PrimeExp> fct;
            merge(fct, d_factorization, FN.d_factorization);
            fct.swap(d_factorization);
            d_N *= FN.d_N;
            return *this;
        }
    
        FactorizedNumber & mulCoprimeLarger(const FactorizedNumber & FN); // multiply with a coprime number whose prime factors are all larger
    
        friend FactorizedNumber GCD(const FactorizedNumber & A, const FactorizedNumber & B)
        {
            FactorizedNumber res;
            GCD(res, A, B);
            return res;
        }
    
        friend void GCD(FactorizedNumber & FN, const FactorizedNumber & A, const FactorizedNumber & B)
        {
            std::list<PrimeExp> fct;
            merge(fct, A.d_factorization, B.d_factorization, modGCD);
            fct.swap(FN.d_factorization);
            FN.d_N = GCD(A.d_N, B.d_N);
        }
    
        friend FactorizedNumber LCM(const FactorizedNumber & A, const FactorizedNumber & B)
        {
            FactorizedNumber res;
            LCM(res, A, B);
            return res;
        }
    
        friend void LCM(FactorizedNumber & FN, const FactorizedNumber & A, const FactorizedNumber & B)
        {
            std::list<PrimeExp> fct;
            merge(fct, A.d_factorization, B.d_factorization, modLCM);
            fct.swap(FN.d_factorization);
            FN.d_N = (A.d_N / GCD(A.d_N, B.d_N)) * B.d_N;
            makeAbs(FN.d_N);
        }
    
        friend void power(FactorizedNumber & res, const FactorizedNumber & FN, long e)
        {
            res.d_N = power(FN.d_N, e);
            res.d_factorization = FN.d_factorization;
            for (std::list<PrimeExp>::iterator i = res.d_factorization.begin(); i != res.d_factorization.end(); ++i)
                i->e *= e;
        }
    
        friend FactorizedNumber power(const FactorizedNumber & FN, long e)
        {
            FactorizedNumber res;
            power(res, FN, e);
            return res;
        }
    
        friend std::ostream & operator << (std::ostream & s, const FactorizedNumber & FN)
        {
            s << FN.d_N << " " << "(";
            if (FN.d_factorization.empty())
                s << "1";
            for (std::list<PrimeExp>::const_iterator i = FN.d_factorization.begin(); i != FN.d_factorization.end(); ++i)
            {
                if (i != FN.d_factorization.begin())
                    s << ", ";
                s << "[" << i->p << "^" << i->e << "]";
            }
            s << ")";
            return s;
        }
    };

    arithmetic::Integer computeNthRoot(const arithmetic::Integer & N, long n);
    // Assume N is non-negative. Returns an integer N' such that (N')^n <= N < (N' + 1)^n.

    void factorize(FactorizedNumber & result, const arithmetic::Integer & N);
    FactorizedNumber factorize(const arithmetic::Integer & N);

    class ChineseRemainder
    {
    private:
        struct CRTInfo
        {
            arithmetic::Integer value;
            arithmetic::Integer basiselement; // an element which is congruent to 1 modulo this prime power, and
            // congruent to 0 modulo the other prime powers
        };
    
        FactorizedNumber d_modulus; // the moduli are the prime powers
        std::list<CRTInfo> d_values; // the values to be CRT'ed
    
        void makeInfo();
    
    public:
        class Iterator
        {
            friend class ChineseRemainder;
        
        private:
            std::list<FactorizedNumber::PrimeExp>::const_iterator it;
            std::list<CRTInfo>::iterator it2;
        
            Iterator(const std::list<FactorizedNumber::PrimeExp>::const_iterator & _it, const std::list<CRTInfo>::iterator & _it2)
                : it(_it), it2(_it2)
            {
            }
        public:
            
            Iterator & operator ++()
            {
                ++it;
                ++it2;
                return *this;
            }
            
            Iterator operator ++(int)
            {
                Iterator cpy = *this;
                ++it;
                ++it2;
                return cpy;
            }
            
            const arithmetic::Integer & prime() const
            {
                return it->p;
            }
            
            const long exponent() const
            {
                return it->e;
            }
            
            void setValue(const arithmetic::Integer & v)
            {
                it2->value = v;
            }
            
            bool operator == (const Iterator & b) const
            {
                return (it == b.it) && (it2 == b.it2);
            }
            
            bool operator != (const Iterator & b) const
            {
                return (it != b.it) || (it2 != b.it2);
            }
        };
    
        ChineseRemainder(const arithmetic::Integer & N)
            : d_modulus(N)
        {
            makeInfo();
        }
    
        ChineseRemainder(const FactorizedNumber & N)
            : d_modulus(N)
        {
            makeInfo();
        }
    
        const FactorizedNumber & modulus() const
        {
            return d_modulus;
        }
    
        Iterator begin()
        {
            return Iterator(d_modulus.factorization().begin(), d_values.begin());
        }
    
        Iterator end()
        {
            return Iterator(d_modulus.factorization().end(), d_values.end());
        }
    
        arithmetic::Integer process(bool signedRes = true) const;
        // Do Chinese Remainder. If signedRes is true, result will be in (-N/2, N/2]. Otherwise, result will be in [0, N).
    };
}

#endif

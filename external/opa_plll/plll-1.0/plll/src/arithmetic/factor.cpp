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

#include <plll/config.hpp>
#include "factor.hpp"
#include "primes.hpp"
#include <iostream>

namespace plll
{
    FactorizedNumber & FactorizedNumber::mulCoprimeLarger(const FactorizedNumber & FN)
    // multiply with a coprime number whose prime factors are all larger
    {
        assert(&FN != this);
        d_N *= FN.d_N;
        d_factorization.insert(d_factorization.end(), FN.d_factorization.begin(), FN.d_factorization.end());
        return *this;
    }

    void FactorizedNumber::merge(std::list<PrimeExp> & res, const std::list<PrimeExp> & A, const std::list<PrimeExp> & B, Mode mode)
    {
        std::list<PrimeExp>::const_iterator a = A.begin(), b = B.begin();
        // Combine lists
        while (a != A.end() && b != B.end())
        {
            if (a->p == b->p)
            {
                long e = a->e;
                if (mode == modMul)
                    e += b->e;
                else if (mode == modGCD)
                    e = std::min(e, b->e);
                else if (mode == modLCM)
                    e = std::max(e, b->e);
                res.push_back(PrimeExp(a->p, e));
                ++a;
                ++b;
            }
            else if (a->p < b->p)
            {
                if (mode != modGCD)
                    res.push_back(*a);
                ++a;
            }
            else
            {
                if (mode != modGCD)
                    res.push_back(*b);
                ++b;
            }
        }
        // Insert leftovers
        if (mode != modGCD)
        {
            res.insert(res.end(), a, A.end());
            res.insert(res.end(), b, B.end());
        }
    }
    
    void FactorizedNumber::optimize(std::list<PrimeExp> & res)
    {
        // Firstly, sort by prime
        res.sort();
        // Secondly, remove double entries
        for (std::list<PrimeExp>::iterator i = res.begin(); ; ++i)
        {
            std::list<PrimeExp>::iterator n = i;
            ++n;
            if (n == res.end())
                break;
            if (i->p == n->p)
            {
                i->e += n->e;
                res.erase(n);
            }
            else
                i = n;
        }
    }

    arithmetic::Integer computeNthRoot(const arithmetic::Integer & N, long n)
    // Assume N is non-negative. Returns an integer N' such that (N')^n <= N < (N' + 1)^n.
    {
        arithmetic::Integer res;
        long bit = (ceilOfLog2(N) + n - 1) / n;
        arithmetic::Integer b;
        setOne(b);
        b <<= bit;
        while (bit)
        {
            // Check whether bit is set in res
            if (power(res + b, n) <= N)
                res += b;
            // Go to next bit
            --bit;
            b >>= 1;
        }
        return res;
    }

    // !!! ??? ... EC-Curve?!
/*
struct ECMCurve
{
    struct Point
    {
        bool isInfinity;
        NTL::ZZ_p x, y;
        
        Point(const NTL::ZZ_p & X, const NTL::ZZ_p & Y)
            : isInfinity(false), x(X), y(Y)
        {
        }
        
        Point()
            : isInfinity(true)
        {
        }
        
        bool operator == (const Point & B) const
        {
            if (isInfinity)
                return B.isInfinity;
            if (B.isInfinity)
                return false;
            return (x == B.x) && (y == B.y);
        }
        
        bool isNegativeOf(const Point & B) const
        {
            if (isInfinity)
                return B.isInfinity;
            if (B.isInfinity)
                return false;
            return (x == B.x) && (y == -B.y);
        }
        
    };
    
    NTL::ZZ_p A, B;
    Point P;

    // internal
    NTL::ZZ_p xd, x;
    NTL::ZZ a, b;
    
    inline Point Double(const Point & P, bool & failure, NTL::ZZ & factor)
    {
        if (P.isInfinity)
            return P;
        if (NTL::IsZero(P.y))
            return Point();
        xd = P.y + P.y;
        XGCD(factor, a, b, rep(xd), NTL::ZZ_p::modulus());
        if (!NTL::IsOne(factor))
        {
            failure = true;
            return P;
        }
        xd = (NTL::to_ZZ_p(3) * sqr(P.x) - P.y) * NTL::to_ZZ_p(a);
        x = sqr(xd) - P.x - P.x;
        return Point(x, P.y + xd * (x - P.x));
    }
    
    inline Point Add(const Point & P, const Point & Q, bool & failure, NTL::ZZ & factor)
    {
        if (P.isInfinity)
            return Q;
        if (Q.isInfinity)
            return P;
        bool equal = false;
        if (P.x == Q.x)
        {
            if (P.y == -Q.y)
                return Point();
            equal = true;
        }
        xd = equal ? (P.y + P.y) : (P.x - Q.x);
        XGCD(factor, a, b, rep(xd), NTL::ZZ_p::modulus());
        if (!NTL::IsOne(factor))
        {
            failure = true;
            return P;
        }
        xd = equal ? ((NTL::to_ZZ_p(3) * sqr(P.x) - P.y) * NTL::to_ZZ_p(a)) : ((P.y - Q.y) * NTL::to_ZZ_p(a));
        x = sqr(xd) - P.x - Q.x;
        return Point(x, P.y + xd * (x - P.x));
    }
    
    Point Mul(const Point & P, unsigned fact, bool & failure, NTL::ZZ & factor);
};

ECMCurve::Point ECMCurve::Mul(const Point & P, unsigned fact, bool & failure, NTL::ZZ & factor)
{
    Point Q = P;
    Point R;
    if (fact & 1)
        R = P;
    fact >>= 1;
    while (fact)
    {
        Q = Double(Q, failure, factor);
        if (failure || Q.isInfinity)
            break;
        if (fact & 1)
            R = Add(R, Q, failure, factor);
        if (failure)
            break;
        fact >>= 1;
    }
    return R;
    
    / *
    if (fact == 1)
        return P;
    Point Q = Double(P, failure, factor);
    if (failure)
        return Q;
    Point R = Mul(Q, fact >> 1, failure, factor);
    if (((fact % 1) == 0) || failure)
        return R;
    return Add(R, P, failure, factor);
    * /
}

bool FactorECM(const NTL::ZZ n, NTL::ZZ & factor, const std::list<int> & smallPrimeList)
// <n> should be coprime to 2 and 3. If it can factor <n>, stores a non-trivial factor in <factor>
// and returns true. Otherwise, returns false.
// Assume that NTL::ZZ_p::init(n) has been called before.
{
    std::list<ECMCurve> curves;
    for (int i = 0; i < 20; ++i)
    {
        // Select a random curve E : y^2 = x^3 + A x + B and a point P = (Px, Py) on E
        ECMCurve curve;
        NTL::ZZ_p Px, Py;
        while (true)
        {
            curve.A = NTL::random_ZZ_p();
            Px = NTL::random_ZZ_p();
            Py = NTL::random_ZZ_p();
            curve.B = sqr(Py) - sqr(Px) * Px - curve.A * Px;
            // Check if curve is smooth
            NTL::ZZ_p Delta = NTL::to_ZZ_p(4) * sqr(curve.A) * curve.A + NTL::to_ZZ_p(27) * sqr(curve.B);
            factor = GCD(rep(Delta), n);
            if (NTL::IsOne(factor))
                break;
            else
                if (!NTL::IsZero(factor))
                    return true;
        };
        curve.P = ECMCurve::Point(Px, Py);
//        std::cout << i << ": trying (" << Px << "," << Py << ") on curve y^2 = x^3 + " << curve.A << " x + " << curve.B << "\n";
        curves.push_back(curve);
    }
    int lp = 0;
    for (std::list<int>::const_iterator i = smallPrimeList.begin(); i != smallPrimeList.end(); ++i)
    {
        if (*i / 10000 != lp)
        {
            lp = *i / 10000;
//            std::cout << lp << " / " << smallPrimeList.back() / 10000 << "\r";
//            std::cout.flush();
        }
        int e = lp < 1 ? 4 : (lp < 2 ? 3 : (lp < 3 ? 2 : 1));
        for (std::list<ECMCurve>::iterator c = curves.begin(); c != curves.end(); )
        {
            for (int j = 0; j < e; ++j)
            {
                unsigned p = *i;
                // compute (Px, Py) = p * (Px, Py);
                bool failure = false;
                c->P = c->Mul(c->P, p, failure, factor);
                if (failure)
                {
                    if (!NTL::IsZero(factor))
                        return true;
                    else
                    {
//                        std::cout << "removing curve (zero)\n";
                        // kick this curve out
                        c = curves.erase(c);
                        continue;
                    }
                }
                if (c->P.isInfinity)
                {
//                    std::cout << "removing curve (identity)\n";
                    // kick this curve out
                    c = curves.erase(c);
                    continue;
                }
            }
            ++c;
        }
    }
    for (std::list<ECMCurve>::iterator c = curves.begin(); c != curves.end(); ++c)
    {
        if (c->P.isInfinity)
            continue;
        // Check if P has `bad coordinates'
        factor = GCD(rep(c->P.x), n);
        if ((!NTL::IsOne(factor)) && (!NTL::IsZero(factor)))
            return true;
        factor = GCD(rep(c->P.y), n);
        if ((!NTL::IsOne(factor)) && (!NTL::IsZero(factor)))
            return true;
    }
    return false;
}
*/

    arithmetic::Integer findFactor(const arithmetic::Integer N)
    // Given: positive integer N, not prime, not a perfect power, not divisible by "small" primes
    // Returns: non-trivial factor of N
    {
        // Create `small prime list'
        static std::list<int> smallPrimeList;
        if (smallPrimeList.begin() == smallPrimeList.end())
        {
            unsigned sievesize = (2 << 26);
            PrimesGenerator pg;
            pg.generateUpTo(sievesize);
            for (unsigned long i = 0; i < pg.count(); ++i)
                smallPrimeList.push_back(pg(i));
        }
        // Try to factor using ECM
        /*
          NTL::ZZ res;
          NTL::ZZ_pBak backup;
          backup.save();
          NTL::ZZ_p::init(N);
          while (!FactorECM(N, res, smallPrimeList))
          ;
          backup.restore();
          return res;
        */
        return arithmetic::Integer(1l);
    }
    
    void factorize(FactorizedNumber & result, const arithmetic::Integer & N)
    {
        assert(!isZero(N));
    
        result.d_N = arithmetic::convert<arithmetic::Integer>(sign(N));
        result.d_factorization.clear();
        arithmetic::Integer n = N;
        makeAbs(n);

        if (isOne(n))
        {
            // Is it one?
        }
        else if (isPrimeMR(n))
        {
            // Is it a prime?
            result.d_N *= n;
            result.d_factorization.push_back(FactorizedNumber::PrimeExp(n, 1));
        }
        else
        {
            // Is it a perfect power?
            unsigned long bound = ceilOfLog2(n);
            unsigned long primebound = std::max((unsigned long)100, bound);
            static PrimesGenerator pg;
            if ((pg.count() == 0) || ((pg.count() != 0) && (pg.last() <= primebound)))
                pg.generateUpTo(primebound);
            for (unsigned long i = 0; i < pg.count(); ++i)
            {
                if (pg(i) > bound)
                    break;
                arithmetic::Integer R = computeNthRoot(n, pg(i));
                if (power(R, pg(i)) == n)
                {
                    factorize(result, R);
                    power(result, result, pg(i));
                    if (sign(N) < 0)
                        result.d_N = -result.d_N;
                    return;
                }
            }
        
            // Apparently, it is not. Check first primes
            bool foundFactor = false;
            for (unsigned long i = 0; i < pg.count(); ++i)
            {
                arithmetic::Integer p = arithmetic::convert<arithmetic::Integer>(pg(i));
                long e = 0;
                while (isZero(n % p))
                {
                    ++e;
                    n /= p;
                }
                if (e)
                {
                    result.d_N *= power(p, e);
                    result.d_factorization.push_back(FactorizedNumber::PrimeExp(arithmetic::convert<arithmetic::Integer>(pg(i)), e));
                    foundFactor = true;
                }
            }
        
            // Found something? Recurse...
            if (foundFactor)
            {
                FactorizedNumber fn;
                factorize(fn, n);
                result *= fn;
                return;
            }

            // Just search for something with some algorithm...
            arithmetic::Integer f = findFactor(n);
            FactorizedNumber fn1, fn2;
            factorize(fn1, f);
            factorize(fn2, n / f);
            result *= fn1;
            result *= fn2;
        }
    }

    FactorizedNumber factorize(const arithmetic::Integer & N)
    {
        FactorizedNumber res;
        factorize(res, N);
        return res;
    }

    void ChineseRemainder::makeInfo()
    {
        for (std::list<FactorizedNumber::PrimeExp>::const_iterator i = d_modulus.factorization().begin(); i != d_modulus.factorization().end(); ++i)
        {
            CRTInfo info;
            setZero(info.value);
            arithmetic::Integer x, y, pp = power(i->p, i->e), rest = d_modulus.number() / pp, d;
            XGCD(d, x, y, pp, rest);
            // d = 1 = x * pp + y * rest
            info.basiselement = y * d_modulus.number() / pp;
            info.basiselement %= d_modulus.number();
            d_values.push_back(info);
        }
    }

    arithmetic::Integer ChineseRemainder::process(bool signedRes) const
    // Do Chinese Remainder. If signedRes is true, result will be in (-N/2, N/2]. Otherwise, result will be in [0, N).
    {
        arithmetic::Integer res;
        for (std::list<CRTInfo>::const_iterator i = d_values.begin(); i != d_values.end(); ++i)
            res += (i->value * i->basiselement) % d_modulus.number();
        res %= d_modulus.number();
        if (sign(res) < 0)
            res += d_modulus.number();
        if (signedRes)
            if (res > (d_modulus.number() >> 1))
                res -= d_modulus.number();
        return res;
    }
}

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

#include <plll/rational.hpp>
#include <cctype>

namespace plll
{
    namespace arithmetic
    {
        void power(Rational & r, const Rational & a, const Integer & e)
        {
            assert(sign(e) >= 0);
            Rational aa(a);
            setOne(r);
            unsigned bitlen = bitLength(e);
            for (unsigned b = 0; b < bitlen; ++b)
            {
                if (bit(e, b))
                    r *= aa;
                if (b < bitlen - 1)
                    square(aa, aa);
            }
        }
        
        std::ostream & operator << (std::ostream & s, const Rational & r)
        {
            s << "(" << r.d_num << "/" << r.d_denom << ")";
            return s;
        }
        
        std::istream & operator >> (std::istream & s, Rational & r)
        {
            char c = ' ';
            while (s && std::isspace(c))
                s >> c;
            if (c != '(')
            {
                s.setstate(std::ios_base::badbit);
                return s;
            }
            s >> r.d_num;
            if (!s)
            {
                setZero(r);
                return s;
            }
            s >> c;
            if (c != '/')
            {
                s.setstate(std::ios_base::badbit);
                return s;
            }
            s >> r.d_denom;
            if (!s || isZero(r.d_denom))
            {
                setZero(r);
                return s;
            }
            r.normalize();
            if (sign(r.d_denom) < 0)
            {
                makeAbs(r.d_denom);
                neg(r.d_num, r.d_num);
            }
            s >> c;
            if (c != ')')
            {
                s.setstate(std::ios_base::badbit);
                return s;
            }
            return s;
        }
    }
}

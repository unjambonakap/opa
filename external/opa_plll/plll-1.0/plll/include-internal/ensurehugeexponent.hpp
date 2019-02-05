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

#ifndef PLLL_INCLUDE_GUARD__ENSURE_HUGE_EXPONENT_HPP
#define PLLL_INCLUDE_GUARD__ENSURE_HUGE_EXPONENT_HPP

namespace plll
{
    namespace arithmetic
    {
        template<class RealTypeContext, bool has_huge> class HugeExponentImpl;
    
        template<class RealTypeContext>
        class HugeExponent
        {
        public:
            typedef typename HugeExponentImpl<RealTypeContext, RealTypeContext::has_huge_exponent>::Type Type;
        };
    
        template<class RealTypeContext>
        class HugeExponentImpl<RealTypeContext, true>
        {
        public:
            typedef typename RealTypeContext::Real Type;
        };
    
        template<class RealTypeContext>
        class HugeExponentAdder;
    
        template<class RealTypeContext>
        class HugeExponentImpl<RealTypeContext, false>
        {
        public:
            typedef HugeExponentAdder<RealTypeContext> Type;
        };
    
        template<class RealTypeContext>
        class HugeExponentAdder;
    
        template<class RealTypeContext>
        inline HugeExponentAdder<RealTypeContext> operator << (const HugeExponentAdder<RealTypeContext> & a, long e);
        template<class RealTypeContext>
        inline HugeExponentAdder<RealTypeContext> operator >> (const HugeExponentAdder<RealTypeContext> & a, long e);
        template<class RealTypeContext>
        inline HugeExponentAdder<RealTypeContext> operator << (const HugeExponentAdder<RealTypeContext> & a, const HugeExponentAdder<RealTypeContext> & e);
        template<class RealTypeContext>
        inline HugeExponentAdder<RealTypeContext> operator >> (const HugeExponentAdder<RealTypeContext> & a, const HugeExponentAdder<RealTypeContext> & e);
        template<class RealTypeContext>
        inline bool operator == (const HugeExponentAdder<RealTypeContext> & a, const HugeExponentAdder<RealTypeContext> & b);
        template<class RealTypeContext>
        inline bool operator != (const HugeExponentAdder<RealTypeContext> & a, const HugeExponentAdder<RealTypeContext> & b);
        template<class RealTypeContext>
        inline bool operator <= (const HugeExponentAdder<RealTypeContext> & a, const HugeExponentAdder<RealTypeContext> & b);
        template<class RealTypeContext>
        inline bool operator >= (const HugeExponentAdder<RealTypeContext> & a, const HugeExponentAdder<RealTypeContext> & b);
        template<class RealTypeContext>
        inline bool operator < (const HugeExponentAdder<RealTypeContext> & a, const HugeExponentAdder<RealTypeContext> & b);
        template<class RealTypeContext>
        inline bool operator > (const HugeExponentAdder<RealTypeContext> & a, const HugeExponentAdder<RealTypeContext> & b);
    
        template<class RealTypeContext>
        inline bool isOne(const HugeExponentAdder<RealTypeContext> &);
        // Tests whether the given value is = 1.
    
        template<class RealTypeContext>
        inline bool isPositive(const HugeExponentAdder<RealTypeContext> &);
        // Tests whether the given value is > 0.
    
        template<class RealTypeContext>
        inline bool isNegative(const HugeExponentAdder<RealTypeContext> &);
        // Tests whether the given value is < 0.
    
        template<class RealTypeContext>
        inline void mul(HugeExponentAdder<RealTypeContext> & r, const HugeExponentAdder<RealTypeContext> & a, const HugeExponentAdder<RealTypeContext> & b);
        template<class RealTypeContext>
        inline void div(HugeExponentAdder<RealTypeContext> & r, const HugeExponentAdder<RealTypeContext> & a, const HugeExponentAdder<RealTypeContext> & b);
        template<class RealTypeContext>
        inline void shl(HugeExponentAdder<RealTypeContext> & r, const HugeExponentAdder<RealTypeContext> & a, long b);
        template<class RealTypeContext>
        inline void shr(HugeExponentAdder<RealTypeContext> & r, const HugeExponentAdder<RealTypeContext> & a, long b);
        template<class RealTypeContext>
        inline void neg(HugeExponentAdder<RealTypeContext> & r, const HugeExponentAdder<RealTypeContext> & a);
        template<class RealTypeContext>
        inline void abs(HugeExponentAdder<RealTypeContext> & r, const HugeExponentAdder<RealTypeContext> & a);
    
        template<class RealTypeContext>
        inline HugeExponentAdder<RealTypeContext> abs(const HugeExponentAdder<RealTypeContext> &);
        template<class RealTypeContext>
        inline void makeAbs(HugeExponentAdder<RealTypeContext> &);
    
        template<class RealTypeContext>
        inline void swap(HugeExponentAdder<RealTypeContext> & a, HugeExponentAdder<RealTypeContext> & b);
    
        template<class RealTypeContext>
        inline void setOne(HugeExponentAdder<RealTypeContext> &);
        // Set to one
    
        template<class RealTypeContext>
        inline int compare(const HugeExponentAdder<RealTypeContext> &, const HugeExponentAdder<RealTypeContext> &);
        // Tests whether the first integer is < (-1), = (0) or > (1) than the second.
        template<class RealTypeContext>
        inline int compareAbsValues(const HugeExponentAdder<RealTypeContext> &, const HugeExponentAdder<RealTypeContext> &);
        // Tests whether the absolute value of the first integer is < (-1), = (0) or > (1) than the
        // absolute value of the second integer.
        template<class RealTypeContext>
        inline int sign(const HugeExponentAdder<RealTypeContext> &);
        // Returns the sign (i.e. -1, 0, 1).
    
        template<class RealTypeContext>
        inline HugeExponentAdder<RealTypeContext> square(const HugeExponentAdder<RealTypeContext> &);
        template<class RealTypeContext>
        inline void square(HugeExponentAdder<RealTypeContext> &, const HugeExponentAdder<RealTypeContext> &);
        // Some basic arithmetic operations
    
        template<class RealTypeContext>
        inline HugeExponentAdder<RealTypeContext> sqrt(const HugeExponentAdder<RealTypeContext> & r);
        template<class RealTypeContext>
        inline void sqrt(HugeExponentAdder<RealTypeContext> & res, const HugeExponentAdder<RealTypeContext> & r);
        template<class RealTypeContext>
        inline HugeExponentAdder<RealTypeContext> power(const HugeExponentAdder<RealTypeContext> &, signed long);
        template<class RealTypeContext>
        inline HugeExponentAdder<RealTypeContext> power(const HugeExponentAdder<RealTypeContext> &, unsigned long);
        template<class RealTypeContext>
        inline HugeExponentAdder<RealTypeContext> power(const HugeExponentAdder<RealTypeContext> &, const Integer &);
        template<class RealTypeContext>
        inline void power(HugeExponentAdder<RealTypeContext> &, const HugeExponentAdder<RealTypeContext> &, signed long);
        template<class RealTypeContext>
        inline void power(HugeExponentAdder<RealTypeContext> &, const HugeExponentAdder<RealTypeContext> &, unsigned long);
        template<class RealTypeContext>
        inline void power(HugeExponentAdder<RealTypeContext> &, const HugeExponentAdder<RealTypeContext> &, const Integer &);
    
        template<class RealTypeContext>
        std::ostream & operator << (std::ostream &, const HugeExponentAdder<RealTypeContext> &);
    
        template<class RealTypeContext>
        class HugeExponentAdder
        {
#ifdef PLLL_INTERNAL_NO_TEMPLATE_FRIENDS
        public:
#else
        private:
#endif
            typename RealTypeContext::Real d_value;
            long d_exponent;
        
            inline void normalize()
            {
                d_exponent += d_value.exponentNormalize();
            }
        
            inline HugeExponentAdder(const typename RealTypeContext::Real & r, long e)
                : d_value(r), d_exponent(e)
            {
            }
        
            inline HugeExponentAdder(const typename RealTypeContext::Real & r, long e, bool)
                : d_value(r), d_exponent(e)
            {
                normalize();
            }
        
        public:
            inline HugeExponentAdder()
                : d_exponent(0)
            {
                setOne(d_value);
            }
        
            explicit inline HugeExponentAdder(const typename RealTypeContext::Real & r)
                : d_value(r), d_exponent(0)
            {
                normalize();
            }
        
            inline HugeExponentAdder(const HugeExponentAdder & r)
                : d_value(r.d_value), d_exponent(r.d_exponent)
            {
            }
        
            inline HugeExponentAdder(const RealTypeContext & rc)
                : d_value(rc), d_exponent(0)
            {
                setOne(d_value);
            }
        
            inline void setContext(const RealTypeContext & rc)
            {
                d_value.setContext(rc);
            }
        
            inline HugeExponentAdder operator - () const
            {
                return HugeExponentAdder(-d_value, d_exponent);
            }
        
            inline HugeExponentAdder operator * (const HugeExponentAdder<RealTypeContext> & b) const
            {
                return HugeExponentAdder(d_value * b.d_value, d_exponent + b.d_exponent, true);
            }
        
            inline HugeExponentAdder operator / (const HugeExponentAdder<RealTypeContext> & b) const
            {
                return HugeExponentAdder(d_value / b.d_value, d_exponent - b.d_exponent, true);
            }
        
#ifndef PLLL_INTERNAL_NO_TEMPLATE_FRIENDS
            friend HugeExponentAdder<RealTypeContext> operator << <>(const HugeExponentAdder<RealTypeContext> & a, long e);
            friend HugeExponentAdder<RealTypeContext> operator >> <>(const HugeExponentAdder<RealTypeContext> & a, long e);
#endif
        
            inline HugeExponentAdder & operator = (const HugeExponentAdder & r)
            {
                d_value = r.d_value;
                d_exponent = r.d_exponent;
                return *this;
            }
        
            inline HugeExponentAdder & operator *= (const HugeExponentAdder & r)
            {
                d_value *= r.d_value;
                d_exponent += r.d_exponent;
                normalize();
                return *this;
            }
        
            inline HugeExponentAdder & operator /= (const HugeExponentAdder & r)
            {
                d_value /= r.d_value;
                d_exponent -= r.d_exponent;
                normalize();
                return *this;
            }
        
            inline HugeExponentAdder & operator <<= (long r)
            {
                shl(*this, *this, r);
                return *this;
            }
        
            inline HugeExponentAdder & operator >>= (long r)
            {
                shr(*this, *this, r);
                return *this;
            }

#ifndef PLLL_INTERNAL_NO_TEMPLATE_FRIENDS
            friend bool operator == <>(const HugeExponentAdder<RealTypeContext> & a, const HugeExponentAdder<RealTypeContext> & b);
            friend bool operator != <>(const HugeExponentAdder<RealTypeContext> & a, const HugeExponentAdder<RealTypeContext> & b);
            friend bool operator <= <>(const HugeExponentAdder<RealTypeContext> & a, const HugeExponentAdder<RealTypeContext> & b);
            friend bool operator >= <>(const HugeExponentAdder<RealTypeContext> & a, const HugeExponentAdder<RealTypeContext> & b);
            friend bool operator < <>(const HugeExponentAdder<RealTypeContext> & a, const HugeExponentAdder<RealTypeContext> & b);
            friend bool operator > <>(const HugeExponentAdder<RealTypeContext> & a, const HugeExponentAdder<RealTypeContext> & b);
            friend bool isOne<>(const HugeExponentAdder<RealTypeContext> & r);
            // Tests whether the given value is = 1.
            friend bool isPositive<>(const HugeExponentAdder<RealTypeContext> & r);
            // Tests whether the given value is > 0.
            friend bool isNegative<>(const HugeExponentAdder<RealTypeContext> & r);
            // Tests whether the given value is < 0.
            friend void mul<>(HugeExponentAdder<RealTypeContext> & r, const HugeExponentAdder<RealTypeContext> & a, const HugeExponentAdder<RealTypeContext> & b);
            friend void div<>(HugeExponentAdder<RealTypeContext> & r, const HugeExponentAdder<RealTypeContext> & a, const HugeExponentAdder<RealTypeContext> & b);
            friend void shl<>(HugeExponentAdder<RealTypeContext> & r, const HugeExponentAdder<RealTypeContext> & a, long b);
            friend void shr<>(HugeExponentAdder<RealTypeContext> & r, const HugeExponentAdder<RealTypeContext> & a, long b);
            friend void neg<>(HugeExponentAdder<RealTypeContext> & r, const HugeExponentAdder<RealTypeContext> & a);
            friend void abs<>(HugeExponentAdder<RealTypeContext> & r, const HugeExponentAdder<RealTypeContext> & a);
            friend HugeExponentAdder<RealTypeContext> abs<>(const HugeExponentAdder<RealTypeContext> & a);
            friend void makeAbs<>(HugeExponentAdder<RealTypeContext> & a);
            friend void setOne<>(HugeExponentAdder<RealTypeContext> & r);
            // Set to one
            friend int compare<>(const HugeExponentAdder<RealTypeContext> & a, const HugeExponentAdder<RealTypeContext> & b);
            // Tests whether the first real is < (-1), = (0) or > (1) than the second.
            friend int compareAbsValues<>(const HugeExponentAdder<RealTypeContext> & a, const HugeExponentAdder<RealTypeContext> & b);
            // Tests whether the absolute value of the first HugeExponentAdder is < (-1), = (0) or > (1) than the
            // absolute value of the second HugeExponentAdder.
            friend int sign<>(const HugeExponentAdder<RealTypeContext> & r);
            // Returns the sign (i.e. -1, 0, 1).
            friend HugeExponentAdder<RealTypeContext> square<>(const HugeExponentAdder<RealTypeContext> & a);
            friend void square<>(HugeExponentAdder<RealTypeContext> & r, const HugeExponentAdder<RealTypeContext> & a);
            friend HugeExponentAdder<RealTypeContext> sqrt<>(const HugeExponentAdder<RealTypeContext> & r);
            friend void sqrt<>(HugeExponentAdder<RealTypeContext> & res, const HugeExponentAdder<RealTypeContext> & r);
            friend HugeExponentAdder<RealTypeContext> power<>(const HugeExponentAdder<RealTypeContext> & a, signed long e);
            friend HugeExponentAdder<RealTypeContext> power<>(const HugeExponentAdder<RealTypeContext> & a, unsigned long e);
            friend HugeExponentAdder<RealTypeContext> power<>(const HugeExponentAdder<RealTypeContext> & a, const Integer & e);
            friend void power<>(HugeExponentAdder<RealTypeContext> & r, const HugeExponentAdder<RealTypeContext> & a, signed long e);
            friend void power<>(HugeExponentAdder<RealTypeContext> & r, const HugeExponentAdder<RealTypeContext> & a, unsigned long e);
            friend void power<>(HugeExponentAdder<RealTypeContext> & r, const HugeExponentAdder<RealTypeContext> & a, const Integer & e);
            friend void swap<>(HugeExponentAdder<RealTypeContext> & a, HugeExponentAdder<RealTypeContext> & b);
            friend std::ostream & operator << <>(std::ostream & s, const HugeExponentAdder<RealTypeContext> & r);
#endif
        
            long getExponent() const // number is +-f * 2^e, where e is the exponent and 0.5 < f <= 1
            {
                return d_exponent;
            }
        
            void setExponent(long e)
            {
                d_exponent = e;
            }
        
            long exponentNormalize() // moves to range (0.5, 1], returns exponent
            {
                long e = d_exponent;
                d_exponent = 0;
                return e;
            }
        };

        template<class RealTypeContext>
        std::ostream & operator << (std::ostream & s, const HugeExponentAdder<RealTypeContext> & r)
        // Output to stream
        {
            s << r.d_value << "*2^";
            if (r.d_exponent < 0)
                s << "(";
            s << r.d_exponent;
            if (r.d_exponent < 0)
                s << ")";
            s << " (" << (r.d_value << r.d_exponent) << ")";
            return s;
        }
        
        template<class RealTypeContext>
        inline HugeExponentAdder<RealTypeContext> operator << (const HugeExponentAdder<RealTypeContext> & a, long e)
        {
            HugeExponentAdder<RealTypeContext> r;
            shl(r, a, e);
            return r;
        }
    
        template<class RealTypeContext>
        inline HugeExponentAdder<RealTypeContext> operator >> (const HugeExponentAdder<RealTypeContext> & a, long e)
        {
            HugeExponentAdder<RealTypeContext> r;
            shr(r, a, e);
            return r;
        }
    
        template<class RealTypeContext>
        inline bool operator == (const HugeExponentAdder<RealTypeContext> & a, const HugeExponentAdder<RealTypeContext> & b)
        {
            return (a.d_value == b.d_value) && (a.d_exponent == b.d_exponent);
        }
        
        template<class RealTypeContext>
        inline bool operator != (const HugeExponentAdder<RealTypeContext> & a, const HugeExponentAdder<RealTypeContext> & b)
        {
            return (a.d_value != b.d_value) || (a.d_exponent != b.d_exponent);
        }
        
        template<class RealTypeContext>
        inline bool operator <= (const HugeExponentAdder<RealTypeContext> & a, const HugeExponentAdder<RealTypeContext> & b)
        {
            return compare(a, b) <= 0;
        }
        
        template<class RealTypeContext>
        inline bool operator >= (const HugeExponentAdder<RealTypeContext> & a, const HugeExponentAdder<RealTypeContext> & b)
        {
            return compare(a, b) >= 0;
        }
        
        template<class RealTypeContext>
        inline bool operator < (const HugeExponentAdder<RealTypeContext> & a, const HugeExponentAdder<RealTypeContext> & b)
        {
            return compare(a, b) < 0;
        }
    
        template<class RealTypeContext>
        inline bool operator > (const HugeExponentAdder<RealTypeContext> & a, const HugeExponentAdder<RealTypeContext> & b)
        {
            return compare(a, b) > 0;
        }
    
        template<class RealTypeContext>
        inline bool isOne(const HugeExponentAdder<RealTypeContext> & r)
        // Tests whether the given value is = 1.
        {
            return (r.d_exponent == 0) && (r.d_value == (RealTypeContext)1);
        }
    
        template<class RealTypeContext>
        inline bool isPositive(const HugeExponentAdder<RealTypeContext> & r)
        // Tests whether the given value is > 0.
        {
            return r.d_value > (RealTypeContext)0;
        }
    
        template<class RealTypeContext>
        inline bool isNegative(const HugeExponentAdder<RealTypeContext> & r)
        // Tests whether the given value is < 0.
        {
            return r.d_value < (RealTypeContext)0;
        }
        
        template<class RealTypeContext>
        inline void mul(HugeExponentAdder<RealTypeContext> & r, const HugeExponentAdder<RealTypeContext> & a, const HugeExponentAdder<RealTypeContext> & b)
        {
            r.d_value = a.d_value * b.d_value;
            r.d_exponent = a.d_exponent + b.d_exponent;
            r.normalize();
        }
        
        template<class RealTypeContext>
        inline void div(HugeExponentAdder<RealTypeContext> & r, const HugeExponentAdder<RealTypeContext> & a, const HugeExponentAdder<RealTypeContext> & b)
        {
            r.d_value = a.d_value / b.d_value;
            r.d_exponent = a.d_exponent - b.d_exponent;
            r.normalize();
        }
        
        template<class RealTypeContext>
        inline void shl(HugeExponentAdder<RealTypeContext> & r, const HugeExponentAdder<RealTypeContext> & a, long b)
        {
            r.d_value = a.d_value;
            r.d_exponent = a.d_exponent + b;
        }
        
        template<class RealTypeContext>
        inline void shr(HugeExponentAdder<RealTypeContext> & r, const HugeExponentAdder<RealTypeContext> & a, long b)
        {
            r.d_value = a.d_value;
            r.d_exponent = a.d_exponent - b;
        }
        
        template<class RealTypeContext>
        inline void neg(HugeExponentAdder<RealTypeContext> & r, const HugeExponentAdder<RealTypeContext> & a)
        {
            r.d_value = -a.d_value;
            r.d_exponent = a.d_exponent;
        }
        
        template<class RealTypeContext>
        inline void abs(HugeExponentAdder<RealTypeContext> & r, const HugeExponentAdder<RealTypeContext> & a)
        {
            r.d_value = abs(a.d_value);
            r.d_exponent = a.d_exponent;
        }
    
        template<class RealTypeContext>
        inline HugeExponentAdder<RealTypeContext> abs(const HugeExponentAdder<RealTypeContext> & a)
        {
            return HugeExponentAdder<RealTypeContext>(abs(a.d_value), a.d_exponent);
        }
        
        template<class RealTypeContext>
        inline void makeAbs(HugeExponentAdder<RealTypeContext> & a)
        {
            makeAbs(a.d_value);
        }
    
        template<class RealTypeContext>
        inline void setOne(HugeExponentAdder<RealTypeContext> & r)
        // Set to one
        {
            setOne(r.d_value);
            r.d_exponent = 0;
        }
    
        template<class RealTypeContext>
        inline int compare(const HugeExponentAdder<RealTypeContext> & a, const HugeExponentAdder<RealTypeContext> & b)
        // Tests whether the first real is < (-1), = (0) or > (1) than the second.
        {
            bool na = isNegative(a.d_value), nb = isNegative(b.d_value);
            if (na != nb)
                return na ? -1 : +1;
            int v = compareAbsValues(a, b);
            return na ? -v : v;
        }
    
        template<class RealTypeContext>
        inline int compareAbsValues(const HugeExponentAdder<RealTypeContext> & a, const HugeExponentAdder<RealTypeContext> & b)
        // Tests whether the absolute value of the first HugeExponentAdder is < (-1), = (0) or > (1) than the
        // absolute value of the second HugeExponentAdder.
        {
            if (a.d_exponent != b.d_exponent)
                return (a.d_exponent < b.d_exponent) ? -1 : 1;
            return compareAbsValues(a.d_value, b.d_value);
        }
    
        template<class RealTypeContext>
        inline int sign(const HugeExponentAdder<RealTypeContext> & r)
        // Returns the sign (i.e. -1, 0, 1).
        {
            return r.d_value < 0 ? -1 : (r.d_value > 0 ? 1 : 0);
        }
    
        template<class RealTypeContext>
        inline HugeExponentAdder<RealTypeContext> square(const HugeExponentAdder<RealTypeContext> & a)
        {
            return HugeExponentAdder<RealTypeContext>(a.d_value * a.d_value, a.d_exponent * 2, true);
        }
    
        template<class RealTypeContext>
        inline void square(HugeExponentAdder<RealTypeContext> & r, const HugeExponentAdder<RealTypeContext> & a)
        {
            mul(r.d_value, a.d_value, a.d_value);
            r.d_exponent = 2 * a.d_exponent;
            r.normalize();
        }
        
        template<class RealTypeContext>
        inline HugeExponentAdder<RealTypeContext> sqrt(const HugeExponentAdder<RealTypeContext> & r)
        {
            long odd = (r.d_exponent & 1) ? 1 : 0;
            return HugeExponentAdder<RealTypeContext>(sqrt(odd ? (r.d_value >> 1) : r.d_value), (r.d_exponent - odd) / 2, true);
        }
    
        template<class RealTypeContext>
        inline void sqrt(const HugeExponentAdder<RealTypeContext> & r, const HugeExponentAdder<RealTypeContext> & a)
        {
            long odd = (a.d_exponent & 1) ? 1 : 0;
            if (odd)
            {
                shr(r.d_value, a.d_value, 1);
                sqrt(r.d_value, r.d_value);
            }
            else
                sqrt(r.d_value, a.d_value);
            r.d_exponent = (r.d_exponent - odd) / 2;
            r.normalize();
        }
    
        template<class RealTypeContext>
        inline HugeExponentAdder<RealTypeContext> power(const HugeExponentAdder<RealTypeContext> & a, signed long e)
        {
            return HugeExponentAdder<RealTypeContext>(power(a.d_value, e), a.d_exponent * e, true);
        }
    
        template<class RealTypeContext>
        inline HugeExponentAdder<RealTypeContext> power(const HugeExponentAdder<RealTypeContext> & a, unsigned long e)
        {
            return HugeExponentAdder<RealTypeContext>(power(a.d_value, e), a.d_exponent * e, true);
        }
        
        template<class RealTypeContext>
        inline HugeExponentAdder<RealTypeContext> power(const HugeExponentAdder<RealTypeContext> & a, const Integer & e)
        {
            return HugeExponentAdder<RealTypeContext>(power(a.d_value, e), a.d_exponent * convert<long>(e), true);
        }
    
        template<class RealTypeContext>
        inline void power(HugeExponentAdder<RealTypeContext> & r, const HugeExponentAdder<RealTypeContext> & a, signed long e)
        {
            power(r.d_value, a.d_value, e);
            r.d_exponent = a.d_exponent * e;
        }
        
        template<class RealTypeContext>
        inline void power(HugeExponentAdder<RealTypeContext> & r, const HugeExponentAdder<RealTypeContext> & a, unsigned long e)
        {
            power(r.d_value, a.d_value, e);
            r.d_exponent = a.d_exponent * e;
        }
        
        template<class RealTypeContext>
        inline void power(HugeExponentAdder<RealTypeContext> & r, const HugeExponentAdder<RealTypeContext> & a, const Integer & e)
        {
            power(r.d_value, a.d_value, e);
            r.d_exponent = a.d_exponent * convert<long>(e);
        }
    
        template<class RealTypeContext>
        inline void swap(HugeExponentAdder<RealTypeContext> & a, HugeExponentAdder<RealTypeContext> & b)
        {
            arithmetic::swap(a.d_value, b.d_value);
            std::swap(a.d_exponent, b.d_exponent);
        }
    }
}

#endif

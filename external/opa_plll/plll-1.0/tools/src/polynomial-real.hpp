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

#ifndef PLLL_INCLUDE_GUARD__POLYNOMIAL_REAL_HPP
#define PLLL_INCLUDE_GUARD__POLYNOMIAL_REAL_HPP

#include <iostream>
#include <vector>
#include <memory>
#include <plll/arithmetic.hpp>

namespace plll
{
    namespace arithmetic
    {
    }
}

template<int n, class T, class Alloc = std::allocator<T> >
class realpoly;
/*
  Represents a polynomial in n indeterminates with coefficients from T. This is implemented
  recursively: realpoly<n, T> is a polynomial in one indeterminate with coefficients in realpoly<n-1,
  T>. The last instance, realpoly<0, T>, is essentially equal to T.
  
  In fact, the complete polynomial type is realpoly<n, T, A>, where A is an allocator for T, for example
  std::allocator<T>; this is also the default parameter. Usually this parameter does not needs to be
  changed.
  
  Let us describe the interface:
  
  * Every polynomial has a predicate isZero(), testing whether the polynomial is the zero
    polynomial.
    
  * Moreover, it has basic arithmetic, i.e. addition, subtraction, multiplication; this is realized
    via operator overloading. Assignment and comparism operators are also provided.
    
  * Polynomials have a degree function degree(), returning the degree as a polynomial in the first
    indeterminate, as well as a function leading(), returning the leading term (which is of type
    const realpoly<n-1, T> &). Note that the degree of the zero polynomial is -1.
    
  * The i-th coefficient can be obtained by writing f[i] if f is of type realpoly<n, T>; then f[i] is of
    type realpoly<n - 1, T>. In case i is larger than the degree of f, a zero polynomial will be
    returned, and in case f is not constant, the internal space for the coefficients will be
    enlarged to have space for the i-th coefficient. Note that one should call f.normalize()
    afterwards if one accessed f[i] for non-const f to ensure that degree and leading coefficients
    will be computed correctly afterwards.
    
  * The polynomial can be evaluated using operator(). If f is of type realpoly<n, T, A>, then f(x) with
    x of type S will return an object of type realpoly_evaluator<n, T, A, S>. This object can be casted
    to realpoly<n-1, S> to obtain f with the first indeterminate evaluated as x. Note that S must be a
    type to which elements of T can be casted. One can also continue evaluating with the
    realpoly_evaluator<n, T, A, S> object, yielding realpoly_evaluator_impl<k, T, ..., A, S> objects, k <
    n. These object hierarchy ensures that polynomial evaluation is efficient, i.e. after (good
    enough) optimiziation of the compiler, is essentially a block of code of (n-k+1) nested
    for-loops.
    
  * Polynomials can be written to std::ostream's using operator<<. There is also a member function
    print(std::ostream &) const. Note that the first indeterminate in realpoly<n, T> is denoted by X_0,
    and the last by X_{n-1}. In case of n == 1, the indeterminate is just called X.
    
  * Finally, there is a method swap() to swap two polynomial's contents. The polynomials have to be
    of the same type.
    
  The coefficients are stored in a std::vector<>; while realpoly<1, T> uses std::vector<T>, realpoly<n, T>
  for n > 1 uses std::vector<realpoly<n-1, T> *>, i.e. pointers to coefficients are stored. This is
  implemented using the "intelligent" vector realpoly_ivector<> template.
  
  Creation can be accomplished using the functions Monomial<T, A>(e_0, ..., e_{n-1}), where A again
  is std::allocator<T> by default. This creates a polynomial of type realpoly<n, T, A> which consists of
  exactly one monomial X_0^{e_0} * ... * X_{n-1}^{e_{n-1}} having coefficient 1. One can create a
  polynomial for example using
      realpoly<2, int> f = Monomial<int>(1, 2) + 3 * Monomial<int>(4, 5);
  to obtain
      f = X_0 X_1^2 + 3 X_0^4 X_1^5.
  
  The library also offers to compute partial derivatives and partial integrals. Given a polynomial f
  of type realpoly<n, T>, one can compute the partial derivative with respect to the j-th variable of f
  by writing derivative<j>(f). Similarly, integral<j>(f) computes an integral with respect to the
  j-th variable of f. Note that T is required to allow multiplication by int's for derivation to
  work, and to allow division by non-zero int's for integration to work.
 */

// A class implementing partial derivatives. This is implemented using a class since C++ does not
// support partial template specialization for template functions. Taking partial derivatives is
// eventually implemented by a function, using this class template.
template<int N, int n, class T, class Alloc = std::allocator<T> >
class real_derivative_computer;

// A class implementing partial integrals. This is implemented using a class since C++ does not
// support partial template specialization for template functions. Taking partial integrals is
// eventually implemented by a function, using this class template.
template<int N, int n, class T, class Alloc = std::allocator<T> >
class real_integral_computer;

template<int n, class T, class Alloc = std::allocator<T>, class S = T>
// Helper template to evaluate a polynomial of type realpoly<n, T, Alloc> in a point of S. This template
// is returned from realpoly<n, T, Alloc>::operator().
// 
// Note that evaluation is more complex, since we want that the compiler optimizes evaluation to a
// set of nested for-loops. The original implementation was done recursively, returning a polynomial
// of type realpoly<n-1, T> evaluated in the first indeterminate, and than the operator() of realpoly<n-1,
// T> was called to evaluate in the next indeterminate, etc. This proved to be very ineffective (not
// very surprisingly). If evaluation is seldomly done, this is ok, but in my case, it was the main
// bottleneck, whence I implemented this more complicated, but also much more efficient solution.
// 
// ...
class realpoly_evaluator;

template<int n, class T, class HL, class Alloc, class S>
// Another helper for polynomial evaluation. This template is returned from realpoly_evaluator<n, T,
// Alloc>::operator().
// 
// The template parameter HL is the type of the "owner" of this template, i.e. either
// realpoly_evaluator<n+1, T, Alloc> or realpoly_evaluator_impl<n+1, T, ..., Alloc>.
class realpoly_evaluator_impl;

template<class T, class HL, class Alloc, class S>
class realpoly_evaluator_impl<1, T, HL, Alloc, S>
{
    template<int nn, class TT, class AA, class SS>
    friend class realpoly_evaluator;
    
    template<int nn, class TT, class HLHL, class AA, class SS>
    friend class realpoly_evaluator_impl;
    
private:
    const HL & d_owner; // The "owner"
    const S & d_evalpoint; // The evaluation point on this level
    
    inline realpoly_evaluator_impl(const HL & owner, const S & evalpoint)
        : d_owner(owner), d_evalpoint(evalpoint)
    {
    }
    
    class eval_fun
    // Evaluates a polynomial of type realpoly<1, T> in d_evalpoint.
    {
        const realpoly_evaluator_impl<1, T, HL, Alloc, S> & d_owner;
        
    public:
        inline eval_fun(const realpoly_evaluator_impl<1, T, HL, Alloc, S> & owner)
            : d_owner(owner)
        {
        }
        
        inline S operator() (const realpoly<1, T, Alloc> & p) const
        // Evaluate
        {
            if (p.d_value.size() > 1)
            {
                S res = (S)(T)p.d_value[p.d_value.size() - 1];
                for (unsigned i = p.d_value.size() - 1; i > 0; --i)
                {
                    res *= d_owner.d_evalpoint;
                    res += (S)(T)p.d_value[i - 1];
                }
                return res;
            }
            else if (p.d_value.size() == 1)
                return (S)(T)p.d_value[0];
            else
            {
                S res;
                setZero(res);
                return res;
            }
        }
    };
    
public:
    inline void evaluate_to(S & res) const
    {
        setZero(res);
        d_owner.evaluate(res, eval_fun(*this));
    }
    
    inline operator S() const
    // Cast to S
    {
        S res;
        evaluate_to(res);
        return res;
    }
    
    inline S operator() () const
    // Explicit evaluate. Essentially calls operator S().
    {
        return (S)(*this);
    }
};

template<int n, class T, class HL, class Alloc, class S>
class realpoly_evaluator_impl
{
    template<int nn, class TT, class AA, class SS>
    friend class realpoly_evaluator;
    
    template<int nn, class TT, class HLHL, class AA, class SS>
    friend class realpoly_evaluator_impl;
    
private:
    const HL & d_owner; // The "owner"
    const S & d_evalpoint; // The evaluation point on this level
    
    inline realpoly_evaluator_impl(const HL & owner, const S & evalpoint)
        : d_owner(owner), d_evalpoint(evalpoint)
    {
    }
    
    template<class SS, class Fun>
    class eval_fun
    // Evaluates a polynomial of type realpoly<n, T> in d_evalpoint; the coefficients of the first
    // indeterminate are evaluated using the given functor of type Fun.
    {
        const realpoly_evaluator_impl<n, T, HL, Alloc, S> & d_owner;
        const Fun & d_evalfun;
        
    public:
        inline eval_fun(const realpoly_evaluator_impl<n, T, HL, Alloc, S> & owner, const Fun & evalfun)
            : d_owner(owner), d_evalfun(evalfun)
        {
        }
        
        inline SS operator() (const realpoly<n, T, Alloc> & p) const
        // Evaluate: the first indeterminate is replaced by d_evalpoint, and the coefficients of the
        // first indeterminate (which are polynomials of type realpoly<n-1, T>) are evaluated using
        // d_evalfun.
        {
            if (p.d_value.size() > 1)
            {
                SS res = d_evalfun(p.d_value[p.d_value.size() - 1]);
                for (unsigned i = p.d_value.size() - 1; i > 0; --i)
                {
                    res *= d_owner.d_evalpoint;
                    res += d_evalfun(p.d_value[i - 1]);
                }
                return res;
            }
            else if (p.d_value.size() == 1)
                return d_evalfun(p.d_value[0]);
            else
                return SS();
        }
    };
    
    template<class SS, class Fun>
    inline void evaluate(SS & res, const Fun & evalfun) const
    // This will be called from a child (i.e. a class of type realpoly_evaluator_impl<n-1, T,
    // realpoly_evaluator<n,T,HL,Alloc,S>, Alloc, S>) to trigger evaluation.
    {
        // We have to pass evaluation on to our owner, but give a new functor which now evaluates
        // polynomials of type realpoly<n, T>.
        d_owner.evaluate(res, eval_fun<SS, Fun>(*this, evalfun));
    }
    
    class eval_fun2
    // Evaluates a polynomial of type realpoly<n, T> in d_evalpoint; the coefficients of the first
    // indeterminate (which are polynomials of type realpoly<n-1, T>) are casted to realpoly<n-1, S>, and
    // the result is of type realpoly<n-1, S> as well.
    {
        const realpoly_evaluator_impl<n, T, HL, Alloc, S> & d_owner;
        
    public:
        inline eval_fun2(const realpoly_evaluator_impl<n, T, HL, Alloc, S> & owner)
            : d_owner(owner)
        {
        }
        
        inline realpoly<n - 1, S, typename Alloc::template rebind<S>::other> operator() (const realpoly<n, T, Alloc> & p) const
        // Evaluate: the first indeterminate is replaced by d_evalpoint, and the coefficients of the
        // first indeterminate (which are polynomials of type realpoly<n-1, T>) are casted to the type
        // realpoly<n-1, S>.
        {
            if (p.d_value.size() > 1)
            {
                realpoly<n - 1, S, typename Alloc::template rebind<S>::other>
                    res = ((realpoly<n - 1, S, typename Alloc::template rebind<S>::other>)p.d_value[p.d_value.size() - 1]);
                for (unsigned i = p.d_value.size() - 1; i > 0; --i)
                {
                    res *= d_owner.d_evalpoint;
                    res += ((realpoly<n - 1, S, typename Alloc::template rebind<S>::other>)p.d_value[i - 1]);
                }
                return res;
            }
            else if (p.d_value.size() == 1)
                return (realpoly<n - 1, S, typename Alloc::template rebind<S>::other>)p.d_value[0];
            else
                return realpoly<n - 1, S, typename Alloc::template rebind<S>::other>(p.get_allocator());
        }
    };
    
public:
    inline void evaluate_to(realpoly<n - 1, S, typename Alloc::template rebind<S>::other> & res) const
    {
        setZero(res);
        d_owner.evaluate(res, eval_fun2(*this));
    }
    
    inline operator realpoly<n - 1, S, typename Alloc::template rebind<S>::other>() const
    // Allows casting to realpoly<n-1, S>.
    {
        realpoly<n - 1, S, typename Alloc::template rebind<S>::other> res; // missing: determine allocator object
        // We need to pass evaluation on to our owner
        d_owner.evaluate(res, eval_fun2(*this));
        return res;
    }
    
    template<class SS>
    inline realpoly_evaluator_impl<n - 1, T, realpoly_evaluator_impl<n, T, HL, Alloc, S>, Alloc, SS> operator() (const SS & x) const
    // Continues evaluation with the next indeterminant.
    {
        return realpoly_evaluator_impl<n - 1, T, realpoly_evaluator_impl<n, T, HL, Alloc, S>, Alloc, SS>(*this, x);
    }
};

template<class T, class Alloc, class S>
class realpoly_evaluator<1, T, Alloc, S>
// The top level polynomial evaluation class, in case n = 1, i.e. does direct evaluation.
{
    friend class realpoly<1, T, Alloc>;

private:
    const realpoly<1, T, Alloc> & d_poly; // The polynomial in question
    const S & d_evalpoint; // The evaluation point
    
    inline realpoly_evaluator(const realpoly<1, T, Alloc> & poly, const S & evalpoint)
        : d_poly(poly), d_evalpoint(evalpoint)
    {
    }
    
public:
    inline void evaluate_to(S & res) const
    {
        if (d_poly.d_value.size() > 1)
        {
            res = (S)(T)d_poly.d_value[d_poly.d_value.size() - 1];
            for (unsigned i = d_poly.d_value.size() - 1; i > 0; --i)
            {
                res *= d_evalpoint;
                res += (S)(T)d_poly.d_value[i - 1];
            }
        }
        else if (d_poly.d_value.size() == 1)
            res = (S)(T)d_poly.d_value[0];
        else
            setZero(res);
    }
    
    inline operator S() const
    // Casting to S is done by evaluation in d_evalpoint.
    {
        if (d_poly.d_value.size() > 1)
        {
            S res = (S)(T)d_poly.d_value[d_poly.d_value.size() - 1];
            for (unsigned i = d_poly.d_value.size() - 1; i > 0; --i)
            {
                res *= d_evalpoint;
                res += (S)(T)d_poly.d_value[i - 1];
            }
            return res;
        }
        else if (d_poly.d_value.size() == 1)
            return (S)(T)d_poly.d_value[0];
        else
        {
            S res;
            setZero(res);
            return res;
        }
    }
    
    inline S operator() () const
    // Evaluates and returns value in S.
    {
        return (S)(*this);
    }
};

template<int n, class T, class Alloc, class S>
class realpoly_evaluator
// The top level polynomial evaluation class, in case n > 1, i.e. in case the coefficients are
// polynomials by themselves.
{
    friend class realpoly<n, T, Alloc>;

    template<int nn, class TT, class HLHL, class AA, class SS>
    friend class realpoly_evaluator_impl;

private:
    const realpoly<n, T, Alloc> & d_poly; // The polynomial in question
    const S & d_evalpoint; // the evaluation point
    
    inline realpoly_evaluator(const realpoly<n, T, Alloc> & poly, const S & evalpoint)
        : d_poly(poly), d_evalpoint(evalpoint)
    {
    }
    
    template<class SS, class Fun>
    inline void evaluate(SS & res, const Fun & evalfun) const
    // Will be called by "child". Evaluates the polynomial into an element of SS (which must not
    // necessarily be S, it can also be realpoly<k, S> for some k < n - 1) using the given functor to
    // evaluate the coefficients, which are of type realpoly<n-1, T>.
    {
        if (d_poly.d_value.size() > 0)
        {
            res = evalfun(d_poly.d_value[d_poly.d_value.size() - 1]);
            for (unsigned i = d_poly.d_value.size() - 1; i > 0; --i)
            {
                res *= d_evalpoint;
                res += evalfun(d_poly.d_value[i - 1]);
            }
        }
        else
            res.d_value.resize(0);
    }
    
public:
    inline void evaluate_to(realpoly<n - 1, S, typename Alloc::template rebind<S>::other> & res) const
    {
        if (d_poly.d_value.size() > 1)
        {
            res = d_poly.d_value[d_poly.d_value.size() - 1];
            for (unsigned i = d_poly.d_value.size() - 1; i > 0; --i)
            {
                res *= d_evalpoint;
                res += realpoly<n - 1, S, typename Alloc::template rebind<S>::other>(d_poly.d_value[i - 1]);
            }
        }
        else if (d_poly.d_value.size() == 1)
            res = d_poly.d_value[0];
        else
            setZero(res);
    }
    
    inline operator realpoly<n - 1, S, typename Alloc::template rebind<S>::other >() const
    // Evaluate to polynomial of type realpoly<n-1, S>.
    {
        if (d_poly.d_value.size() > 1)
        {
            realpoly<n - 1, S, typename Alloc::template rebind<S>::other>
                res = realpoly<n - 1, S, typename Alloc::template rebind<S>::other>(d_poly.d_value[d_poly.d_value.size() - 1]);
            for (unsigned i = d_poly.d_value.size() - 1; i > 0; --i)
            {
                res *= d_evalpoint;
                res += realpoly<n - 1, S, typename Alloc::template rebind<S>::other>(d_poly.d_value[i - 1]);
            }
            return res;
        }
        else if (d_poly.d_value.size() == 1)
            return realpoly<n - 1, S, typename Alloc::template rebind<S>::other>(d_poly.d_value[0]);
        else
            return realpoly<n - 1, S, typename Alloc::template rebind<S>::other>(d_poly.get_allocator());
    }
    
    template<class SS>
    inline realpoly_evaluator_impl<n - 1, T, realpoly_evaluator<n, T, Alloc, S>, Alloc, SS> operator() (const SS & x) const
    // Continue evaluation to lower level.
    {
        return realpoly_evaluator_impl<n - 1, T, realpoly_evaluator<n, T, Alloc, S>, Alloc, SS>(*this, x);
    }
};

template<class T, class Alloc>
class realpoly<0, T, Alloc>
// Stores a polynomial of degree 0, i.e. a scalar of type T. We assume that the type T is not "too"
// complex, otherwise this class will be partially not very effective.
{
private:
    Alloc d_allocator;
    T d_value;

public:
    inline realpoly()
    {
        plll::arithmetic::setZero(d_value);
    }
    
    inline realpoly(const T & v, const Alloc & allocator = Alloc())
        : d_allocator(allocator), d_value(v)
    {
    }
    
    inline realpoly(const Alloc & allocator)
        : d_allocator(allocator), d_value(T())
    {
    }
    
    inline bool isZero() const
    {
        return plll::arithmetic::isZero(d_value);
    }
    
    inline void setZero()
    {
        plll::arithmetic::setZero(d_value);
    }
    
    inline int degree() const
    // Returns the degree of this polynomial. If this is the zero polynomial, the degree is -1.
    {
        return isZero() ? -1 : 0;
    }
    
    inline int total_degree() const
    // Returns the total degree of this polynomial. If this is the zero polynomial, the degree is -1.
    {
        return isZero() ? -1 : 0;
    }
    
    inline operator const T & () const
    {
        return d_value;
    }
    
    inline operator T & ()
    {
        return d_value;
    }
    
    inline realpoly & operator = (const T & v)
    {
        d_value = v;
        return *this;
    }

    inline T operator() () const // evaluate
    {
        return d_value;
    }
    
    inline realpoly operator * (const T & v) const
    {
        return realpoly(d_value * v);
    }
    
    inline realpoly operator / (const T & v) const
    {
        return realpoly(d_value / v);
    }
    
    inline realpoly operator + (const T & v) const
    {
        return realpoly(d_value + v);
    }
    
    inline realpoly operator - (const T & v) const
    {
        return realpoly(d_value - v);
    }
    
    inline realpoly operator - () const
    {
        return realpoly(-d_value);
    }
    
    inline realpoly & operator *= (const T & v)
    {
        d_value *= v;
        return *this;
    }
    
    inline realpoly & operator /= (const T & v)
    {
        d_value /= v;
        return *this;
    }
    
    inline realpoly & operator += (const T & v)
    {
        d_value += v;
        return *this;
    }
    
    inline realpoly & operator -= (const T & v)
    {
        d_value -= v;
        return *this;
    }
    
    inline bool operator == (const T & v) const
    {
        return d_value == v;
    }
    
    inline bool operator != (const T & v) const
    {
        return d_value != v;
    }
    
    void print(std::ostream & s, int N = 0) const
    {
        s << d_value;
    }
    
    inline void swap(realpoly & p)
    // Swaps two polynomials.
    {
        std::swap(d_value, p.d_value);
    }
    
    Alloc get_allocator() const
    {
        return d_allocator;
    }
};

// Next, we want to define the storage class realpoly_ivector<T, Alloc, usePointers>. It behaves like a
// subset of std::vector<T>'s capabilities (i.e. access elements, get size, set size, get last
// element, swap with other storage of same type), but uses std::vector<T*> in case usePointers is
// true.
// 
// The advantage of this approach is that if T is a more complex object, reallocation done with
// resize() can be very costly.

template<class T, class Alloc = std::allocator<T>, bool usePointers = false>
// The version just using std::vector<T>.
class realpoly_ivector
{
private:
    std::vector<T, Alloc> d_vec;
    
public:
    realpoly_ivector(const Alloc & allocator = Alloc())
        : d_vec(allocator)
    {
    }
    
    realpoly_ivector(unsigned size, const Alloc & allocator = Alloc())
        : d_vec(size, T(), allocator)
    {
    }
    
    realpoly_ivector(unsigned size, const T & entry, const Alloc & allocator = Alloc())
        : d_vec(size, entry, allocator)
    {
    }
    
    unsigned size() const
    {
        return d_vec.size();
    }
    
    void resize(unsigned size, const T & entry = T())
    {
        d_vec.resize(size, entry);
    }
    
    const T & operator[] (unsigned i) const
    {
        return d_vec[i];
    }
    
    T & operator[] (unsigned i)
    {
        return d_vec[i];
    }
    
    const T & back() const
    {
        return d_vec.back();
    }
    
    T & back()
    {
        return d_vec.back();
    }
    
    void swap(realpoly_ivector & v)
    {
        d_vec.swap(v.d_vec);
    }
    
    Alloc get_allocator() const
    {
        return d_vec.get_allocator();
    }
};

template<class T, class Alloc>
// The version using std::vector<T*>, but behaving like std::vector<T>.
class realpoly_ivector<T, Alloc, true>
{
private:
    Alloc d_allocator;
    std::vector<typename Alloc::pointer, typename Alloc::template rebind<typename Alloc::pointer>::other> d_vec;
    
    void create(unsigned begin, unsigned end, typename Alloc::const_reference entry)
    {
        for (unsigned i = begin; i < end; ++i)
        {
            d_vec[i] = d_allocator.allocate(sizeof(T));
            d_allocator.construct(d_vec[i], entry);
        }
    }
    
    void free(unsigned begin, unsigned end)
    {
        for (unsigned i = begin; i < end; ++i)
        {
            d_allocator.destroy(d_vec[i]);
            d_allocator.deallocate(d_vec[i], sizeof(T));
        }
    }
    
    template<class A>
    void copy_from(const std::vector<typename Alloc::pointer, A> & source)
    {
        for (unsigned i = 0; i < d_vec.size(); ++i)
        {
            d_vec[i] = d_allocator.allocate(sizeof(T));
            d_allocator.construct(d_vec[i], *source[i]);
        }
    }
    
public:
    realpoly_ivector(const Alloc & allocator = Alloc())
        : d_allocator(allocator), d_vec(allocator)
    {
    }
    
    realpoly_ivector(unsigned size)
        : d_allocator(Alloc()), d_vec(size)
    {
        create(0, size, T());
    }
    
    realpoly_ivector(unsigned size, const T & entry, const Alloc & allocator = Alloc())
        : d_allocator(allocator), d_vec(size)
    {
        create(0, size, entry);
    }
    
    realpoly_ivector(const realpoly_ivector & v)
        : d_vec(v.size())
    {
        copy_from(v.d_vec);
    }
    
    ~realpoly_ivector()
    {
        free(0, d_vec.size());
    }

    realpoly_ivector & operator = (const realpoly_ivector & v)
    {
        if (&v != this)
        {
            free(0, d_vec.size());
            d_vec.resize(v.size());
            copy_from(v.d_vec);
        }
        return *this;
    }
    
    unsigned size() const
    {
        return d_vec.size();
    }
    
    void resize(unsigned size, const T & entry)
    {
        unsigned oldsize = d_vec.size();
        if (oldsize > size)
            free(size, oldsize);
        d_vec.resize(size);
        if (oldsize < size)
            create(oldsize, size, entry);
    }
    
    void resize(unsigned size)
    {
        unsigned oldsize = d_vec.size();
        if (oldsize > size)
            free(size, oldsize);
        d_vec.resize(size);
        if (oldsize < size)
            create(oldsize, size, T());
    }
    
    const T & operator[] (unsigned i) const
    {
        return *d_vec[i];
    }
    
    T & operator[] (unsigned i)
    {
        return *d_vec[i];
    }
    
    const T & back() const
    {
        return *d_vec.back();
    }
    
    T & back()
    {
        return *d_vec.back();
    }
    
    void swap(realpoly_ivector & v)
    {
        d_vec.swap(v.d_vec);
    }
    
    Alloc get_allocator() const
    {
        return d_vec.get_allocator();
    }
};

template<int n, class T, class Alloc>
class make_univariate_computer;

template<int n, class T, class Alloc>
class realpoly
// The main polynomial class. Stores a polynomial in n indeterminates, with n > 0.
{
    template<int NN, int nn, class TT, class AA>
    friend class real_derivative_computer;
    
    template<int NN, int nn, class TT, class AA>
    friend class real_integral_computer;
    
    template<int nn, class TT, class AA, class SS>
    friend class realpoly_evaluator;
    
    template<int nn, class TT, class HLHL, class AA, class SS>
    friend class realpoly_evaluator_impl;
    
    template<int nn, class TT, class AA>
    friend class make_univariate_computer;
    
private:
    // We want to make sure that arithmetic::initArithmeticThreadAllocators() is called as seldom as
    // possible. In fact, we want it to be called before d_zero is initialized only.
    class MemoryInit
    {
    public:
        MemoryInit()
        {
            plll::arithmetic::initArithmeticThreadAllocators();
        }
    };
    
    class ZeroContainer
    {
    private:
        MemoryInit d_mem_init;
    public:
        realpoly<n - 1, T, Alloc> d_zero;
    };
    
    static ZeroContainer d_zero;
    realpoly_ivector<realpoly<n - 1, T, Alloc>, typename Alloc::template rebind<realpoly<n - 1, T, Alloc> >::other, (n>1)> d_value;
    
    inline realpoly(bool, unsigned s, const Alloc & allocator)
        : d_value(s) // initialize array of length s
    {
    }
    
public:
    inline void normalize()
    // Adjusts the size of d_value that the leading term and degree can be computed trivially. This
    // must be called only after calls to the non-const operator[], in which the degree of the
    // polynomial has potentially been changed.
    {
        unsigned dp1 = d_value.size();
        while (dp1)
        {
            if (d_value[dp1 - 1].isZero())
                --dp1;
            else
                break;
        }
        if (dp1 < d_value.size())
            d_value.resize(dp1);
    }
    
    inline realpoly(const Alloc & allocator = Alloc())
    // Constructs a zero polynomial
        : d_value(allocator)
    {
    }
    
    inline realpoly(const T & v, const Alloc & allocator = Alloc())
    // Constructs a constant polynomial with constant term v.
        : d_value(1, realpoly<n - 1, T, Alloc>(v), allocator)
    {
    }
    
    inline realpoly(const realpoly<n - 1, T, Alloc> & pdm1, const Alloc & allocator = Alloc())
    // Constructs a polynomial of type realpoly<n, T> which is initialized with a polynomial of type realpoly<n-1, T>.
        : d_value(1, pdm1)
    {
    }
    
    template<class S, class A>
    inline realpoly(const realpoly<n, S, A> & p, const Alloc & allocator = Alloc())
    // Casts a polynomial of type realpoly<n, S> to a polynomial of type realpoly<n, T>.
        : d_value(p.degree() + 1, realpoly<n - 1, T, Alloc>(), allocator)
    {
        for (unsigned i = 0; i < d_value.size(); ++i)
            d_value[i] = p[i];
        normalize();
    }
    
    template<class S, class A>
    inline realpoly & operator = (const realpoly<n, S, A> & p)
    // Copies a polynomial of type realpoly<n, S> to this polynomial (of type realpoly<n, T>).
    {
        d_value.resize(p.degree() + 1);
        for (unsigned i = 0; i < d_value.size(); ++i)
            d_value[i] = p[i];
        normalize();
        return *this;
    }
    
    inline int degree() const
    // Returns the degree of this polynomial. If this is the zero polynomial, the degree is -1.
    {
        return d_value.size() - 1;
    }
    
    inline int total_degree() const
    // Returns the total degree of this polynomial. If this is the zero polynomial, the degree is -1.
    {
        if (n == 1)
            return d_value.size() - 1;
        else
        {
            // case n > 1: the coefficients are (multivariate) polynomials themselves
            int result = -1;
            for (unsigned i = 0; i < d_value.size(); ++i)
            {
                int r = d_value[i].total_degree();
                if (r >= 0)
                    if (r + (int)i > result)
                        result = r + i;
            }
            return result;
        }
    }
    
    inline const realpoly<n - 1, T, Alloc> & leading() const
    // Returns the leading term (of type realpoly<n-1, T>) of the first indeterminate. Returns 0 (of
    // type realpoly<n-1, T>) in case of the zero polynomial.
    {
        return d_value.size() ? d_value.back() : d_zero.d_zero;
    }
    
    inline bool isZero() const
    // Tests whether this polynomial is the zero polynomial.
    {
        return d_value.size() == 0;
    }
    
    inline void setZero()
    {
        d_value.resize(0);
    }
    
    inline realpoly<n - 1, T, Alloc> & operator[] (unsigned i)
    // Returns a reference to the i-th coefficient. If i > degree(), the array d_value is
    // enlarged. Afterwards, one should better call normalize() to make sure future operations are
    // correct.
    {
        if (i >= d_value.size())
            d_value.resize(i + 1);
        return d_value[i];
    }
    
    inline const realpoly<n - 1, T, Alloc> & operator[] (unsigned i) const
    // Returns a reference to the i-th coefficient, or zero if i > degree().
    {
        return i < d_value.size() ? d_value[i] : d_zero.d_zero;
    }

    inline realpoly_evaluator<n, T, Alloc, T> operator() (const T & x) const
    // Evaluate in x.
    {
        return realpoly_evaluator<n, T, Alloc, T>(*this, x);
    }

    template<class S>
    inline realpoly_evaluator<n, T, Alloc, S> operator() (const S & x) const
    // Evaluate in x of type S.
    {
        return realpoly_evaluator<n, T, Alloc, S>(*this, x);
    }
    
    inline realpoly operator * (const T & v) const
    // Multiply by constant.
    {
        realpoly r(*this);
        for (unsigned i = 0; i < d_value.size(); ++i)
            r[i] *= v;
        return r;
    }
    
    inline realpoly operator / (const T & v) const
    // Division by constant.
    {
        realpoly r(*this);
        for (unsigned i = 0; i < d_value.size(); ++i)
            r[i] /= v;
        return r;
    }
    
    inline realpoly & operator *= (const realpoly & p)
    // Multiplication by a polynomial.
    {
        return *this = *this * p;
    }
    
    inline realpoly & operator *= (const T & v)
    // Multiplication by a constant.
    {
        realpoly r(*this);
        for (unsigned i = 0; i < d_value.size(); ++i)
            d_value[i] *= v;
        return *this;
    }
    
    inline realpoly & operator /= (const T & v)
    // Division by a constant.
    {
        for (unsigned i = 0; i < d_value.size(); ++i)
            d_value[i] /= v;
        return *this;
    }
    
    friend inline realpoly operator * (const T & v, const realpoly & p)
    // Multiplication by a constant from the left.
    {
        realpoly r(p);
        for (unsigned i = 0; i < p.d_value.size(); ++i)
            r[i] *= v;
        return r;
    }
    
    inline realpoly operator - () const
    // Returns the additive inverse of the polynomial.
    {
        realpoly r(true, d_value.size(), get_allocator());
        for (unsigned i = 0; i < d_value.size(); ++i)
            r[i] = -d_value[i];
        return r;
    }
    
    inline realpoly operator + (const realpoly & q) const
    // Computes the sum of two polynomials.
    {
        realpoly r(*this);
        if (r.d_value.size() < q.d_value.size())
            r.d_value.resize(q.d_value.size());
        for (unsigned i = 0; i < q.d_value.size(); ++i)
            r[i] += q[i];
        r.normalize();
        return r;
    }
    
    inline realpoly operator - (const realpoly & q) const
    // Computes the difference of two polynomials.
    {
        realpoly r(*this);
        if (r.d_value.size() < q.d_value.size())
            r.d_value.resize(q.d_value.size());
        for (unsigned i = 0; i < q.d_value.size(); ++i)
            r[i] -= q[i];
        r.normalize();
        return r;
    }
    
    inline realpoly & operator += (const realpoly & q)
    // Adds q to this polynomial.
    {
        if (d_value.size() < q.d_value.size())
            d_value.resize(q.d_value.size());
        for (unsigned i = 0; i < q.d_value.size(); ++i)
            d_value[i] += q[i];
        normalize();
        return *this;
    }
    
    inline realpoly & operator -= (const realpoly & q)
    // Subtracts q from this polynomial.
    {
        if (d_value.size() < q.d_value.size())
            d_value.resize(q.d_value.size());
        for (unsigned i = 0; i < q.d_value.size(); ++i)
            d_value[i] -= q[i];
        normalize();
        return *this;
    }
    
    inline realpoly operator + (const realpoly<n - 1, T, Alloc> & q) const
    // Computes the sum of this polynomial with one of degree one less.
    {
        realpoly r(*this);
        if (r.d_value.size() < 1)
            r.d_value.resize(1);
        r[0] += q;
        r.normalize();
        return r;
    }
    
    inline realpoly operator - (const realpoly<n - 1, T, Alloc> & q) const
    // Computes the difference of this polynomial with one of degree one less.
    {
        realpoly r(*this);
        if (r.d_value.size() < 1)
            r.d_value.resize(1);
        r[0] -= q;
        r.normalize();
        return r;
    }
    
    inline realpoly operator + (const T & v) const
    // Computes the sum with a constant.
    {
        realpoly r(*this);
        if (r.d_value.size() < 1)
            r.d_value.resize(1);
        r[0] += v;
        r.normalize();
        return r;
    }
    
    inline realpoly operator - (const T & v) const
    // Computes the difference to a constant.
    {
        realpoly r(*this);
        if (r.d_value.size() < 1)
            r.d_value.resize(1);
        r[0] -= v;
        r.normalize();
        return r;
    }
    
    inline realpoly & operator += (const realpoly<n - 1, T, Alloc> & q)
    // Adds a polynomial of degree n-1.
    {
        if (d_value.size() < 1)
            d_value.resize(1);
        d_value[0] += q;
        normalize();
        return *this;
    }
    
    inline realpoly & operator -= (const realpoly<n - 1, T, Alloc> & q)
    // Subtracts a polynomial of degree n-1.
    {
        if (d_value.size() < 1)
            d_value.resize(1);
        d_value[0] -= q;
        normalize();
        return *this;
    }
    
    inline realpoly & operator += (const T & v)
    // Adds a constant.
    {
        if (d_value.size() < 1)
            d_value.resize(1);
        d_value[0] += v;
        normalize();
        return *this;
    }
    
    inline realpoly & operator -= (const T & v)
    // Subtracts a constant.
    {
        if (d_value.size() < 1)
            d_value.resize(1);
        d_value[0] -= v;
        normalize();
        return *this;
    }
    
    inline bool operator == (const realpoly & q) const
    // Compares two polynomials.
    {
        if (d_value.size() != q.d_value.size())
            return false;
        for (unsigned i = 0; i < d_value.size(); ++i)
            if (d_value[i] != q[i])
                return false;
        return true;
    }
    
    inline bool operator != (const realpoly & q) const
    // Compares two polynomials for inequality.
    {
        return !(*this == q);
    }
    
    inline bool operator == (const T & v) const
    // Compares a polynomial to a constant.
    {
        if (isZero(v) && (d_value.size() == 0))
            return true;
        if (d_value.size() != 1)
            return false;
        return d_value[0] == v;
    }
    
    inline bool operator != (const T & v) const
    // Compares a polynomial to a constant for inequality.
    {
        return !(*this == v);
    }
    
    void print(std::ostream & s, int N = n) const
    // Prints this polynomial to the stream s. N gives the number of variables; this is needed for
    // recursive printing.
    {
        if (isZero())
        {
            T v;
            plll::arithmetic::setZero(v);
            s << v;
        }
        else
        {
            unsigned nonzero = 0;
            for (unsigned i = 0; i < d_value.size(); ++i)
                if (!d_value[i].isZero())
                    ++nonzero;
            if (nonzero > 1)
                s << "(";
            bool first = true;
            for (unsigned i = 0; i < d_value.size(); ++i)
                if (!d_value[i].isZero())
                {
                    if (first)
                        first = false;
                    else
                        s << " + ";
                    d_value[i].print(s, N);
                    if (i > 0)
                    {
                        s << " ";
                        s << "X";
                        if ((n != N) || (n != 1))
                            s << "_" << N - n;
                        if (i > 1)
                            s << "^" << i;
                    }
                    
                }
            if (nonzero > 1)
                s << ")";
        }
    }
    
    friend inline std::ostream & operator << (std::ostream & s, const realpoly & p)
    // Outputs p to s using print().
    {
        p.print(s);
        return s;
    }
    
    inline realpoly operator * (const realpoly & p) const
    // Multiplies the polynomial with another polynomial.
    {
        int d = p.degree() + degree();
        if (d < -1)
            d = -1;
        realpoly r(true, d + 1, get_allocator());
        if (!isZero() && !p.isZero())
            for (unsigned i = 0; i < r.d_value.size(); ++i)
                for (unsigned j = 0; j < d_value.size(); ++j)
                    if (i < j + p.d_value.size())
                        r[i] += d_value[j] * p[i - j];
        r.normalize();
        return r;
    }
    
    inline void swap(realpoly & p)
    // Swaps two polynomials.
    {
        d_value.swap(p.d_value);
    }
    
    Alloc get_allocator() const
    {
        return d_value.get_allocator();
    }
};

template<int n, class T, class Alloc>
inline void setZero(realpoly<n, T, Alloc> & p)
{
    p.setZero();
}

template<class T, class Alloc>
inline realpoly<1, T, Alloc> RealMonomial(unsigned e)
// Creates a monomial in one indeterminate.
{
    realpoly<1, T, Alloc> p;
    setOne(p[e]);
    return p;
}

template<class T, class Alloc>
inline realpoly<2, T, Alloc> RealMonomial(unsigned e, unsigned f)
// Creates a monomial in two indeterminates.
{
    realpoly<2, T, Alloc> p;
    setOne(p[e][f]);
    return p;
}

template<class T, class Alloc>
inline realpoly<3, T, Alloc> RealMonomial(unsigned e, unsigned f, unsigned g)
// Creates a monomial in three indeterminates.
{
    realpoly<3, T, Alloc> p;
    setOne(p[e][f][g]);
    return p;
}

template<class T, class Alloc>
inline realpoly<4, T, Alloc> RealMonomial(unsigned e, unsigned f, unsigned g, unsigned h)
// Creates a monomial in four indeterminates.
{
    realpoly<4, T, Alloc> p;
    setOne(p[e][f][g][h]);
    return p;
}

template<class T>
inline realpoly<1, T, std::allocator<T> > RealMonomial(unsigned e)
// Creates a monomial in one indeterminate.
{
    realpoly<1, T, std::allocator<T> > p;
    setOne(p[e]);
    return p;
}

template<class T>
inline realpoly<2, T, std::allocator<T> > RealMonomial(unsigned e, unsigned f)
// Creates a monomial in two indeterminates.
{
    realpoly<2, T, std::allocator<T> > p;
    setOne(p[e][f]);
    return p;
}

template<class T>
inline realpoly<3, T, std::allocator<T> > RealMonomial(unsigned e, unsigned f, unsigned g)
// Creates a monomial in three indeterminates.
{
    realpoly<3, T, std::allocator<T> > p;
    setOne(p[e][f][g]);
    return p;
}

template<class T>
inline realpoly<4, T, std::allocator<T> > RealMonomial(unsigned e, unsigned f, unsigned g, unsigned h)
// Creates a monomial in four indeterminates.
{
    realpoly<4, T, std::allocator<T> > p;
    setOne(p[e][f][g][h]);
    return p;
}

template<int n, class T, class Alloc>
class real_derivative_computer<0, n, T, Alloc>
// Computes the derivative of a given polynomial with respect to the first indeterminate. The
// coefficients of the first indeterminate are left untouched.
{
public:
    static inline void computeDerivative(const realpoly<n, T, Alloc> & src, realpoly<n, T, Alloc> & dest)
    {
        dest.d_value.resize(src.degree() >= 0 ? src.degree() : 0);
        for (int i = 1; i <= src.degree(); ++i)
        {
            dest[i - 1] = src[i];
            dest[i - 1] *= (T)(long)i;
        }
        dest.normalize();
    }
};

template<int N, int n, class T, class Alloc>
class real_derivative_computer
// Computes the derivative of a given polynomial with respect to the N-th indeterminate, N > 0. Uses
// real_derivative_computer<N-1, n-1, T, Alloc> to partially derive the coefficients of the first
// indeterminate.
{
public:
    static inline void computeDerivative(const realpoly<n, T, Alloc> & src, realpoly<n, T, Alloc> & dest)
    {
        dest.d_value.resize(src.degree() + 1);
        for (int i = 0; i <= src.degree(); ++i)
            real_derivative_computer<N - 1, n - 1, T, Alloc>::computeDerivative(src[i], dest[i]);
        dest.normalize();
    }
};

template<class T, class Alloc>
class real_derivative_computer<0, 0, T, Alloc>
// Provides an error message in case N = n.
{
public:
    class ERROR_N_must_be_strictly_less_than_n_in_derivative_template;
    
    static inline void computeDerivative(const realpoly<0, T, Alloc> &, realpoly<0, T, Alloc> &)
    {
        ERROR_N_must_be_strictly_less_than_n_in_derivative_template();
    }
};

template<int N, class T, class Alloc>
class real_derivative_computer<N, 0, T, Alloc>
// Provides an error message in case N > n.
{
public:
    class ERROR_N_must_be_strictly_less_than_n_in_derivative_template;
    
    static inline void computeDerivative(const realpoly<0, T, Alloc> &, realpoly<0, T, Alloc> &)
    {
        ERROR_N_must_be_strictly_less_than_n_in_derivative_template();
    }
};

template<int N, int n, class T, class Alloc>
inline realpoly<n, T, Alloc> derivative(const realpoly<n, T, Alloc> & p)
// Computes the partial derivative of p with respect to the N-th indeterminate. We assume that 0 <=
// N < n.
{
    realpoly<n, T, Alloc> res(p.get_allocator());
    real_derivative_computer<N, n, T, Alloc>::computeDerivative(p, res);
    return res;
}

template<int N, int n, class T>
inline realpoly<n, T, std::allocator<T> > derivative(const realpoly<n, T, std::allocator<T> > & p)
// Computes the partial derivative of p with respect to the N-th indeterminate. We assume that 0 <=
// N < n.
{
    realpoly<n, T, std::allocator<T> > res(p.get_allocator());
    real_derivative_computer<N, n, T, std::allocator<T> >::computeDerivative(p, res);
    return res;
}

template<int N, int n, class T, class Alloc>
inline void derivative(realpoly<n, T, Alloc> & res, const realpoly<n, T, Alloc> & p)
// Computes the partial derivative of p with respect to the N-th indeterminate. We assume that 0 <=
// N < n.
{
    real_derivative_computer<N, n, T, Alloc>::computeDerivative(p, res);
}

template<int N, int n, class T>
inline void derivative(realpoly<n, T, std::allocator<T> > & res, const realpoly<n, T, std::allocator<T> > & p)
// Computes the partial derivative of p with respect to the N-th indeterminate. We assume that 0 <=
// N < n.
{
    real_derivative_computer<N, n, T, std::allocator<T> >::computeDerivative(p, res);
}

template<int n, class T, class Alloc>
class real_integral_computer<0, n, T, Alloc>
// Computes the integral of a given polynomial with respect to the first indeterminate. The
// coefficients of the first indeterminate are left untouched.
{
public:
    static inline void computeIntegral(const realpoly<n, T, Alloc> & src, realpoly<n, T, Alloc> & dest)
    {
        dest.d_value.resize(src.degree() >= 0 ? src.degree() + 2 : 0);
        dest[0] = dest.d_zero.d_zero;
        for (int i = 0; i <= src.degree(); ++i)
        {
            dest[i + 1] = src[i];
            dest[i + 1] /= (T)(long)(i + 1);
        }
        dest.normalize();
    }
};

template<int N, int n, class T, class Alloc>
class real_integral_computer
// Computes the integral of a given polynomial with respect to the N-th indeterminate, N > 0. Uses
// real_integral_computer<N-1, n-1, T, Alloc> to partially derive the coefficients of the first
// indeterminate.
{
public:
    static inline void computeIntegral(const realpoly<n, T, Alloc> & src, realpoly<n, T, Alloc> & dest)
    {
        dest.d_value.resize(src.degree() + 1);
        for (int i = 0; i <= src.degree(); ++i)
            real_integral_computer<N - 1, n - 1, T, Alloc>::computeIntegral(src[i], dest[i]);
        dest.normalize();
    }
};

template<class T, class Alloc>
class real_integral_computer<0, 0, T, Alloc>
// Provides an error message in case N = n.
{
public:
    class ERROR_N_must_be_strictly_less_than_n_in_integral_template;
    
    static inline void computeIntegral(const realpoly<0, T, Alloc> &, realpoly<0, T, Alloc> &)
    {
        ERROR_N_must_be_strictly_less_than_n_in_integral_template();
    }
};

template<int N, class T, class Alloc>
class real_integral_computer<N, 0, T, Alloc>
// Provides an error message in case N > n.
{
public:
    class ERROR_N_must_be_strictly_less_than_n_in_integral_template;
    
    static inline void computeIntegral(const realpoly<0, T, Alloc> &, realpoly<0, T, Alloc> &)
    {
        ERROR_N_must_be_strictly_less_than_n_in_integral_template();
    }
};

template<int N, int n, class T, class Alloc>
inline realpoly<n, T, Alloc> integral(const realpoly<n, T, Alloc> & p)
// Computes the partial integral of p with respect to the N-th indeterminate. We assume that 0 <=
// N < n.
{
    realpoly<n, T, Alloc> res(p.get_allocator());
    real_integral_computer<N, n, T, Alloc>::computeIntegral(p, res);
    return res;
}

template<int N, int n, class T>
inline realpoly<n, T, std::allocator<T> > integral(const realpoly<n, T, std::allocator<T> > & p)
// Computes the partial integral of p with respect to the N-th indeterminate. We assume that 0 <=
// N < n.
{
    realpoly<n, T, std::allocator<T> > res(p.get_allocator());
    real_integral_computer<N, n, T, std::allocator<T> >::computeIntegral(p, res);
    return res;
}

template<int N, int n, class T, class Alloc>
inline void integral(realpoly<n, T, Alloc> & res, const realpoly<n, T, Alloc> & p)
// Computes the partial integral of p with respect to the N-th indeterminate. We assume that 0 <=
// N < n.
{
    real_integral_computer<N, n, T, Alloc>::computeIntegral(p, res);
}

template<int N, int n, class T>
inline void integral(realpoly<n, T, std::allocator<T> > & res, const realpoly<n, T, std::allocator<T> > & p)
// Computes the partial integral of p with respect to the N-th indeterminate. We assume that 0 <=
// N < n.
{
    real_integral_computer<N, n, T, std::allocator<T> >::computeIntegral(p, res);
}

template<int n, class T, class Alloc>
typename realpoly<n, T, Alloc>::ZeroContainer realpoly<n, T, Alloc>::d_zero; // Declare the zero coefficient.

template<class T, class Alloc>
class make_univariate_computer<1, T, Alloc>
{
    template<int nn, class TT, class AA>
    friend class make_univariate_computer;
    
private:
    static void helper(realpoly<1, T, Alloc> & dest, const realpoly<1, T, Alloc> & src, unsigned offset)
    {
        for (unsigned i = 0; i < src.d_value.size(); ++i)
            dest.d_value[i + offset] += src.d_value[i];
    }
    
public:
    static void make_univariate(realpoly<1, T, Alloc> & dest, const realpoly<1, T, Alloc> & src)
    {
        dest = src;
    }
};

template<int n, class T, class Alloc>
class make_univariate_computer
{
    template<int nn, class TT, class AA>
    friend class make_univariate_computer;
    
private:
    static void helper(realpoly<1, T, Alloc> & dest, const realpoly<n, T, Alloc> & src, unsigned offset)
    {
        for (unsigned i = 0; i < src.d_value.size(); ++i)
            make_univariate_computer<n - 1, T, Alloc>::helper(dest, src.d_value[i], offset + i);
    }
    
public:
    static void make_univariate(realpoly<1, T, Alloc> & dest, const realpoly<n, T, Alloc> & src)
    {
        int total = src.total_degree();
        if (total < 0)
        {
            dest.setZero();
            return;
        }
        dest.d_value.resize(total + 1);
        for (unsigned i = 0; i < dest.d_value.size(); ++i)
            setZero(dest.d_value[i]);
        helper(dest, src, 0);
        dest.normalize();
    }
};

template<int n, class T, class Alloc>
void make_univariate(realpoly<1, T, Alloc> & dest, const realpoly<n, T, Alloc> & src)
{
    make_univariate_computer<n, T, Alloc>::make_univariate(dest, src);
}

#endif

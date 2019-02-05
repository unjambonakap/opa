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

#ifndef PLLL_INCLUDE_GUARD__DOCUMENTATION_HPP
#define PLLL_INCLUDE_GUARD__DOCUMENTATION_HPP

/**
   \file
   \brief Doxygen Documentation.
   
   This header provides extensive <a href="https://en.wikipedia.org/wiki/Doxygen">doxygen</a>
   documentation for `plll`. No code or definitions are contained in this file.
*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// INTRODUCTION ////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
   \mainpage Introduction
   
   All lattices considered are integral lattices in \f$\mathbb{Z}^n\f$, given by a generating system
   \f$b_1, \ldots, b_k\f$. Note that \f$k\f$ can be in any relation to \f$n\f$ since we do not
   assume that \f$b_1, \ldots, b_k\f$ forms a basis of the lattice.
   
   \section introlatticered Introduction to lattice reduction
   
   Let \f$b_1, \ldots, b_k\f$ be \f$k\f$ linearly independent vectors in \f$\mathbb{R}^n\f$, where
   \f$n \ge k\f$. Then the _lattice of rank_ \f$k\f$ _generated_ by \f$b_1, \ldots, b_k\f$ is the
   subgroup \f[ \Lambda(b_1, \ldots, b_k) := \biggl\{ \sum_{i=1}^k \lambda_i b_i \;\biggm|\;
   \lambda_1, \dots, \lambda_k \in \mathbb{Z} \biggr\} \subseteq \mathbb{R}^n. \f] Any subgroup of
   \f$\mathbb{R}^n\f$ of this form is called a _lattice_, and the corresponding number \f$k\f$ its
   _rank_. The dimension \f$n\f$ of the ambient space is called the _dimension_ of the lattice, and
   the vectors \f$b_1, \ldots, b_k\f$ are called a _basis_ of the lattice. In case the lattice is
   contained in \f$\mathbb{Z}^n\f$, we say that the lattice is _integral_.
   
   An example for a two-dimensional lattice of rank 2 is given in the following picture:
   
   \image html lattice.png "A two-dimensional lattice"
   \image latex lattice.eps "A two-dimensional lattice" width=10cm
   
   This lattice can be given implicitly as \f[ \Lambda = \{ (x, y) \in \mathbb{Z}^2 \mid x + y
   \equiv 0 \pmod{2} \} \f] The large dot in the center represents the zero vector. This lattice can
   be generated in different ways. A simple basis is \f[ b_1 = \begin{pmatrix} 3 \\ 1 \end{pmatrix},
   \quad b_2 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}, \f] which is shown in the following picture:
   
   \image html lattice-basis1.png "A simple basis of the above lattice"
   \image latex lattice-basis1.eps "A simple basis of the above lattice" width=10cm
   
   But also other bases are possible. In fact, as soon as \f$k \ge 2\f$, there are infinitely many
   bases for every lattice of rank \f$k\f$. The following picture shows a more complicated basis of
   the same lattice as before:
   
   \image html lattice-basis2.png "A more complicated basis of the above lattice"
   \image latex lattice-basis2.eps "A more complicated basis of the above lattice" width=10cm
   
   The aim of _lattice reduction_ is to start with a (complicated) basis and try to find a better,
   "nicer" basis. One goal is to find as short as possible basis vectors. The vectors in the first
   basis \f$b_1, b_2\f$ are a relatively good choice, but not optimal: a better basis would be \f[
   c_1 = \begin{pmatrix} 1 \\ -1 \end{pmatrix}, \quad c_2 = \begin{pmatrix} 1 \\ 1 \end{pmatrix},
   \f] as demonstrated by the following picture:
   
   \image html lattice-basis3.png "An 'optimal' basis for the above lattice"
   \image latex lattice-basis3.eps "An 'optimal' basis for the above lattice" width=10cm
   
   Note that both in this basis as in the first one, a shortest non-zero vector is part of the
   basis, namely \f$b_2 = c_2 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}\f$ and also \f$c_1\f$ in the
   last basis: their Euclidean length is \f$\sqrt{1^2 + 1^2} = \sqrt{2} \approx 1.4142\f$.
   
   Since finding such shortest vectors is hard (for \f$k \to \infty\f$), as M. Ajtai showed in 1998
   \cite ajtai-SVPisNPhard, most lattice basis reduction algorithms try to find good approximations
   of a shortest vector, where the approximation quality usually increases with running time, but is
   often far away from optimal. The first polynomial time lattice basis reduction algorithm, 1982's
   celebrated LLL algorithm \cite lenstra-lenstra-lovasz-LLL, has a exponential approximation bound
   (which Schnorr showed is optimal). For some problems, such bounds already suffice, while for
   others, better approximations are needed. Since 1982, many new algorithms have been designed to
   find better approximations. The `plll` library tries to make many of them available.
   
   For a more detailed description of lattice reduction algorithms, see \ref latticereduction.
   
   \section usage-example Usage
   A simple example of how to use the `plll` library is the following:
   \include example-01.cpp
   It applies LLL with standard parameters (that is, reduction parameter 0.99) to a basis read from
   standard input, and writes the resulting basis to standard output.
   
   A more thorough discussion of this example can be found in \ref example-01. More examples can be 
   found in \ref examples. For a detailed explanation of the interface of the `plll::LatticeReduction` 
   class, see its description. For a detailed description of how to work with matrices, including how 
   to input and output them, see \ref matrixvector. To know more about how to work with the arithmetic 
   functions provided by `plll`, see \ref arithctxts. And for a detailed description of the supported 
   lattice reduction algorithms, see \ref desc-algs.
   
   \section attrib Attributions
   
   All code except the memory allocator was written by Felix Fontein.
   
   All algorithms implemented in this library are taken from published papes. Some were slightly 
   modified. See \ref desc-algs for descriptions of the algorithms, attributions to their inventors
   and references to the relevant publications.
   
   The memory allocator in the files plll/src/utilities/dlalloc.c and plll/src/utilities/dlalloc.h
   is Doug Lea's dlmalloc, which was released to the public domain. For more information, check out 
   plll/src/utilities/dlalloc.c or 
   <a href="ftp://gee.cs.oswego.edu/pub/misc/malloc.c">ftp://gee.cs.oswego.edu/pub/misc/malloc.c</a>.
   
   \section license License
   
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
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// ARITHMETIC CONTEXTS /////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
   \page arithctxts Arithmetic Contexts
   
   \section arithctxts_intro Introduction
   
   With aritrary precision floating point types, the main problem is how to store the default
   precision for new variables. Low-level libraries such as MPFR, but also high level libraries such
   as NTL, use a global (static) variable to store the precision, while allowing to change the
   precision of variables during initalization (MPFR) or later (MPFR and NTL).
   
   In the context of multithreading, using one global variable for all threads is not
   acceptable. Thus, one has to carry around the required precision per thread (or logical part of
   the program), and initialize floating point variables using this precision. On the other hand,
   for fixed precision floating point types such as the native float, double and long double, such a
   precision is not needed.
   
   The idea in the `plll` library is to carry around contextes, which are lightweight objects
   containing no (for native types) or very little (for arbitrary precision floating point types)
   information which are used to initialize new variables in a uniform way: if `c` is an instance of a
   context type `Context`, then a variable `v` of type `Context::Type` is initialized via
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
       Context::Type v(c);
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   with the correct precision given by `c` in the case it carries a precision, or `c` is ignored in
   case of native types or arbitrary precision integer types. For other types later added, a context
   could in theory also hold other information; code written using contexts should be agnostic about
   this extra information.
   
   
   \section arithctxts_contexts Contexts
   
   The basic information captured by every context is a typedef `Type`, which yields the type
   described by this context, together with an `enum` providing certain properties:
   
   - `is_cputype`: `true` for native types supported by CPU (or emulated natively on CPUs, such as
     `long double` on most UltraSPARC CPUs).
     
   - `is_realtype`: `true` for types representing (most or all) rational numbers, such as floating
     point types.
     
   - `is_inttype`: `true` for types representing (a range of or all) integer types. Only signed
     integer types are allowed.
     
   - `is_exact`: `true` for types which do use approximations (like floating point types). For
     example, for GMP and CPU integers, this is always `true`.
     
   - `has_infinity`: `true` for types having objects representing \f$+\infty\f$ and
     \f$-\infty\f$.
     
   - `has_uniform_rng`: `true` if the type provides a uniform random number generator. Such a
     context exhibits a type `UniformRNG`, whose interface will be described in more detail below in
     \ref arithctxts_urng. Note that the rational number context from `rational.hpp` is an example
     of a context not providing a uniform random number generator.
   
   Integer types also have a typedef `Integer`, being the same as `Type`. Besides the properties
   listed above, they also provide:
   
   - `is_modulo`: `true` for types which are "wrap around", such as native integer arithmetic which
     in fact computes in rings such as \f$\mathbb{Z}/2^n\mathbb{Z}\f$.
   
   Real types also have a typedef `Real`, being the same as `Type`. Besides the properties listed
   above, they also provide:
   
   - `is_variable_precision`: `true` if the precision of the current type can be modified. (See
     below for the interface.)
     
   - `has_squareroot`: `true` if (approximated) square roots are supported. This will be `false` for
     rational numbers, but `true` for most (or even all) floating point types. For more information,
     see \ref arithctxts_sqrtfp.
     
   - `has_full_power`: `true` if elements can be raised to fractional powers (where the powers are
     coming from the same context). This will be `false` for rational numbers, but `true` for most
     (or even all) floating point types. For more information, see \ref arithctxts_sqrtfp.
     
   - `has_special_fns`: `true` if special functions such as `lgamma` etc. are supported for this
     type. A list of special functions with explanations can be found below in \ref
     arithctxts_specfns.
     
   - `has_huge_exponent`: `true` for floating point types using at least the range of a `signed int`
     for the exponent. If the exponent range is (much) smaller, this has to be set to `false`.
     
   - `has_constants`: `true` if special constants such as \f$\varepsilon\f$, \f$\pi\f$, \f$e\f$ and
     \f$\log 2\f$ are provided by the context. See below in \ref arithctxts_consts for more
     information.
     
   - `has_trigonometric`: `true` if trigonometric functions such as `sin` and `cos` are
     provided. See below in \ref arithctxts_trig for more information.
   
   Integer contexts do not provide any functions. Real contexts, on the other hand, provide a few
   standard functions:
   
   - `void setRealPrecision(unsigned long)`: sets the precision of the context to (at least) the
     given value, if possible. If this is not possible, the context is free should chose a precision
     as high as possible. In case the precision cannot be changed (i.e. `is_variable_precision` is
     `false`), this call can be ignored.
     
   - `unsigned long getRealPrecision() const`: retrieves the current precision of the context. For
     types not having a precision, the return value should be a large number such as
     `std::numeric_limits<unsigned long>::max()`.
     
   - `unsigned long getMinRealPrecision() const`: retrieves the minimal possible precision. For
     types not having a precision, the return value should be a large number such as
     `std::numeric_limits<unsigned long>::max()`.
     
   - `unsigned long getMaxRealPrecision() const`: retrieves the maximal possible precision. For
     types not having a precision, the return value should be a large number such as
     `std::numeric_limits<unsigned long>::max()`. In any case, the value returned must always be
     larger than the one returned by `getMinRealPrecision()`.
   
   
   \section arithctxts_ops Operands
   
   Context types are expected to overload any reasonable operator, such as the arithmetic operators
   `+`, `-`, `*`, `/`, `%`, `<<`, `>>`, `++`, `--`, the arithmetic assignment operators `+=`, `-=`,
   `*=`, `/=`, `%=`, `<<=`, `>>=`, and the comparison operators `==`, `!=`, `<`, `>`, `<=`, `>=`.
   
   Additionally, excepted some operator-like functions:
   
   - `Context::Type abs(const Context::Type &)`: takes the absolute value of the operand and returns
     the result.
   
   - `Context::Type square(const Context::Type &)`: squares the operand and returns the result.
   
   For real types, there are variants which also accept a context:
   
   - `Context::Type abs(const Context::Type &, const Context &)`: takes the absolute value of the
     first operand and returns the result. Uses precision etc. from the context (second argument).
   
   - `Context::Type square(const Context::Type &, const Context &)`: squares the first operand and
     returns the result. Uses precision etc. from the context (second argument).
   
   
   \section arithctxts_top Three-operand forms
   
   Besides the above operators, context types are also expected to have _three-operand forms_ for
   most of these operators (for binary operations, and two-operand forms for unary operators). For
   all types, the following three-operand forms are expected:
   
   - `void add(Context::Type &, const Context::Type &, const Context::Type &)`: adds the second and
     third operand and stores the result in the first.
   
   - `void addmul(Context::Type &, const Context::Type &, const Context::Type &)` multiplies the
     second and third operand and adds the result to the first operand.
     
   - `void sub(Context::Type &, const Context::Type &, const Context::Type &)`: subtracts the third
     from the second operand and stores the result in the first.
     
   - `void submul(Context::Type &, const Context::Type &, const Context::Type &)` multiplies the
     second and third operand and subtracts the result from the first operand.
     
   - `void mul(Context::Type &, const Context::Type &, const Context::Type &)`: multiplies the
     second with the third operand and stores the result in the first.
   
   - `void div(Context::Type &, const Context::Type &, const Context::Type &)`: divides the second
     by the third operand and stores the quotient in the first.
     
   - `void mod(Context::Type &, const Context::Type &, const Context::Type &)`: divides the second
     by the third operand and stores the remainder in the first.
     
   - `void shl(Context::Type &, const Context::Type &, const Context::Type &)`: multiplies the
     second operand by 2 to the power of the third operand and stores the result in the first.
     
   - `void shl(Context::Type &, const Context::Type &, long)`: multiplies the second operand by 2 to
     the power of the third operand and stores the result in the first.
     
   - `void shl(Context::Type &, const Context::Type &, const arithmetic::Integer &)`: multiplies the
     second operand by 2 to the power of the third operand and stores the result in the first.
     
   - `void shr(Context::Type &, const Context::Type &, const Context::Type &)`: multiplies the
     second operand by 2 to the power of the third operand and stores the result in the first.
     
   - `void shr(Context::Type &, const Context::Type &, long)`: multiplies the second operand by 2 to
     the power of the third operand and stores the result in the first.
     
   - Context::Type &, const arithmetic::Integer &)`: multiplies the second operand by 2 to the power
     of the third operand and stores the result in the first.
   
   The following two-operand forms exist for uniary operators:
   
   - `void neg(Context::Type &, const Context::Type &)`: negates the second operand and stores the
     result in the first.
   
   - `void abs(Context::Type &, const Context::Type &)`: takes the absolute value of the second
     operand and stores the result in the first.
   
   - `void square(Context::Type &, const Context::Type &)`: squares the second operand and stores
     the result in the first.
   
   - `void increment(Context::Type &, const Context::Type &)`: increments the second operand by one
     and stores the result in the first.
   
   - `void decrement(Context::Type &, const Context::Type &)`: decrements the second operand by one
     and stores the result in the first.
   
   The following is an abbrevation:
   
   - `void makeAbs(Context::Type &)`: makes the operand non-negative by flipping its sign if
     necessary.
   
   For comparisons, there should exist the following functions:
   
   - `int compare(const Context::Type &, const Context::Type &)`: compares the first to the second
     operand. Returns a negative number if the first is smaller, a positive number if the first is
     larger, and zero if they are equal.
   
   - `int compareAbsValues(const Context::Type &, const Context::Type &)`: compares the absolute
     value of the first to the absolute value of the second operand. Returns a negative number if
     the first absolute value is smaller, a positive number if the first absolute value is larger,
     and zero if the absolute values are equal.
   
   For integer types, the following bitwise three-operand and two-operand forms exist:
   
   - `void band(Context::Integer &, const Context::Integer &, const Context::Integer &)`: returns
     the bitwise and of the second and third operand in the first operand.
   
   - `void bor(Context::Integer &, const Context::Integer &, const Context::Integer &)`: returns
     the bitwise or of the second and third operand in the first operand.
     
   - `void bxor(Context::Integer &, const Context::Integer &, const Context::Integer &)`: returns
     the bitwise exclusive or of the second and third operand in the first operand.
     
   - `void bneg(Context::Integer &, const Context::Integer &)`: returns the bitwise negation of the
     second operand in the first operand.
   
   
   \section arithctxts_preds Predicates
   
   For all types, the following predicates exist:
   
   - `bool isZero(const Context::Type &)`: returns `true` if and only if the argument is zero.
   
   - `bool isOne(const Context::Type &)`: returns `true` if and only if the argument is one.
   
   - `bool isPositive(const Context::Type &)`: returns `true` if and only if the argument is
     positive.
   
   - `bool isNegative(const Context::Type &)`: returns `true` if and only if the argument is
     negative.
   
   - `bool isNonPositive(const Context::Type &)`: returns `true` if and only if the argument is
     negative or zero.
   
   - `bool isNonNegative(const Context::Type &)`: returns `true` if and only if the argument is
     positive or zero.
   
   - `int sign(const Context::Type &)`: returns a negative number if the argument is negative, a
     positive number if it is positive, and zero if it is zero.
   
   For integer types, the following additional predicates exist:
   
   - `bool isPMOne(const Context::Integer &)`: returns `true` if and only if the argument is one or
     minus one.
   
   - `bool isPMTwo(const Context::Integer &)`: returns `true` if and only if the argument is one or
     minus two.
   
   
   \section arithctxts_ints Integer Functions
   
   The following functions exists only for integer types:
   
   - Bit access
     
     - `int bit(const Context::Integer &, long n)`: returns the value of the `n`-th bit of the
       absolute value first argument.
     
     - `void setbit(const Context::Integer &, long n, bool value = true)`: sets the value of the
       `n`-th bit of the first argument to `value`.
     
     - `long bitLength(const Context::Integer & x)`: returns `n` such that \f$2^{n-1} \le |x| <
       2^n\f$.

   - Approximated Base-2 Logarithm
   
     - `long approxLog2(const Context::Integer & x)`: returns \f$\approx \log_2 |x|\f$.
     
     - `long ceilOfLog2(const Context::Integer & x)`: returns \f$\lceil \log_2 |x| \rceil\f$.
     
     - `long floorOfLog2(const Context::Integer & x)`: returns \f$\lfloor \log_2 |x| \rfloor\f$.
   
   - Rounded division
   
     - `void floorDiv(Context::Integer &, const Context::Integer & a, const Context::Integer & b)`:
       Computes \f$\lfloor \tfrac{a}{b} \rfloor\f$ and stores the result in the first operand.
     
     - `Context::Integer floorDiv(const Context::Integer & a, const Context::Integer & b)`: Computes
       \f$\lfloor \tfrac{a}{b} \rfloor\f$ and returns the result.
       
     - `void ceilDiv(Context::Integer &, const Context::Integer & a, const Context::Integer & b)`:
       Computes \f$\lceil \tfrac{a}{b} \rceil\f$ and stores the result in the first operand.
       
     - `Context::Integer ceilDiv(const Context::Integer & a, const Context::Integer & b)`: Computes
       \f$\lceil \tfrac{a}{b} \rceil\f$ and returns the result.
     
     - `void roundDiv(Context::Integer &, const Context::Integer & a, const Context::Integer & b)`:
       Computes \f$\lfloor \tfrac{a}{b} \rceil\f$ (rounding to the next integer) and stores the
       result in the first operand.
     
     - `Context::Integer roundDiv(const Context::Integer & a, const Context::Integer & b)`: Computes
        \f$\lfloor \tfrac{a}{b} \rceil\f$ (rounding to the next integer) and returns the result.
       
   - Euclidean division
   
     - `void euclideanDivision(Context::Integer & q, Context::Integer & r, const Context::Integer &
       a, const Context::Integer & b)`: Computes an Euclidean Division of `a` by `b`, i.e. computes
       `q` and `r` such that `a == q * b + r` and that \f$0 \le |r| \le |b|\f$ and \f$b \cdot r \ge
       0\f$.
     
     - `void euclideanDivisionPos(Context::Integer & q, Context::Integer & r, const Context::Integer
       & a, const Context::Integer & b)`: Computes an Euclidean Division of `a` by `b`,
       i.e. computes `q` and `r` such that `a == q * b + r` and that \f$0 \le r \le |b|\f$.

   - Greatest Common Divisors and Least Common Multiples
     
     - `void GCD(Context::Integer &, const Context::Integer & x, const Context::Integer & y)`:
       Computes the non-negative Greatest Common Divisior of `x` and `y` and stores it in the first
       argument.
     
     - `Context::Integer GCD(const Context::Integer & x, const Context::Integer & y)`: Computes the
       non-negative Greatest Common Divisior of `x` and `y` and returns it.
     
     - `void XGCD(Context::Integer & r, Context::Integer & a, Context::Integer & b, const
       Context::Integer & x, const Context::Integer & y)`: Computes the non-negative extended
       Greatest Common Divisior `r` of `x` and `y`, i.e. computes also `a` and `b` such that `r == x
       * a + b * y`.
     
     - `void LCM(Context::Integer &, const Context::Integer & x, const Context::Integer & y)`:
       Computes the non-negative Least Common Multiple `r` of `x` and `y` and stores it in the first
       operand.
     
     - `Context::Integer LCM(const Context::Integer & x, const Context::Integer & y)`: Computes the
       non-negative Least Common Multiple of `x` and `y` and returns it.
   
   
   \section arithctxts_urng Uniform Random Number Generator
   
   If the `has_uniform_rng` property of a context type `Context` is `true`, then
   `Context::UniformRNG` describes a lightweight object which takes as input a
   `plll::arithmetic::RandomNumberGenerator` object (see below), and which provides basic uniform
   random number generation. The following interface is always provided:
   
   - `UniformRNG(RandomNumberGenerator &)`: creates a new `Context::UniformRNG` object based on the
     given random number generator instance.
  
   For integer types, `Context::UniformRNG` provides the following interface:
   
   - `void random(Context::Integer & result, const Context::Integer & bound)`: creates a random
     number \f$result\f$ satisfying \f$0 \le result < bound\f$.
     
   - `Context::Integer random(const Context::Integer & bound, const Context &)`: creates and returns
     a random number \f$x\f$ satisfying \f$0 \le x < bound\f$.
     
   - `void randomBits(Context::Integer & result, unsigned long bits)`: creates `bits` random bits
     into `result`. That is, it generates a random integer \f$result\f$ satisfying \f$0 \le result <
     2^{bits}\f$.
     
   - `Context::Integer randomBits(unsigned long bits, const Context &)`: creates and returns `bits`
     random bits. That is, it generates a random integer \f$x\f$ satisfying \f$0 \le x <
     2^{bits}\f$.
     
   - `void randomLen(Context::Integer & result, unsigned long bits)`: creates a random integer in
     `result` of precisely `bits` bits. That is, it generates a random integer \f$result\f$
     satisfying \f$2^{bits-1} \le result < 2^{bits}\f$.
     
   - `Context::Integer randomLen(unsigned long bits, const Context &)`: creates a random integer of
     precisely `bits` bits. That is, it generates a random integer \f$x\f$ satisfying \f$2^{bits-1}
     \le x < 2^{bits}\f$.
   
   For real types, on the other hand, the following interface is provided:
   
   - `void randomUniform(Context::Real &)`: creates a uniform random number in \f$[0, 1)\f$.
   
   - `Context::Real randomUniform(const Context &)`: creates a uniform random number in \f$[0, 1)\f$
     using the given context.
   
   
   \section arithctxts_sqrtfp Square Roots and Full Powers
   
   Certain real types provide support for square roots and arbitrary powers. In case of square
   roots, this is indicated by the `has_squareroot` property of the context. If it is `true`, the
   following functions will be provided:
   
   - `void sqrt(Context::Real &, const Context::Real &)`: computes the square root of the second
     argument and stores it in the first argument.
   
   - `Context::Real sqrt(const Context::Real &)`: computes and returns the square root of the
     argument. Note that the return type can be different when expression templates are employed;
     the return type must be implicitly castable to the type `Context::Real`.
   
   - `Context::Real sqrt(const Context::Real &, const Context &)`: computes and returns the square
     root of the argument with the precision given by the context. Note that the return type can be
     different when expression templates are employed; the return type must be implicitly castable
     to the type `Context::Real`.
   
   In any case, the following functions are provided to compute powers -- both for real types and
   integer types:
   
   - `void square(Context::Type &, const Context::Type &)`: squares the second argument, and returns
     the result in the first argument.
   
   - `Context::Type square(const Context::Type &)`: squares the first argument, and returns the
     result.
   
   - `void power(Context::Type &, const Context::Type &, long)`: raises the second argument to the
     third argument's power, and returns the result in the first argument.
   
   - `Context::Type power(const Context::Type &, long)`: raises the first argument to the second
     argument's power, and returns the result.
   
   - `void power(Context::Type &, const Context::Type &, const plll::arithmetic::Integer &)`: raises
     the second argument to the third argument's power, and returns the result in the first
     argument.
   
   - `Context::Type power(const Context::Type &, const plll::arithmetic::Integer &)`: raises the
     first argument to the second argument's power, and returns the result.
   
   - `Context::Type operator << (const Context::Type &, long)`: multiplies the first argument by 2
     to the power of the second argument, and returns the result.
   
   - `Context::Type operator >> (const Context::Type &, long)`: divides the first argument by 2 to
     the power of the second argument, and returns the result. In the case of integers, rounds
     towards zero.
   
   - `Context::Type & operator <<= (Context::Type &, long)`: multiplies the first argument by 2 to
     the power of the second argument, and stores the result in the first argument.
   
   - `Context::Type & operator >>= (Context::Type &, long)`: divides the first argument by 2 to the
     power of the second argument, and stores the result in the first argument. In the case of
     integers, rounds towards zero.
   
   - `void shl(Context::Type &, const Context::Type &, long)`: multiplies the second argument by 2
     to the power of the third argument and stores the result in the first argument.
   
   - `void shr(Context::Type &, const Context::Type &, long)`: divides the second argument by 2 to
     the power of the third argument and stores the result (rounded towards zero) in the first
     argument.
   
   For real types, the functions returning a `Context::Real` also accept a `const Context &` as the
   third parameter whose precision is used for the result.
   
   For integer types, the following functions are always provided:
   
   - `void power(Context::Integer &, const Context::Integer &, const Context::Integer &)`: raises
     the second argument to the third argument's power, and returns the result in the first
     argument.
   
   - `Context::Integer power(const Context::Integer &, const Context::Integer &)`: raises the first
     argument to the second argument's power, and returns the result.`
   
   For real types, the following functions are _only_ provided if `Context::has_full_power` is set
   to `true`:
   
   - `void power(Context::Real &, const Context::Real &, const Context::Real &)`: raises the second
     argument to the third argument's power, and returns the result in the first argument.

   - `Context::Real power(const Context::Real &, const Context::Real &)`: raises the first argument
     to the second argument's power, and returns the result.`

   - `Context::Real power(const Context::Real &, const Context::Real &, const Context &)`: raises
     the first argument to the second argument's power with the precision given by the context, and
     returns the result.`
   
   - `Context::Real operator << (const Context::Real &, const Context::Real)`: multiplies the first
     argument by 2 to the power of the second argument, and returns the result.
   
   - `Context::Real operator >> (const Context::Real &, const Context::Real)`: divides the first
     argument by 2 to the power of the second argument, and returns the result.
     
   - `Context::Real & operator <<= (Context::Real &, const Context::Real)`: multiplies the first
     argument by 2 to the power of the second argument, and stores the result in the first argument.
   
   - `Context::Real & operator >>= (Context::Real &, const Context::Real)`: divides the first
     argument by 2 to the power of the second argument, and stores the result in the first
     argument.
   
   - `void shl(Context::Real &, const Context::Real &, unsigned const Context::Real)`: multiplies
     the second argument by 2 to the power of the third argument and stores the result in the first
     argument.
   
   - `void shr(Context::Real &, const Context::Real &, unsigned const Context::Real)`: divides the
     second argument by 2 to the power of the third argument and stores the result in the first
     argument.
   
   
   \section arithctxts_trig Trigonometric Functions
   
   If a real type's context `Context` has `Context::has_trigonometric == true`, the following
   trigonometric functions are provided for that type:
   
   - Sine:
     - `void sin(Context::Real &, const Context::Real &)`: computes the sine \f$\sin x\f$ of the
       second argument \f$x\f$ and stores the result in the first argument.
     - `Context::Real sin(const Context::Real &)`: computes the sine \f$\sin x\f$ of the first
       argument \f$x\f$ and returns the result.
     - `Context::Real sin(const Context::Real &, const Context &)`: computes the sine \f$\sin x\f$
       of the first argument \f$x\f$ with the precision given by the context and returns the result.

   - Cosine:
     - `void cos(Context::Real &, const Context::Real &)`: computes the cosine \f$\cos x\f$ of the
       second argument \f$x\f$ and stores the result in the first argument.
     - `Context::Real cos(const Context::Real &)`: computes the cosine \f$\cos x\f$ of the first
       argument \f$x\f$ and returns the result.
     - `Context::Real cos(const Context::Real &, const Context &)`: computes the cosine \f$\cos x\f$
       of the first argument \f$x\f$ with the precision given by the context and returns the result.
   
   - Tangent:
     - `void tan(Context::Real &, const Context::Real &)`: computes the tangent \f$\tan x\f$ of the
       second argument \f$x\f$ and stores the result in the first argument.
     - `Context::Real tan(const Context::Real &)`: computes the tangent \f$\tan x\f$ of the first
       argument \f$x\f$ and returns the result.
     - `Context::Real tan(const Context::Real &, const Context &)`: computes the tangent \f$\tan
       x\f$ of the first argument \f$x\f$ with the precision given by the context and returns the
       result.
   
   In case these functions are provided, also the following inverse trigonometric functions are
   provided for that type:
   
   - Arcsine:
     - `void asin(Context::Real &, const Context::Real &)`: computes the arcsine \f$\arcsin x \in
       [-\tfrac{\pi}{2}, \tfrac{\pi}{2}]\f$ of the second argument \f$x\f$ and stores the result in
       the first argument.
     - `Context::Real asin(const Context::Real &)`: computes the arcsine \f$\arcsin x \in
       [-\tfrac{\pi}{2}, \tfrac{\pi}{2}]\f$ of the first argument \f$x\f$ and returns the result.
     - `Context::Real asin(const Context::Real &, const Context &)`: computes the arcsine \f$\arcsin x
       \in [-\tfrac{\pi}{2}, \tfrac{\pi}{2}]\f$ of the first argument \f$x\f$ with the precision given
       by the context and returns the result.
   
   - Arccosine:
     - `void acos(Context::Real &, const Context::Real &)`: computes the arccosine \f$\arccos x \in
       [0, \pi]\f$ of the second argument \f$x\f$ and stores the result in the first argument.
     - `Context::Real acos(const Context::Real &)`: computes the arccosine \f$\arccos x \in [0,
       \pi]\f$ of the first argument \f$x\f$ and returns the result.
     - `Context::Real acos(const Context::Real &, const Context &)`: computes the arccosine \f$\arccos
       x \in [0, \pi]\f$ of the first argument \f$x\f$ with the precision given by the context and
       returns the result.
   
   - Arctangent:
     - `void atan(Context::Real &, const Context::Real &)`: computes the arctangent \f$\arctan x \in
       [-\tfrac{\pi}{2}, \tfrac{\pi}{2}]\f$ of the second argument \f$x\f$ and stores the result in
       the first argument.
     - `Context::Real atan(const Context::Real &)`: computes the arctangent \f$\arctan x \in
       [-\tfrac{\pi}{2}, \tfrac{\pi}{2}]\f$ of the first argument \f$x\f$ and returns the result.
     - `Context::Real atan(const Context::Real &, const Context &)`: computes the arctangent
       \f$\arctan x \in [-\tfrac{\pi}{2}, \tfrac{\pi}{2}]\f$ of the first argument \f$x\f$ with the
       precision given by the context and returns the result.
     - `void atan2(Context::Real &, const Context::Real & y, const Context::Real & x)`: computes the
       arctangent \f$\arctan \frac{y}{x} \in [-\pi, \pi]\f$ of the second argument \f$y\f$ and the
       third argument \f$x\f$ and stores the result in the first argument. The signs of \f$x\f$ and
       \f$y\f$ are used to determine the quadrant of the result.
     - `Context::Real atan2(const Context::Real & y, const Context::Real & x)`: computes the
       arctangent \f$\arctan \frac{y}{x} \in [-\pi, \pi]\f$ of the first argument \f$y\f$ and the
       second argument \f$x\f$ and returns the result. The signs of \f$x\f$ and \f$y\f$ are used to
       determine the quadrant of the result.
     - `Context::Real atan2(const Context::Real & y, const Context::Real & x, const Context &)`:
       computes the arctangent \f$\arctan \frac{y}{x} \in [-\pi, \pi]\f$ of the first argument
       \f$y\f$ and the second argument \f$x\f$ with the precision given by the context and returns
       the result. The signs of \f$x\f$ and \f$y\f$ are used to determine the quadrant of the
       result.
   
   
   \section arithctxts_specfns Special Functions
   
   If a real type's context `Context` has `Context::has_special_fns == true`, the following special
   functions are provided for that type:
   
   - Exponential function:
     - `void exp(Context::Real &, const Context::Real &)`: computes the exponential function \f$\exp
       x\f$ of the second argument \f$x\f$ and stores the result in the first argument. Note that
       \f$\exp x = \sum_{n=0}^\infty \frac{x^n}{n!}\f$.
     - `Context::Real exp(const Context::Real &)`: computes the exponential function \f$\exp x\f$ of
       the first argument \f$x\f$ and returns the result. Note that \f$\exp x = \sum_{n=0}^\infty
       \frac{x^n}{n!}\f$.
     - `Context::Real exp(const Context::Real &, const Context &)`: computes the exponential
       function \f$\exp x\f$ of the first argument \f$x\f$ with the precision given by the context
       and returns the result. Note that \f$\exp x = \sum_{n=0}^\infty \frac{x^n}{n!}\f$.
   
   - Natural logarithm:
     - `void log(Context::Real &, const Context::Real &)`: computes the natural logarithm \f$\log
       x\f$ of the second argument \f$x\f$ and stores the result in the first argument. Note that
       \f$\exp \log x = x\f$ for all \f$x > 0\f$.
     - `Context::Real log(const Context::Real &)`: computes the natural logarithm \f$\log x\f$ of
       the first argument \f$x\f$ and returns the result. Note that \f$\exp \log x = x\f$ for all
       \f$x > 0\f$.
     - `Context::Real log(const Context::Real &, const Context &)`: computes the natural logarithm
       \f$\log x\f$ of the first argument \f$x\f$ with the precision given by the context and
       returns the result. Note that \f$\exp \log x = x\f$ for all \f$x > 0\f$.
   
   - Gamma function:
     - `void gamma(Context::Real &, const Context::Real &)`: computes the Gamma function
       \f$\Gamma(x)\f$ of the second argument \f$x\f$ and stores the result in the first
       argument. Note that \f$\Gamma(x) = \int_0^\infty t^{x-1} e^{-t} \; dt\f$ for \f$x > 0\f$, and
       that \f$\Gamma(x)\f$ has poles at 0 and the negative integers.
     - `Context::Real gamma(const Context::Real &)`: computes the Gamma function \f$\Gamma(x)\f$ of
       the first argument \f$x\f$ and returns the result. Note that \f$\Gamma(x) = \int_0^\infty
       t^{x-1} e^{-t} \; dt\f$ for \f$x > 0\f$, and that \f$\Gamma(x)\f$ has poles at 0 and the
       negative integers.
     - `Context::Real gamma(const Context::Real &, const Context &)`: computes the Gamma function
       \f$\Gamma(x)\f$ of the first argument \f$x\f$ with the precision given by the context and
       returns the result. Note that \f$\Gamma(x) = \int_0^\infty t^{x-1} e^{-t} \; dt\f$ for \f$x >
       0\f$, and that \f$\Gamma(x)\f$ has poles at 0 and the negative integers.

   - Natural logarithm of the absolute value of the Gamma function:
     - `void lgamma(Context::Real &, const Context::Real &)`: computes the logarithm \f$\log
       |\Gamma(x)|\f$ of the absolute value of the Gamma function of the second argument \f$x\f$ and
       stores the result in the first argument.
     - `Context::Real lgamma(const Context::Real &)`: computes the logarithm \f$\log |\Gamma(x)|\f$
       of the absolute value of the Gamma function of the first argument \f$x\f$ and returns the
       result.
     - `Context::Real lgamma(const Context::Real &, const Context &)`: computes the logarithm
       \f$\log |\Gamma(x)|\f$ of the absolute value of the Gamma function of the first argument
       \f$x\f$ with the precision given by the context and returns the result.
     - `void lgamma(Context::Real &, int &, const Context::Real &)`: computes the logarithm \f$\log
       |\Gamma(x)|\f$ of the absolute value of the Gamma function of the third argument \f$x\f$ and
       stores the result in the first argument. The sign of \f$\Gamma(x)\f$ is stored in the second
       variable by 1 for \f$\Gamma(x) \ge 0\f$ and -1 for \f$\Gamma(x) < 0\f$.
     - `Context::Real lgamma(int &, const Context::Real &)`: computes the logarithm \f$\log
       |\Gamma(x)|\f$ of the absolute value of the Gamma function of the second argument \f$x\f$ and
       returns the result. The sign of \f$\Gamma(x)\f$ is stored in the second variable by 1 for
       \f$\Gamma(x) \ge 0\f$ and -1 for \f$\Gamma(x) < 0\f$.
     - `Context::Real lgamma(int &, const Context::Real &, const Context &)`: computes the logarithm
       \f$\log |\Gamma(x)|\f$ of the absolute value of the Gamma function of the second argument
       \f$x\f$ with the precision given by the context and returns the result. The sign of
       \f$\Gamma(x)\f$ is stored in the second variable by 1 for \f$\Gamma(x) \ge 0\f$ and -1 for
       \f$\Gamma(x) < 0\f$.
   
   
   \section arithctxts_consts Constants
   
   If a real type's context `Context` has `Context::has_constants == true`, the following constants
   are provided for that type:

   - \f$\varepsilon\f$: the smallest positive number such that \f$1 + \varepsilon \neq 1\f$ (this
     only makes sense for approximate floating point types). This is also called the _machine
     precision_ constant.

   - \f$\pi\f$: the circle number, such that a circle of radius \f$r\f$ has area \f$\pi r^2\f$ and
     circumference \f$2 \pi r\f$.

   - \f$e\f$: the Euler number \f$e = \sum_{n=0}^\infty \frac{1}{n!} = \exp 1\f$.
   
   - \f$\log 2\f$: the natural logarithm of 2.
   
   These variables can be accessed by the following member functions of a real `Context`:
   
   - `Context::Real getEpsilon() const`: returns the machine precision constant \f$\varepsilon\f$.
   
   - `void getEpsilon(Context::Real &) const`: stores the machine precision constant
     \f$\varepsilon\f$ in the argument.
   
   - `Context::Real getPi() const`: returns \f$\pi\f$.
   
   - `void getPi(Context::Real &) const`: stores \f$\pi\f$ in the argument.
   
   - `Context::Real getEuler() const`: returns \f$e = \exp 1\f$.
   
   - `void getEuler(Context::Real &) const`: stores \f$e = \exp 1\f$ in the argument.
   
   - `Context::Real getLog2() const`: returns \f$\log 2\f$.
   
   - `void getLog2(Context::Real &) const`: stores \f$\log 2\f$ in the argument.
   
   
   \section arithctxts_traits Traits
   
   If `T` is a type, then `plll::arithmetic::traits::type_traits<T>` retrieves some basic
   information and arithmetic properties of `T`. The traits are mainly used by the conversion
   functions (see \ref arithctxts_convs), but can also be used independently of these.
   
   Traits are only defined for native types such as `signed int`, `unsigned long` and `double`, for
   string types such as `const char *` and `std::string`, and for arithmetic types given by contexts
   as above.
   
   The following properties are defined in all defined `type_traits<>` instances:
   
   - `is_number`: `true` for arithmetic (number) types.
   
   - `is_realtype`: `true` for arithmetic types which are real types (such as floating point types
     or rational numbers). Always `false` if `is_number` is `false`.
   
   - `is_inttype`: `true` for arithmetic types which are integer types. Always `false` if
     `is_number` is `false`.
   
   - `is_string`: `true` for string types. Always `false` if `is_number` is `true`.
      
      \todo Implement support for other character types than `char` (like `wchar_t`, and also new C++11 strings)
   
   - `is_cpp_string`: `true` for `std::string`. Always `false` if `is_string` is `false`.
   
   - `is_c_string`: `true` for `char *` and `const char *`. Always `false` if `is_string` is
     `false`.
   
   For number types, i.e. if `is_number` is `true`, the following properties are defined:
   
   - `is_cputype`: `true` for native types supported by CPU (or emulated natively on CPUs, such as
     `long double` on most UltraSPARC CPUs).
   
   - `is_exact`: `true` for types which do use approximations (like floating point types). For
     example, for GMP and CPU integers, this is always `true`.
   
   - `has_infinity`: `true` for types having objects representing \f$+\infty\f$ and
     \f$-\infty\f$.
   
   - `has_uniform_rng`: `true` if the type provides a uniform random number generator. Such a
     context exhibits a type `UniformRNG`, whose interface will be described in more detail above in
     \ref arithctxts_urng. Note that the rational number context from `rational.hpp` is an example
     of a context not providing a uniform random number generator.
   
   - `has_context`: `true` if this type has an associated context as described in \ref
     arithctxts_contexts. If `has_context` is `true`, a typedef `Context` yields the context type.
   
   - `is_native`: `true` only for native types `int`, `long`, `long long`, `float`, `double`, `long
     double` and the unsigned counterparts of the integer types.
   
   For number types, also the following types are defined:
   
   - `PromoteType`: another type to which automatic promotes are a good idea. For example, `int` and
     `unsigned int` have `PromoteType` defined as `long` and `unsigned long`, respectively. For most
     types, `PromoteType` equals `T`. There must always be implicit conversions from `T` to
     `PromoteType`, which are required to be very efficient.
   
   - `ConstReferenceType`: the preferred way of passing an equivalent of a const reference. This is
     often equal to `const T &`, but might also be `T` for native types and encapsulated native
     types for which a reference yields unnecessary overhead (because copying is very cheap).
   
   For integer types, the following property is available in addition to the above general
   properties:
   
   - `is_modulo`: `true` for types which are "wrap around", such as native integer arithmetic which
     in fact computes in rings such as \f$\mathbb{Z}/2^n\mathbb{Z}\f$.
   
   For integer types, all other properties defined for real types -- `is_variable_precision`,
   `has_squareroot`, `has_full_power`, `has_special_fns`, `has_huge_exponents`, `has_constants`,
   `has_trigonometric` -- are set to `false`. Vice versa, for real types, `is_modulo` is set to
   `false`.
   
   For real types, `is_modulo` is always `false`. Real types have the following properties:
   
   - `is_variable_precision`: `true` if the precision of the type can be modified.
     
   - `has_squareroot`: `true` if (approximated) square roots are supported. This will be `false` for
     rational numbers, but `true` for most (or even all) floating point types. For more information,
     see \ref arithctxts_sqrtfp.
     
   - `has_full_power`: `true` if elements can be raised to fractional powers (where the powers are
     coming from the same context). This will be `false` for rational numbers, but `true` for most
     (or even all) floating point types. For more information, see \ref arithctxts_sqrtfp.
     
   - `has_special_fns`: `true` if special functions such as `lgamma` etc. are supported for this
     type. A list of special functions with explanations can be found above in \ref
     arithctxts_specfns.
     
   - `has_huge_exponent`: `true` for floating point types using at least the range of a `signed int`
     for the exponent. If the exponent range is (much) smaller, this has to be set to `false`.
     
   - `has_constants`: `true` if special constants such as \f$\varepsilon\f$, \f$\pi\f$, \f$e\f$ and
     \f$\log 2\f$ are provided by the context. See above in \ref arithctxts_consts for more
     information.
     
   - `has_trigonometric`: `true` if trigonometric functions such as `sin` and `cos` are
     provided. See above in \ref arithctxts_trig for more information.
   
   
   \section arithctxts_convs Conversions
   
   \subsection arithctxts_convs_w_context Conversions with Contexts
   
   Assume we are given a context `Context` with instantiation `c`, a variable `v` of type
   `Context::Type`, and another variable `w` of another type. To store into `v` the converted value
   of `w`, we can use one of the two following calls:
   
   - `convert(v, w, c)`,
   - `c = convert(w, c)`.
   
   If `w2` is another variable of the same type as `w`, and if `Context::is_realtype` is `true`, we
   can also convert the fraction `w/w2` by calling:
   
   - `convert_fraction(v, w, w2, c)`,
   - `v = convert_fraction(w, w2, c)`.
   
   In case `w` is a real type and `v` is an integer type (i.e. `Context::is_inttype` is `true`), the
   following calls allow to convert the floor of `w`, the ceiling of `w` or a rounded value of `w`
   to `v`:
   
   - `convert_floor(v, w, c)`,
   - `v = convert_floor(w, c)`,
   - `convert_ceil(v, w, c)`,
   - `v = convert_ceil(w, c)`,
   - `convert_round(v, w, c)`,
   - `v = convert_round(w, c)`.
   
   In case we want to know whether `w` was rounded up or down to achieve `v` in the last two calls,
   we can use the following syntax:
   
   - `convert_round(v, w, rounded_up, c)`,
   - `v = convert_round(w, rounded_up, c)`;
   
   here, `rounded_up` must be a non-`const` variable of type `bool`, which is set to `true` in case
   `w` was rounded up, and set to `false` in case `w` was rounded down. Note that this is not always
   accurate, as the rounding is often done before the conversion and the conversion might modify the
   value due to precision restrictions of the destination type.
   
   \subsection arithctxts_convs_native Conversions to Native Types
   
   Note that while the above functions also work with native types for `w` and `w2`, they cannot be
   used to convert to a native type (since they don't have contexts). To convert `w` to, say, an
   `int`, we have to use special forms:
   
   - `v = convert<T>(w)`,
   - `v = convert_floor<T>(w)`,
   - `v = convert_ceil<T>(w)`,
   - `v = convert_round<T>(w)`,
   - `v = convert_round<T>(w, rounded_up)`.
   
   Here, `T` can be any of `int`, `long`, `long long`, `float`, `double`, `long double` and their
   unsigned counterparts.
   
   As a special case, `T` can also be `plll::arithmetic::Integer`, which allows to create arbitrary
   precision integers using the same syntax.
   
   \subsection arithctxts_convs_string Conversions to and from Strings
   
   To (try to) convert a string to a type, with or without context, one can also use the above
   conversion functions. If `s` is a string, i.e. of types `char *`, `const char *` or
   `std::string`, it can be converted to a variable `v` of type `Context::Type` by
   
   - `convert(v, s, c)`,
   - `v = convert(s, c)`,

   or to a variable `n` of native type `T` (or `plll::arithmetic::Integer`) by
   
   - `n = convert<T>(s)`.
   
   To convert a variable to a string, note that only conversion to `std::string` is allowed to avoid
   the many pitfalls with C strings (`char *` and `const char *`). For this, also the native
   conversion can be used:
   
   - `s = convert<std::string>(v)`;
   
   here, `v` can be a type with context, or a native type without context.
   
   Another way to convert numbers to strings is to use the `arithmetic::StringContext` and its
   derived types, `arithmetic::HexStringContext` and `arithmetic::OctalStringContext`. With them,
   one can use the usual `convert()` function: if `v` is a number and `s` a `std::string` variable,
   
   - `s = convert(v, s, arithmetic::StringContext())`,
   - `s = convert(v, arithmetic::StringContext())`.
   
   Instead of using the default constructor `arithmetic::StringContext::StringContext()`, which
   generates a `arithmetic::StringContext` object for base 10 (decimal) conversion, one can also
   specify a base via the `arithmetic::StringContext::StringContext(unsigned)` context, or use
   `arithmetic::HexStringContext` to fix base 16 (hexadecimal) or `arithmetic::OctalStringContext`
   to fix base 8 (octal).
   
   Note that types do not have to support all bases. Floating point types might support only base
   10, and integer types might only support bases 8, 10 and 16. Only `plll::arithmetic::Integer` is
   guaranteed to support bases 2 to 36.
   
   
   \section arithctxts_resvals Result Types of Arithmetic Operations
   
   In C++11, it is simple to find out the result of a operation such as `a + b` using <a
   href="http://en.cppreference.com/w/cpp/language/decltype">`decltype()`</a>. Unfortunately, in
   C++03 and C++98, this is not available. Additionally, when using expression templates, a finer
   distinction between intermediate result of `a + b`, which might be something like `addition<Type,
   Type>(a, b)`, and the finaly type of `a + b`, which will be something like `Type`, can be quite
   important.
   
   The utility templates `plll::arithmetic::binary_operation<>` and
   `plll::arithmetic::unary_operation<>` provide such information. Each instance provides some
   properties as well as two types: `ResultType`, the final result type, and `IntermetiateType`, the
   intermediate result type.
   
   The two properties are:

   - `supported`: `true` if information about these types is available. Usually `false` never
     appears since in that case, the template is usually not defined.

   - `intermediate_expression`: `true` if an intermediate expression is returned which is not of the
     final type. More precisely, `true` if and only if `ResultType` is not equal to
     `IntermediateType`.

   To query the result of a binary expression such as `a + b`, one has to use the
   `plll::arithmetic::binary_operation` template:

       plll::arithmetic::binary_operation<Type_of_A, Type_of_B, operation>

   Here, `operation` is one of the classes `plll::arithmetic::op::addition` for addition,
   `plll::arithmetic::op::subtraction` for subtraction, `plll::arithmetic::op::multiplication` for
   multiplication, `plll::arithmetic::op::division` for division, and `plll::arithmetic::op::modulo`
   for modulo.
   
   To query the result of a unary expression such as `-a`, one has to use the
   `plll::arithmetic::unary_operation` template:

       plll::arithmetic::unary_operation<Type_of_A, operation>

   Here, `operation` must be `plll::arithmetic::op::negation` for negation.
   
   A possible use-case is a proxy function which should return the product of two elements, like
   coefficients of a vector or a matrix. We know that the coefficients have types `A` and `B`. Then
   the return type of the function should better be `plll::arithmetic::binary_operation<A, B,
   plll::arithmetic::op::multiplication>::IntermediateType` to avoid too early evaluation. Extensive
   examples can be found in the matrix and vector classes; see \ref matrixvector.
   
   
   \section arithctxts_types Available Contexts
   
   \subsection arithctxts_types_main Main Arithmetic
   The following contexts are always made available in the `plll` library:
   
   - `plll::arithmetic::IntegerContext`: a context providing arbitrary-precision integer arithmetic;
   - `plll::arithmetic::RealContext`: a context providing arbitrary-precision floating point
     arithmetic;
   - `plll::arithmetic::RationalContext` (in `rational.hpp`): a context providing
     arbitrary-precision <a href="https://en.wikipedia.org/wiki/Rational_number">rational</a>
     arithmetic;
   - `plll::arithmetic::NIntContext<T>` (in `arithmetic-nint.hpp`): an integer context for the
     native CPU integer types `T == int`, `T == long` and `T == long long`.
   
   Note that the corresponding types for the first three contexts, `plll::arithmetic::Integer`,
   `plll::arithmetic::Real` and `plll::arithmetic::Rational`, all employ <a
   href="https://en.wikipedia.org/wiki/Expression_templates">expression templates</a> to increase
   efficiency. Therefore, writing
   
       v = a * b;
   
   instead of
   
       mul(v, a, b);
   
   will make no difference. In fact, both lines should compile to exactly the same result. This is
   also the case for add-mul and sub-mul statements such as
   
       v += a * b
       v -= a * b
   
   which should compile to
   
       addmul(v, a, b);
       submul(v, a, b);
   
   For `plll::arithmetic::Integer` and `plll::arithmetic::Real`, these instructions have
   implementations tailored to this situation which will often appear during vector and matrix
   operations, such as computing dot products, squared norms and matrix-matrix and matrix-vector
   products.
   
   \subsection arithctxts_types_internal Internal Floating Point Arithmetic
   Besides the above basic contexts, more contexts are used internally which are (currently) not
   exposed to the end user:

   - `plll::arithmetic::NFPContext<T>`: a real context for the native floating point types `T ==
     float`, `T == double` and `T == long double`;

   - `plll::arithmetic::DDQDContext<T>`: a real context for the double double and quad double types
     `dd_real` and `qd_real` from the `qd` <a
     href="http://crd-legacy.lbl.gov/~dhbailey/mpdist/">quad double library</a>.
   
   Note that the double double and quad double types are not available on every platform. Currently,
   to use these contexts outside the library, one has to include the corresponding files from the
   `plll/include-internal` path inside the `plll` distribution.
   
   \subsection arithctxts_types_combine_gmpmpfr Combining with GMP and MPFR
   
   Since `plll::arithmetic::Integer` is based on <a
   href="https://en.wikipedia.org/wiki/GNU_Multiple_Precision_Arithmetic_Library">GMP</a>'s `mpz_t`
   and `plll::arithmetic::Real` is based on <a href="https://en.wikipedia.org/wiki/MPFR">MPFR</a>'s
   `mpfr_t`, it is possible to manipulate such objects directly with GMP or MPFR functions, or
   transfer data between `plll::arithmetic::Integer` and `plll::arithmetic::Real` objects and their
   GMP and MPFR counterparts. To make this possible, both `plll::arithmetic::Integer` and
   `plll::arithmetic::Real` provide methods `getInternal()`, which provide a reference to the
   underlying `mpz_t` (for integers) respectively `mpfr_t` (for floating point numbers) object.
   
   Use these functions with care, as most of `plll` is agnostic to the underlying library used for
   arbitrary-precision integers and floating point numbers, and thus the underlying libraries (now
   GMP and MPFR) might be exchanged with other libraries at some time.
   
   Also note that `plll` changes the memory allocators for both GMP and MPFR when calling
   `initArithmeticThreadAllocators()`. This should be done once at the very beginning of the
   program, before any `plll::arithmetic::Integer` or `plll::arithmetic::Real` object is created,
   and once at the start of every new thread. This will enable thread-local allocators for GMP and
   MPFR objects, which increases the performance of `plll` in multithreaded environments
   dramatically, especially if many threads are used concurrently (for example, during parallel
   enumeration).
 */
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// MATRICES AND VECTORS ////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
   \page matrixvector Matrices and Vectors
   
   \section matrixvector_intro Introduction
   
   The matrix and vector library of `plll` supplies flexible matrix and vector types which can be
   used with any types, and which are efficient when used with \ref arithctxts.
   
   Similarly to the C++ template library <a
   href="https://en.wikipedia.org/wiki/Eigen_%28C%2B%2B_library%29">Eigen</a>, its aim is to provide
   efficient arithmetic while allowing to write operations in a natural way. For example, writing `v
   = a + b + c` for vectors `v`, `a`, `b` and `c` should be as efficient as directly writing
   
       for (unsigned i = 0; i < v.size(); ++i)
           v[i] = a[i] + b[i] + c[i];
   
   This is achieved using <a href="https://en.wikipedia.org/wiki/Template_metaprogramming">template
   meta programming</a>.
   
   
   \section matrixvector_defs The Basic Templates
   
   The most basic template is the `plll::linalg::base_matrix<>` template, which can describe
   matrices and vectors of any kind. In fact, all other templates are more or less typedefs (or
   wrappers, as templated typedefs are not available in C++98 and C++03) build on
   `plll::linalg::base_matrix<>`.
   
   The `plll::linalg::base_matrix<>` template accepts up to five template parameters:
   
        template<typename T,
                 int Rows = plll::linalg::Flexible,
                 int Cols = plll::linalg::Flexible,
                 typename StorageTraits = plll::linalg::storage_traits<T>,
                 bool MathObject = false>
        class plll::linalg::base_matrix;

   The first template parameter, `T`, specifies the coefficient type. The next two template
   parameters, `Rows` and `Cols`, specify whether the number of rows and columns of the matrix are
   already determined on compile-time -- in that case, they can be provided here by a non-negative
   integer -- or whether they are flexible and can be chosen and changed during runtime -- in that
   case, use `plll::linalg::Flexible`.
   
   The fourth template parameter, `StorageTraits`, allows to modify the default storage behavior. By
   default, memory is allocated on the heap using the default allocator, but with a different
   `StorageTraits` class, one could also use a custom allocator, or store the data somewhere
   else. For more information on how such a `StorageTraits` object works, see the example
   implementations in `matrix-mem.hpp`.
   
   The fifth (and last) template parameter, `MathObject`, allows to enable or disable "math
   mode". Only matrices with `MathObject` set to `true` can be added, subtracted or
   multiplied. Also, if `MathObject` is set to `true`, the library assumes that a function `void
   setZero(T &)` is avaiable for type `T` to set instances of type `T` to the value zero.
   
   The following templates are more or less specializations of `plll::linalg::base_matrix<>`:
   
   - The template
     
         template<typename T, int Rows = Flexible, int Cols = Flexible, typename StorageTraits = storage_traits<T> >
         class plll::linalg::math_matrix

     is a short form for `plll::linalg::base_matrix<T, Rows, Cols, StorageTraits, true>`.
   
   - The template
     
         template<typename T, int Cols = Flexible, typename StorageTraits = storage_traits<T>, bool MathObject = false>
         class plll::linalg::base_rowvector

     is a short form for `plll::linalg::base_matrix<T, 1, Cols, StorageTraits, MathObject>`.
   
   - The template
     
         template<typename T, int Rows = Flexible, typename StorageTraits = storage_traits<T>, bool MathObject = false>
         class plll::linalg::base_colvector
     
     is a short form for `plll::linalg::base_matrix<T, Rows, 1, StorageTraits, MathObject>`.
   
   - The template
     
         template<typename T, int Cols = Flexible, typename StorageTraits = storage_traits<T> >
         class plll::linalg::math_rowvector
     
     is a short form for `plll::linalg::base_rowvector<T, Cols, StorageTraits, true>`.
     
   - The template
     
         template<typename T, int Rows = Flexible, typename StorageTraits = storage_traits<T> >
         class plll::linalg::math_colvector
     
     is a short form for `plll::linalg::base_colvector<T, Rows, StorageTraits, true>`.
   
   The templates `plll::linalg::math_matrix<>`, `plll::linalg::math_rowvector<>` and
   `plll::linalg::math_colvector<>` provide math objects, while the templates
   `plll::linalg::base_matrix<>`, `plll::linalg::base_rowvector<>` and
   `plll::linalg::base_colvector<>` provide by default non-math objects.
   
   \warning Note that when compiled in C++98/03 mode, these "aliases" are in fact derived template
            classes. Therefore, if a function requires a `plll::linalg::math_rowvector<>` and
            something else is passed which should be of the same form (for example of type
            `plll::linalg::base_rowvector<T, Cols, StorageTraits, true>`), a
            `plll::linalg::math_rowvector<T, Cols, StorageTraits>` temporary object will be
            copy-constructed from the given one and will be used as the argument.
            
            When compiled in C++11 mode, both types are identical and no temporary creation and copy
            construction is necessary.
   
   
   \section matrixvector_ctors Constructors
   
   - Create an empty matrix with each dimension being either zero, or the (at compile time)
     pre-determined dimension. All coefficients will be default constructed.
     - `plll::linalg::base_matrix::base_matrix()`
   
   - Create a copy of the given matrix or matrix expression.
     - `plll::linalg::base_matrix::base_matrix(const base_matrix<T, Rows, Cols, StorageTraits,
       MathObject> &)`
     - `template<int Rs, int Cs, bool MO> plll::linalg::base_matrix::base_matrix(const base_matrix<T,
       Rs, Cs, StorageTraits, MO> &)`
     - `template<typename S, int Rs, int Cs, typename ST, bool MO>
       plll::linalg::base_matrix::base_matrix(const base_matrix<S, Rs, Cs, ST, MO> &)`
     - `template<typename MatrixType> plll::linalg::base_matrix::base_matrix(const MatrixType &)`
     - `template<template<typename> class Op, typename Data> plll::linalg::base_matrix::base_matrix<T,
       Rows, Cols, StorageTraits, MathObject> & operator = (const implementation::expressions::expr<Op, Data> &)`
   
   - Create a vector with the given number of entries. These constructors are only avaiable if at
     least one of the dimensions was fixed at compile time to 1. The entries are
     default-constructed.
     
     - `plll::linalg::base_matrix::base_matrix(size_type)`
     - `plll::linalg::base_matrix::base_matrix(int)`
   
   - Create a vector with the given number of entries. These constructors are only avaiable if at
     least one of the dimensions was fixed at compile time to 1. The entries are copy-constructed
     from the given initializer object.
        
     - `template<typename S, bool def> plll::linalg::base_matrix::base_matrix(size_type, const
       implementation::Initialize_Impl<S, def> &)`

   - Create a matrix with the given number of rows and columns. The entries are default-constructed.
   
     - `plll::linalg::base_matrix::base_matrix(size_type rows, size_type cols)`
   
   - Create a matrix with the given number of rows and columns. The entries are copy-constructed
     from the given initializer object.
     
     - `template<typename S, bool def> plll::linalg::base_matrix::base_matrix(size_type rows,
       size_type cols, const implementation::Initialize_Impl<S, def> &)`
   
   Note that an initializer object `obj` can be specified by writing
   `plll::linalg::Initialize(obj)`. A default constructed object of type `S` is given by
   `plll::linalg::Initialize<S>()`.
   
   Note that also the vector templates `plll::linalg::base_rowvector<>`,
   `plll::linalg::base_colvector<>`, `plll::linalg::math_rowvector<>` and
   `plll::linalg::math_colvector<>` provide all of the above listed constructors as well.
   
   
   \section matrixvector_ops Operations
   
   \subsection matrixvector_ops_general General Operations
   
   The following global operations are available for all vectors and matrices:
   
   - `void plll::linalg::swap(MatrixTypeA & A, MatrixTypeB & B)`: swaps the contents of matrices `A`
     and `B`, if compile time dimension specifiers do not obstruct this and if the coefficient types
     can be swapped.

   - `void plll::linalg::assign(MatrixTypeA & A, const MatrixTypeB & B)`: copies the content of `B`
     into the matrix `A`. Note that `A` is resized if necessary. Requires that the coefficient type
     of `MatrixTypeB` can be implicitly converted to the coefficient type of `MatrixTypeB`.

   - `void plll::linalg::transpose(MatrixTypeA & A, const MatrixTypeB & B)`: copies the transposed
     content of `B` into the matrix `A`. Note that `A` is resized if necessary. Requires that the
     coefficient type of `MatrixTypeB` can be implicitly converted to the coefficient type of
     `MatrixTypeB`.

     Note that currently, transposing a matrix via `plll::linalg::transpose(A, A)` is _not_ allowed
     and will result in undefined behavior. Most likely, half of the entries of the matrix are lost,
     and the other half duplicated.
   
   - Comparisons using `operator ==` and `operator !=` are available for all matrices having
     coefficent types which can be tested for equality or inequality. These operators test whether
     the matrices have the same format and, if that is the case, whether the entries compare to be
     equal.
   
   - Comparisons using `operator <`, `operator >`, `operator <=` and `operator >=` are available for
     all matrices having coefficent types which can be tested for equality or inequality.
     
     For these comparison operators, a lexicographic approach is used. First, the size of the
     matrices are lexicographically compared: if the first matrix has \f$r_1\f$ rows and \f$c_1\f$
     columns and the second matrix has \f$r_2\f$ rows and \f$c_2\f$ columns, then we obtain the result
     \f[\begin{cases}
       r_1 < r_2: & \text{first matrix is smaller} \\
       r_1 = r_2 \wedge c_1 < c_2: & \text{first matrix is smaller} \\
       c_1 = r_2 \wedge c_1 = c_2: & \text{content has to be checked} \\
       r_1 = r_2 \wedge c_1 > c_2: & \text{first matrix is larger} \\
       r_1 > r_2: & \text{first matrix is larger} \\
     \end{cases}\f].
     In case \f$r_1 = r_2 \wedge c_1 = c_2\f$, i.e. when the matrices are of the same format, the
     content is checked as follows. If the first matrix is \f[\begin{pmatrix} a_{11} & a_{12} &
     \dots \\ a_{21} & a_{22} & \dots \\ \vdots & \vdots & \ddots \end{pmatrix}\f] and the second
     \f[\begin{pmatrix} b_{11} & b_{12} & \dots \\ b_{21} & b_{22} & \dots \\ \vdots & \vdots &
     \ddots \end{pmatrix},\f] then the vectors \f$(a_{11}, a_{12}, \dots, a_{21}, a_{22}, \dots,
     \dots)\f$ and \f$(b_{11}, b_{12}, \dots, b_{21}, b_{22}, \dots, \dots)\f$ are compared
     lexicographically.
     
   \subsection matrixvector_ops_math Math Operations
   
   The following global operations are avaiable for math matrices and vectors:

   - Computing the sum \f$\sum_{i,j} v_{ij}\f$ of all entries \f$(v_{ij})_{ij}\f$ of \f$A\f$:
     - `template<typename RetType> void sum(RetType & result, const MatrixTypeA & A)`
     - `CoefficientType sum(const MatrixTypeA & A)`
   
   - Computing the dot product \f$\sum_{i,j} a_{ij} b_{ij}\f$ of \f$A = (a_{ij})_{ij}\f$ and \f$B =
     (b_{ij})_{ij}\f$:
     - `template<typename RetType> void dot(RetType & result, const MatrixTypeA & A, const MatrixTypeB
       & B)`
     - `CoefficientType dot(const MatrixTypeA & A, const MatrixTypeB & B)`
   
   - Computing the squared norm \f$\sum_{i,j} v_{ij}^2\f$ of all entries \f$(v_{ij})_{ij}\f$ of
     \f$A\f$:
     - `template<typename RetType> void normSq(RetType & result, const MatrixTypeA & A)`
     - `CoefficientType normSq(const MatrixTypeA & A)`
   
   - Addition of matrices:
     - `add(ResultType &, const MatrixTypeA & A, const MatrixTypeB & B)`
     - `ResultType operator + (const MatrixTypeA & A, const MatrixTypeB & B)`
     - `MatrixTypeA & operator += (MatrixTypeA & A, const MatrixTypeB & B)`
     
   - Subtraction:
     - `sub(ResultType &, const MatrixTypeA & A, const MatrixTypeB & B)`
     - `ResultType operator - (const MatrixTypeA & A, const MatrixTypeB & B)`
     - `MatrixTypeA & operator -= (MatrixTypeA & A, const MatrixTypeB & B)`
      
   - Negation:
     - `neg(ResultType &, const MatrixTypeA & A)`
     - `ResultType operator - (const MatrixTypeA & A)`
   
   - Multiplication:
     - `mul(ResultType &, const MatrixTypeA & A, const MatrixTypeB & B)`: Note that there, `A` or
       `B` can also be scalars.
     - `ResultType operator * (const MatrixTypeA & A, const MatrixTypeB & B)`: Note that there, `A`
       or `B` can also be scalars.
     - `MatrixTypeA & operator *= (MatrixTypeA & A, const MatrixTypeB & B)`: Note that there, `B`
       can also be a scalar.
   
   - Componentwise multiplication, i.e. given matrices \f$A = (a_{ij})_{ij}\f$ and \f$B =
     (a_{ij})_{ij}\f$, computes the matrix \f$R = (a_{ij} b_{ij})_{ij}\f$:
     - `void componentwise_mul(ResultType &, const MatrixTypeA & A, const MatrixTypeB & B)`
     - `ResultType componentwise_mul(const MatrixTypeA & A, const MatrixTypeB & B)`
        
   - Combined addition/subtraction and multiplication:
     - `void addmul(MatrixTypeR & R, const MatrixTypeA & A, const MatrixTypeB & B)`: adds the result
       of `A * B` to `R`. Note that there, `A` or `B` can also be scalars.
     - `void submul(MatrixTypeR & R, const MatrixTypeA & A, const MatrixTypeB & B)`: subtracts the
       result of `A * B` from `R`. Note that there, `A` or `B` can also be scalars.
   
   - Division:
     - `div(ResultType &, const MatrixTypeA & A, const ScalarTypeB & B)`
     - `ResultType operator * (const MatrixTypeA & A, const ScalarTypeB & B)`
     - `MatrixTypeA & operator *= (MatrixTypeA & A, const ScalarTypeB & B)`
   
   - Componentwise division, i.e. given matrices \f$A = (a_{ij})_{ij}\f$ and \f$B =
     (a_{ij})_{ij}\f$, computes the matrix \f$R = (\frac{a_{ij}}{b_{ij}})_{ij}\f$:
     - `void componentwise_div(ResultType &, const MatrixTypeA & A, const MatrixTypeB & B)`
     - `ResultType componentwise_div(const MatrixTypeA & A, const MatrixTypeB & B)`
        
   - Modulo:
     - `mod(ResultType &, const MatrixTypeA & A, const ScalarTypeB & B)`
     - `ResultType operator % (const MatrixTypeA & A, const ScalarTypeB & B)`
     - `MatrixTypeA & operator %= (MatrixTypeA & A, const ScalarTypeB & B)`
   
   - Componentwise modulo, i.e. given matrices \f$A = (a_{ij})_{ij}\f$ and \f$B =
     (a_{ij})_{ij}\f$, computes the matrix \f$R = (a_{ij} \bmod b_{ij})_{ij}\f$:
     - `void componentwise_mod(ResultType &, const MatrixTypeA & A, const MatrixTypeB & B)`
     - `ResultType componentwise_mod(const MatrixTypeA & A, const MatrixTypeB & B)`
   
   
   \section matrixvector_io Input and Output
   
   The usual ways of reading and writing objects in C++ are also provided for matrices and
   vectors. Namely, one can read them from an input stream using `operator >>`, and write them to an
   output stream using `operator <<`. Errors are reported by setting the stream's
   `std::ios_base::badbit` flag.
   
   The output formats are:
   
   - `mat<2, 3>[[1, 2, 3], [4, 5, 6]]`: the \f$2 \times 3\f$ matrix \f$\begin{pmatrix} 1 & 2 & 3 \\
     4 & 5 & 6 \end{pmatrix}\f$
   
   - `vec<3>[1, 2, 3]`: the column vector \f$\begin{pmatrix} 1 \\ 2 \\ 3 \end{pmatrix}\f$
   
   In case the given matrix has precisely one column (determined at compile time), the second format
   is used. If this is not the case, always the first format is used.
   
   To input matrices and vectors, they can be provided in both the above formats, as well as the
   following format:
   
   - `[[1, 2, 3], [4, 5, 6]]`: the \f$2 \times 3\f$ matrix \f$\begin{pmatrix} 1 & 2 & 3 \\ 4 & 5 & 6
     \end{pmatrix}\f$
   
   This allows to exchange information with other libraries such as `NTL` and `fplll`, which can
   read and write matrices in this format. Note that also the <a
   href="http://latticechallenge.org/svp-challenge/">SVP challenge</a> matrices are given in this
   format.
   
   To allow more options while reading or writing matrices, the following functions are provided:
   
   - `void print_matrix(std::ostream &, const MatrixType & mat, bool forceGeneric = false)`: prints
     the matrix `mat` on the given output stream. In case `forceGeneric` is `true`, the output will
     always be of the form `mat<Rows, Cols>[...]`, while in case `forceGeneric` is `false`, column
     vectors will be printed as `vec<Rows>[...]`.
   
   - `bool scan_matrix(std::istream &, MatrixType & mat)`: reads the matrix `mat` from the input
     stream. In case of errors, returns `false`. Elements are created by first default-constructing
     them, and then using their `operator >>` to read their value in.
   
   - `bool scan_matrix(std::istream &, MatrixType & mat, const S & default_object)`: reads the
     matrix `mat` from the input stream. In case of errors, returns `false`. Elements are created by
     first copy-constructing them from `default_object`, and then using their `operator >>` to read
     their value in. Note that `default_object` could be a context object as described in \ref
     arithctxts_contexts to create entries of a given precision.
 */
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// LATTICE REDUCTION ///////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
   \page latticereduction Lattice Reduction
   
   \section howitworks How lattice reduction works
   
   If we are given a basis \f$b_1, \dots, b_k\f$ of a lattice \f$\Lambda\f$, we have \f[ \det
   \Lambda \le \prod_{i=1}^k \|b_i\|, \f] with equality if and only if the \f$b_i\f$'s are pairwise
   orthogonal. (Here, \f$\|\cdot\|\f$ denotes the Euclidean norm on \f$\mathbb{R}^n\f$.) The less
   orthogonal the vectors are, the larger the product on the right-hand side is. One notion of
   lattice reduction is to minimize the difference between the product on the right-hand side and
   the determinant \f$\det \Lambda\f$: then the vectors \f$b_i\f$ must be relatively short.
   
   This motivates the idea of *size reduction*: for this, one first computes the Gram-Schmidt
   orthogonalization of \f$b_1, \dots, b_k\f$ by iteratively, \f$i = 1, \dots, k\f$, computing: \f[
   \mu_{ij} := \frac{\langle \hat{b}_i, \hat{b}_j \rangle}{\langle \hat{b}_i, \hat{b}_i \rangle}
   \quad \text{for } j = 1, \dots, i - 1, \quad \text{and} \quad \hat{b}_i := b_i - \sum_{j=1}^{i-1}
   \mu_{ij} \hat{b}_j. \f] The resulting vectors \f$\hat{b}_1, \dots, \hat{b}_k\f$ are pairwise
   orthogonal and the subvector space spanned by \f$b_1, \dots, b_i\f$ is also generated by
   \f$\hat{b}_1, \dots, \hat{b}_i\f$, \f$0 \le i \le k\f$.

   Essentially, \f$\hat{b}_i\f$ is the orthogonal projection of \f$b_i\f$ onto the orthogonal
   complement of the span of \f$b_1, \dots, b_{i-1}\f$. Denote this projection by \f$\pi_i\f$; then
   \f$\pi_i(b_i) = \hat{b}_i\f$, and more generally \f[ \pi_j(b_i) = \sum_{\ell=j+1}^i \mu_{i\ell}
   \hat{b}_\ell, \qquad 0 \le j \le i. \f] Since the \f$\hat{b}_i\f$'s are pairwise orthogonal, we
   get \f[ \|\pi_j(b_i)\|^2 = \sum_{\ell=j+1}^i \mu_{i\ell}^2 \|\hat{b}_\ell\|^2. \f]
   
   While the vectors \f$\hat{b}_1, \dots, \hat{b}_k\f$ are pairwise orthogonal and generate the same
   subvector space as \f$b_1, \dots, b_k\f$, they do in general _not_ generate the same lattice. In
   fact, even if \f$b_i \in \mathbb{Z}^n\f$ for every \f$i\f$, it could be that \f$\hat{b}_i \in
   \mathbb{Q}^n \setminus \mathbb{Z}^n\f$ for every \f$i > 1\f$.
   
   To obtain a basis of the same lattice which is somewhat more orthogonal than the original basis,
   one can replace \f$b_2\f$ by \f$b_2 - \lfloor \mu_{21} \rceil b_2\f$; if one computes
   \f$\mu_{21}\f$ of the resulting new basis, one obtains \f$-\tfrac{1}{2} \le \mu_{21} \le
   \tfrac{1}{2} \f$. Iteratively repeating this, one obtains a basis \f$\tilde{b}_1, \dots,
   \tilde{b}_k \f$ of the same lattice which satisfies \f$-\tfrac{1}{2} \le \mu_{ij} \le
   \tfrac{1}{2} \f$ for every \f$i, j \f$ with \f$1 \le j < i \le k\f$. Such a basis is called *size
   reduced*, and the process of reaching such a basis is called *size reduction*.
   
   Essentially all lattice reduction algorithms can be described by alternating the process of size
   reduction with another operation.
   
   \subsection howitworks-lll Basic reduction: LLL and friends
   
   For a basis to be LLL-reduced (with reduction parameter \f$\alpha \in (1/4, 1]\f$), one requires
   it to be size reduced and satisfy the Lovász condition: \f[ \alpha \cdot \|\pi_i(b_i)\|^2 \le
   \|\pi_i(b_{i+1})\|^2 \quad \text{for all } i = 1, \dots, k - 1. \f] For \f$\alpha = 1\f$, this
   condition ensures that \f$\pi_i(b_i)\f$ is a shortest vector in the lattice spanned by
   \f$\pi_i(b_i)\f$ and \f$\pi_i(b_{i+1})\f$.
   
   The LLL reduction algorithm is performed by finding the smallest \f$i\f$ for a sized reduced
   basis \f$b_1, \dots, b_k\f$ for which this inequality is not satisfied, and switching \f$b_i\f$
   and \f$b_{i+1}\f$ and repeating by first size reducing and then finding the then minimal
   \f$i\f$. If no such \f$i\f$ is found, the algorithm is done. In case \f$\alpha < 1\f$,
   A. Lenstra, H. Lenstra and L. Lovász \cite lenstra-lenstra-lovasz-LLL showed that the algorithm
   terminates in time polynomially in \f$k\f$ and \f$\max\{ \|b_i\| \mid 1 \le i \le k \}\f$
   (assuming \f$b_i \in \mathbb{Z}^n\f$ for every \f$i\f$). For \f$\alpha = 1\f$, it is not even
   known whether the algorithm terminates in general.
   
   Note that a \f$\alpha\f$-LLL reduced basis satisfies \f[ \|b_1\| \le (\alpha - 1/4)^{-(k-1)/2}
   \cdot \lambda_1(\Lambda) \quad \text{and} \quad \|b_1\| \le (\alpha - 1/4)^{-(n - 1)/4} \cdot
   (\det \Lambda)^{1/n}. \f] Here, \f[ \lambda_1(\Lambda) := \min\{ \|v\| \mid v \in \Lambda
   \setminus \{ 0 \} \} \f] denotes the length of a shortest non-zero vector of \f$\Lambda\f$.
   
   Note that if \f$k = 2\f$ and \f$\alpha = 1\f$, the algorithm is guaranteed to terminate and yield
   a shortest vector of the lattice. In that case, the algorithm is equal to an algorithm already
   known go Gauss, called *Gauss reduction*.
   
   There exist variations of the LLL algorithm. They usually relax/replace the Lovász condition,
   such as by the Siegel condition \f[ \tfrac{3}{4} \cdot \alpha \cdot \|\pi_k(b_k)\|^2 >
   \|\pi_{k+1}(b_{k+1})\|^2, \f], or by adding additional constraints, such as in the case of Deep
   Insertions \cite schnorr-euchner-BKZ.
   
   \subsection howitworks-svphkz Strong reduction: SVP and HKZ bases
   
   While for \f$k = 2\f$ and \f$\alpha = 1\f$, LLL reduction produces a shortest vector, this is not
   necessarily true as soon as \f$k > 2\f$ or \f$\alpha < 1\f$. In fact, Schnorr showed that the
   bound \f$\|b_1\| \le (\alpha - 1/4)^{-(n - 1)/4} \cdot (\det \Lambda)^{1/n}\f$ is sharp.
   
   Nonetheless, a shortest vector exists for the lattice. We call a basis \f$b_1, \dots, b_k\f$ of
   \f$\Lambda\f$ a *SVP basis* if \f$b_1\f$ is a shortest vector of \f$\Lambda\f$; in other words,
   it is a SVP basis iff \f$\|b_1\| = \lambda_1(\Lambda)\f$. Here, SVP stands for *Shorest Vector
   Problem", which denotes the problem of finding a vector \f$v \in \Lambda\f$ with \f$\|v\| =
   \lambda_1(\Lambda)\f$.
   
   A SVP basis can be computed practically, but the fastest algorithms are of single exponential
   complexity in \f$k\f$; in fact, it was proven by Ajtai that computing \f$v \in \Lambda\f$ with
   \f$\|v\| = \lambda_1(\Lambda)\f$ is NP hard (under randomized reductions)
   \cite ajtai-SVPisNPhard. Assuming that there are no subexponential or even polynomial algorithms
   to solve NP hard problems, such algorithms are essentially optimal from an asymptotic point of
   view.
   
   In case \f$\|b_1\|^2 \le \beta \cdot \lambda_1(\Lambda)^2\f$ for some \f$\beta \ge 1\f$, \f$b_1,
   \dots, b_k\f$ is called a *\f$\beta\f$-SVP basis*. Most of the time, such bases are computed by
   computing a shortest vector, comparing its length to \f$\|b_1\|\f$, and replacing \f$b_1\f$ by
   that vector only if it is shorter by a factor of at least \f$\beta\f$. The notion of
   \f$\beta\f$-SVP bases is mostly used during other reduction algorithms.
   
   Since (\f$\beta\f$-)SVP bases are only defined by a property on their first vector, all other
   vectors can be quite large and "ugly". Therefore, one is interested in also reducing them. A good
   notion is the one of a *\f$\beta\f$-Hermite-Korkine-Zolotarev (HKZ) basis*, where one requires a
   basis \f$b_1, \dots, b_k\f$ to be size reduced and to satisfy \f[ \|\pi_i(b_i)\|^2 \le \beta
   \cdot \lambda_1(\Lambda(\pi_i(b_i), \dots, \pi_i(b_k)))^2 \quad \text{for all } i \in \{ 1,
   \dots, k \}. \f] Again, if \f$\beta = 1\f$, \f$b_1, \dots, b_k\f$ is just called a *HKZ basis*
   (instead of 1-HKZ basis). A HKZ basis can be computed by iteratively computing SVP bases and
   doing size reduction. The required running time is up to a polynomial factor identical to the one
   for computing a single SVP basis for lattices of rank \f$k\f$.
   
   \subsection howitworks-bkz Intermediate reduction: BKZ and friends
   
   While the exponential time SVP/HKZ basis reduction algorithms compute a shortest vector for the
   whole lattice, the polynomial time LLL algorithm only computes shortest vectors for
   two-dimensional orthogonally projected lattices \f$\Lambda(\pi_i(b_i), \pi_i(b_{i+1}))\f$, \f$1
   \le i < k\f$. Schnorr suggested \cite schnorr-hierarchyLBRA to interpolate between these two
   extrema by considering orthogonally projected lattices of larger rank.

   Fixing a block size \f$b \in \{ 2, \dots, k \}\f$, a basis is *\f$\alpha\f$-BKZ reduced* with
   *blocksize \f$b\f$* if it is size-reduced and if \f[ \alpha \cdot \|\pi_i(b_i)\|^2 \le
   \lambda_1(\Lambda(\pi_i(b_i), \dots, \pi_i(b_{\min\{ i+b-1, k \}})))^2 \quad \text{for all } i
   \in \{ 1, \dots, k \}. \f] Schnorr and Euchner showed \cite schnorr-euchner-BKZ that a 1-BKZ
   reduced basis with blocksize \f$b\f$ satisfies \f[ \|b_1\|^2 \le b^{(1+\log b) (k-1)/(b-1)} \cdot
   \lambda_1(\Lambda)^2 \f] in case \f$b-1\f$ divides \f$k-1\f$.
   
   BKZ is currently the most practical lattice reduction algorithm to find short vectors in all
   dimensions. Playing with the block size \f$b\f$ and the reduction parameter \f$\alpha\f$ as well
   as introducing early aborts (i.e. stopping the lattice reduction after a specified time or when a
   short enough vector is found) allow to adjust the algorithm to many practical situations.
   
   
   \section desc-algs Supported algorithms
   
   In this section, we will discuss the algorithms provided by this library.
   
   
   \subsection desc-algs-lll LLL-like reduction
   
   `plll` implements several variants of the LLL algorithm. They differ only by their choice of the
   swap condition. Deep Insertions are implemented orthogonally to this and can be enabled for LLL-
   and BKZ-type algorithms; this is described in \ref desc-algconf-di.
   
   \subsubsection desc-algs-lll-classic Classic LLL
   
   The original condition by A. Lenstra, H. Lenstra and L. Lovász \cite lenstra-lenstra-lovasz-LLL,
   called the *Lovász condition* (see \ref howitworks-lll), which compares the projected lengths of
   two adjacent vectors, where the vectors are projected onto the orthogonal complement of all
   previous vectors.
   
   In terms of the Gram-Schmidt orthogonalization \f$\pi_1(b_1), \dots, \pi_k(b_k)\f$ of the basis
   \f$b_1, \dots, b_k\f$, this condition can be expressed as \f$\alpha \cdot \|\pi_k(b_k)\|^2 >
   \|\pi_k(b_{k+1})\|^2\f$.
   
   Classic LLL reduction can be used by calling `plll::LatticeReduction::lll()` with the LLL mode
   `plll::LatticeReduction::LLL_Classic`.
   
   \subsubsection desc-algs-lll-unproj Unprojected LLL
   
   This variant of LLL uses a simpler condition, which compares the _unprojected_ lengths of two
   adjacent basis vector. This essentially attempts to sort the vectors by increasing (unprojected)
   length while at the same time size-reducing the basis.
   
   In terms of the basis \f$b_1, \dots, b_k\f$, this condition can be expressed as \f$\alpha \cdot
   \|b_k\|^2 > \|b_{k+1}\|^2\f$.
   
   Unprojected LLL reduction can be used by calling `plll::LatticeReduction::lll()` with the LLL mode
   `plll::LatticeReduction::LLL_Unprojected`.
   
   \subsubsection desc-algs-lll-siegel Siegel LLL
   
   This variant of LLL uses a simpler condition, which allows to prove the same bounds on the output
   quality, called the *Siegel condition*. Note that the Lovász condition implies the Siegel
   condition.
   
   In terms of the Gram-Schmidt orthogonalization \f$\pi_1(b_1), \dots, \pi_k(b_k)\f$ of the basis
   \f$b_1, \dots, b_k\f$, this condition can be expressed as \f$\frac{3}{4} \cdot \alpha \cdot
   \|\pi_k(b_k)\|^2 > \|\pi_{k+1}(b_{k+1})\|^2\f$.
   
   Siegel LLL reduction can be used by calling `plll::LatticeReduction::lll()` with the LLL mode
   `plll::LatticeReduction::LLL_Siegel`.
   
   
   \subsection desc-algs-bkz BKZ-like reduction
   
   The `plll` library implements many different BKZ variants. They all have in common that they
   opreate by applying SVP or HKZ basis computations to projected blocks (or their duals).
   
   \subsubsection dest-algs-bkz-schnorreuchner Schnorr-Euchner BKZ
   
   The classical Schnorr-Euchner BKZ algorithm, as described by C.-P. Schnorr and M. Euchner in
   Section 6 of
   \cite schnorr-euchner-BKZ.
   
   Schnorr-Euchner BKZ reduction can be used by calling `plll::LatticeReduction::bkz()` with the BKZ
   mode `plll::LatticeReduction::BKZ_SchnorrEuchner`.
   
   \subsubsection dest-algs-bkz-simplified Simplified BKZ
   
   A simplified version of the BKZ algorithm, as described by G. Hanrot, X. Pujol and D. Stehle in
   \cite hanrot-pujol-stehle-BKZdynamic.
   
   Simplified BKZ reduction can be used by calling `plll::LatticeReduction::bkz()` with the BKZ mode
   `plll::LatticeReduction::BKZ_Simplified`.
   
   \subsubsection dest-algs-bkz-terminating Terminating BKZ
   
   A "Terminating BKZ" variant of the simplified version of the BKZ algorithm, as described by
   G. Hanrot, X. Pujol and D. Stehle in \cite hanrot-pujol-stehle-BKZdynamic. There are two
   versions. The first computes Hermite-Korkine-Zolotarev (HKZ) bases for every block, and the
   second Shortest Vector (SVP) bases for every block.
   
   Terminating BKZ reduction can be used by calling `plll::LatticeReduction::bkz()` with the BKZ
   mode `plll::LatticeReduction::BKZ_HanrotPujolStehleHKZ` respectively
   `plll::LatticeReduction::BKZ_HanrotPujolStehleSVP`.
   
   \subsubsection dest-algs-bkz-primaldual Primal-Dual BKZ
   
   The primal-dual BKZ algorithm by H. Koy, as described in Section 3 of
   \cite schnorr-bkzrevisited.
   
   This implementation is a mixture of three slightly varying descriptions of the algorithms:
   
   - the description by H. Koy in the slides from 2004, available <a
     href="http://www.mi.informatik.uni-frankfurt.de/research/papers/primdual.ps">here</a>;
   
   - the description by C.-P. Schnorr in Section 3 of
     \cite schnorr-bkzrevisited;
   
   - the description by C.-P. Schnorr in the lecture notes <a
     href="http://www.mi.informatik.uni-frankfurt.de/teaching/lecture_notes/gitterAKTUELL.pdf">"Gitter
     und Kryptographie"</a>.
   
   The basic idea of the algorithm is to first run LLL on the whole basis and then HKZ on the first
   block is from
   \cite schnorr-bkzrevisited. That each round starts by HKZ-reducing the block \f$\ell+1\f$ is
   shared by
   \cite schnorr-bkzrevisited and the lecture notes, while in Koy's slides, only an SVP basis is
   computed.
      
   As the next step, in Koy's slides a DSVP basis is computed for block \f$\ell\f$, while in the
   lecture notes a dual HKZ basis is computed, and
   \cite schnorr-bkzrevisited computes a dual HKZ basis but only applies the transformation
   conditionally. In the three sources, the aim is to maximize the norm of the last GS vector in the
   \f$\ell\f$-th block.
      
   In
   \cite schnorr-bkzrevisited and the lecture notes, a LLL step is applied to both blocks \f$\ell\f$
   and \f$\ell+1\f$; in case a swap appeared between the two blocks, the changes are taken. (In
   \cite schnorr-bkzrevisited, the HKZ reduction of the dual of block \f$\ell\f$ is only applied to
   the basis in this case. In the lecture notes, the HKZ reduction of the dual is always applied.)
   In Koy's slides, one checks the Lovász condition for the two vectors where the blocks meet, and
   if it is not satisfied LLL reduction is applied to both blocks.
      
   Finally, in
   \cite schnorr-bkzrevisited, no size reduction is applied. In the lecture notes, size reduction is
   applied only at the very end of the algorithm (and the algorithm never increases \f$\ell\f$). In
   Koy's slides, size reduction is not mentioned.
      
   This implementation is close to Koy's slides, with an additional swap for the adjacent vectors
   before running LLL. Additionally, we run size reduction at the end as in the lecture notes.
   
   Primal-Dual BKZ reduction can be used by calling `plll::LatticeReduction::bkz()` with the BKZ
   mode `plll::LatticeReduction::BKZ_PrimalDual`.
   
   \subsubsection dest-algs-bkz-slide Slide Reduction
   
   The slide reduction algorithm, as described by N. Gama and P. Q. Nguyen in
   \cite gama-nguyen-mordellsinequality.
   
   Slide reduction can be used by calling `plll::LatticeReduction::bkz()` with the BKZ mode
   `plll::LatticeReduction::BKZ_SlideReduction`.
   
   \subsubsection dest-algs-bkz-improvslide Improved Slide Reduction
   
   A improved and accelerated version of Slide Reduction, as described in
   \cite schnorr-accslide by C.-P. Schnorr (Section 3).
   
   Improved slide reduction can be used by calling `plll::LatticeReduction::bkz()` with the BKZ mode
   `plll::LatticeReduction::BKZ_ImprovedSlideReduction`.
   
   Note that in the paper, Schnorr describes two variations:
   
   - One variation uses a larger dual SVP enumeration and a single primal SVP enumeration. That one
     can be used with the BKZ mode `plll::LatticeReduction::BKZ_ImprovedSlideReduction2`.

   - The other variation uses a single larger primal SVP enumeration and one dual SVP
     enumeration. That one can be used with the BKZ mode
     `plll::LatticeReduction::BKZ_ImprovedSlideReduction3`.
   
   \subsubsection dest-algs-bkz-semiblock2k Semi-Block-2k-Reduction
   
   The semi block \f$2k\f$-reduction, as described in
   \cite schnorr-bkzrevisited by C.-P. Schnorr (Section 2).
   
   This implementation is a mixture of four varying descriptions of the algorithms:
   
   - C.-P. Schnorr: "A Hierarchy of Polynomial Time Lattice Basis Reduction"
     \cite schnorr-hierarchyLBRA;
   
   - C.-P. Schnorr: "Blockwise Lattice Basis Reduction Revisited"
     \cite schnorr-bkzrevisited;
   
   - C.-P. Schnorr: "Gitter und Kryptographie"; lecture notes, available <a
     href="http://www.mi.informatik.uni-frankfurt.de/teaching/lecture_notes/gitterAKTUELL.pdf">here</a>;
   
   - C.-P. Schnorr: "Progress on LLL and Lattice Reduction"; Chapter 4 in
     \cite nguyen-vallee-LLLbook.
   
   In
   \cite schnorr-hierarchyLBRA, first all blocks are HKZ reduced. Then, one selects the smallest i
   such that the determinant or the Siegel condition are violated. In the latter case, an LLL step
   is applied to the two vectors, and both adjacent blocks are HKZ reduced. In the former case, one
   applies HKZ to the combined block (of double size). This is repeated until the two conditions are
   always satisfies.
   
   In
   \cite schnorr-bkzrevisited, one first computes a HKZ basis of the first block, then applies LLL
   to the whole lattice. Then, for each \f$\ell\f$, one first HKZ-reduces block \f$\ell+1\f$. Then one
   tests whether swapping the adjacent vectors from the two blocks shortens the GS norm of the last
   vector of the \f$\ell\f$-th block by at least a factor of \f$\alpha\f$. Finally, the double block
   is always HKZ reduced, though the changes are not directly applied. (It is not totally clear to
   me if the double block HKZ reduction is always applied, or only if the swap occured.)  Then, the
   new determinant for block \f$\ell\f$ (after HKZ-reduction of double block) is computed and
   compared to the old; if it decreased by a factor of at least \f$\sqrt{\alpha}\f$, the changes
   from the double block HKZ reduction are applied and \f$\ell\f$ is decreased. Otherwise, the
   double-block HKZ-reduction changes are ignored and \f$\ell\f$ is increased.
   
   In the lecture notes, one first computes an LLL basis of the whole lattice, and then a HKZ basis
   of the first block. Then, for each \f$\ell\f$, one first HKZ-reduces block \f$\ell+1\f$. Then, one
   applies HKZ reduction to the double block \f$\ell,\ell+1\f$, but does not directly apply the
   changes. As in
   \cite schnorr-bkzrevisited, they are only applied if the determinant of block \f$\ell\f$ decreases
   by a certain factor (which can be different from the LLL \f$\alpha\f$). At the very end of the
   algorithm, size reduction is applied to the basis.
   
   In
   \cite nguyen-vallee-LLLbook, one first computes an LLL basis of the whole lattice, and then a HKZ
   basis of the first block. Then, for each \f$\ell\f$ (starting with \f$\ell = 2\f$), one first
   HKZ-reduces block \f$\ell+1\f$. Then, one applies LLL reduction to the double block
   \f$\ell,\ell+1\f$, but does not directly apply the changes. If an LLL swap connects the two
   blocks, the LLL reduction is applied, \f$\ell\f$ decreased, and one restarts with this
   \f$\ell\f$. In case no connecting swap appears, one throws away the changes and instead
   HKZ-reduces the double block \f$\ell,\ell+1\f$. As in
   \cite schnorr-bkzrevisited and the lecture notes, one only applies the changes from the HKZ
   reduction in case the determinant decreases by a certain factor; in that case, \f$\ell\f$ is
   decreased, otherwise increased. In case \f$\ell\f$ is 1, one restarts the whole algorithm
   (beginning with LLL-reducing the basis); otherwise one continues with the new value of
   \f$\ell\f$.
   
   Semi block \f$2k\f$-reduction can be used by calling `plll::LatticeReduction::bkz()` with the BKZ
   mode `plll::LatticeReduction::BKZ_SemiBlock2k`.
   
   \subsubsection dest-algs-bkz-sr Sampling Reduction
   
   The Sampling Reduction, as described by J. A. Buchmann and C. Ludwig in
   \cite buchmann-ludwig-samplingreduction.
   
   Sampling reduction can be used by calling `plll::LatticeReduction::bkz()` with the BKZ
   mode `plll::LatticeReduction::BKZ_SamplingReduction`.
   
   
   \subsection desc-algs-svp SVP solvers
   
   Up to dimension 30, all SVP problems are solved by enumeration, except if a specific SVP solver
   is requested (and then only on the highest reduction hierarchy level). The enumeration code there
   is a simplified version of the Kannan-Schnorr-Euchner enumeration below and does not call any
   callback function.
   
   The default choice for enumeration is the parallel Kannan-Schnorr-Euchner enumeration if more
   than one core is available, and the usual Kannan-Schnorr-Euchner enumeration otherwise. In
   practice, it is recommended to use parallel enumeration only for higher dimensions, say 50 and
   above.
   
   \subsubsection desc-algs-svp-kse Kannan-Schnorr-Euchner enumeration
   
   A variant of the Kannan-Schnorr-Euchner enumeration algorithm \cite schnorr-euchner-BKZ with
   various improvements; see, for example, Appendix B of
   \cite gama-nguyen-regev-extremepruning by N. Gama, P. Q. Nguyen and O. Regev.
   
   A parallelized variant of the Kannan-Schnorr-Euchner enumeration algorithm with various
   improvements. The parallelization uses multiple cores. There is almost no communication between
   cores, except if a shorter vector is found by one, or if one core runs out of work and asks
   others to split their workload.
   
   Parallel or non-parallel Kannan-Schnorr-Euchner enumeration can be enabled by calling
   `plll::LatticeReduction::setSVPMode()` with the SVP mode
   `plll::LatticeReduction::SVP_ParallelKannanSchnorrEuchner` respectively
   `plll::LatticeReduction::SVP_KannanSchnorrEuchner`.
   
   \subsubsection desc-algs-svp-schnorr Schnorr's new SVP solver
   
   A (single threaded) version of Schnorr's new SVP solver, as described in
   \cite schnorr-factoringcvp by C.-P. Schnorr.
   
   Schnorr's new SVP solver can be enabled by calling `plll::LatticeReduction::setSVPMode()` with
   the SVP mode `plll::LatticeReduction::SVP_SchnorrFast`.
   
   This implementation is still experimental.
   
   \warning This algorithm is not (yet) working well in high dimensions, and might be slow in medium
            dimensions.
   
   \subsubsection desc-algs-svp-voronoi Voronoi Cell computation
   
   A deterministic SVP solver which first computes the voronoi cell of the lattice. While being
   deterministic and asymptotically faster than enumeration, in practice this algorithm is only
   usable for very low dimensions, say at most 10, where enumeration is still much faster. This
   algorithm was described by D. Micciancio and P. Voulgaris in
   \cite micciancio-voulgaris-voronoicell.
   
   The Voronoi Cell SVP solver can be enabled by calling `plll::LatticeReduction::setSVPMode()` with
   the SVP mode `plll::LatticeReduction::SVP_VoronoiCellSVP`.
   
   \warning This algorithm is only practical in _very low dimensions_!
   
   \subsubsection desc-algs-svp-list List Sieve
   
   The probabilistic list sieve SVP solver. It was described by D. Micciancio and P. Voulgaris in
   Section 3.1 of
   \cite micciancio-voulgaris-fastersvp.
   
   The List Sieve SVP solver can be enabled by calling `plll::LatticeReduction::setSVPMode()` with
   the SVP mode `plll::LatticeReduction::SVP_ListSieve`.
   
   This implementation is still experimental.
   
   \warning This algorithm is not (yet) working well in high dimensions due to large memory
            consumption, and might be quite slow in medium dimensions.
   
   \subsubsection desc-algs-svp-lsb List Sieve Birthday
   
   The probabilistic list sieve (with birthday paradox exploitation) SVP solver. It was described by
   X. Pujol and D. Stehle in
   \cite pujol-stehle-birthdaysieve.
   
   The List Sieve Birthday SVP solver can be enabled by calling
   `plll::LatticeReduction::setSVPMode()` with the SVP mode
   `plll::LatticeReduction::SVP_ListSieveBirthday`.
   
   This implementation is still experimental.
   
   \warning This algorithm is not (yet) working well in high dimensions due to large memory
            consumption, and might be quite slow in medium dimensions.
   
   \subsubsection desc-algs-svp-gauss Gauss Sieve
   
   The probabilistic Gauss sieve SVP solver. It was described by D. Micciancio and P. Voulgaris in
   Section 3.2 of
   \cite micciancio-voulgaris-fastersvp. Our implementation is similar to <a
   href="http://cseweb.ucsd.edu/~pvoulgar/impl.html">one by P. Voulgaris</a>.
   
   The Gauss sieve can be enabled by calling `plll::LatticeReduction::setSVPMode()` with
   the SVP mode `plll::LatticeReduction::SVP_GaussSieve`.
   
   \warning This algorithm is not (yet) working well in high dimensions due to large memory
            consumption.
   
   
   \section desc-algconf Algorithm configuration
   
   There are essentially two algorithmic configurations available for the lattice reduction
   algorithms:
   
   - different kind of Deep Insertions;
   
   - Simulated Annealing.
   
   `plll` supports Schnorr-Euchner Deep Insertions as well as potential minimizing Deep
   Insertions. They are supported (directly and indirectly) by essentially all lattice reduction
   algorithms provided by `plll`.
   
   The Simulated Annealing support is restricted to LLL and BKZ (classical and simplified). It is
   very experimental and should at the moment only be used for further experiments.
   
   \subsection desc-algconf-di Deep Insertions
   
   Classical LLL has two operations: size reduction, and swapping adjacent basis vectors to move a
   shorter vector more to the front. In
   \cite Schnorr-Euchner-BKZ, C.-P. Schnorr and M. Euchner not only describe BKZ, but also a
   modification to classic LLL which allows to insert basis vectors much more in the front if they
   fit somewhere earlier as well: the so-called Deep Insertions variant of LLL.
   
   Another variant of Deep Insertions are potentially minimizing Deep Insertions, introduced by
   F. Fontein, M. Schneider and U. Wagner \cite fontein-schneider-wagner-minpot-wcc
   \cite fontein-schneider-wagner-potlll. This led to two algorithms PotLLL and PotLLL2 which are
   LLL with some kind of Deep Insertions which, opposed to the Schnorr-Euchner Deep Insertions,
   allow to still prove polynomial running time for the algorithm. For Schnorr-Euchner Deep
   Insertions, it is not known whether the algorithm has polynomial running time or not, although
   C.-P. Schnorr and M. Euchner claim this for restricted variants (without proof) in their original
   paper (Comment 2 at the end of Section 3 of
   \cite Schnorr-Euchner-BKZ).
   
   There are four different Deep Insertion methods, which can be set using the
   `plll::LatticeReduction::setDeepInsertionMethod()` function as the first parameter:
   
   - `plll::LatticeReduction::DI_None`: disables Deep Insertions; this is the default;
   
   - `plll::LatticeReduction::DI_Classic`: uses Schnorr-Euchner Deep Insertions;
   
   - `plll::LatticeReduction::DI_MinimizePotential1`: uses potential minimizing Deep Insertions as
     described in
     \cite fontein-schneider-wagner-minpot-wcc and
     \cite fontein-schneider-wagner-potlll;
   
   - `plll::LatticeReduction::DI_MinimizePotential2`: uses a slightly different method for potential
     minimizing Deep Insertions as described in
     \cite fontein-schneider-wagner-potlll.
   
   The second argument of `plll::LatticeReduction::setDeepInsertionMethod()` determines the Deep
   Insertions mode, which decides whether to do Deep Insertions before or/and after size reduction:
   
   - `plll::LatticeReduction::DIM_BeforeSR`: does Deep Insertions before size reduction;
   
   - `plll::LatticeReduction::DIM_AfterSR`: does Deep Insertions after size reduction;
   
   - `plll::LatticeReduction::DIM_Both`: does Deep Insertions both before and after size reduction;
     this is the default.
   
   Note that for `plll::LatticeReduction::DIM_Both`, the Deep Insertion after the size reduction is
   only called if the size reduction did modify the basis. Also note that the mode can be set
   independently from the method via `plll::LatticeReduction::setDeepInsertionMode()`.
   
   Finally, it is possible to configure the range where vectors can be inserted by calling
   `plll::LatticeReduction::setDeepInsertionChoice()`:
   
   - `plll::LatticeReduction::DIC_First`: if the second argument to
     `plll::LatticeReduction::setDeepInsertionChoice()` is `t`, then the current vector will only be
     tried to be inserted among the first `t` vectors of the base;
   
   - `plll::LatticeReduction::DIC_Block`: if the second argument to
     `plll::LatticeReduction::setDeepInsertionChoice()` is `t`, then the current vector will only be
     tried to be inserted among the (at most) `t` vectors before the current basis vector;
   
   - `plll::LatticeReduction::DIC_FirstBlock`: this is a combination of
     `plll::LatticeReduction::DIC_First` and `plll::LatticeReduction::DIC_Block`, i.e. the current
     basis vector can be inserted both at the beginning and in the block before the current vector;
   
   - `plll::LatticeReduction::DIC_All`: the current basis vector can be inserted at any position
     before itself; this is the default.
   
   Note that the default setting, ``plll::LatticeReduction::DIC_All`, seems to yield exponential
   running time. Therefore, one usually chooses one of the other modes.
   
   Currently it is not known which choices of Deep Insertions are particularly good; some research
   indicates that Deep Insertions can be as good as BKZ in some cases
   \cite gama-nguyen-predictinglatticereduction. We are currently running experiments to provide
   guidelines which parameters to choose.
   
   \subsection desc-algconf-anneal Simulated Annealing
   
   \todo Write documentation for Simulated Annealing.
   
   \section desc-genconf General configuration
   
   \subsection desc-genconf-arith Arithmetic

   There are two arithmetics which can be configured: real arithmetic and integer arithmetic.
   
   The real arithmetic is used for computing Gram-Schmidt orthogonalizations or duals of projected
   lattices. It can be set via `plll::LatticeReduction::setArithmetic()` and queried via
   `plll::LatticeReduction::getArithmetic()`. The following arithmetics are available:

   - `plll::LatticeReduction::A_Rational`: rational arithmetic with arbitrary precision
     integers. Will always return correct results, but is quite slow.
   
   - `plll::LatticeReduction::A_Real`: arbitrary precision floating point arithmetic. With a
     numerically stable Gram-Schmidt orthogonalization, it should always return quite correct
     results and be quite fast.
   
   - `plll::LatticeReduction::A_Double`: uses native `double` floating point arithmetic. Should be
     fine up to certain dimensions and very fast. Might result in hanging for too high versions.
   
   - `plll::LatticeReduction::A_LongDouble`: uses native `long double` floating point
     arithmetic. Should be fine up to certain dimensions and very fast. Might result in hanging for
     too high versions.
   
   - `plll::LatticeReduction::A_DoubleDouble`: uses `double double` floating point arithmetic which
     is obtained by representing a floating point numbers as the _sum_ of two `double` values. Has
     higher precision then a `double` variable, but is slower and has no larger exponent range. Might
     not be available on every system.
   
   - `plll::LatticeReduction::A_QuadDouble`: uses `quad double` floating point arithmetic which is
     obtained by representing a floating point numbers as the _sum_ of four `double` values. Has
     higher precision then `double`, `long double` and `double double` variables, but is slower and
     has no larger exponent range than `double` variables. Might not be available on every system.
     
     \warning While conversions between `double double` and `quad double` and some native types (`int`, 
           `float` and `double`) are fully supported (by `libqd`), conversions to and from `long double`, 
           `long int`, `long long` and the GMP and MPFR arbitrary precision types are not properly 
           implemented (yet). So use on your own risk.
     
     \todo Fully implement conversions between `double double`, `quad double` and `long double`, 
           `long int`, `long long` and the GMP and MPFR arbitrary precision types.
   
   Note that for arbitrary precision floating point arithmetic, i.e. for
   `plll::LatticeReduction::A_Real`, one can enforce a minimal precision by calling
   `plll::LatticeReduction::ensurePrecision()`.
   
   On the other hand, the integer arithmetic is used to represent the (integral) lattice's
   entries. It can be set via `plll::LatticeReduction::setIntegers()` and queried via
   `plll::LatticeReduction::getIntegers()`. There are two main different integer arithmetics:
   
   - `plll::LatticeReduction::I_ArbitraryPrecision`: uses arbitrary precision integers. They can
     represent any integer which fits into memory, and no overflow will occur. Slow, but always
     correct.
   
   - `plll::LatticeReduction::I_LongInt`: uses native `long int` integers. They can represent only a
     limited range of integers, and overflow might occur. Much faster, but only suitable for lattice
     bases with small entries.
   
   There is a third option, `plll::LatticeReduction::I_Auto`, which tries to balance between the two
   above choices. It constantly monitors the size of the basis coefficients and switches between the
   two above arithmetics depending on the coefficient sizes.
   
   The default choice is `plll::LatticeReduction::I_Auto`. Note that it can be slower than
   `plll::LatticeReduction::I_ArbitraryPrecision`!
   
   \subsection desc-genconf-gs Gram-Schmidt orthogonalization
   
   `plll` provides a choice of several Gram-Schmidt orthogonalizations. The choice can be set with
   `plll::LatticeReduction::setGramSchmidt()` and queried with
   `plll::LatticeReduction::getGramSchmidt()`.
   
   - `plll::LatticeReduction::G_Classic`: classical Gram-Schmidt orthogonalization. Uses formulae to
                                          update Gram-Schmidt coefficients on swaps, transformations
                                          etc., which can introduce additional error with floating
                                          point arithmetic.
   
   - `plll::LatticeReduction::G_ClassicInteger`: Integer-based Gram-Schmidt
                                                 orthogonalization. Internally uses arbitrary
                                                 precision integers respectively rational
                                                 numbers. Slow, but always accurate.
   
   - `plll::LatticeReduction::G_Givens`: Uses Givens rotations to compute Gram-Schmidt
                                         orthogonalization. Only available for floating-point
                                         arithmetic (i.e. every arithmetic except
                                         `plll::LatticeReduction::A_Rational`). Uses formulae to
                                         update Gram-Schmidt coefficients on swaps, transformations
                                         etc., which can introduce additional error with floating
                                         point arithmetic.
                                         
                                         Note that this is only available for floating point
                                         arithmetic, and not for rational arithmetic.
   
   - `plll::LatticeReduction::G_NumStable`: The default Gram-Schmidt orthogonalization. Uses
                                            numerically stable Gram-Schmidt orthogonalization as
                                            described by P. Q. Nguyen and D. Stehle in
                                            \cite nguyen-stehle-fplll-revisited. This arithmetic is
                                            more resilient against approximation errors than other
                                            approaches. Instead of updating Gram-Schmidt
                                            coefficients for swaps, transformations etc., these are
                                            recomputed to ensure maximal precision.
   
   Note that some of the Gram-Schmidt orthogonalizations support spontaneous restarts (in case
   precision is too bad). This can be toggled via `plll::LatticeReduction::setGramSchmidtRestart()`
   and queried via `plll::LatticeReduction::getGramSchmidtRestart()`. This is only supported for
   `plll::LatticeReduction::G_Classic` and `plll::LatticeReduction::G_Givens` and currently still
   experimental.
   
   \subsection desc-genconf-trans Recording transformation matrices
   
   It is possible to record all transformations done during lattice reduction in form of a
   transformation matrix. Assume that \f$b_1, \dots, b_k\f$ is the input system, and \f$\hat{b}_1,
   \dots, \hat{b}_\ell\f$ the output system.
   
   For notational reasons, let us define two matrices. Let \f$B\f$ be the matrix with \f$k\f$ rows
   having \f$b_i\f$ as its \f$i\f$-th row, and \f$\hat{B}\f$ the matrix with \f$\ell\f$ rows having
   \f$\hat{b}_j\f$ as its \f$j\f$-th row.
   
   There are two transformation recording modes, which can be enabled by calling
   `plll::LatticeReduction::enableTransform()`:
   
   - `plll::LatticeReduction::T_Normal`: computes a transformation matrix \f$T_1\f$ such that
     \f$\hat{B} = T_1 \cdot B\f$;
   
   - `plll::LatticeReduction::T_Inverse`: computes a transformation matrix \f$T_1\f$ such that
     \f$B = T_2 \cdot \hat{B}\f$.
   
   In case \f$k = \ell\f$, both transformation matrices are invertible and \f$T_2 = T_1^{-1}\f$.
   
   Note that \f$B\f$ represents the matrix at the point when
   `plll::LatticeReduction::enableTransform()` was called. A subsequent call, even with the same
   mode, resets the current transformation matrix to the identity matrix.
   
   Transformation recording can be disabled by calling
   `plll::LatticeReduction::disableTransform()`. Its status can be queried by calling
   `plll::LatticeReduction::isTransformationRecorded()` and
   `plll::LatticeReduction::getTransformationMode()`.
   
   In case `plll::LatticeReduction::isTransformationRecorded()` returns `true`, the function
   `plll::LatticeReduction::getTransformation()` returns a non-`NULL` pointer to a matrix containing
   the transformation. In case `plll::LatticeReduction::isTransformationRecorded()` returns `false`,
   it returns `NULL`.
   
   \subsection desc-genconf-callback Callbacks
   
   There exist several <a
   href="https://en.wikipedia.org/wiki/Callback_%28computer_programming%29">callback</a> hooks in
   the `plll` library:
   
   - The most important hook is a callback function which is called in regular intervals from during
     the lattice reduction code. The default interval is approximately every 5 minutes (it can take
     longer when no lattice reduction code is executed, but for example an SVP solver).
     
     The callback function can be set with `plll::LatticeReduction::setCallbackFunction()`. It
     accepts either one argument of type `plll::LatticeReduction::CallbackFunction`, one argument of
     type `plll::LatticeReduction::CallbackFunction_LI`, or two arguments of type
     `plll::LatticeReduction::CallbackFunction` and
     `plll::LatticeReduction::CallbackFunction_LI`. These are `boost::function<>` objects accepting
     a `const` reference to the current lattice basis and should return `void`. Please refer to the
     documentation of the function object types for more details on the function arguments.
     
     Both currently set function objects can be obtained by calling
     `plll::LatticeReduction::getCallbackFunction()`.
     
     The interval can be set by calling `plll::LatticeReduction::setCallbackInterval()`; its
     argument is in seconds (as a `double` floating point number), and the default value is
     `60.0*5.0`, i.e. 5 minutes. This is the minimal waiting time between two calls of the callback
     function.
     
     Note that the callback function can throw an exception of type
     `plll::LatticeReduction::stop_reduction` to interrupt the reduction process. In that case, the
     reduction process will be stopped as fast as possible and control is handed back to the caller
     of the lattice reduction library.
   
   - The next most important hook is a callback function which will be called any time the code
     finds a newest shortest vector. For this, the (absolute) lengths of basis vectors will be
     compared all the time, which might add a small performance penalty. Usually it is negligible.
     
     Such a callback function can be set with `plll::LatticeReduction::setMinCallbackFunction()`. It
     accepts either one argument of type `plll::LatticeReduction::MinCallbackFunction`, one argument
     of type `plll::LatticeReduction::MinCallbackFunction_LI`, or two arguments of type
     `plll::LatticeReduction::MinCallbackFunction` and
     `plll::LatticeReduction::MinCallbackFunction_LI`. These are `boost::function<>` objects
     accepting a `const` reference to the current lattice basis, an `unsigned` index into the
     lattice which indicates the currently shortest touched vector, and a `const arithmetic::Integer
     &` specifying its squared Euclidean norm. The return type should be `void`. Please refer to the
     documentation of the function object types for more details on the function arguments.
     
     Both currently set function objects can be obtained by calling
     `plll::LatticeReduction::getMinCallbackFunction()`.
     
     Note that this callback function will not be called during enumerations (only at the end of the
     enumeration) or during processing of dual lattices.
     
     Note that also this callback function can throw an exception of type
     `plll::LatticeReduction::stop_reduction` to interrupt the reduction process. In that case, the
     reduction process will be stopped as fast as possible and control is handed back to the caller
     of the lattice reduction library.
   
   - The third kind of hook are enumeration callback functions: they are called during enumeration
     (or other SVP solvers) when during this enumeration a new shortest vector (in the projected
     sublattice the enumeration is working on) is detected.
     
     Such a callback function can be set with
     `plll::LatticeReduction::setEnumCallbackFunction()`. It accepts either one argument of type
     `plll::LatticeReduction::EnumCallbackFunction`, one argument of type
     `plll::LatticeReduction::EnumCallbackFunction_LI`, or two arguments of type
     `plll::LatticeReduction::EnumCallbackFunction` and
     `plll::LatticeReduction::EnumCallbackFunction_LI`. These are `boost::function<>` objects
     accepting a `const` reference to the current lattice basis, an `int` index specifying the first
     basis vector involved in the linear combination, as well as a row vector with integer
     coefficients which specifies the linear combination of the shortest vector. The return type
     should be `void`. Please refer to the documentation of the function object types for more
     details on the function arguments.
     
     Both currently set function objects can be obtained by calling
     `plll::LatticeReduction::getEnumCallbackFunction()`.
     
     If the vector is `vec` and the `int` index is `p`, and the basis vectors (i.e. the rows of the
     matrix) are denoted by `b[0]`, `b[1]`, etc., then the vector found is \f[ \sum_{i=0}^{\text{\tt
     vec.size()} - 1} \text{\tt vec[}i\text{\tt ]} \cdot \text{\tt b[p} + i \text{]}. \f]
     
     Note that this callback function will be only called during enumerations of primal lattices,
     not during enumeration of dual lattices.
     
     Note that also this callback function can throw an exception of type
     `plll::LatticeReduction::stop_reduction` to interrupt the reduction process. In that case, both
     enumeration and the reduction process will be stopped as fast as possible and control is handed
     back to the caller of the lattice reduction library.
     
     But it can also throw an exception of type `plll::LatticeReduction::stop_enumeration` to stop
     just enumeration. Then the so far shortest vector will be returned.
     
     Note that in case of multi-threaded enumeration (\ref desc-algs-svp-kse), it can happen that
     new shortest vectors are found during enumeration while the handler is already running in
     another thread. In such cases, the handler has to pay attention to this problem and make sure
     that no data races or (fatal) errors occur.
   
   If for one of the three cases above, two callback functions are given, the better fitting one is
   selected depending on the currently used integer arithmetic (see \ref desc-genconf-arith). In
   case only one callback function is given, it is used for both integer arithmetics.
   
   \warning In case only one callback function is specified in one case, and the incorrect integer
            arithmetic is used, the arguments have to be converted for _every_ function call. This
            can slow down execution a lot!
   
   \subsection desc-genconf-verbose Verbose output
   
   During operation, the library generates a lot of more and less informative messages. By default,
   the most important of these messages are written to `std::cerr`. By calling
   `plll::LatticeReduction::setVerbose()`, one can change the verbosity level (first argument) and
   optionally set a verbose callback function which can output the messages somewhere else, or
   simply store or ignore them.
   
   The verbosity level can attain the following values:
   
   - `plll::LatticeReduction::VOL_None`: nothing is reported;
   - `plll::LatticeReduction::VOL_Warnings`: only warnings and errors are reported;
   - `plll::LatticeReduction::VOL_Informative`: informative messages, warnings and errors are
                                                reported;
   - `plll::LatticeReduction::VOL_Full`: *everything* is reported, including a lot of unnecessary
                                         and potentially annoying messages.
   
   The current verbosity level can be queried by calling
   `plll::LatticeReduction::getVerboseOutputLevel()`, and the current verbose callback function can
   be queried by calling `plll::LatticeReduction::getVerboseFunction()`.
   
   Verbose callback functions accept a verbose level (for the current message) of type
   `plll::LatticeReduction::VerboseLevel` as the first argument, and a `const std::string &` with
   the message itself as the second argument. They return `void`.
   
   The verbose level can have the following values:
   
   - `plll::LatticeReduction::VL_Error`: the message is a (fatal) error;
   - `plll::LatticeReduction::VL_Warning`: the message is a (non-fatal) warning;
   - `plll::LatticeReduction::VL_Information`: the message is purely informational;
   - `plll::LatticeReduction::VL_Chatter`: the message can be safely ignored. It can be helpful for
                                           debugging purposes or to see more about how the
                                           algorithms work.
   
   \subsection desc-genconf-range Reduction Range
   
   All reduction algorithms and SVP solvers can be restricted to a certain projected sublattice. If
   \f$b_1, \dots, b_k\f$ is the current generating system and a range of \f$[begin, end]\f$ is
   given, the lattice generated by \f[ \pi_{begin}(b_{begin}), \pi_{begin}(b_{begin+1}), \dots,
   \pi_{begin}(b_{end}) \f] is operated on.
   
   The current range can be queried by calling `plll::LatticeReduction::getRange()`, and it can be
   set by calling `plll::LatticeReduction::setRange()`. If only one argument is given, this is taken
   as `begin`, while `end` is set to the maximum; otherwise, the input is of the form `(begin,
   end)`.
   
   Note that ranges can be modified if zero vectors are found (due to linear dependencies) and
   eliminated, or when new vectors are inserted (SVP solving without turning the result into a
   basis).
   
   \subsection desc-genconf-multithreading Multi-Threading
   
   Parts of the `plll` library support active multi-threading. At the moment, this is limited to
   parallel enumeration (see \ref desc-algs-svp-kse). To control the maximal number of threads for
   parallel enumeration, one can set this number by calling
   `plll::LatticeReduction::setMaximalCoreUsage()`. Setting it to 0 lets the system decide, which
   usually results in as many threads as the machine has (logical) cores.
   
   The current number of parallel threads can be queried by calling
   `plll::LatticeReduction::getMaximalCoreUsage()`; the default value is 0.
*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// LATTICE REDUCTION ///////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
   \page examples Examples
   
   Note that most features are implemented in the `plll` command line program, which can be found in
   `tools/src/plll.cpp`. Therefore, to see how a certain feature is used which is not presented in this
   section, you can always check out `tools/src/plll.cpp`.
   
   \section example-01 A simple example
   
   This is the same sample as in \ref usage-example. It can be found in `examples/src/example-01.cpp`:
   \include example-01.cpp
   The example first creates an empty (i.e. \f$0 \times 0\f$) matrix with arbitrary precision integer entries,
   called `A`, and then reads a matrix from standard input (`std::cin`) and stores the result in `A`.
   
   Then, an instance called `lr` of the lattice reduction class is created. It is initialized with the lattice
   generated by the rows of the matrix `A`. By calling `lr.lll()`, the \ref desc-algs-lll-classic algorithm is
   applied to this lattice basis with the default parameter \f$\alpha = 0.99\f$.
   
   Finally, the program reads the reduced lattice basis by `lr.getLattice()`---the result is a matrix whose
   rows yield the basis---and writes the matrix to standard output (`std::cout`).
   
   To demonstrate how this example works on real-world data, let us consider the lattice in 
   \ref introlatticered. The first basis given there is \f[ b_1 = \begin{pmatrix} 3 \\ 1 \end{pmatrix},
   \quad b_2 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}, \f], resulting in the following picture:

   \image html lattice-basis1.png "The input basis"
   \image latex lattice-basis1.eps "The input basis" width=10cm
   
   To apply the program to this basis, we have to construct a matrix containing these vectors as its rows:
   \f[ A = \begin{pmatrix} 3 & 1 \\ 1 & 1 \end{pmatrix} \f] To enter this matrix into our program, we have 
   to write it as follows: `[[3, 1], [1, 1]]`
   
   If we enter this at the beginning of our program, we get the following output:
   \verbatim
   [[3, 1], [1, 1]]
   mat<2,2>[[1, 1], [1, -1]]
   \endverbatim
   (The first line is the echo of the input.) Thus, as the output, we obtain the nice lattice basis in 
   \ref introlatticered, which yields the following picture:
   
   \image html lattice-basis3.png "The 'optimal' output basis"
   \image latex lattice-basis3.eps "The 'optimal' output basis" width=10cm
   
   
   \section example-02 Finding a very short vector
   
   While finding a shortest vector is NP hard \cite ajtai-SVPisNPhard, it is possible to approximate a
   short vector by using lattice reduction algorithms. An algorithm very helpful in practice is the
   \ref dest-algs-bkz-schnorreuchner algorithm. Unfortunately, in particular for higher block sizes,
   the algorithm often already obtained the shortest vector it will return long before it terminates.
   For that reason, it is wishful to be able to always see what is the shortest vector found so far.
   `plll` provides a callback hook for this reason which is called when the reduction algorithm stumbles
   over a vector with Euclidean norm shorter than any vector it encountered before while processing
   the current basis (\sa desc-genconf-callback).
   
   The following example shows how to do this in practice. It can be found in 
   `examples/src/example-02.cpp`:
   \include example-02.cpp
   It is similar to the previous example (\ref example-01), with the main differences that:
   
   - it does not output the resulting basis;
   
   - it does not use \ref desc-algs-lll-classic, but \ref dest-algs-bkz-schnorreuchner (with 
     \f$\alpha = 0.99\f$ and a blocksize of 40);
   
   - it sets a minimal vector callback function using `plll::LatticeReduction::setMinCallbackFunction`.
   
   For the latter, it provides two functions: one for arbitrary precision integers 
   (`plll::arithmetic::Integer`) and one for CPU integers (`plll::arithmetic::NInt<long>`). Since by
   default, `plll` determines which integer arithmetic to use on the fly (\sa desc-genconf-arith), it
   is useful to provide both alternatives so that `plll` is not forced to convert the internal
   representation from one format to the other to be able to call the callback function.
   
   To illustrate how this example works, let us use the SVP challenge basis for dimension 50 (with seed 0); 
   it can be found <a href="http://www.latticechallenge.org/svp-challenge/download/challenges/svpchallengedim50seed0.txt">here</a>.
   By starting the program and pasting the basis' textual representation into it (or by using pipes or
   redirecting standard input), we let the program process this lattice basis.
   
   The library internally first applies LLL, which converts this basis (given in Hermite Normal Form)
   into something more compact. The first output line is
   \verbatim
   The currently shortest vector is mat<1,50>[[-9770253083323343916965275912883206062922561406584472239
   11084680829795919819354762585959322324750268253379484023871452642838548115289419328942293214818, 1, 
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]] with squared norm 954578453121893086923031172722240141
   5470291413191066642863537393601930079752899447459843620501562783729668907328211198000701122696066760
   0771389257710458775685171721625238187868083121323492585428385068254902186858647242938612587927038568
   6748545805753278169971353652684403255259051086907300041494773125 (Integer)
   \endverbatim
   (line breaks inserted for readability), which equals the second basis vector reduced by a multiple of
   the first basis vector. The squared norm here is quite large. At some point, `plll` decides to switch
   from arbitrary precision arithmetic to using CPU integers. This happens at this point:
   \verbatim
   The currently shortest vector is mat<1,50>[[-541, 416, 284, 24, 233, 164, -150, 71, -367, -288, -170,
   71, -49, 125, 260, 143, -168, -30, 465, 49, 418, -349, -309, 134, 468, 109, -195, -218, -96, -104, 
   -153, -608, -193, 657, 43, -155, -275, -234, -58, 896, -917, -23, 328, 356, 570, 155, 47, -606, -1,
   0]] with squared norm 5626057
   Changing lattice reduction interface.
   (starting 4 threads)
   The currently shortest vector is mat<1,50>[[364, 142, 213, 304, 786, -923, -71, 123, -7, -81, 323,
   -436, 63, 209, 344, 279, 10, -287, -324, 267, -17, 327, -901, 239, 328, 190, 92, 664, 348, -254, -5,
   340, -833, -157, -477, 561, -82, -41, 203, 271, 612, -458, -130, -26, 383, 16, -252, -75, 0, 0]] with
   squared norm 6656166 (NInt<long>)
   \endverbatim
   The vector returned before restarting the algorithm with CPU integer arithmetic already found a vector
   of Euclidean norm close to 2372, which is much less than the input vectors' norms. After restarting,
   the algorithm starts again with tracking shortest vectors, and thus first reports something much
   longer, namely a vector of Euclidean norm close to 2580. But shortly after that, it already finds a
   vector of Euclidean norm close to 2220, which beats the previously found shortest vector. Finally, it
   stops with a vector of Euclidean norm close to 1893:
   \verbatim
   The currently shortest vector is mat<1,50>[[-13, -124, -146, 277, -107, -180, 673, -311, -167, 47,
   200, 395, 167, -25, -136, -392, 117, -165, 147, -515, 185, 637, 343, 8, 247, 44, -220, -146, 52, 135,
   -347, -369, -332, -102, 469, -285, 1, 167, 397, 84, -97, -138, -135, 218, 567, 141, 72, 21, 312,
   -41]] with squared norm 3584092 (NInt<long>)
   \endverbatim
   It turns out that the shortest vectors in this lattice also have Euclidean norm close to 1893 (see
   <a href="http://www.latticechallenge.org/svp-challenge/halloffame.php">the hall of fame</a> for the
   SVP challenges), whence we can be quite sure that this is one of the shortest vectors. To verify this
   ourselves, we would have to do enumeration, which can be done via calling `plll::LatticeReduction::svp`.
   
   Note that in particular in higher dimensions, the last output will come a long time before the algorithm
   terminates. Thus for an application where you need a short enough vector, and given a vector you can
   decide on the fly whether it was short enough, you could use this callback mechanism to perform this
   computation for every new shortest vector found. Note that from the callback functions, you can
   interrupt computation by throwing an exception of type `plll::LatticeReduction::stop_reduction`. If this
   exception isn't thrown, the reduction will continue until the reduction algorithm terminates.

*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif

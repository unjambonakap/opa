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

inline void setZero(int & i)
{
    i = 0;
}

inline void setZero(long & i)
{
    i = 0;
}

inline void setZero(long long & i)
{
    i = 0;
}

inline void setZero(double & i)
{
    i = 0;
}

inline int square(int i)
{
    return i * i;
}

inline long square(long i)
{
    return i * i;
}

inline long long square(long long i)
{
    return i * i;
}

inline double square(double i)
{
    return i * i;
}

#include <plll/matrix.hpp>
#include <plll/arithmetic.hpp>
#include <plll/arguments.hpp>
#include "old_code/vecmat.hpp"
#include <timer.hpp>
#include <profiling.hpp>
#include <vector>
#include <algorithm>

using namespace plll;

//////////////////////////////////////////////////////////////////////////////////////////////////
/// Helper routines
//////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T>
void touch(T &)
{
}

unsigned seed = 1;

inline unsigned urandom()
{
    return seed = (seed * 213421984 + 1) % 213421987;
}

inline signed srandom()
{
    signed v = urandom();
    return (urandom() & 1) ? -v : v;
}

inline void generate_random_integer(arithmetic::Integer & dest)
{
    unsigned u1 = urandom(), u2 = urandom();
    mpz_set_ui(dest.getInternal(), u2);
    mpz_mul_2exp(dest.getInternal(), dest.getInternal(), 32);
    mpz_add_ui(dest.getInternal(), dest.getInternal(), u1 >> 1);
    if (u1 & 1)
        mpz_neg(dest.getInternal(), dest.getInternal());
}

inline void generate_random_integer(int & dest)
{
    unsigned u = urandom();
    dest = u >> 1;
    if (u & 1)
        dest = -dest;
}

inline void generate_random_integer(long & dest)
{
    unsigned u = urandom();
    dest = u >> 1;
    if (u & 1)
        dest = -dest;
}

inline void generate_random_integer(long long & dest)
{
    unsigned u = urandom();
    dest = u >> 1;
    if (u & 1)
        dest = -dest;
}

inline void generate_random_integer(double & dest)
{
    unsigned u = urandom();
    long d = u >> 1;
    if (u & 1)
        d = -d;
    dest = ldexp(d, -10);
}

template<class matrix>
void fill_matrix(matrix & m)
{
    for (unsigned i = 0; i < m.rows(); ++i)
        for (unsigned j = 0; j < m.cols(); ++j)
            generate_random_integer(m(i, j));
}

template<class vector>
void fill_vector(vector & m)
{
    for (unsigned j = 0; j < m.size(); ++j)
        generate_random_integer(m[j]);
}

//////////////////////////////////////////////////////////////////////////////////////////////////
/// Test code: vector filling
//////////////////////////////////////////////////////////////////////////////////////////////////

template<class Int, class row_vector, class col_vector, class matrix>
class test_fill_vector
{
public:
    void operator() ()
    {
        row_vector v;
        v.resize(200);
        for (unsigned j = 0; j < 30000; ++j)
            fill_vector(v);
    }
    
    void direct()
    {
        row_vector v;
        v.resize(200);
        for (unsigned j = 0; j < 30000; ++j)
            fill_vector(v);
    }
};

template<class Int, class row_vector, class col_vector, class matrix>
class test_fill_matrix
{
public:
    void operator() ()
    {
        matrix A;
        A.resize(200, 200);
        for (unsigned j = 0; j < 150; ++j)
            fill_matrix(A);
    }
    
    void direct()
    {
        matrix A;
        A.resize(200, 200);
        for (unsigned j = 0; j < 150; ++j)
            fill_matrix(A);
    }
};

//////////////////////////////////////////////////////////////////////////////////////////////////
/// Test code: vector add-assignment
//////////////////////////////////////////////////////////////////////////////////////////////////

template<class Int, class row_vector, class col_vector, class matrix>
class test_inplace_add_vector_vector
{
public:
    void operator() ()
    {
        row_vector v, w;
        v.resize(200);
        w.resize(200);
        fill_vector(v);
        fill_vector(w);
        for (unsigned j = 0; j < 30000; ++j)
            v += w;
    }
    
    void direct()
    {
        row_vector v, w;
        v.resize(200);
        w.resize(200);
        fill_vector(v);
        fill_vector(w);
        for (unsigned j = 0; j < 30000; ++j)
            for (unsigned i = 0; i < v.size(); ++i)
                v[i] += w[i];
    }
};

template<class Int, class row_vector, class col_vector, class matrix>
class test_inplace_add_vector_matrixrow
{
public:
    void operator() ()
    {
        row_vector v;
        v.resize(200);
        matrix A;
        A.resize(200, 200);
        fill_vector(v);
        fill_matrix(A);
        for (unsigned j = 0; j < 300; ++j)
            for (unsigned i = 0; i < 100; ++i)
                v += A.row(i);
    }
    
    void direct()
    {
        row_vector v;
        v.resize(200);
        matrix A;
        A.resize(200, 200);
        fill_vector(v);
        fill_matrix(A);
        for (unsigned j = 0; j < 300; ++j)
            for (unsigned i = 0; i < 100; ++i)
                for (unsigned k = 0; k < v.size(); ++k)
                    v[k] += A(i, k);
    }
};

template<class Int, class row_vector, class col_vector, class matrix>
class test_inplace_add_matrixrow_vector
{
public:
    void operator() ()
    {
        matrix A;
        A.resize(200, 200);
        row_vector v;
        v.resize(200);
        fill_matrix(A);
        fill_vector(v);
        for (unsigned j = 0; j < 300; ++j)
            for (unsigned i = 0; i < 100; ++i)
                A.row(i) += v;
    }
    
    void direct()
    {
        matrix A;
        A.resize(200, 200);
        row_vector v;
        v.resize(200);
        fill_matrix(A);
        fill_vector(v);
        for (unsigned j = 0; j < 300; ++j)
            for (unsigned i = 0; i < 100; ++i)
                for (unsigned k = 0; k < v.size(); ++k)
                    A(i, k) += v[k];
    }
};

template<class Int, class row_vector, class col_vector, class matrix>
class test_inplace_add_matrixrow_matrixrow
{
public:
    void operator() ()
    {
        matrix A;
        A.resize(200, 200);
        fill_matrix(A);
        for (unsigned j = 0; j < 300; ++j)
            for (unsigned i = 0; i < 100; ++i)
                A.row(i) += A.row(0);
    }
    
    void direct()
    {
        matrix A;
        A.resize(200, 200);
        fill_matrix(A);
        for (unsigned j = 0; j < 300; ++j)
            for (unsigned i = 0; i < 100; ++i)
                for (unsigned k = 0; k < A.cols(); ++k)
                    A(i, k) += A(0, k);
    }
};

//////////////////////////////////////////////////////////////////////////////////////////////////
/// Test code: vector add-assignment with scalar multiple
//////////////////////////////////////////////////////////////////////////////////////////////////

template<class Int, class row_vector, class col_vector, class matrix>
class test_inplace_add_vector_vector_mul
{
public:
    void operator() ()
    {
        Int m(15);
        row_vector v, w;
        v.resize(200);
        w.resize(200);
        fill_vector(v);
        fill_vector(w);
        for (unsigned j = 0; j < 30000; ++j)
            v += m * w;
    }
    
    void direct()
    {
        Int m(15);
        row_vector v, w;
        v.resize(200);
        w.resize(200);
        fill_vector(v);
        fill_vector(w);
        for (unsigned j = 0; j < 30000; ++j)
            for (unsigned k = 0; k < v.size(); ++k)
                v[k] += m * w[k];
    }
};

template<class Int, class row_vector, class col_vector, class matrix>
class test_inplace_add_vector_matrixrow_mul
{
  public:
    void operator() ()
    {
        Int m(15);
        row_vector v;
        v.resize(200);
        matrix A;
        A.resize(200, 200);
        fill_vector(v);
        fill_matrix(A);
        for (unsigned j = 0; j < 300; ++j)
            for (unsigned i = 0; i < 100; ++i)
                v += m * A.row(i);
    }
    
    void direct()
    {
        Int m(15);
        row_vector v;
        v.resize(200);
        matrix A;
        A.resize(200, 200);
        fill_vector(v);
        fill_matrix(A);
        for (unsigned j = 0; j < 300; ++j)
            for (unsigned i = 0; i < 100; ++i)
                for (unsigned k = 0; k < v.size(); ++k)
                    v[k] += m * A(i, k);
    }
};

template<class Int, class row_vector, class col_vector, class matrix>
class test_inplace_add_matrixrow_vector_mul
{
public:
    void operator() ()
    {
        Int m(15);
        matrix A;
        A.resize(200, 200);
        row_vector v;
        v.resize(200);
        fill_matrix(A);
        fill_vector(v);
        for (unsigned j = 0; j < 300; ++j)
            for (unsigned i = 0; i < 100; ++i)
                A.row(i) += m * v;
    }
    
    void direct()
    {
        Int m(15);
        matrix A;
        A.resize(200, 200);
        row_vector v;
        v.resize(200);
        fill_matrix(A);
        fill_vector(v);
        for (unsigned j = 0; j < 300; ++j)
            for (unsigned i = 0; i < 100; ++i)
                for (unsigned k = 0; k < v.size(); ++k)
                    A(i, k) += m * v[k];
    }
};

template<class Int, class row_vector, class col_vector, class matrix>
class test_inplace_add_matrixrow_matrixrow_mul
{
public:
    void operator() ()
    {
        Int m(15);
        matrix A;
        A.resize(200, 200);
        fill_matrix(A);
        for (unsigned j = 0; j < 300; ++j)
            for (unsigned i = 0; i < 100; ++i)
                A.row(i) += m * A.row(0);
    }
    
    void direct()
    {
        Int m(15);
        matrix A;
        A.resize(200, 200);
        fill_matrix(A);
        for (unsigned j = 0; j < 300; ++j)
            for (unsigned i = 0; i < 100; ++i)
                for (unsigned k = 0; k < A.cols(); ++k)
                    A(i, k) += m * A(0, k);
    }
};

//////////////////////////////////////////////////////////////////////////////////////////////////
/// Test code: dot products
//////////////////////////////////////////////////////////////////////////////////////////////////

template<class Int, class row_vector, class col_vector, class matrix>
class test_dot_vector
{
public:
    void operator() ()
    {
        Int x;
        row_vector v, w;
        v.resize(200);
        w.resize(200);
        fill_vector(v);
        fill_vector(w);
        for (unsigned j = 0; j < 30000; ++j)
            x = dot(v, w);
        touch(x);
    }
    
    void direct()
    {
        Int x;
        row_vector v, w;
        v.resize(200);
        w.resize(200);
        fill_vector(v);
        fill_vector(w);
        for (unsigned j = 0; j < 30000; ++j)
        {
            x = v[0] * w[0];
            for (unsigned k = 1; k < v.size(); ++k)
                x += v[k] * w[k];
        }
        touch(x);
    }
};

template<class Int, class row_vector, class col_vector, class matrix>
class test_dot_matrixrow
{
public:
    void operator() ()
    {
        Int x;
        matrix A;
        A.resize(200, 200);
        fill_matrix(A);
        for (unsigned j = 0; j < 300; ++j)
            for (unsigned i = 0; i < 100; ++i)
                x = dot(A.row(0), A.row(i));
        touch(x);
    }
    
    void direct()
    {
        Int x;
        matrix A;
        A.resize(200, 200);
        fill_matrix(A);
        for (unsigned j = 0; j < 300; ++j)
            for (unsigned i = 0; i < 100; ++i)
            {
                x = A(0, 0) * A(i, 0);
                for (unsigned k = 1; k < A.cols(); ++k)
                    x += A(0, k) * A(i, k);
            }
        touch(x);
    }
};

template<class Int, class row_vector, class col_vector, class matrix>
class test_dot_vector_TO
{
public:
    void operator() ()
    {
        Int x;
        row_vector v, w;
        v.resize(200);
        w.resize(200);
        fill_vector(v);
        fill_vector(w);
        for (unsigned j = 0; j < 30000; ++j)
            dot(x, v, w);
        touch(x);
    }
    
    void direct()
    {
        Int x;
        row_vector v, w;
        v.resize(200);
        w.resize(200);
        fill_vector(v);
        fill_vector(w);
        for (unsigned j = 0; j < 30000; ++j)
        {
            x = v[0] * w[0];
            for (unsigned k = 1; k < v.size(); ++k)
                x += v[k] * w[k];
        }
        touch(x);
    }
};

template<class Int, class row_vector, class col_vector, class matrix>
class test_dot_matrixrow_TO
{
public:
    void operator() ()
    {
        Int x;
        matrix A;
        A.resize(200, 200);
        fill_matrix(A);
        for (unsigned j = 0; j < 300; ++j)
            for (unsigned i = 0; i < 100; ++i)
                dot(x, A.row(0), A.row(i));
        touch(x);
    }
    
    void direct()
    {
        Int x;
        matrix A;
        A.resize(200, 200);
        fill_matrix(A);
        for (unsigned j = 0; j < 300; ++j)
            for (unsigned i = 0; i < 100; ++i)
            {
                x = A(0, 0) * A(i, 0);
                for (unsigned k = 1; k < A.cols(); ++k)
                    x += A(0, k) * A(i, k);
            }
        touch(x);
    }
};

template<class Int, class row_vector, class col_vector, class matrix>
class test_normSq_vector
{
public:
    void operator() ()
    {
        Int x;
        row_vector v;
        v.resize(200);
        fill_vector(v);
        for (unsigned j = 0; j < 30000; ++j)
            x = normSq(v);
        touch(x);
    }
    
    void direct()
    {
        Int x;
        row_vector v;
        v.resize(200);
        fill_vector(v);
        for (unsigned j = 0; j < 30000; ++j)
        {
            x = v[0] * v[0];
            for (unsigned k = 1; k < v.size(); ++k)
                x += v[k] * v[k];
        }
        touch(x);
    }
};

template<class Int, class row_vector, class col_vector, class matrix>
class test_normSq_matrixrow
{
public:
    void operator() ()
    {
        Int x;
        matrix A;
        A.resize(200, 200);
        fill_matrix(A);
        for (unsigned j = 0; j < 300; ++j)
            for (unsigned i = 0; i < 100; ++i)
                x = normSq(A.row(i));
        touch(x);
    }
    
    void direct()
    {
        Int x;
        matrix A;
        A.resize(200, 200);
        fill_matrix(A);
        for (unsigned j = 0; j < 300; ++j)
            for (unsigned i = 0; i < 100; ++i)
            {
                x = A(i, 0) * A(i, 0);
                for (unsigned k = 1; k < A.cols(); ++k)
                    x += A(i, k) * A(i, k);
            }
        touch(x);
    }
};

template<class Int, class row_vector, class col_vector, class matrix>
class test_normSq_vector_TO
{
public:
    void operator() ()
    {
        Int x;
        row_vector v;
        v.resize(200);
        fill_vector(v);
        for (unsigned j = 0; j < 30000; ++j)
            normSq(x, v);
        touch(x);
    }
    
    void direct()
    {
        Int x;
        row_vector v;
        v.resize(200);
        fill_vector(v);
        for (unsigned j = 0; j < 30000; ++j)
        {
            x = v[0] * v[0];
            for (unsigned k = 1; k < v.size(); ++k)
                x += v[k] * v[k];
        }
        touch(x);
    }
};

template<class Int, class row_vector, class col_vector, class matrix>
class test_normSq_matrixrow_TO
{
public:
    void operator() ()
    {
        Int x;
        matrix A;
        A.resize(200, 200);
        fill_matrix(A);
        for (unsigned j = 0; j < 300; ++j)
            for (unsigned i = 0; i < 100; ++i)
                normSq(x, A.row(i));
        touch(x);
    }
    
    void direct()
    {
        Int x;
        matrix A;
        A.resize(200, 200);
        fill_matrix(A);
        for (unsigned j = 0; j < 300; ++j)
            for (unsigned i = 0; i < 100; ++i)
            {
                x = A(i, 0) * A(i, 0);
                for (unsigned k = 1; k < A.cols(); ++k)
                    x += A(i, k) * A(i, k);
            }
        touch(x);
    }
};
//////////////////////////////////////////////////////////////////////////////////////////////////
/// Test code: matrix-vector and matrix-matrix multiplication
//////////////////////////////////////////////////////////////////////////////////////////////////

template<class Int, class row_vector, class col_vector, class matrix>
class test_matrix_vector_multiplication
{
public:
    void operator() ()
    {
        col_vector v, w;
        v.resize(200);
        w.resize(200);
        matrix A;
        A.resize(200, 200);
        fill_vector(v);
        fill_matrix(A);
        for (unsigned j = 0; j < 500; ++j)
            v = A * w;
    }
    
    void direct()
    {
        col_vector v, w;
        v.resize(200);
        w.resize(200);
        matrix A;
        A.resize(200, 200);
        fill_vector(v);
        fill_matrix(A);
        for (unsigned j = 0; j < 500; ++j)
            for (unsigned i = 0; i < A.rows(); ++i)
            {
                v[i] = A(i, 0) * w[0];
                for (unsigned k = 1; k < A.cols(); ++k)
                    v[i] += A(i, k) * w[k];
            }
    }
};

template<class Int, class row_vector, class col_vector, class matrix>
class test_vector_matrix_multiplication
{
public:
    void operator() ()
    {
        row_vector v, w;
        v.resize(200);
        w.resize(200);
        matrix A;
        A.resize(200, 200);
        fill_vector(v);
        fill_matrix(A);
        for (unsigned j = 0; j < 500; ++j)
            v = w * A;
    }
    
    void direct()
    {
        row_vector v, w;
        v.resize(200);
        w.resize(200);
        matrix A;
        A.resize(200, 200);
        fill_vector(v);
        fill_matrix(A);
        for (unsigned j = 0; j < 500; ++j)
            for (unsigned i = 0; i < A.rows(); ++i)
            {
                v[i] = w[0] * A(0, i);
                for (unsigned k = 1; k < A.cols(); ++k)
                    v[i] += w[k] * A(k, i);
            }
    }
};

template<class Int, class row_vector, class col_vector, class matrix>
class test_matrix_matrixcol_multiplication
{
public:
    void operator() ()
    {
        matrix A, B;
        A.resize(200, 200);
        B.resize(200, 200);
        fill_matrix(A);
        fill_matrix(B);
        for (unsigned j = 0; j < 200; ++j)
            B.col(0) = A * B.col(j);
    }
    
    void direct()
    {
        matrix A, B;
        A.resize(200, 200);
        B.resize(200, 200);
        fill_matrix(A);
        fill_matrix(B);
        for (unsigned j = 0; j < 200; ++j)
            for (unsigned i = 0; i < A.rows(); ++i)
            {
                B(i, 0) = A(i, 0) * B(0, j);
                for (unsigned k = 1; k < A.cols(); ++k)
                    B(i, 0) += A(i, k) * B(k, j);
            }
    }
};

template<class Int, class row_vector, class col_vector, class matrix>
class test_matrixrow_matrix_multiplication
{
public:
    void operator() ()
    {
        matrix A, B;
        A.resize(200, 200);
        B.resize(200, 200);
        fill_matrix(A);
        fill_matrix(B);
        for (unsigned j = 0; j < 200; ++j)
            B.row(0) = B.row(j) * A;
    }
    
    void direct()
    {
        matrix A, B;
        A.resize(200, 200);
        B.resize(200, 200);
        fill_matrix(A);
        fill_matrix(B);
        for (unsigned j = 0; j < 200; ++j)
            for (unsigned i = 0; i < A.cols(); ++i)
            {
                B(0, i) = B(j, 0) * A(0, i);
                for (unsigned k = 1; k < A.rows(); ++k)
                    B(0, i) += B(j, k) * A(k, i);
            }
    }
};

template<class Int, class row_vector, class col_vector, class matrix>
class test_matrix_matrix_multiplication
{
public:
    void operator() ()
    {
        matrix A, B, C;
        A.resize(200, 200);
        B.resize(200, 200);
        C.resize(200, 200);
        fill_matrix(A);
        fill_matrix(B);
        for (unsigned ell = 0; ell < 1; ++ell)
            C = A * B;
    }
    
    void direct()
    {
        matrix A, B, C;
        A.resize(200, 200);
        B.resize(200, 200);
        C.resize(200, 200);
        fill_matrix(A);
        fill_matrix(B);
        for (unsigned ell = 0; ell < 1; ++ell)
            for (unsigned i = 0; i < A.rows(); ++i)
                for (unsigned j = 0; j < B.cols(); ++j)
                {
                    C(i, j) = A(i, 0) * B(0, j);
                    for (unsigned k = 1; k < A.cols(); ++k)
                        C(i, j) += A(i, k) * B(k, j);
                }
    }
};

//////////////////////////////////////////////////////////////////////////////////////////////////
/// The test management
//////////////////////////////////////////////////////////////////////////////////////////////////

unsigned repetitions = 1;
unsigned tries = 4;

enum { field_width = 70 };

template<template<typename, typename, typename, typename> class F, typename Int>
void run_test(const std::string & title, int ignore)
{
    // Create strings for timings
    std::vector<char> title_old, title_new, title_oldd, title_newd;
    std::string title_aligned = title;
    while (title_aligned.size() < field_width)
        title_aligned += ' ';
    std::ostringstream ss;
    ss << title << ": old";
    std::string s = ss.str();
    while (s.size() < field_width)
        s += ' ';
    title_old.resize(s.size() + 1, 0);
    std::copy(s.begin(), s.end(), title_old.begin());
    ss << " direct";
    s = ss.str();
    while (s.size() < field_width)
        s += ' ';
    title_oldd.resize(s.size() + 1, 0);
    std::copy(s.begin(), s.end(), title_oldd.begin());
    ss.str("");
    ss << title << ": new";
    s = ss.str();
    while (s.size() < field_width)
        s += ' ';
    title_new.resize(s.size() + 1, 0);
    std::copy(s.begin(), s.end(), title_new.begin());
    ss << " direct";
    s = ss.str();
    while (s.size() < field_width)
        s += ' ';
    title_newd.resize(s.size() + 1, 0);
    std::copy(s.begin(), s.end(), title_newd.begin());
    TimedDataCollector td_old(&title_old.front());
    TimedDataCollector td_oldd(&title_oldd.front());
    TimedDataCollector td_new(&title_new.front());
    TimedDataCollector td_newd(&title_newd.front());
    // Instantiate tests
    F<Int, linalgold::math_vector<Int>, linalgold::math_vector<Int>, linalgold::math_matrix<Int> > old_test;
    F<Int, linalg::math_rowvector<Int>, linalg::math_colvector<Int>, linalg::math_matrix<Int> > new_test;
    // Do tests
    for (unsigned r = 0; r < repetitions; ++r)
    {
        old_test();
        {
            Profiler<TimedDataCollector> p(td_old);
            for (unsigned t = 0; t < tries; ++t)
            {
                if (t) p.restart();
                old_test();
            }
        }
        new_test();
        {
            Profiler<TimedDataCollector> p(td_new);
            for (unsigned t = 0; t < tries; ++t)
            {
                if (t) p.restart();
                new_test();
            }
        }
        old_test.direct();
        {
            Profiler<TimedDataCollector> p(td_oldd);
            for (unsigned t = 0; t < tries; ++t)
            {
                if (t) p.restart();
                old_test.direct();
            }
        }
        new_test.direct();
        {
            Profiler<TimedDataCollector> p(td_newd);
            for (unsigned t = 0; t < tries; ++t)
            {
                if (t) p.restart();
                new_test.direct();
            }
        }
    }
    std::cout << title_aligned << ": "
              << (long double)td_new.getTimer().elapsed() / (long double)td_old.getTimer().elapsed() * 100.0 << "% (new-old), "
              << (long double)td_new.getTimer().elapsed() / (long double)td_newd.getTimer().elapsed() * 100.0 << "% (new-direct), "
              << (long double)td_old.getTimer().elapsed() / (long double)td_oldd.getTimer().elapsed() * 100.0 << "% (old-direct) "
              << "\n";
}

int waste_time()
{
    int i = 0;
    enum { N = 500 };
    for (unsigned j = 0; j < N; ++j)
    {
        i += 1;
        for (unsigned k = 0; k < N; ++k)
        {
            i *= 5;
            for (unsigned ell = 0; ell < N; ++ell)
            {
                i *= 3;
                i %= 123127;
            }
        }
    }
    return i;
}

template<typename Int>
void do_tests(int res)
{
    run_test<test_matrix_vector_multiplication, Int>("matrix-vector multiplication", res);
    run_test<test_vector_matrix_multiplication, Int>("vector-matrix multiplication", res);
    run_test<test_matrix_matrixcol_multiplication, Int>("matrix-vector (vector == matrix column) multiplication", res);
    run_test<test_matrixrow_matrix_multiplication, Int>("vector-matrix (vector == matrix row) multiplication", res);
    run_test<test_matrix_matrix_multiplication, Int>("matrix-matrix multiplication", res);
    run_test<test_fill_vector, Int>("vector filling", res);
    run_test<test_fill_matrix, Int>("matrix filling", res);
    run_test<test_inplace_add_vector_vector, Int>("operator += (vector, vector)", res);
    run_test<test_inplace_add_vector_matrixrow, Int>("operator += (vector, matrix.row())", res);
    run_test<test_inplace_add_matrixrow_vector, Int>("operator += (matrix.row(), vector)", res);
    run_test<test_inplace_add_matrixrow_matrixrow, Int>("operator += (matrix.row(), matrix.row())", res);
    run_test<test_inplace_add_vector_vector_mul, Int>("operator += (vector, vector * scalar)", res);
    run_test<test_inplace_add_vector_matrixrow_mul, Int>("operator += (vector, matrix.row() * scalar)", res);
    run_test<test_inplace_add_matrixrow_vector_mul, Int>("operator += (matrix.row(), vector * scalar)", res);
    run_test<test_inplace_add_matrixrow_matrixrow_mul, Int>("operator += (matrix.row(), matrix.row() * scalar)", res);
    run_test<test_dot_vector, Int>("dot product (vectors)", res);
    run_test<test_dot_matrixrow, Int>("dot product (matrix rows)", res);
    run_test<test_dot_vector_TO, Int>("dot product (vectors), three-operand variant", res);
    run_test<test_dot_matrixrow_TO, Int>("dot product (matrix rows), three-operand variant", res);
    run_test<test_normSq_vector, Int>("squared norm (vector)", res);
    run_test<test_normSq_matrixrow, Int>("squared norm (matrix row)", res);
    run_test<test_normSq_vector_TO, Int>("squared norm (vector), two-operand variant", res);
    run_test<test_normSq_matrixrow_TO, Int>("squared norm (matrix row), two-operand variant", res);
}

int main(int argc, char **argv)
{
    arithmetic::initArithmeticThreadAllocators();
    
    enum Type { T_Integer, T_int, T_long, T_longlong, T_double };
    
    helper::ArgumentParser args(argc, argv);
    const helper::ArgumentParser::Value * v;
    if ((v = args.getValue("repetitions")) == NULL)
        v = args.getValue("r");
    if (v)
    {
        if (v->isInteger())
            if (v->getInteger() > 0)
            {
                repetitions = v->getInteger();
                std::cout << "Doing " << repetitions << " repetitions.\n";
            }
    }

    if ((v = args.getValue("tries")) == NULL)
        v = args.getValue("tr");
    if (v)
    {
        if (v->isInteger())
            if (v->getInteger() > 0)
            {
                tries = v->getInteger();
                std::cout << "Doing " << tries << " tries (per repetition).\n";
            }
    }

    Type type = T_Integer;
    if ((v = args.getValue("type")) == NULL)
        v = args.getValue("t");
    if (v)
    {
        if ((v->getText() == "bi") || (v->getText() == "bigint"))
            type = T_Integer;
        if ((v->getText() == "i") || (v->getText() == "int"))
            type = T_int;
        if ((v->getText() == "l") || (v->getText() == "long"))
            type = T_long;
        if ((v->getText() == "ll") || (v->getText() == "longlong"))
            type = T_longlong;
        if ((v->getText() == "d") || (v->getText() == "double"))
            type = T_double;
    }
    
    CheckTime ct;
    int res = waste_time();
    std::cout << "WANT: < 100%\n";
    switch (type)
    {
    case T_Integer:  std::cout << "Using arithmetic::Integer\n"; do_tests<arithmetic::Integer>(res); break;
    case T_int:      std::cout << "Using int\n";                 do_tests<int>(res);                 break;
    case T_long:     std::cout << "Using long\n";                do_tests<long>(res);                break;
    case T_longlong: std::cout << "Using long long\n";           do_tests<long long>(res);           break;
    case T_double:   std::cout << "Using double\n";              do_tests<double>(res);              break;
    }
}

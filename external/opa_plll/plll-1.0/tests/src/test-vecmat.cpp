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

void setZero(int & x)
{
    x = 0;
}

void setZero(double & x)
{
    x = 0;
}

void setZero(long & x)
{
    x = 0;
}

int square(int x)
{
    return x * x;
}

#include <plll/matrix.hpp>

using namespace plll;

void compileTest()
// To see whether all constructs compile
{
    class TestClass
    {
        void testBaseVector()
        {
            linalg::base_rowvector<int> v;
            linalg::base_rowvector<int> w(1);
            linalg::base_rowvector<int> x(2, linalg::Initialize(0));
            linalg::base_rowvector<int> y(x);
            linalg::base_rowvector<int> z;
            linalg::base_rowvector<long> q = linalg::base_rowvector<long>(z);
            assign(v, w);
            linalg::swap(v, w);
            v.swap(w);
            v = w;
            v = q;
            q = v;
            assert(v.size() == 0);
            v.resize(2);
            assert(v.data()[0] == 0);
            v[1] = x[0];
            v(0, 0) = x(1, 0) + 1;
            assert(v != w);
            assert(v == v);
        }

        void testMathVector()
        {
            linalg::math_rowvector<int> v;
            linalg::math_rowvector<int> w(1);
            linalg::math_rowvector<int> x(2, linalg::Initialize(0));
            linalg::math_rowvector<int> y(x);
            linalg::math_rowvector<int> z;
            v = w;
            w %= 2;
            w /= 3;
            w *= 4;
            z = v;
            v += z;
            v -= z;
            v = -z;
            assert(v < x);
            assert(v <= x);
            assert(!(v > x));
            assert(!(v >= x));
            assert(v != x);
            assert(v == v);
            z = v;
            z = v % 1;
            z = v / 2;
            z = v * 3;
            z = 4 * v;
            z = x + y;
            z = x - y;
            assert(normSq(z) == 0);
            z = v;
            mod(z, v, 2);
            div(z, v, 3);
            mul(z, v, 4);
            z = x;
            add(z, x, y);
            sub(z, x, y);
            addmul(z, x, 2);
            submul(z, y, 3);
            neg(x, y);
            int t;
            normSq(t, x);
            t = dot(x, y);
            dot(t, x, y);
        }

        void testBaseMatrix()
        {
            linalg::base_matrix<int> A;
            linalg::base_matrix<int> B(1, 2);
            linalg::base_matrix<int> C(3, 4, linalg::Initialize(5));
            linalg::base_matrix<int> D(B);
            linalg::base_matrix<long> E(D); // explicit?
            assign(A, B);
            linalg::swap(A, C);
            A.swap(C);
            B = D;
            E = C;
            assert(A.rows() == 1);
            assert(B.cols() == 2);
            D.resize(5, 6);
            assert(A(0, 1) == 0);
            assert(A.data()[0] == 0);
            assert(A.data(0) == 0);
            assert(E == C);
            assert(A == B);
            assert(E != A);
            C = A.transpose();
        }

        void testMathMatrix()
        {
            linalg::math_matrix<int> A;
            linalg::math_matrix<int> B(2, 3);
            linalg::math_matrix<int> C(4, 5, linalg::Initialize(6));
            linalg::math_matrix<int> D(B);
            linalg::math_matrix<long> E(C);
            linalg::math_matrix<double> F;
            A = D;
            F = A;
            A %= 5;
            A /= 6;
            A *= 7;
            A += F;
            A -= D;
            B = -C;
            assert(A < B);
            assert(A <= B);
            assert(!(A > B));
            assert(!(A >= B));
            assert(A != B);
            assert(A == A);
            A = B % 2;
            A = C / 3;
            A = D * 4;
            A = 5 * D;
            C = C + E;
            F = E;
            F = F - E;
            A = B.transpose();
            C.resize(B.cols(), 5);
            B *= C;
            F.resize(B.cols(), B.cols());
            B *= F;
            A = B;
            mod(A, B, 1);
            A = C;
            div(A, C, 2);
            D.resize(A.cols(), A.cols());
            mul(A, D, 3);
            C.resize(B.rows(), B.cols());
            A.resize(B.rows(), B.cols());
            add(A, B, C);
            sub(A, B, C);
            neg(A, B);
            A.resize(B.cols(), B.rows());
            transpose(A, B);
            E = A;
            F = A;
            mod(A, E, 1);
            div(A, E, 2);
            mul(A, E, 3);
            add(A, E, F);
            sub(A, E, F);
            neg(A, E);
            A.resize(F.cols(), F.rows());
            transpose(A, F);
            C.resize(B.cols(), 2);
            A.resize(B.rows(), C.cols());
            mul(A, B, C);
        }

        void testMathMatrixRow()
        {
            // ...
        }

        void testMathMatrixCol()
        {
            // ...
        }

        void testMathMatrixInteraction()
        // Interaction with linalg::math_rowvector, Row, Column, CRow, CColumn
        {
            // ...
        }

    public:
        TestClass()
        {
            testBaseVector();
            testMathVector();
            testBaseMatrix();
            testMathMatrix();
            testMathMatrixRow();
            testMathMatrixCol();
            testMathMatrixInteraction();
        }
    };
    
    TestClass A;
}

int main()
{
    compileTest();
}

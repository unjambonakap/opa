#ifndef PLLL_INCLUDE_GUARD__VECMAT_IMPL_CPP
#define PLLL_INCLUDE_GUARD__VECMAT_IMPL_CPP

namespace plll
{
    namespace linalgold
    {
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Vectors
        
        template<class T, class StorageTraits, class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
        inline void assign(base_vector<T, StorageTraits> & r, const V & v)
        {
            r.assign(v);
        }
        
        template<class T, class StorageTraits>
        inline void swap(base_vector<T, StorageTraits> & A, base_vector<T, StorageTraits> & B)
        {
            A.swap(B);
        }
        
        template<class T, class StorageTraits, class V>
        inline void swap(base_vector<T, StorageTraits> & v1, V & v2)
        {
            v1.swap(v2);
        }
        
        template<class T, class ST, class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
        inline void mod(math_vector<T, ST> & r, const V & a, const T & b)
        {
            assert(a.size() == r.size());
            for (typename base_vector<T, ST>::size_type i = 0; i < a.size(); ++i)
                r[i] = a[i] % b;
        }

        template<class T, class ST, class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
        inline void div(math_vector<T, ST> & r, const V & a, const T & b)
        {
            assert(a.size() == r.size());
            for (typename base_vector<T, ST>::size_type i = 0; i < a.size(); ++i)
                r[i] = a[i] / b;
        }

        template<class T, class ST, class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
        inline void mul(math_vector<T, ST> & r, const V & a, const T & b)
        {
            assert(a.size() == r.size());
            for (typename base_vector<T, ST>::size_type i = 0; i < a.size(); ++i)
                r[i] = a[i] * b;
        }

        template<class T, class ST, class V2, class V3> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
        inline void add(math_vector<T, ST> & r, const V2 & a, const V3 & b)
        {
            assert(a.size() == b.size());
            assert(a.size() == r.size());
            for (typename base_vector<T, ST>::size_type i = 0; i < a.size(); ++i)
                r[i] = a[i] + b[i];
        }
    
        template<class T, class ST, class V2, class V3> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
        inline void sub(math_vector<T, ST> & r, const V2 & a, const V3 & b)
        {
            assert(a.size() == b.size());
            assert(a.size() == r.size());
            for (typename base_vector<T, ST>::size_type i = 0; i < a.size(); ++i)
                r[i] = a[i] - b[i];
        }

        template<class T, class ST, class V2, class TT> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
        inline void addmul(math_vector<T, ST> & r, const V2 & a, const TT & b)
        {
            assert(a.size() == r.size());
            for (typename base_vector<T, ST>::size_type i = 0; i < a.size(); ++i)
                r[i] += a[i] * b;
        }

        template<class T, class ST, class V2, class TT> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
        inline void submul(math_vector<T, ST> & r, const V2 & a, const TT & b)
        {
            assert(a.size() == r.size());
            for (typename base_vector<T, ST>::size_type i = 0; i < a.size(); ++i)
                r[i] -= a[i] * b;
        }

        template<class T, class ST, class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
        inline void neg(math_vector<T, ST> & r, const V & a)
        {
            assert(a.size() == r.size());
            for (typename base_vector<T, ST>::size_type i = 0; i < a.size(); ++i)
                r[i] = -a[i];
        }

        template<class T, class ST>
        inline void normSq(T & r, const math_vector<T, ST> & a)
        {
            setZero(r);
            for (typename base_vector<T, ST>::size_type i = 0; i < a.size(); ++i)
                r += a[i] * a[i];
        }

        template<class T, class ST, class V>
        inline void dot(T & r, const math_vector<T, ST> & a, const V & b)
        {
            assert(a.size() == b.size());
            setZero(r);
            for (typename base_vector<T, ST>::size_type i = 0; i < a.size(); ++i)
                r += a[i] * b[i];
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Matrices

        template<class T, class StorageTraits, class S, class A>
        inline void assign(base_matrix<T, StorageTraits> & r, const base_matrix<S, A> & m)
        {
            r.assign(m);
        }
        
        template<class T, class StorageTraits>
        inline void swap(base_matrix<T, StorageTraits> & A, base_matrix<T, StorageTraits> & B)
        {
            A.swap(B);
        }
        
        template<class T, class ST>
        inline math_matrix<T, ST> operator * (const T & s, const math_matrix<T, ST> & m)
        {
            math_matrix<T, ST> ret(m.rows(), m.cols());
            mul(ret, m, s);
            return ret;
        }

        template<class T, class ST, class T2, class ST2>
        inline math_matrix<T, ST> operator * (const math_matrix<T, ST> & m1, const math_matrix<T2, ST2> & m2)
        {
            math_matrix<T, ST> ret(m1.rows(), m2.cols());
            mul(ret, m1, m2);
            return ret;
        }

        template<class T, class ST, class T2, class ST2>
        inline math_vector<T, ST> operator * (const math_matrix<T, ST> & A, const math_vector<T2, ST2> & v)
        {
            math_vector<T, ST> ret(A.rows());
            mul(ret, A, v);
            return ret;
        }

        template<class T, class ST, class T2, class ST2>
        inline math_vector<T, ST> operator * (const math_matrix<T, ST> & A, const MatrixRow<T2, ST2> & v)
        {
            math_vector<T, ST> ret(A.rows());
            mul(ret, A, v);
            return ret;
        }

        template<class T, class ST, class T2, class ST2>
        inline math_vector<T, ST> operator * (const math_matrix<T, ST> & A, const MatrixCRow<T2, ST2> & v)
        {
            math_vector<T, ST> ret(A.rows());
            mul(ret, A, v);
            return ret;
        }

        template<class T, class ST, class T2, class ST2>
        inline math_vector<T, ST> operator * (const math_matrix<T, ST> & A, const MatrixColumn<T2, ST2> & v)
        {
            math_vector<T, ST> ret(A.rows());
            mul(ret, A, v);
            return ret;
        }

        template<class T, class ST, class T2, class ST2>
        inline math_vector<T, ST> operator * (const math_matrix<T, ST> & A, const MatrixCColumn<T2, ST2> & v)
        {
            math_vector<T, ST> ret(A.rows());
            mul(ret, A, v);
            return ret;
        }

        template<class T, class ST, class T2, class ST2>
        inline math_vector<T, ST> operator * (const math_vector<T2, ST2> & v, const math_matrix<T, ST> & A)
        {
            math_vector<T, ST> ret(A.cols());
            mul(ret, v, A);
            return ret;
        }

        template<class T, class ST, class T2, class ST2>
        inline math_vector<T, ST> operator * (const MatrixRow<T2, ST2> & v, const math_matrix<T, ST> & A)
        {
            math_vector<T, ST> ret(A.cols());
            mul(ret, v, A);
            return ret;
        }

        template<class T, class ST, class T2, class ST2>
        inline math_vector<T, ST> operator * (const MatrixCRow<T2, ST2> & v, const math_matrix<T, ST> & A)
        {
            math_vector<T, ST> ret(A.cols());
            mul(ret, v, A);
            return ret;
        }

        template<class T, class ST, class T2, class ST2>
        inline math_vector<T, ST> operator * (const MatrixColumn<T2, ST2> & v, const math_matrix<T, ST> & A)
        {
            math_vector<T, ST> ret(A.cols());
            mul(ret, v, A);
            return ret;
        }

        template<class T, class ST, class T2, class ST2>
        inline math_vector<T, ST> operator * (const MatrixCColumn<T2, ST2> & v, const math_matrix<T, ST> & A)
        {
            math_vector<T, ST> ret(A.cols());
            mul(ret, v, A);
            return ret;
        }

// Temporary-less operations

        template<class T, class ST, class T2, class ST2>
        inline void mod(math_matrix<T, ST> & r, const math_matrix<T2, ST2> & a, const T & b)
        {
            assert(a.rows() == r.rows());
            assert(a.cols() == r.cols());
            typename base_matrix<T, ST>::size_type c = r.rows() * r.cols();
            for (typename base_matrix<T, ST>::size_type i = 0; i < c; ++i)
                r.d_data[i] = a.d_data[i] % b;
        }

        template<class T, class ST, class T2, class ST2>
        inline void div(math_matrix<T, ST> & r, const math_matrix<T2, ST2> & a, const T & b)
        {
            assert(a.rows() == r.rows());
            assert(a.cols() == r.cols());
            typename base_matrix<T, ST>::size_type c = r.rows() * r.cols();
            for (typename base_matrix<T, ST>::size_type i = 0; i < c; ++i)
                r.d_data[i] = a.d_data[i] / b;
        }

        template<class T, class ST, class T2, class ST2>
        inline void mul(math_matrix<T, ST> & r, const math_matrix<T2, ST2> & a, const T & b)
        {
            assert(a.rows() == r.rows());
            assert(a.cols() == r.cols());
            typename base_matrix<T, ST>::size_type c = r.rows() * r.cols();
            for (typename base_matrix<T, ST>::size_type i = 0; i < c; ++i)
                r.d_data[i] = a.d_data[i] * b;
        }

        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        inline void add(math_matrix<T, ST> & r, const math_matrix<T2, ST2> & a, const math_matrix<T3, ST3> & b)
        {
            assert(a.rows() == b.rows());
            assert(a.cols() == b.cols());
            assert(a.rows() == r.rows());
            assert(a.cols() == r.cols());
            typename base_matrix<T, ST>::size_type c = r.rows() * r.cols();
            for (typename base_matrix<T, ST>::size_type i = 0; i < c; ++i)
                r.d_data[i] = a.d_data[i] + b.d_data[i];
        }

        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        inline void sub(math_matrix<T, ST> & r, const math_matrix<T2, ST2> & a, const math_matrix<T3, ST3> & b)
        {
            assert(a.rows() == b.rows());
            assert(a.cols() == b.cols());
            assert(a.rows() == r.rows());
            assert(a.cols() == r.cols());
            typename base_matrix<T, ST>::size_type c = r.rows() * r.cols();
            for (typename base_matrix<T, ST>::size_type i = 0; i < c; ++i)
                r.d_data[i] = a.d_data[i] - b.d_data[i];
        }

        template<class T, class ST, class T2, class ST2>
        inline void neg(math_matrix<T, ST> & r, const math_matrix<T2, ST2> & a)
        {
            assert(a.rows() == r.rows());
            assert(a.cols() == r.cols());
            typename base_matrix<T, ST>::size_type c = r.rows() * r.cols();
            for (typename base_matrix<T, ST>::size_type i = 0; i < c; ++i)
                r.d_data[i] = -a.d_data[i];
        }

        template<class T, class ST, class T2, class ST2>
        inline void transpose(math_matrix<T, ST> & r, const math_matrix<T2, ST2> & a)
        {
            assert(a.rows() == r.cols());
            assert(a.cols() == r.rows());
            typename base_matrix<T, ST>::size_type idx = 0;
            for (unsigned i = 0; i < a.cols(); ++i)
            {
                typename base_matrix<T, ST>::size_type idx2 = i;
                for (unsigned j = 0; j < a.rows(); ++j, ++idx, idx2 += a.cols())
                    r.d_data[idx] = a.d_data[idx2];
            }
        }

        template<class T, class ST>
        inline void transpose(math_matrix<T, ST> & r, const math_matrix<T, ST> & a)
        {
            assert(a.rows() == r.cols());
            assert(a.cols() == r.rows());
            if (&r == &a)
            {
                using std::swap;
                // Do transpose in-place
                for (unsigned i = 0; i < r.cols(); ++i)
                    for (unsigned j = 0; j < i; ++j)
                        swap(r(i, j), r(j, i));
            }
            else
            {
                typename base_matrix<T, ST>::size_type idx = 0;
                for (unsigned i = 0; i < a.cols(); ++i)
                {
                    typename base_matrix<T, ST>::size_type idx2 = i;
                    for (unsigned j = 0; j < a.rows(); ++j, ++idx, idx2 += a.cols())
                        r.d_data[idx] = a.d_data[idx2];
                }
            }
        }

// Temporary-less matrix-vector and matrix-matrix multiplications

        template<class T, class ST, class T2, class ST2, class V>
        inline void mul2(math_vector<T, ST> & r, const math_matrix<T2, ST2> & a, const V & b)
        // Promise: &r is not equal to &b
        {
            assert(a.rows() == r.size());
            assert(a.cols() == b.size());
            typename base_matrix<T, ST>::size_type idx = 0;
            for (typename base_matrix<T, ST>::size_type i = 0; i < a.rows(); ++i)
            {
                r[i] = a.d_data[idx++] * b[0];
                for (typename base_matrix<T, ST>::size_type j = 1; j < a.cols(); ++j, ++idx)
                    r[i] += a.d_data[idx] * b[j];
            }
        }

        template<class T, class ST, class T2, class ST2, class V>
        inline void mul2(math_vector<T, ST> & r, const V & a, const math_matrix<T2, ST2> & b)
        // Promise: &r is not equal to &a
        {
            assert(b.rows() == a.size());
            assert(b.cols() == r.size());
            for (typename base_matrix<T, ST>::size_type i = 0; i < b.cols(); ++i)
            {
                typename base_matrix<T, ST>::size_type idx = i;
                r[i] = a[0] * b.d_data[idx];
                idx += b.cols();
                for (typename base_matrix<T, ST>::size_type j = 1; j < b.rows(); ++j, idx += b.cols())
                    r[i] += a[j] * b.d_data[idx];
            }
        }

        template<class T, class ST, class T2, class ST2, class V>
        inline void mul2(const MatrixRow<T2, ST2> & r, const math_matrix<T, ST> & a, const V & b)
        // Promise: &r is not equal to &b or to a row/column of a
        {
            assert(a.rows() == r.size());
            assert(a.cols() == b.size());
            typename base_matrix<T, ST>::size_type idx = 0;
            for (typename base_matrix<T, ST>::size_type i = 0; i < a.rows(); ++i)
            {
                r[i] = a.d_data[idx++] * b[0];
                for (typename base_matrix<T, ST>::size_type j = 1; j < a.cols(); ++j, ++idx)
                    r[i] += a.d_data[idx] * b[j];
            }
        }

        template<class T, class ST, class T2, class ST2, class V>
        inline void mul2(const MatrixRow<T2, ST2> & r, const V & a, const math_matrix<T, ST> & b)
        // Promise: &r is not equal to &a or to a row/column of b
        {
            assert(b.rows() == a.size());
            assert(b.cols() == r.size());
            for (typename base_matrix<T, ST>::size_type i = 0; i < b.cols(); ++i)
            {
                typename base_matrix<T, ST>::size_type idx = i;
                r[i] = a[0] * b.d_data[idx];
                idx += b.cols();
                for (typename base_matrix<T, ST>::size_type j = 1; j < b.rows(); ++j, idx += b.cols())
                    r[i] += a[j] * b.d_data[idx];
            }
        }

        template<class T, class ST, class T2, class ST2, class V>
        inline void mul2(const MatrixColumn<T2, ST2> & r, const math_matrix<T, ST> & a, const V & b)
        // Promise: &r is not equal to &b or to a row/column of a
        {
            assert(a.rows() == r.size());
            assert(a.cols() == b.size());
            typename base_matrix<T, ST>::size_type idx = 0;
            for (typename base_matrix<T, ST>::size_type i = 0; i < a.rows(); ++i)
            {
                r[i] = a.d_data[idx++] * b[0];
                for (typename base_matrix<T, ST>::size_type j = 1; j < a.cols(); ++j, ++idx)
                    r[i] += a.d_data[idx] * b[j];
            }
        }

        template<class T, class ST, class T2, class ST2, class V>
        inline void mul2(const MatrixColumn<T2, ST2> & r, const V & a, const math_matrix<T, ST> & b)
        // Promise: &r is not equal to &a or to a row/column of b
        {
            assert(b.rows() == a.size());
            assert(b.cols() == r.size());
            for (typename base_matrix<T, ST>::size_type i = 0; i < b.cols(); ++i)
            {
                typename base_matrix<T, ST>::size_type idx = i;
                r[i] = a[0] * b.d_data[idx];
                idx += b.cols();
                for (typename base_matrix<T, ST>::size_type j = 1; j < b.rows(); ++j, idx += b.cols())
                    r[i] += a[j] * b.d_data[idx];
            }
        }

        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        inline void mul2(math_matrix<T, ST> & r, const math_matrix<T2, ST2> & a, const base_matrix<T3, ST3> & b)
        // Promise: &r is neither &a nor &b
        {
            assert(a.rows() == r.rows());
            assert(a.cols() == b.rows());
            assert(b.cols() == r.cols());
            typename base_matrix<T, ST>::size_type oi = 0, ii1 = 0;
            for (typename base_matrix<T, ST>::size_type i = 0; i < a.rows(); ++i, ii1 += a.cols())
                for (typename base_matrix<T, ST>::size_type j = 0; j < b.cols(); ++j, ++oi)
                {
                    typename base_matrix<T, ST>::size_type ii2 = j;
                    r.d_data[oi] = a.d_data[ii1] * b.data(ii2);
                    ii2 += b.cols();
                    for (typename base_matrix<T, ST>::size_type k = 1; k < a.cols(); ++k, ii2 += b.cols())
                        r.d_data[oi] += a.d_data[ii1 + k] * b.data(ii2);
                }
        }

// First argument: math_vector<>; matrix times vector
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        inline void mul(math_vector<T, ST> & r, const math_matrix<T3, ST3> & a, const math_vector<T2, ST2> & b)
        {
            if ((void *)&r == (void *)&b)
            {
                using std::swap;
                // Use temporary!
                math_vector<T, ST> tmp(r.size());
                mul2(tmp, a, b);
                swap(r, tmp);
            }
            else
                // Do it directly
                mul2(r, a, b);
        }

        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        inline void mul(math_vector<T, ST> & r, const math_matrix<T3, ST3> & a, const MatrixRow<T2, ST2> & b)
        { mul2(r, a, b); }
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        inline void mul(math_vector<T, ST> & r, const math_matrix<T3, ST3> & a, const MatrixCRow<T2, ST2> & b)
        { mul2(r, a, b); }
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        inline void mul(math_vector<T, ST> & r, const math_matrix<T3, ST3> & a, const MatrixColumn<T2, ST2> & b)
        { mul2(r, a, b); }
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        inline void mul(math_vector<T, ST> & r, const math_matrix<T3, ST3> & a, const MatrixCColumn<T2, ST2> & b)
        { mul2(r, a, b); }

// First argument: math_vector<>; vector times matrix
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        inline void mul(math_vector<T, ST> & r, const math_vector<T2, ST2> & a, const math_matrix<T3, ST3> & b)
        {
            if ((void *)&r == (void *)&a)
            {
                using std::swap;
                // Use temporary!
                math_vector<T, ST> tmp(r.size());
                mul2(tmp, a, b);
                swap(r, tmp);
            }
            else
                // Do it directly
                mul2(r, a, b);
        }

        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        inline void mul(math_vector<T, ST> & r, const MatrixRow<T2, ST2> & a, const math_matrix<T3, ST3> & b)
        { mul2(r, a, b); }
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        inline void mul(math_vector<T, ST> & r, const MatrixCRow<T2, ST2> & a, const math_matrix<T3, ST3> & b)
        { mul2(r, a, b); }
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        inline void mul(math_vector<T, ST> & r, const MatrixColumn<T2, ST2> & a, const math_matrix<T3, ST3> & b)
        { mul2(r, a, b); }
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        inline void mul(math_vector<T, ST> & r, const MatrixCColumn<T2, ST2> & a, const math_matrix<T3, ST3> & b)
        { mul2(r, a, b); }

// First argument: MatrixRow<>; matrix times vector
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        inline void mul(const MatrixRow<T2, ST2> & r, const math_matrix<T, ST> & a, const math_vector<T3, ST3> & b)
        {
            if ((void *)&r.matrix() == (void *)&a)
            {
                // Use temporary!
                math_vector<T, ST> tmp(r.size());
                mul2(tmp, a, b);
                assign(r, tmp);
            }
            else
                // Do it directly
                mul2(r, a, b);
        }

        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        inline void mul(const MatrixRow<T2, ST2> & r, const math_matrix<T, ST> & a, const MatrixRow<T3, ST3> & b)
        {
            if ((((void *)&r.matrix() == (void *)&b.matrix()) && (r.row() == b.row())) || ((void *)&r.matrix() == (void *)&a))
            {
                // Use temporary!
                math_vector<T, ST> tmp(r.size());
                mul2(tmp, a, b);
                assign(r, tmp);
            }
            else
                // Do it directly
                mul2(r, a, b);
        }

        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        inline void mul(const MatrixRow<T2, ST2> & r, const math_matrix<T, ST> & a, const MatrixCRow<T3, ST3> & b)
        {
            if ((((void *)&r.matrix() == (void *)&b.matrix()) && (r.row() == b.row())) || ((void *)&r.matrix() == (void *)&a))
            {
                // Use temporary!
                math_vector<T, ST> tmp(r.size());
                mul2(tmp, a, b);
                assign(r, tmp);
            }
            else
                // Do it directly
                mul2(r, a, b);
        }

        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        inline void mul(const MatrixRow<T2, ST2> & r, const math_matrix<T, ST> & a, const MatrixColumn<T3, ST3> & b)
        {
            if (((void *)&r.matrix() == (void *)&b.matrix()) || ((void *)&r.matrix() == (void *)&a))
            {
                // Use temporary!
                math_vector<T, ST> tmp(r.size());
                mul2(tmp, a, b);
                assign(r, tmp);
            }
            else
                // Do it directly
                mul2(r, a, b);
        }

        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        inline void mul(const MatrixRow<T2, ST2> & r, const math_matrix<T, ST> & a, const MatrixCColumn<T3, ST3> & b)
        {
            if (((void *)&r.matrix() == (void *)&b.matrix()) || ((void *)&r.matrix() == (void *)&a))
            {
                // Use temporary!
                math_vector<T, ST> tmp(r.size());
                mul2(tmp, a, b);
                assign(r, tmp);
            }
            else
                // Do it directly
                mul2(r, a, b);
        }

// First argument: MatrixRow<>; vector times matrix
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        inline void mul(const MatrixRow<T2, ST2> & r, const math_vector<T3, ST3> & a, const math_matrix<T, ST> & b)
        {
            if ((void *)&r.matrix() == (void *)&b)
            {
                // Use temporary!
                math_vector<T, ST> tmp(r.size());
                mul2(tmp, a, b);
                assign(r, tmp);
            }
            else
                // Do it directly
                mul2(r, a, b);
        }

        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        inline void mul(const MatrixRow<T2, ST2> & r, const MatrixRow<T3, ST3> & a, const math_matrix<T, ST> & b)
        {
            if ((((void *)&r.matrix() == (void *)&a.matrix()) && (r.row() == a.row())) || ((void *)&r.matrix() == (void *)&b))
            {
                // Use temporary!
                math_vector<T, ST> tmp(r.size());
                mul2(tmp, a, b);
                assign(r, tmp);
            }
            else
                // Do it directly
                mul2(r, a, b);
        }

        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        inline void mul(const MatrixRow<T2, ST2> & r, const MatrixCRow<T3, ST3> & a, const math_matrix<T, ST> & b)
        {
            if ((((void *)&r.matrix() == (void *)&a.matrix()) && (r.row() == a.row())) || ((void *)&r.matrix() == (void *)&b))
            {
                // Use temporary!
                math_vector<T, ST> tmp(r.size());
                mul2(tmp, a, b);
                assign(r, tmp);
            }
            else
                // Do it directly
                mul2(r, a, b);
        }

        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        inline void mul(const MatrixRow<T2, ST2> & r, const MatrixColumn<T3, ST3> & a, const math_matrix<T, ST> & b)
        {
            if (((void *)&r.matrix() == (void *)&a.matrix()) || ((void *)&r.matrix() == (void *)&b))
            {
                // Use temporary!
                math_vector<T, ST> tmp(r.size());
                mul2(tmp, a, b);
                assign(r, tmp);
            }
            else
                // Do it directly
                mul2(r, a, b);
        }

        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        inline void mul(const MatrixRow<T2, ST2> & r, const MatrixCColumn<T3, ST3> & a, const math_matrix<T, ST> & b)
        {
            if (((void *)&r.matrix() == (void *)&a.matrix()) || ((void *)&r.matrix() == (void *)&b))
            {
                // Use temporary!
                math_vector<T, ST> tmp(r.size());
                mul2(tmp, a, b);
                assign(r, tmp);
            }
            else
                // Do it directly
                mul2(r, a, b);
        }

// First argument: MatrixColumn<>; matrix times vector
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        inline void mul(const MatrixColumn<T2, ST2> & r, const math_matrix<T, ST> & a, const math_vector<T3, ST3> & b)
        {
            if ((void *)&r.matrix() == (void *)&a)
            {
                // Use temporary!
                math_vector<T, ST> tmp(r.size());
                mul2(tmp, a, b);
                assign(r, tmp);
            }
            else
                // Do it directly
                mul2(r, a, b);
        }

        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        inline void mul(const MatrixColumn<T2, ST2> & r, const math_matrix<T, ST> & a, const MatrixRow<T3, ST3> & b)
        {
            if (((void *)&r.matrix() == (void *)&b.matrix()) || ((void *)&r.matrix() == (void *)&a))
            {
                // Use temporary!
                math_vector<T, ST> tmp(r.size());
                mul2(tmp, a, b);
                assign(r, tmp);
            }
            else
                // Do it directly
                mul2(r, a, b);
        }

        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        inline void mul(const MatrixColumn<T2, ST2> & r, const math_matrix<T, ST> & a, const MatrixCRow<T3, ST3> & b)
        {
            if (((void *)&r.matrix() == (void *)&b.matrix()) || ((void *)&r.matrix() == (void *)&a))
            {
                // Use temporary!
                math_vector<T, ST> tmp(r.size());
                mul2(tmp, a, b);
                assign(r, tmp);
            }
            else
                // Do it directly
                mul2(r, a, b);
        }

        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        inline void mul(const MatrixColumn<T2, ST2> & r, const math_matrix<T, ST> & a, const MatrixColumn<T3, ST3> & b)
        {
            if ((((void *)&r.matrix() == (void *)&b.matrix()) && (r.column() == b.column())) || ((void *)&r.matrix() == (void *)&a))
            {
                // Use temporary!
                math_vector<T, ST> tmp(r.size());
                mul2(tmp, a, b);
                assign(r, tmp);
            }
            else
                // Do it directly
                mul2(r, a, b);
        }

        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        inline void mul(const MatrixColumn<T2, ST2> & r, const math_matrix<T, ST> & a, const MatrixCColumn<T3, ST3> & b)
        {
            if ((((void *)&r.matrix() == (void *)&b.matrix()) && (r.column() == b.column())) || ((void *)&r.matrix() == (void *)&a))
            {
                // Use temporary!
                math_vector<T, ST> tmp(r.size());
                mul2(tmp, a, b);
                assign(r, tmp);
            }
            else
                // Do it directly
                mul2(r, a, b);
        }

// First argument: MatrixColumn<>; vector times matrix
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        inline void mul(const MatrixColumn<T2, ST2> & r, const math_vector<T3, ST3> & a, const math_matrix<T, ST> & b)
        {
            if ((void *)&r.matrix() == (void *)&b)
            {
                // Use temporary!
                math_vector<T, ST> tmp(r.size());
                mul2(tmp, a, b);
                assign(r, tmp);
            }
            else
                // Do it directly
                mul2(r, a, b);
        }

        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        inline void mul(const MatrixColumn<T2, ST2> & r, const MatrixRow<T3, ST3> & a, const math_matrix<T, ST> & b)
        {
            if (((void *)&r.matrix() == (void *)&a.matrix()) || ((void *)&r.matrix() == (void *)&b))
            {
                // Use temporary!
                math_vector<T, ST> tmp(r.size());
                mul2(tmp, a, b);
                assign(r, tmp);
            }
            else
                // Do it directly
                mul2(r, a, b);
        }

        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        inline void mul(const MatrixColumn<T2, ST2> & r, const MatrixCRow<T3, ST3> & a, const math_matrix<T, ST> & b)
        {
            if (((void *)&r.matrix() == (void *)&a.matrix()) || ((void *)&r.matrix() == (void *)&b))
            {
                // Use temporary!
                math_vector<T, ST> tmp(r.size());
                mul2(tmp, a, b);
                assign(r, tmp);
            }
            else
                // Do it directly
                mul2(r, a, b);
        }

        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        inline void mul(const MatrixColumn<T2, ST2> & r, const MatrixColumn<T3, ST3> & a, const math_matrix<T, ST> & b)
        {
            if ((((void *)&r.matrix() == (void *)&a.matrix()) && (r.column() == a.column())) || ((void *)&r.matrix() == (void *)&b))
            {
                // Use temporary!
                math_vector<T, ST> tmp(r.size());
                mul2(tmp, a, b);
                assign(r, tmp);
            }
            else
                // Do it directly
                mul2(r, a, b);
        }

        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        inline void mul(const MatrixColumn<T2, ST2> & r, const MatrixCColumn<T3, ST3> & a, const math_matrix<T, ST> & b)
        {
            if ((((void *)&r.matrix() == (void *)&a.matrix()) && (r.column() == a.column())) || ((void *)&r.matrix() == (void *)&b))
            {
                // Use temporary!
                math_vector<T, ST> tmp(r.size());
                mul2(tmp, a, b);
                assign(r, tmp);
            }
            else
                // Do it directly
                mul2(r, a, b);
        }

// First argument: math_matrix<>; matrix times matrix
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        inline void mul(math_matrix<T, ST> & r, const math_matrix<T2, ST2> & a, const math_matrix<T3, ST3> & b)
        {
            if (((void *)&r == (void *)&a) || ((void *)&r == (void *)&b))
            {
                using std::swap;
                // Use temporary!
                math_matrix<T, ST> tmp(r.rows(), r.cols());
                mul2(tmp, a, b);
                swap(r, tmp);
            }
            else
                // Do it directly
                mul2(r, a, b);
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Matrices rows & columns

        template<class T, class ST, class V>
        inline void assign(const MatrixRow<T, ST> & R, const V & v)
        {
            assert(R.size() == v.size());
            for (typename base_matrix<T, ST>::size_type j = 0; j < R.size(); ++j)
                R[j] = v[j];
        }

        template<class T, class ST, class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
        inline void mod(const MatrixRow<T, ST> & r, const V & a, const T & b)
        {
            assert(a.size() == r.size());
            for (typename base_vector<T, ST>::size_type i = 0; i < a.size(); ++i)
                r[i] = a[i] % b;
        }

        template<class T, class ST, class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
        inline void div(const MatrixRow<T, ST> & r, const V & a, const T & b)
        {
            assert(a.size() == r.size());
            for (typename base_vector<T, ST>::size_type i = 0; i < a.size(); ++i)
                r[i] = a[i] / b;
        }

        template<class T, class ST, class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
        inline void mul(const MatrixRow<T, ST> & r, const V & a, const T & b)
        {
            assert(a.size() == r.size());
            for (typename base_vector<T, ST>::size_type i = 0; i < a.size(); ++i)
                r[i] = a[i] * b;
        }

        template<class T, class ST, class V2, class V3> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
        inline void add(const MatrixRow<T, ST> & r, const V2 & a, const V3 & b)
        {
            assert(a.size() == b.size());
            assert(a.size() == r.size());
            for (typename base_vector<T, ST>::size_type i = 0; i < a.size(); ++i)
                r[i] = a[i] + b[i];
        }

        template<class T, class ST, class V2, class V3> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
        inline void sub(const MatrixRow<T, ST> & r, const V2 & a, const V3 & b)
        {
            assert(a.size() == b.size());
            assert(a.size() == r.size());
            for (typename base_vector<T, ST>::size_type i = 0; i < a.size(); ++i)
                r[i] = a[i] - b[i];
        }

        template<class T, class ST, class V2, class TT> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
        inline void addmul(const MatrixRow<T, ST> & r, const V2 & a, const TT & b)
        {
            assert(a.size() == r.size());
            for (typename base_vector<T, ST>::size_type i = 0; i < a.size(); ++i)
                r[i] += a[i] * b;
        }

        template<class T, class ST, class V2, class TT> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
        inline void submul(const MatrixRow<T, ST> & r, const V2 & a, const TT & b)
        {
            assert(a.size() == r.size());
            for (typename base_vector<T, ST>::size_type i = 0; i < a.size(); ++i)
                r[i] -= a[i] * b;
        }

        template<class T, class ST, class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
        inline void neg(const MatrixRow<T, ST> & r, const V & a)
        {
            assert(a.size() == r.size());
            for (typename base_vector<T, ST>::size_type i = 0; i < a.size(); ++i)
                r[i] = -a[i];
        }
        
        template<class T, class ST, class V>
        inline void swap(const MatrixRow<T, ST> & v1, V & v2)
        {
            assert(v1.size() == v2.size());
            for (unsigned i = 0; i < v1.size(); ++i)
                swap(v1[i], v2[i]);
        }
        
        template<class T, class ST, class V>
        inline void swap(const MatrixRow<T, ST> & v1, const V & v2)
        {
            assert(v1.size() == v2.size());
            for (unsigned i = 0; i < v1.size(); ++i)
                swap(v1[i], v2[i]);
        }
        
        template<class T, class ST>
        inline void normSq(T & r, const MatrixRow<T, ST> & a)
        {
            setZero(r);
            for (typename base_vector<T, ST>::size_type i = 0; i < a.size(); ++i)
                r += a[i] * a[i];
        }

        template<class T, class ST, class V>
        inline void dot(T & r, const MatrixRow<T, ST> & a, const V & b)
        {
            assert(a.size() == b.size());
            setZero(r);
            for (typename base_vector<T, ST>::size_type i = 0; i < a.size(); ++i)
                r += a[i] * b[i];
        }

        template<class T, class ST, class V>
        inline void assign(const MatrixColumn<T, ST> & C, const V & v)
        {
            assert(C.size() == v.size());
            for (typename base_matrix<T, ST>::size_type j = 0; j < C.size(); ++j)
                C[j] = v[j];
        }

        template<class T, class ST, class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
        inline void mod(const MatrixColumn<T, ST> & r, const V & a, const T & b)
        {
            assert(a.size() == r.size());
            for (typename base_vector<T, ST>::size_type i = 0; i < a.size(); ++i)
                r[i] = a[i] % b;
        }

        template<class T, class ST, class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
        inline void div(const MatrixColumn<T, ST> & r, const V & a, const T & b)
        {
            assert(a.size() == r.size());
            for (typename base_vector<T, ST>::size_type i = 0; i < a.size(); ++i)
                r[i] = a[i] / b;
        }

        template<class T, class ST, class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
        inline void mul(const MatrixColumn<T, ST> & r, const V & a, const T & b)
        {
            assert(a.size() == r.size());
            for (typename base_vector<T, ST>::size_type i = 0; i < a.size(); ++i)
                r[i] = a[i] * b;
        }

        template<class T, class ST, class V2, class V3> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
        inline void add(const MatrixColumn<T, ST> & r, const V2 & a, const V3 & b)
        {
            assert(a.size() == b.size());
            assert(a.size() == r.size());
            for (typename base_vector<T, ST>::size_type i = 0; i < a.size(); ++i)
                r[i] = a[i] + b[i];
        }

        template<class T, class ST, class V2, class V3> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
        inline void sub(const MatrixColumn<T, ST> & r, const V2 & a, const V3 & b)
        {
            assert(a.size() == b.size());
            assert(a.size() == r.size());
            for (typename base_vector<T, ST>::size_type i = 0; i < a.size(); ++i)
                r[i] = a[i] - b[i];
        }

        template<class T, class ST, class V2, class TT> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
        inline void addmul(const MatrixColumn<T, ST> & r, const V2 & a, const TT & b)
        {
            assert(a.size() == r.size());
            for (typename base_vector<T, ST>::size_type i = 0; i < a.size(); ++i)
                r[i] += a[i] * b;
        }

        template<class T, class ST, class V2, class TT> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
        inline void submul(const MatrixColumn<T, ST> & r, const V2 & a, const TT & b)
        {
            assert(a.size() == r.size());
            for (typename base_vector<T, ST>::size_type i = 0; i < a.size(); ++i)
                r[i] -= a[i] * b;
        }

        template<class T, class ST, class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
        inline void neg(const MatrixColumn<T, ST> & r, const V & a)
        {
            assert(a.size() == r.size());
            for (typename base_vector<T, ST>::size_type i = 0; i < a.size(); ++i)
                r[i] = -a[i];
        }
        
        template<class T, class ST, class V>
        inline void swap(const MatrixColumn<T, ST> & v1, V & v2)
        {
            assert(v1.size() == v2.size());
            for (unsigned i = 0; i < v1.size(); ++i)
                swap(v1[i], v2[i]);
        }
        
        template<class T, class ST, class V>
        inline void swap(const MatrixColumn<T, ST> & v1, const V & v2)
        {
            assert(v1.size() == v2.size());
            for (unsigned i = 0; i < v1.size(); ++i)
                swap(v1[i], v2[i]);
        }
        
        template<class T, class ST>
        inline void normSq(T & r, const MatrixColumn<T, ST> & a)
        {
            setZero(r);
            for (typename base_vector<T, ST>::size_type i = 0; i < a.size(); ++i)
                r += a[i] * a[i];
        }

        template<class T, class ST, class V>
        inline void dot(T & r, const MatrixColumn<T, ST> & a, const V & b)
        {
            assert(a.size() == b.size());
            setZero(r);
            for (typename base_vector<T, ST>::size_type i = 0; i < a.size(); ++i)
                r += a[i] * b[i];
        }

        // Constant accessors for rows and columns (constant and non-constant); these behave similar to
        // a math_vector<>.

        template<class T, class ST>
        inline void normSq(T & r, const MatrixCRow<T, ST> & a)
        {
            setZero(r);
            for (typename base_vector<T, ST>::size_type i = 0; i < a.size(); ++i)
                r += a[i] * a[i];
        }

        template<class T, class ST, class V>
        inline void dot(T & r, const MatrixCRow<T, ST> & a, const V & b)
        {
            assert(a.size() == b.size());
            setZero(r);
            for (typename base_vector<T, ST>::size_type i = 0; i < a.size(); ++i)
                r += a[i] * b[i];
        }

        template<class T, class ST>
        inline void normSq(T & r, const MatrixCColumn<T, ST> & a)
        {
            setZero(r);
            for (typename base_vector<T, ST>::size_type i = 0; i < a.size(); ++i)
                r += a[i] * a[i];
        }

        template<class T, class ST, class V>
        inline void dot(T & r, const MatrixCColumn<T, ST> & a, const V & b)
        {
            assert(a.size() == b.size());
            setZero(r);
            for (typename base_vector<T, ST>::size_type i = 0; i < a.size(); ++i)
                r += a[i] * b[i];
        }

        template<class T, class ST, class V>
        inline T dot(const math_vector<T, ST> & a, const V & b)
        {
            T r;
            dot(r, a, b);
            return r;
        }
        
        template<class T, class ST, class V>
        inline T dot(const MatrixRow<T, ST> & a, const V & b)
        {
            T r;
            dot(r, a, b);
            return r;
        }
        
        template<class T, class ST, class V>
        inline T dot(const MatrixCRow<T, ST> & a, const V & b)
        {
            T r;
            dot(r, a, b);
            return r;
        }
        
        template<class T, class ST, class V>
        inline T dot(const MatrixColumn<T, ST> & a, const V & b)
        {
            T r;
            dot(r, a, b);
            return r;
        }
        
        template<class T, class ST, class V>
        inline T dot(const MatrixCColumn<T, ST> & a, const V & b)
        {
            T r;
            dot(r, a, b);
            return r;
        }
        
        template<class O1, class D1, class T1, class V>
        inline T1 dot(const VecExpression<O1, D1, T1> & a, const V & b)
        {
            T1 r;
            setZero(r);
            for (unsigned i = 0; i < a.size(); ++i)
                r += a.evalComponent(i) * b[i];
            return r;
        }

        template<class T, class ST>
        inline T normSq(const math_vector<T, ST> & a)
        {
            T r;
            normSq(r, a);
            return r;
        }
        
        template<class T, class ST>
        inline T normSq(const MatrixRow<T, ST> & a)
        {
            T r;
            normSq(r, a);
            return r;
        }
        
        template<class T, class ST>
        inline T normSq(const MatrixCRow<T, ST> & a)
        {
            T r;
            normSq(r, a);
            return r;
        }
        
        template<class T, class ST>
        inline T normSq(const MatrixColumn<T, ST> & a)
        {
            T r;
            normSq(r, a);
            return r;
        }
        
        template<class T, class ST>
        inline T normSq(const MatrixCColumn<T, ST> & a)
        {
            T r;
            normSq(r, a);
            return r;
        }
        
        template<class O1, class D1, class T1>
        inline T1 normSq(const VecExpression<O1, D1, T1> & a)
        {
            T1 r, t;
            setZero(r);
            for (unsigned i = 0; i < a.size(); ++i)
            {
                t = a.evalComponent(i);
                r += t * t;
            }
            return r;
        }
    }
}

#endif

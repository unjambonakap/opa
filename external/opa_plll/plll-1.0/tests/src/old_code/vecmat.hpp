#ifndef PLLL_INCLUDE_GUARD__VECMAT_HPP
#define PLLL_INCLUDE_GUARD__VECMAT_HPP

#include <cassert>
#include <algorithm> // for std::swap()

#include "vecmat-mem.hpp"

namespace plll
{
    namespace linalgold
    {
        // Forward Declarations
        
        template<class Container, class ScT>
        class VecWrapper;
        template<class Op, class Data, class ScT>
        class VecExpression;

        template<class T, class StorageTraits>
        class base_vector;

        template<class T, class ST, class V>
        void assign(base_vector<T, ST> & r, const V & v);
        
        template<class T, class ST>
        class math_vector;

        template<class T, class ST, class V>
        void mod(math_vector<T, ST> & r, const V & a, const T & b);
        template<class T, class ST, class V>
        void div(math_vector<T, ST> & r, const V & a, const T & b);
        template<class T, class ST, class V>
        void mul(math_vector<T, ST> & r, const V & a, const T & b);
        template<class T, class ST, class V2, class V3>
        void add(math_vector<T, ST> & r, const V2 & a, const V3 & b);
        template<class T, class ST, class V2, class V3>
        void sub(math_vector<T, ST> & r, const V2 & a, const V3 & b);
        template<class T, class ST, class V2, class TT>
        void addmul(math_vector<T, ST> & r, const V2 & a, const TT & b);
        template<class T, class ST, class V2, class TT>
        void submul(math_vector<T, ST> & r, const V2 & a, const TT & b);
        template<class T, class ST, class V>
        void neg(math_vector<T, ST> & r, const V & a);
        template<class T, class ST>
        void normSq(T & r, const math_vector<T, ST> & a);
        template<class T, class ST, class V>
        void dot(T & r, const math_vector<T, ST> & a, const V & b);

        template<class T, class StorageTraits>
        class base_matrix;

        template<class T, class ST, class S, class A>
        void assign(base_matrix<T, ST> & r, const base_matrix<S, A> & m);
        
        template<class T, class ST, class V>
        void swap(base_vector<T, ST> & v1, V & v2);
        template<class T, class ST>
        void swap(base_vector<T, ST> & A, base_vector<T, ST> & B);
        
        template<class T, class ST>
        class math_matrix;

        template<class T, class ST>
        class MatrixRow;
        template<class T, class ST>
        class MatrixColumn;
        template<class T, class ST>
        class MatrixCRow;
        template<class T, class ST>
        class MatrixCColumn;

        template<class T, class ST, class T2, class ST2>
        void mod(math_matrix<T, ST> & r, const math_matrix<T2, ST2> & a, const T & b);
        template<class T, class ST, class T2, class ST2>
        void div(math_matrix<T, ST> & r, const math_matrix<T2, ST2> & a, const T & b);
        template<class T, class ST, class T2, class ST2>
        void mul(math_matrix<T, ST> & r, const math_matrix<T2, ST2> & a, const T & b);
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        void add(math_matrix<T, ST> & r, const math_matrix<T2, ST2> & a, const math_matrix<T3, ST3> & b);
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        void sub(math_matrix<T, ST> & r, const math_matrix<T2, ST2> & a, const math_matrix<T3, ST3> & b);
        template<class T, class ST, class T2, class ST2>
        void neg(math_matrix<T, ST> & r, const math_matrix<T2, ST2> & a);
        template<class T, class ST, class T2, class ST2>
        void transpose(math_matrix<T, ST> & r, const math_matrix<T2, ST2> & a);
        template<class T, class ST>
        void transpose(math_matrix<T, ST> & r, const math_matrix<T, ST> & a);
        template<class T, class ST, class T2, class ST2, class V>
        void mul2(math_vector<T, ST> & r, const math_matrix<T2, ST2> & a, const V & b);
        template<class T, class ST, class T2, class ST2, class V>
        void mul2(math_vector<T, ST> & r, const V & a, const math_matrix<T2, ST2> & b);
        template<class T, class ST, class T2, class ST2, class V>
        void mul2(const MatrixRow<T2, ST2> & r, const math_matrix<T, ST> & a, const V & b);
        template<class T, class ST, class T2, class ST2, class V>
        void mul2(const MatrixRow<T2, ST2> & r, const V & a, const math_matrix<T, ST> & b);
        template<class T, class ST, class T2, class ST2, class V>
        void mul2(const MatrixColumn<T2, ST2> & r, const math_matrix<T, ST> & a, const V & b);
        template<class T, class ST, class T2, class ST2, class V>
        void mul2(const MatrixColumn<T2, ST2> & r, const V & a, const math_matrix<T, ST> & b);
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        void mul2(math_matrix<T, ST> & r, const math_matrix<T2, ST2> & a, const base_matrix<T3, ST3> & b);
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        void mul(math_vector<T, ST> & r, const math_matrix<T3, ST3> & a, const math_vector<T2, ST2> & b);
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        void mul(math_vector<T, ST> & r, const math_matrix<T3, ST3> & a, const MatrixRow<T2, ST2> & b);
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        void mul(math_vector<T, ST> & r, const math_matrix<T3, ST3> & a, const MatrixCRow<T2, ST2> & b);
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        void mul(math_vector<T, ST> & r, const math_matrix<T3, ST3> & a, const MatrixColumn<T2, ST2> & b);
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        void mul(math_vector<T, ST> & r, const math_matrix<T3, ST3> & a, const MatrixCColumn<T2, ST2> & b);
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        void mul(math_vector<T, ST> & r, const math_vector<T2, ST2> & a, const math_matrix<T3, ST3> & b);
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        void mul(math_vector<T, ST> & r, const MatrixRow<T2, ST2> & a, const math_matrix<T3, ST3> & b);
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        void mul(math_vector<T, ST> & r, const MatrixCRow<T2, ST2> & a, const math_matrix<T3, ST3> & b);
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        void mul(math_vector<T, ST> & r, const MatrixColumn<T2, ST2> & a, const math_matrix<T3, ST3> & b);
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        void mul(math_vector<T, ST> & r, const MatrixCColumn<T2, ST2> & a, const math_matrix<T3, ST3> & b);
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        void mul(const MatrixRow<T2, ST2> & r, const math_matrix<T, ST> & a, const math_vector<T3, ST3> & b);
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        void mul(const MatrixRow<T2, ST2> & r, const math_matrix<T, ST> & a, const MatrixRow<T3, ST3> & b);
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        void mul(const MatrixRow<T2, ST2> & r, const math_matrix<T, ST> & a, const MatrixCRow<T3, ST3> & b);
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        void mul(const MatrixRow<T2, ST2> & r, const math_matrix<T, ST> & a, const MatrixColumn<T3, ST3> & b);
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        void mul(const MatrixRow<T2, ST2> & r, const math_matrix<T, ST> & a, const MatrixCColumn<T3, ST3> & b);
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        void mul(const MatrixRow<T2, ST2> & r, const math_vector<T3, ST3> & a, const math_matrix<T, ST> & b);
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        void mul(const MatrixRow<T2, ST2> & r, const MatrixRow<T3, ST3> & a, const math_matrix<T, ST> & b);
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        void mul(const MatrixRow<T2, ST2> & r, const MatrixCRow<T3, ST3> & a, const math_matrix<T, ST> & b);
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        void mul(const MatrixRow<T2, ST2> & r, const MatrixColumn<T3, ST3> & a, const math_matrix<T, ST> & b);
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        void mul(const MatrixRow<T2, ST2> & r, const MatrixCColumn<T3, ST3> & a, const math_matrix<T, ST> & b);
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        void mul(const MatrixColumn<T2, ST2> & r, const math_matrix<T, ST> & a, const math_vector<T3, ST3> & b);
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        void mul(const MatrixColumn<T2, ST2> & r, const math_matrix<T, ST> & a, const MatrixRow<T3, ST3> & b);
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        void mul(const MatrixColumn<T2, ST2> & r, const math_matrix<T, ST> & a, const MatrixCRow<T3, ST3> & b);
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        void mul(const MatrixColumn<T2, ST2> & r, const math_matrix<T, ST> & a, const MatrixColumn<T3, ST3> & b);
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        void mul(const MatrixColumn<T2, ST2> & r, const math_matrix<T, ST> & a, const MatrixCColumn<T3, ST3> & b);
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        void mul(const MatrixColumn<T2, ST2> & r, const math_vector<T3, ST3> & a, const math_matrix<T, ST> & b);
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        void mul(const MatrixColumn<T2, ST2> & r, const MatrixRow<T3, ST3> & a, const math_matrix<T, ST> & b);
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        void mul(const MatrixColumn<T2, ST2> & r, const MatrixCRow<T3, ST3> & a, const math_matrix<T, ST> & b);
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        void mul(const MatrixColumn<T2, ST2> & r, const MatrixColumn<T3, ST3> & a, const math_matrix<T, ST> & b);
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        void mul(const MatrixColumn<T2, ST2> & r, const MatrixCColumn<T3, ST3> & a, const math_matrix<T, ST> & b);
        template<class T, class ST, class T2, class ST2, class T3, class ST3>
        void mul(math_matrix<T, ST> & r, const math_matrix<T2, ST2> & a, const math_matrix<T3, ST3> & b);
        
        template<class T, class ST>
        void swap(base_matrix<T, ST> & A, base_matrix<T, ST> & B);
        
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Vectors

        template<class T, class StorageTraits = storage_traits<T> >
        class base_vector
        {
        public:
            typedef unsigned int size_type;
    
            enum { has_direct_access = true };
    
#ifdef PLLL_INTERNAL_NO_TEMPLATE_FRIENDS
        public:
#else
        protected:
#endif
            size_type d_size;
            typename StorageTraits::pointer_type d_data;

            inline void assign(const base_vector<T, StorageTraits> & v)
            {
                if (reinterpret_cast<typename StorageTraits::pointer_type>(v.data()) != d_data)
                {
                    StorageTraits::free(d_data, d_size);
                    d_size = v.size();
                    d_data = StorageTraits::clone(v.data(), v.size());
                }
            }
    
            template<class V> // this can be something else which behaves like a vector
            inline void assign(const V & v)
            {
                resize(v.size());
                for (size_type i = 0; i < v.size(); ++i)
                    d_data[i] = v[i];
            }
    
            template<class O1, class D1, class T1>
            inline void assign(const VecExpression<O1, D1, T1> & e)
            {
                resize(e.size());
                if (e.usesVecNF(*this))
                {
                    T tmp;
                    for (size_type i = 0; i < d_size; ++i)
                    {
                        e.evalComponentTo(tmp, i);
                        d_data[i] = tmp;
                    }
                }
                else
                {
                    for (size_type i = 0; i < d_size; ++i)
                        e.evalComponentTo(d_data[i], i);
                }
            }
    
        public:
            inline base_vector()
                : d_size(0), d_data(NULL)
            {
            }
    
            inline base_vector(size_type size)
                : d_size(size), d_data(StorageTraits::alloc(size))
            {
            }
    
            template<class TT>
            inline base_vector(size_type size, const TT & obj)
                : d_size(size), d_data(StorageTraits::alloc(obj, size))
            {
            }
    
            inline base_vector(const base_vector<T, StorageTraits> & v)
                : d_size(v.size()), d_data(StorageTraits::clone(v.data(), v.size()))
            {
            }

            template<class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            inline explicit base_vector(const V & v)
                : d_size(v.size()), d_data(StorageTraits::clone(v.data(), v.size()))
            {
            }
    
            template<class O1, class D1, class T1>
            inline explicit base_vector(const VecExpression<O1, D1, T1> & e)
                : d_size(e.size()), d_data(StorageTraits::alloc(e.size()))
            {
                if (e.usesVecNF(*this))
                {
                    T tmp;
                    for (size_type i = 0; i < d_size; ++i)
                    {
                        e.evalComponentTo(tmp, i);
                        d_data[i] = tmp;
                    }
                }
                else
                {
                    for (size_type i = 0; i < d_size; ++i)
                        e.evalComponentTo(d_data[i], i);
                }
            }
    
            inline base_vector & operator = (const base_vector<T, StorageTraits> & v)
            {
                StorageTraits::check(d_data, d_size);
                assign(v);
                return *this;
            }
    
            template<class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            inline base_vector & operator = (const V & v)
            {
                StorageTraits::check(d_data, d_size);
                assign(v);
                return *this;
            }
    
            template<class O1, class D1, class T1>
            inline base_vector & operator = (const VecExpression<O1, D1, T1> & e)
            {
                StorageTraits::check(d_data, d_size);
                assign(e);
                return *this;
            }
    
            inline ~base_vector()
            {
                StorageTraits::free(d_data, d_size);
            }
    
            inline size_type size() const
            {
                return d_size;
            }
    
            void resize(size_type size)
            {
                if (size == d_size)
                    return;
                typename StorageTraits::pointer_type data = StorageTraits::clone(d_data, size, d_size > size ? size : d_size);
                StorageTraits::free(d_data, d_size);
                d_data = data;
                d_size = size;
            }
    
            template<class S>
            void resize(size_type size, const S & copyobj)
            {
                if (size == d_size)
                    return;
                typename StorageTraits::pointer_type data = StorageTraits::clone(d_data, size, d_size > size ? size : d_size, copyobj);
                StorageTraits::free(d_data, d_size);
                d_data = data;
                d_size = size;
            }
    
            inline typename StorageTraits::pointer_type data() const
            {
                StorageTraits::check(d_data, d_size);
                return d_data;
            }
    
            inline typename StorageTraits::constref_type operator [] (size_type i) const
            {
                StorageTraits::check(d_data, d_size);
                return d_data[i];
            }
    
            inline typename StorageTraits::ref_type operator [] (size_type i)
            {
                StorageTraits::check(d_data, d_size);
                return d_data[i];
            }
    
            inline typename StorageTraits::constref_type operator () (size_type i) const
            {
                StorageTraits::check(d_data, d_size);
                assert((0 <= i) && (i < d_size));
                return d_data[i];
            }
    
            inline typename StorageTraits::ref_type operator () (size_type i)
            {
                StorageTraits::check(d_data, d_size);
                assert((0 <= i) && (i < d_size));
                return d_data[i];
            }
    
            template<class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            inline bool operator == (const V & v) const
            {
                StorageTraits::check(d_data, d_size);
                if (v.size() != d_size)
                    return false;
                for (size_type i = 0; i < d_size; ++i)
                    if (d_data[i] != v.d_data[i])
                        return false;
                return true;
            }
    
            template<class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            inline bool operator != (const V & v) const
            {
                StorageTraits::check(d_data, d_size);
                if (v.size() != d_size)
                    return true;
                for (size_type i = 0; i < d_size; ++i)
                    if (d_data[i] != v.d_data[i])
                        return true;
                return false;
            }
    
            inline void swap(base_vector<T, StorageTraits> & B)
            {
                std::swap(d_size, B.d_size);
                std::swap(d_data, B.d_data);
            }
    
            template<class V>
            inline void swap(V & v)
            {
                assert(size() == v.size());
                for (unsigned i = 0; i < size(); ++i)
                    d_data[i].swap(v[i]);
            }

#ifndef PLLL_INTERNAL_NO_TEMPLATE_FRIENDS
            template<class TT, class STT, class V>
            friend void assign(base_vector<TT, STT> & r, const V & v);
#endif
        };

        template<class T, class ST = storage_traits<T> >
        class math_vector : public base_vector<T, ST>
        {
        public:
            enum { has_direct_access = true };
    
            inline math_vector()
            {
            }
    
            inline math_vector(typename base_vector<T, ST>::size_type size)
                : base_vector<T, ST>(size)
            {
            }
    
            inline math_vector(typename base_vector<T, ST>::size_type size, typename ST::constref_type obj)
                : base_vector<T, ST>(size, obj)
            {
            }
    
            inline math_vector(const math_vector<T, ST> & v)
                : base_vector<T, ST>(v)
            {
            }
    
            template<class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            inline explicit math_vector(const V & v)
                : base_vector<T, ST>(v)
            {
            }
    
            template<class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            inline math_vector & operator = (const V & v)
            {
                linalgold::assign(*this, v);
                return *this;
            }
    
            inline math_vector<T, ST> & operator %= (const T & s)
            {
                mod(*this, *this, s);
                return *this;
            }
    
            inline math_vector<T, ST> & operator /= (const T & s)
            {
                div(*this, *this, s);
                return *this;
            }
    
            inline math_vector<T, ST> & operator *= (const T & s)
            {
                mul(*this, *this, s);
                return *this;
            }
    
            template<class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            inline math_vector<T, ST> & operator += (const V & v)
            {
                add(*this, *this, v);
                return *this;
            }
    
            template<class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            inline math_vector<T, ST> & operator -= (const V & v)
            {
                sub(*this, *this, v);
                return *this;
            }
    
            template<class O1, class D1, class T1>
            inline math_vector<T, ST> & operator += (const VecExpression<O1, D1, T1> & e)
            {
                assert(this->d_size == e.size());
                if (e.usesVecNF(*this))
                {
                    T tmp;
                    for (unsigned i = 0; i < this->d_size; ++i)
                    {
                        e.evalComponentTo(tmp, i);
                        this->d_data[i] += tmp;
                    }
                }
                else
                {
                    for (unsigned i = 0; i < this->d_size; ++i)
                        e.evalComponentAddTo(this->d_data[i], i);
                }
                return *this;
            }
    
            template<class O1, class D1, class T1>
            inline math_vector<T, ST> & operator -= (const VecExpression<O1, D1, T1> & e)
            {
                assert(this->d_size == e.size());
                if (e.usesVecNF(*this))
                {
                    T tmp;
                    for (unsigned i = 0; i < this->d_size; ++i)
                    {
                        e.evalComponentTo(tmp, i);
                        this->d_data[i] -= tmp;
                    }
                }
                else
                {
                    for (unsigned i = 0; i < this->d_size; ++i)
                        e.evalComponentSubFrom(this->d_data[i], i);
                }
                return *this;
            }
    
            template<class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            inline bool operator < (const V & v) const
            {
                if (v.size() != this->d_size)
                    return false;
                for (typename base_vector<T, ST>::size_type i = 0; i < this->d_size; ++i)
                {
                    if ((*this)[i] < v[i])
                        return true;
                    if ((*this)[i] > v[i])
                        return false;
                }
                return false;
            }
    
            template<class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            inline bool operator <= (const V & v) const
            {
                if (v.size() != this->d_size)
                    return false;
                for (typename base_vector<T, ST>::size_type i = 0; i < this->d_size; ++i)
                {
                    if ((*this)[i] < v[i])
                        return true;
                    if ((*this)[i] > v[i])
                        return false;
                }
                return true;
            }
    
            template<class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            inline bool operator > (const V & v) const
            {
                return !(*this <= v);
            }
    
            template<class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            inline bool operator >= (const V & v) const
            {
                return !(*this < v);
            }
    
            inline T normSq() const
            {
                T res;
                linalgold::normSq(res, *this);
                return res;
            }

#ifndef PLLL_INTERNAL_NO_TEMPLATE_FRIENDS
            template<class TT, class STT, class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            friend void mod(math_vector<TT, STT> & r, const V & a, const TT & b);
            template<class TT, class STT, class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            friend void div(math_vector<TT, STT> & r, const V & a, const TT & b);
            template<class TT, class STT, class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            friend void mul(math_vector<TT, STT> & r, const V & a, const TT & b);
            template<class TT, class STT, class V2, class V3> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            friend void add(math_vector<TT, STT> & r, const V2 & a, const V3 & b);
            template<class TT, class STT, class V2, class V3> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            friend void sub(math_vector<TT, STT> & r, const V2 & a, const V3 & b);
            template<class TT, class STT, class V2, class TTT> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            friend void addmul(math_vector<TT, STT> & r, const V2 & a, const TTT & b);
            template<class TT, class STT, class V2, class TTT> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            friend void submul(math_vector<TT, STT> & r, const V2 & a, const TTT & b);
            template<class TT, class STT, class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            friend void neg(math_vector<TT, STT> & r, const V & a);
            friend void linalgold::normSq<>(T & r, const math_vector<T, ST> & a);
            template<class TT, class STT, class V>
            friend void dot(TT & r, const math_vector<TT, STT> & a, const V & b);
#endif
        };

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Matrices

        template<class T, class StorageTraits = storage_traits<T> >
        class base_matrix
        {
        public:
            typedef unsigned int size_type;
    
            enum { has_direct_access = true };
    
#ifdef PLLL_INTERNAL_NO_TEMPLATE_FRIENDS
        public:
#else
        protected:
#endif
            size_type d_rows, d_cols;
            typename StorageTraits::pointer_type d_data;
    
            template<class S, class A>
            inline void assign(const base_matrix<S, A> & m)
            {
                if (reinterpret_cast<typename StorageTraits::pointer_type>(m.data()) != d_data)
                {
                    StorageTraits::free(d_data, d_rows * d_cols);
                    d_rows = m.rows();
                    d_cols = m.cols();
                    d_data = StorageTraits::clone(m.data(), d_rows * d_cols);
                }
            }
    
        public:
#ifndef PLLL_INTERNAL_NO_TEMPLATE_FRIENDS
            template<class TT, class STT, class S, class A>
            friend void linalgold::assign(base_matrix<TT, STT> & r, const base_matrix<S, A> & m);
#endif
    
            inline void swap(base_matrix<T, StorageTraits> & B)
            {
                std::swap(d_rows, B.d_rows);
                std::swap(d_cols, B.d_cols);
                std::swap(d_data, B.d_data);
            }

#ifndef PLLL_INTERNAL_NO_TEMPLATE_FRIENDS
            friend void linalgold::swap<>(base_matrix<T, StorageTraits> & A, base_matrix<T, StorageTraits> & B);
#endif
    
            inline base_matrix()
                : d_rows(0), d_cols(0), d_data(NULL)
            {
            }
    
            inline base_matrix(size_type rows, size_type cols)
                : d_rows(rows), d_cols(cols), d_data(StorageTraits::alloc(d_rows * d_cols))
            {
            }
    
            inline base_matrix(size_type rows, size_type cols, typename StorageTraits::constref_type obj)
                : d_rows(rows), d_cols(cols), d_data(StorageTraits::alloc(obj, d_rows * d_cols))
            {
            }
    
            inline base_matrix(const base_matrix<T, StorageTraits> & m)
                : d_rows(m.rows()), d_cols(m.cols()), d_data(StorageTraits::clone(m.data(), m.rows() * m.cols()))
            {
            }
    
            template<class S, class A>
            inline explicit base_matrix(const base_matrix<S, A> & m)
                : d_rows(m.rows()), d_cols(m.cols()), d_data(StorageTraits::clone(m.data(), m.rows() * m.cols()))
            {
            }
    
            inline base_matrix & operator = (const base_matrix<T, StorageTraits> & m)
            {
                StorageTraits::check(d_data, d_rows * d_cols);
                assign(m);
                return *this;
            }
    
            template<class S, class A>
            inline base_matrix & operator = (const base_matrix<S, A> & m)
            {
                StorageTraits::check(d_data, d_rows * d_cols);
                assign(m);
                return *this;
            }
    
            inline ~base_matrix()
            {
                StorageTraits::free(d_data, d_rows * d_cols);
            }
    
            inline size_type rows() const
            {
                return d_rows;
            }
    
            inline size_type cols() const
            {
                return d_cols;
            }
    
            void resize(size_type rows, size_type cols)
            {
                if ((rows == d_rows) && (cols == d_cols))
                    return;
                StorageTraits::check(d_data, d_rows * d_cols);
                typename StorageTraits::pointer_type data = StorageTraits::alloc_dontconstruct(rows * cols);
                // Copy existing and construct the rest
                for (size_type i = 0; (i < rows) && (i < d_rows); ++i)
                {
                    for (size_type j = 0; (j < cols) && (j < d_cols); ++j)
                        StorageTraits::copy_construct(data[i * cols + j], d_data[i * d_cols + j]);
                    if (d_cols < cols)
                        for (size_type j = d_cols; j < cols; ++j)
                            StorageTraits::construct(data[i * cols + j]);
                }
                // Construct the rest
                if (d_rows < rows)
                    for (size_type i = d_rows; i < rows; ++i)
                        for (size_type j = 0; j < cols; ++j)
                            StorageTraits::construct(data[i * cols + j]);
                StorageTraits::free(d_data, d_rows * d_cols);
                d_data = data;
                d_rows = rows;
                d_cols = cols;
                StorageTraits::check(d_data, d_rows * d_cols);
            }
    
            template<class S>
            void resize(size_type rows, size_type cols, const S & initobj)
            {
                if ((rows == d_rows) && (cols == d_cols))
                    return;
                StorageTraits::check(d_data, d_rows * d_cols);
                typename StorageTraits::pointer_type data = StorageTraits::alloc_dontconstruct(rows * cols);
                // Copy existing and construct the rest
                for (size_type i = 0; (i < rows) && (i < d_rows); ++i)
                {
                    for (size_type j = 0; (j < cols) && (j < d_cols); ++j)
                        StorageTraits::copy_construct(data[i * cols + j], d_data[i * d_cols + j]);
                    if (d_cols < cols)
                        for (size_type j = d_cols; j < cols; ++j)
                            StorageTraits::copy_construct(data[i * cols + j], initobj);
                }
                // Construct the rest
                if (d_rows < rows)
                    for (size_type i = d_rows; i < rows; ++i)
                        for (size_type j = 0; j < cols; ++j)
                            StorageTraits::copy_construct(data[i * cols + j], initobj);
                StorageTraits::free(d_data, d_rows * d_cols);
                d_data = data;
                d_rows = rows;
                d_cols = cols;
                StorageTraits::check(d_data, d_rows * d_cols);
            }
    
            inline typename StorageTraits::constref_type operator () (size_type i, size_type j) const
            {
                StorageTraits::check(d_data, d_rows * d_cols);
                assert((0 <= i) && (i < d_rows));
                assert((0 <= j) && (j < d_cols));
                return d_data[i * d_cols + j];
            }
    
            inline typename StorageTraits::ref_type operator () (size_type i, size_type j)
            {
                StorageTraits::check(d_data, d_rows * d_cols);
                assert((0 <= i) && (i < d_rows));
                assert((0 <= j) && (j < d_cols));
                return d_data[i * d_cols + j];
            }
    
            inline typename StorageTraits::pointer_type data() const
            {
                StorageTraits::check(d_data, d_rows * d_cols);
                return d_data;
            }
    
            inline typename StorageTraits::constref_type data(size_type i) const
            {
                StorageTraits::check(d_data, d_rows * d_cols);
                return d_data[i];
            }
    
            inline typename StorageTraits::ref_type data(size_type i)
            {
                StorageTraits::check(d_data, d_rows * d_cols);
                return d_data[i];
            }
    
            template<class S, class A>
            inline bool operator == (const base_matrix<S, A> & m) const
            {
                StorageTraits::check(d_data, d_rows * d_cols);
                if ((m.rows() != d_rows) || (m.cols() != d_cols))
                    return false;
                size_type c = d_rows * d_cols;
                for (size_type i = 0; i < c; ++i)
                    if (d_data[i] != m.data(i))
                        return false;
                return true;
            }
    
            template<class S, class A>
            inline bool operator != (const base_matrix<S, A> & m) const
            {
                return !(*this == m);
            }
    
            inline base_matrix<T, StorageTraits> transpose() const
            {
                base_matrix<T, StorageTraits> res(d_cols, d_rows);
                for (unsigned i = 0; i < d_cols; ++i)
                    for (unsigned j = 0; j < d_rows; ++j)
                        res(i, j) = (*this)(j, i);
                return res;
            }
        };

        template<class T, class ST = storage_traits<T> >
        class math_matrix : public base_matrix<T, ST>
        {
        public:
            enum { has_direct_access = true };
    
            // Accessors for rows and columns (constant and non-constant); these behave similar to a
            // math_vector<>.
    
            typedef MatrixRow<T, ST> Row;
            typedef MatrixColumn<T, ST> Column;
            typedef MatrixCRow<T, ST> CRow;
            typedef MatrixCColumn<T, ST> CColumn;
    
            // The matrix class
    
            inline math_matrix()
            {
            }
    
            inline math_matrix(typename base_matrix<T, ST>::size_type rows, typename base_matrix<T, ST>::size_type cols)
                : base_matrix<T, ST>(rows, cols)
            {
            }
    
            inline math_matrix(typename base_matrix<T, ST>::size_type rows, typename base_matrix<T, ST>::size_type cols, typename ST::constref_type obj)
                : base_matrix<T, ST>(rows, cols, obj)
            {
            }
    
            inline math_matrix(const math_matrix<T, ST> & m)
                : base_matrix<T, ST>(m)
            {
            }
    
            inline Row row(typename base_matrix<T, ST>::size_type i)
            {
                return Row(*this, i);
            }
    
            inline CRow row(typename base_matrix<T, ST>::size_type i) const
            {
                return CRow(*this, i);
            }
    
            inline Column col(typename base_matrix<T, ST>::size_type i)
            {
                return Column(*this, i);
            }
    
            inline CColumn col(typename base_matrix<T, ST>::size_type i) const
            {
                return CColumn(*this, i);
            }
    
            template<class S, class A>
            inline explicit math_matrix(const base_matrix<S, A> & m)
                : base_matrix<T, ST>(m)
            {
            }
    
            template<class S, class A>
            inline math_matrix<T, ST> & operator = (const base_matrix<S, A> & m)
            {
                linalgold::assign(*this, m);
                return *this;
            }
    
            inline math_matrix<T, ST> & operator %= (const T & s)
            {
                mod(*this, *this, s);
                return *this;
            }
    
            inline math_matrix<T, ST> & operator /= (const T & s)
            {
                div(*this, *this, s);
                return *this;
            }
    
            inline math_matrix<T, ST> & operator *= (const T & s)
            {
                mul(*this, *this, s);
                return *this;
            }
    
            template<class T2, class A2>
            inline math_matrix<T, ST> & operator += (const math_matrix<T2, A2> & v)
            {
                add(*this, *this, v);
                return *this;
            }
    
            template<class T2, class A2>
            inline math_matrix<T, ST> & operator -= (const math_matrix<T2, A2> & v)
            {
                sub(*this, *this, v);
                return *this;
            }
    
            template<class S, class A>
            inline bool operator < (const math_matrix<S, A> & v) const
            {
                if ((v.rows() != this->d_rows) || (v.cols() != this->d_cols))
                    return false;
                typename base_matrix<T, ST>::size_type c = this->d_rows * this->d_cols;
                for (typename base_matrix<T, ST>::size_type i = 0; i < c; ++i)
                {
                    if (this->d_data[i] < v.data()[i])
                        return true;
                    if (this->d_data[i] > v.data()[i])
                        return false;
                }
                return false;
            }
    
            template<class S, class A>
            inline bool operator <= (const math_matrix<S, A> & v) const
            {
                if ((v.rows() != this->d_rows) || (v.cols() != this->d_cols))
                    return false;
                typename base_matrix<T, ST>::size_type c = this->d_rows * this->d_cols;
                for (typename base_matrix<T, ST>::size_type i = 0; i < c; ++i)
                {
                    if (this->d_data[i] < v.data()[i])
                        return true;
                    if (this->d_data[i] > v.data()[i])
                        return false;
                }
                return true;
            }
    
            template<class S, class A>
            inline bool operator > (const math_matrix<S, A> & v) const
            {
                return !(*this <= v);
            }
    
            template<class S, class A>
            inline bool operator >= (const math_matrix<S, A> & v) const
            {
                return !(*this < v);
            }
    
            inline math_matrix<T, ST> operator % (const T & s) const
            {
                math_matrix<T, ST> ret(this->d_rows, this->d_cols);
                mod(ret, *this, s);
                return ret;
            }
    
            inline math_matrix<T, ST> operator / (const T & s) const
            {
                math_matrix<T, ST> ret(this->d_rows, this->d_cols);
                div(ret, *this, s);
                return ret;
            }
    
            inline math_matrix<T, ST> operator * (const T & s) const
            {
                math_matrix<T, ST> ret(this->d_rows, this->d_cols);
                mul(ret, *this, s);
                return ret;
            }

#ifndef PLLL_INTERNAL_NO_TEMPLATE_FRIENDS
            template<class TT, class STT>
            friend math_matrix<TT, STT> operator * (const TT & s, const math_matrix<TT, STT> & m);
#endif
    
            template<class T2, class ST2>
            inline math_matrix<T, ST> operator + (const math_matrix<T2, ST2> & m2) const
            {
                math_matrix<T, ST> ret(this->d_rows, this->d_cols);
                add(ret, *this, m2);
                return ret;
            }
    
            template<class T2, class ST2>
            inline math_matrix<T, ST> operator - (const math_matrix<T2, ST2> & m2) const
            {
                math_matrix<T, ST> ret(this->d_rows, this->d_cols);
                sub(ret, *this, m2);
                return ret;
            }
    
            inline math_matrix<T, ST> operator - () const
            {
                math_matrix<T, ST> ret(this->d_rows, this->d_cols);
                neg(ret, *this);
                return ret;
            }
    
            inline math_matrix<T, ST> transpose() const
            {
                math_matrix<T, ST> res(this->d_cols, this->d_rows);
                linalgold::transpose<T, ST, T, ST>(res, *this); // we know that &res != this, hence we call the more "general"
                // function which assumes that source and destination are not the same.
                return res;
            }
    
            // Matrix-vector and matrix-matrix multiplications
    
            template<class T2, class ST2>
            inline math_matrix<T, ST> & operator *= (const math_matrix<T2, ST2> & m2)
            {
                mul(*this, *this, m2);
                return *this;
            }

#ifndef PLLL_INTERNAL_NO_TEMPLATE_FRIENDS
            template<class TT, class STT, class T2, class ST2>
            friend math_matrix<TT, STT> operator * (const math_matrix<TT, STT> & m1, const math_matrix<T2, ST2> & m2);
    
            template<class TT, class STT, class T2, class ST2>
            friend math_vector<TT, STT> operator * (const math_matrix<TT, STT> & A, const math_vector<T2, ST2> & v);
    
            template<class TT, class STT, class T2, class ST2>
            friend math_vector<TT, STT> operator * (const math_matrix<TT, STT> & A, const MatrixRow<T2, ST2> & v);
    
            template<class TT, class STT, class T2, class ST2>
            friend math_vector<TT, STT> operator * (const math_matrix<TT, STT> & A, const MatrixCRow<T2, ST2> & v);
    
            template<class TT, class STT, class T2, class ST2>
            friend math_vector<TT, STT> operator * (const math_matrix<TT, STT> & A, const MatrixColumn<T2, ST2> & v);
    
            template<class TT, class STT, class T2, class ST2>
            friend math_vector<TT, STT> operator * (const math_matrix<TT, STT> & A, const MatrixCColumn<T2, ST2> & v);
    
            template<class TT, class STT, class T2, class ST2>
            friend math_vector<TT, STT> operator * (const math_vector<T2, ST2> & v, const math_matrix<TT, STT> & A);
    
            template<class TT, class STT, class T2, class ST2>
            friend math_vector<TT, STT> operator * (const MatrixRow<T2, ST2> & v, const math_matrix<TT, STT> & A);
    
            template<class TT, class STT, class T2, class ST2>
            friend math_vector<TT, STT> operator * (const MatrixCRow<T2, ST2> & v, const math_matrix<TT, STT> & A);
    
            template<class TT, class STT, class T2, class ST2>
            friend math_vector<TT, STT> operator * (const MatrixColumn<T2, ST2> & v, const math_matrix<TT, STT> & A);
    
            template<class TT, class STT, class T2, class ST2>
            friend math_vector<TT, STT> operator * (const MatrixCColumn<T2, ST2> & v, const math_matrix<TT, STT> & A);
    
            // Temporary-less operations

            template<class TT, class STT, class T2, class ST2>
            friend void mod(math_matrix<TT, STT> & r, const math_matrix<T2, ST2> & a, const TT & b);
            template<class TT, class STT, class T2, class ST2>
            friend void div(math_matrix<TT, STT> & r, const math_matrix<T2, ST2> & a, const TT & b);
            template<class TT, class STT, class T2, class ST2>
            friend void mul(math_matrix<TT, STT> & r, const math_matrix<T2, ST2> & a, const TT & b);
            template<class TT, class STT, class T2, class ST2, class T3, class ST3>
            friend void add(math_matrix<TT, STT> & r, const math_matrix<T2, ST2> & a, const math_matrix<T3, ST3> & b);
            template<class TT, class STT, class T2, class ST2, class T3, class ST3>
            friend void sub(math_matrix<TT, STT> & r, const math_matrix<T2, ST2> & a, const math_matrix<T3, ST3> & b);
            template<class TT, class STT, class T2, class ST2>
            friend void neg(math_matrix<TT, STT> & r, const math_matrix<T2, ST2> & a);
            template<class TT, class STT, class T2, class ST2>
            friend void linalgold::transpose(math_matrix<TT, STT> & r, const math_matrix<T2, ST2> & a);
            template<class TT, class STT>
            friend void linalgold::transpose(math_matrix<TT, STT> & r, const math_matrix<TT, STT> & a);
    
            // Temporary-less matrix-vector and matrix-matrix multiplications
    
            template<class TT, class STT, class T2, class ST2, class V>
            friend void mul2(math_vector<TT, STT> & r, const math_matrix<T2, ST2> & a, const V & b);
    
            template<class TT, class STT, class T2, class ST2, class V>
            friend void mul2(math_vector<TT, STT> & r, const V & a, const math_matrix<T2, ST2> & b);
    
            template<class TT, class STT, class T2, class ST2, class V>
            friend void mul2(const MatrixRow<T2, ST2> & r, const math_matrix<TT, STT> & a, const V & b);
    
            template<class TT, class STT, class T2, class ST2, class V>
            friend void mul2(const MatrixRow<T2, ST2> & r, const V & a, const math_matrix<TT, STT> & b);
    
            template<class TT, class STT, class T2, class ST2, class V>
            friend void mul2(const MatrixColumn<T2, ST2> & r, const math_matrix<TT, STT> & a, const V & b);
    
            template<class TT, class STT, class T2, class ST2, class V>
            friend void mul2(const MatrixColumn<T2, ST2> & r, const V & a, const math_matrix<TT, STT> & b);
    
            template<class TT, class STT, class T2, class ST2, class T3, class ST3>
            friend void mul2(math_matrix<TT, STT> & r, const math_matrix<T2, ST2> & a, const base_matrix<T3, ST3> & b);

            // First argument: math_vector<>; matrix times vector
            template<class TT, class STT, class T2, class ST2, class T3, class ST3>
            friend void mul(math_vector<TT, STT> & r, const math_matrix<T3, ST3> & a, const math_vector<T2, ST2> & b);
            template<class TT, class STT, class T2, class ST2, class T3, class ST3>
            friend void mul(math_vector<TT, STT> & r, const math_matrix<T3, ST3> & a, const MatrixRow<T2, ST2> & b);
            template<class TT, class STT, class T2, class ST2, class T3, class ST3>
            friend void mul(math_vector<TT, STT> & r, const math_matrix<T3, ST3> & a, const MatrixCRow<T2, ST2> & b);
            template<class TT, class STT, class T2, class ST2, class T3, class ST3>
            friend void mul(math_vector<TT, STT> & r, const math_matrix<T3, ST3> & a, const MatrixColumn<T2, ST2> & b);
            template<class TT, class STT, class T2, class ST2, class T3, class ST3>
            friend void mul(math_vector<TT, STT> & r, const math_matrix<T3, ST3> & a, const MatrixCColumn<T2, ST2> & b);
    
            // First argument: math_vector<>; vector times matrix
            template<class TT, class STT, class T2, class ST2, class T3, class ST3>
            friend void mul(math_vector<TT, STT> & r, const math_vector<T2, ST2> & a, const math_matrix<T3, ST3> & b);
            template<class TT, class STT, class T2, class ST2, class T3, class ST3>
            friend void mul(math_vector<TT, STT> & r, const MatrixRow<T2, ST2> & a, const math_matrix<T3, ST3> & b);
            template<class TT, class STT, class T2, class ST2, class T3, class ST3>
            friend void mul(math_vector<TT, STT> & r, const MatrixCRow<T2, ST2> & a, const math_matrix<T3, ST3> & b);
            template<class TT, class STT, class T2, class ST2, class T3, class ST3>
            friend void mul(math_vector<TT, STT> & r, const MatrixColumn<T2, ST2> & a, const math_matrix<T3, ST3> & b);
            template<class TT, class STT, class T2, class ST2, class T3, class ST3>
            friend void mul(math_vector<TT, STT> & r, const MatrixCColumn<T2, ST2> & a, const math_matrix<T3, ST3> & b);
    
            // First argument: MatrixRow<>; matrix times vector
            template<class TT, class STT, class T2, class ST2, class T3, class ST3>
            friend void mul(const MatrixRow<T2, ST2> & r, const math_matrix<TT, STT> & a, const math_vector<T3, ST3> & b);
            template<class TT, class STT, class T2, class ST2, class T3, class ST3>
            friend void mul(const MatrixRow<T2, ST2> & r, const math_matrix<TT, STT> & a, const MatrixRow<T3, ST3> & b);
            template<class TT, class STT, class T2, class ST2, class T3, class ST3>
            friend void mul(const MatrixRow<T2, ST2> & r, const math_matrix<TT, STT> & a, const MatrixCRow<T3, ST3> & b);
            template<class TT, class STT, class T2, class ST2, class T3, class ST3>
            friend void mul(const MatrixRow<T2, ST2> & r, const math_matrix<TT, STT> & a, const MatrixColumn<T3, ST3> & b);
            template<class TT, class STT, class T2, class ST2, class T3, class ST3>
            friend void mul(const MatrixRow<T2, ST2> & r, const math_matrix<TT, STT> & a, const MatrixCColumn<T3, ST3> & b);
    
            // First argument: MatrixRow<>; vector times matrix
            template<class TT, class STT, class T2, class ST2, class T3, class ST3>
            friend void mul(const MatrixRow<T2, ST2> & r, const math_vector<T3, ST3> & a, const math_matrix<TT, STT> & b);
            template<class TT, class STT, class T2, class ST2, class T3, class ST3>
            friend void mul(const MatrixRow<T2, ST2> & r, const MatrixRow<T3, ST3> & a, const math_matrix<TT, STT> & b);
            template<class TT, class STT, class T2, class ST2, class T3, class ST3>
            friend void mul(const MatrixRow<T2, ST2> & r, const MatrixCRow<T3, ST3> & a, const math_matrix<TT, STT> & b);
            template<class TT, class STT, class T2, class ST2, class T3, class ST3>
            friend void mul(const MatrixRow<T2, ST2> & r, const MatrixColumn<T3, ST3> & a, const math_matrix<TT, STT> & b);
            template<class TT, class STT, class T2, class ST2, class T3, class ST3>
            friend void mul(const MatrixRow<T2, ST2> & r, const MatrixCColumn<T3, ST3> & a, const math_matrix<TT, STT> & b);
    
            // First argument: MatrixColumn<>; matrix times vector
            template<class TT, class STT, class T2, class ST2, class T3, class ST3>
            friend void mul(const MatrixColumn<T2, ST2> & r, const math_matrix<TT, STT> & a, const math_vector<T3, ST3> & b);
            template<class TT, class STT, class T2, class ST2, class T3, class ST3>
            friend void mul(const MatrixColumn<T2, ST2> & r, const math_matrix<TT, STT> & a, const MatrixRow<T3, ST3> & b);
            template<class TT, class STT, class T2, class ST2, class T3, class ST3>
            friend void mul(const MatrixColumn<T2, ST2> & r, const math_matrix<TT, STT> & a, const MatrixCRow<T3, ST3> & b);
            template<class TT, class STT, class T2, class ST2, class T3, class ST3>
            friend void mul(const MatrixColumn<T2, ST2> & r, const math_matrix<TT, STT> & a, const MatrixColumn<T3, ST3> & b);
            template<class TT, class STT, class T2, class ST2, class T3, class ST3>
            friend void mul(const MatrixColumn<T2, ST2> & r, const math_matrix<TT, STT> & a, const MatrixCColumn<T3, ST3> & b);
    
            // First argument: MatrixColumn<>; vector times matrix
            template<class TT, class STT, class T2, class ST2, class T3, class ST3>
            friend void mul(const MatrixColumn<T2, ST2> & r, const math_vector<T3, ST3> & a, const math_matrix<TT, STT> & b);
            template<class TT, class STT, class T2, class ST2, class T3, class ST3>
            friend void mul(const MatrixColumn<T2, ST2> & r, const MatrixRow<T3, ST3> & a, const math_matrix<TT, STT> & b);
            template<class TT, class STT, class T2, class ST2, class T3, class ST3>
            friend void mul(const MatrixColumn<T2, ST2> & r, const MatrixCRow<T3, ST3> & a, const math_matrix<TT, STT> & b);
            template<class TT, class STT, class T2, class ST2, class T3, class ST3>
            friend void mul(const MatrixColumn<T2, ST2> & r, const MatrixColumn<T3, ST3> & a, const math_matrix<TT, STT> & b);
            template<class TT, class STT, class T2, class ST2, class T3, class ST3>
            friend void mul(const MatrixColumn<T2, ST2> & r, const MatrixCColumn<T3, ST3> & a, const math_matrix<TT, STT> & b);
    
            // First argument: math_matrix<>; matrix times matrix
            template<class TT, class STT, class T2, class ST2, class T3, class ST3>
            friend void mul(math_matrix<TT, STT> & r, const math_matrix<T2, ST2> & a, const math_matrix<T3, ST3> & b);
#endif
        };
        
        template<class T, class ST, class V>
        void swap(const MatrixRow<T, ST> & v1, V & v2);
        template<class T, class ST, class V>
        void swap(const MatrixRow<T, ST> & v1, const V & v2);

        template<class T, class ST, class V>
        void assign(const MatrixRow<T, ST> & R, const V & v);
        template<class T, class ST, class V>
        void mod(const MatrixRow<T, ST> & r, const V & a, const T & b);
        template<class T, class ST, class V>
        void div(const MatrixRow<T, ST> & r, const V & a, const T & b);
        template<class T, class ST, class V>
        void mul(const MatrixRow<T, ST> & r, const V & a, const T & b);
        template<class T, class ST, class V2, class V3>
        void add(const MatrixRow<T, ST> & r, const V2 & a, const V3 & b);
        template<class T, class ST, class V2, class V3>
        void sub(const MatrixRow<T, ST> & r, const V2 & a, const V3 & b);
        template<class T, class ST, class V2, class TT>
        void addmul(const MatrixRow<T, ST> & r, const V2 & a, const TT & b);
        template<class T, class ST, class V2, class TT>
        void submul(const MatrixRow<T, ST> & r, const V2 & a, const TT & b);
        template<class T, class ST, class V>
        void neg(const MatrixRow<T, ST> & r, const V & a);
        
        template<class T, class ST>
        void normSq(T & r, const MatrixRow<T, ST> & a);
        template<class T, class ST, class V>
        void dot(T & r, const MatrixRow<T, ST> & a, const V & b);

        template<class T, class ST, class V>
        void assign(const MatrixColumn<T, ST> & C, const V & v);
        template<class T, class ST, class V>
        void mod(const MatrixColumn<T, ST> & r, const V & a, const T & b);
        template<class T, class ST, class V>
        void div(const MatrixColumn<T, ST> & r, const V & a, const T & b);
        template<class T, class ST, class V>
        void mul(const MatrixColumn<T, ST> & r, const V & a, const T & b);
        template<class T, class ST, class V2, class V3>
        void add(const MatrixColumn<T, ST> & r, const V2 & a, const V3 & b);
        template<class T, class ST, class V2, class V3>
        void sub(const MatrixColumn<T, ST> & r, const V2 & a, const V3 & b);
        template<class T, class ST, class V2, class TT>
        void addmul(const MatrixColumn<T, ST> & r, const V2 & a, const TT & b);
        template<class T, class ST, class V2, class TT>
        void submul(const MatrixColumn<T, ST> & r, const V2 & a, const TT & b);
        template<class T, class ST, class V>
        void neg(const MatrixColumn<T, ST> & r, const V & a);
        
        template<class T, class ST, class V>
        void swap(const MatrixColumn<T, ST> & v1, V & v2);
        template<class T, class ST, class V>
        void swap(const MatrixColumn<T, ST> & v1, const V & v2);
        
        template<class T, class ST>
        void normSq(T & r, const MatrixColumn<T, ST> & a);
        template<class T, class ST, class V>
        void dot(T & r, const MatrixColumn<T, ST> & a, const V & b);

        template<class T, class ST>
        void normSq(T & r, const MatrixCRow<T, ST> & a);
        template<class T, class ST, class V>
        void dot(T & r, const MatrixCRow<T, ST> & a, const V & b);

        template<class T, class ST>
        void normSq(T & r, const MatrixCColumn<T, ST> & a);
        template<class T, class ST, class V>
        void dot(T & r, const MatrixCColumn<T, ST> & a, const V & b);

        template<class T, class ST>
        class MatrixRow
        {
            friend class math_matrix<T, ST>;

#ifdef PLLL_INTERNAL_NO_TEMPLATE_FRIENDS
        public:
#else
        private:
#endif
            math_matrix<T, ST> & d_m;
            typename base_matrix<T, ST>::size_type d_i;
    
            inline MatrixRow(math_matrix<T, ST> & m, typename base_matrix<T, ST>::size_type i)
                : d_m(m), d_i(i)
            {
            }
    
        public:
            enum { has_direct_access = true };
    
            inline operator MatrixRow() const
            {
                return MatrixRow(d_m, d_i);
            }
    
            inline math_matrix<T, ST> & matrix() const
            {
                return d_m;
            }
    
            inline typename ST::pointer_type data() const
            {
                return d_m.data() + d_i * d_m.cols();
            }
    
            inline typename base_matrix<T, ST>::size_type row() const
            {
                return d_i;
            }
    
            inline typename base_matrix<T, ST>::size_type size() const
            {
                return d_m.cols();
            }
    
            inline typename ST::ref_type operator [] (typename base_matrix<T, ST>::size_type j) const
            {
                return d_m(d_i, j);
            }

#ifndef PLLL_INTERNAL_NO_TEMPLATE_FRIENDS
            template<class TT, class STT, class V>
            friend void assign(const MatrixRow<TT, STT> & R, const V & v);
#endif
    
            template<class V>
            inline void assign(const V & v) const
            {
                linalgold::assign(*this, v);
            }
    
            template<class O1, class D1, class T1>
            inline void assign(const VecExpression<O1, D1, T1> & e)
            {
                assert(size() == e.size());
                if (e.usesVecNF(*this))
                {
                    T tmp;
                    for (typename base_matrix<T, ST>::size_type i = 0; i < size(); ++i)
                    {
                        e.evalComponentTo(tmp, i);
                        (*this)[i] = tmp;
                    }
                }
                else
                {
                    for (typename base_matrix<T, ST>::size_type i = 0; i < size(); ++i)
                        e.evalComponentTo((*this)[i], i);
                }
            }
    
            inline const MatrixRow & operator = (const MatrixRow & v) const
            {
                assign(v);
                return *this;
            }
    
            template<class V>
            inline const MatrixRow & operator = (const V & v) const
            {
                assign(v);
                return *this;
            }
    
            inline const MatrixRow & operator %= (const T & s) const
            {
                mod(*this, *this, s);
                return *this;
            }
        
            inline const MatrixRow & operator /= (const T & s) const
            {
                div(*this, *this, s);
                return *this;
            }
        
            inline const MatrixRow & operator *= (const T & s) const
            {
                mul(*this, *this, s);
                return *this;
            }
        
            template<class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            inline const MatrixRow & operator += (const V & v) const
            {
                add(*this, *this, v);
                return *this;
            }
        
            template<class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            inline const MatrixRow & operator -= (const V & v) const
            {
                sub(*this, *this, v);
                return *this;
            }
    
            template<class O1, class D1, class T1>
            inline const MatrixRow & operator += (const VecExpression<O1, D1, T1> & e) const
            {
                assert(size() == e.size());
                if (e.usesVecNF(*this))
                {
                    T tmp;
                    for (unsigned i = 0; i < size(); ++i)
                    {
                        e.evalComponentTo(tmp, i);
                        (*this)[i] += tmp;
                    }
                }
                else
                {
                    for (unsigned i = 0; i < size(); ++i)
                        e.evalComponentAddTo((*this)[i], i);
                }
                return *this;
            }
    
            template<class O1, class D1, class T1>
            inline const MatrixRow & operator -= (const VecExpression<O1, D1, T1> & e) const
            {
                assert(size() == e.size());
                if (e.usesVecNF(*this))
                {
                    T tmp;
                    for (unsigned i = 0; i < size(); ++i)
                    {
                        e.evalComponentTo(tmp, i);
                        (*this)[i] -= tmp;
                    }
                }
                else
                {
                    for (unsigned i = 0; i < size(); ++i)
                        e.evalComponentSubFrom((*this)[i], i);
                }
                return *this;
            }
        
            template<class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            inline bool operator < (const V & v) const
            {
                if (v.size() != this->size())
                    return false;
                for (typename base_vector<T, ST>::size_type i = 0; i < this->size(); ++i)
                {
                    if ((*this)[i] < v[i])
                        return true;
                    if ((*this)[i] > v[i])
                        return false;
                }
                return false;
            }
        
            template<class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            inline bool operator <= (const V & v) const
            {
                if (v.size() != this->size())
                    return false;
                for (typename base_vector<T, ST>::size_type i = 0; i < this->size(); ++i)
                {
                    if ((*this)[i] < v[i])
                        return true;
                    if ((*this)[i] > v[i])
                        return false;
                }
                return true;
            }
        
            template<class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            inline bool operator > (const V & v) const
            {
                return !(*this <= v);
            }
        
            template<class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            inline bool operator >= (const V & v) const
            {
                return !(*this < v);
            }
        
            inline T normSq() const
            {
                T res;
                linalgold::normSq<T, ST>(res, *this);
                return res;
            }

#ifndef PLLL_INTERNAL_NO_TEMPLATE_FRIENDS
            template<class TT, class STT, class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            friend void mod(const MatrixRow<TT, STT> & r, const V & a, const TT & b);
            template<class TT, class STT, class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            friend void div(const MatrixRow<TT, STT> & r, const V & a, const TT & b);
            template<class TT, class STT, class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            friend void mul(const MatrixRow<TT, STT> & r, const V & a, const TT & b);
            template<class TT, class STT, class V2, class V3> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            friend void add(const MatrixRow<TT, STT> & r, const V2 & a, const V3 & b);
            template<class TT, class STT, class V2, class V3> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            friend void sub(const MatrixRow<TT, STT> & r, const V2 & a, const V3 & b);
            template<class TT, class STT, class V2, class TTT> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            friend void addmul(const MatrixRow<TT, STT> & r, const V2 & a, const TTT & b);
            template<class TT, class STT, class V2, class TTT> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            friend void submul(const MatrixRow<TT, STT> & r, const V2 & a, const TTT & b);
            template<class TT, class STT, class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            friend void neg(const MatrixRow<TT, STT> & r, const V & a);
            template<class TT, class STT, class V>
            friend void linalgold::swap(const MatrixRow<TT, STT> & v1, V & v2);
            template<class TT, class STT, class V>
            friend void linalgold::swap(const MatrixRow<TT, STT> & v1, const V & v2);
            template<class TT, class STT>
            friend void linalgold::normSq(TT & r, const MatrixRow<TT, STT> & a);
            template<class TT, class STT, class V>
            friend void dot(TT & r, const MatrixRow<TT, STT> & a, const V & b);
#endif
        };
    
        template<class T, class ST>
        class MatrixColumn
        {
            friend class math_matrix<T, ST>;

#ifdef PLLL_INTERNAL_NO_TEMPLATE_FRIENDS
        public:
#else
        private:
#endif
            math_matrix<T, ST> & d_m;
            typename base_matrix<T, ST>::size_type d_i;
        
            inline MatrixColumn(math_matrix<T, ST> & m, typename base_matrix<T, ST>::size_type i)
                : d_m(m), d_i(i)
            {
            }

        public:
            enum { has_direct_access = false };
    
            inline operator MatrixColumn() const
            {
                return MatrixColumn(d_m, d_i);
            }
        
            inline math_matrix<T, ST> & matrix() const
            {
                return d_m;
            }
        
            inline typename base_matrix<T, ST>::size_type column() const
            {
                return d_i;
            }
        
            inline typename base_matrix<T, ST>::size_type size() const
            {
                return d_m.rows();
            }
        
            inline typename ST::ref_type operator [] (typename base_matrix<T, ST>::size_type j) const
            {
                return d_m(j, d_i);
            }

#ifndef PLLL_INTERNAL_NO_TEMPLATE_FRIENDS
            template<class TT, class STT, class V>
            friend void assign(const MatrixColumn<TT, STT> & C, const V & v);
#endif
    
            template<class V>
            inline void assign(const V & v) const
            {
                linalgold::assign(*this, v);
            }
    
            template<class O1, class D1, class T1>
            inline void assign(const VecExpression<O1, D1, T1> & e)
            {
                assert(size() == e.size());
                if (e.usesVecNF(*this))
                {
                    T tmp;
                    for (typename base_matrix<T, ST>::size_type i = 0; i < size(); ++i)
                    {
                        e.evalComponentTo(tmp, i);
                        (*this)[i] = tmp;
                    }
                }
                else
                {
                    for (typename base_matrix<T, ST>::size_type i = 0; i < size(); ++i)
                        e.evalComponentTo((*this)[i], i);
                }
            }
    
            inline const MatrixColumn & operator = (const MatrixColumn & v) const
            {
                assign(v);
                return *this;
            }
    
            template<class V>
            inline const MatrixColumn & operator = (const V & v) const
            {
                assign(v);
                return *this;
            }

            inline const MatrixColumn & operator %= (const T & s) const
            {
                mod(*this, *this, s);
                return *this;
            }
        
            inline const MatrixColumn & operator /= (const T & s) const
            {
                div(*this, *this, s);
                return *this;
            }
        
            inline const MatrixColumn & operator *= (const T & s) const
            {
                mul(*this, *this, s);
                return *this;
            }
        
            template<class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            inline const MatrixColumn & operator += (const V & v) const
            {
                add(*this, *this, v);
                return *this;
            }
        
            template<class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            inline const MatrixColumn & operator -= (const V & v) const
            {
                sub(*this, *this, v);
                return *this;
            }
    
            template<class O1, class D1, class T1>
            inline const MatrixColumn & operator += (const VecExpression<O1, D1, T1> & e) const
            {
                assert(size() == e.size());
                if (e.usesVecNF(*this))
                {
                    T tmp;
                    for (unsigned i = 0; i < size(); ++i)
                    {
                        e.evalComponentTo(tmp, i);
                        (*this)[i] += tmp;
                    }
                }
                else
                {
                    for (unsigned i = 0; i < size(); ++i)
                        e.evalComponentAddTo((*this)[i], i);
                }
                return *this;
            }
    
            template<class O1, class D1, class T1>
            inline const MatrixColumn & operator -= (const VecExpression<O1, D1, T1> & e) const
            {
                assert(size() == e.size());
                if (e.usesVecNF(*this))
                {
                    T tmp;
                    for (unsigned i = 0; i < size(); ++i)
                    {
                        e.evalComponentTo(tmp, i);
                        (*this)[i] -= tmp;
                    }
                }
                else
                {
                    for (unsigned i = 0; i < size(); ++i)
                        e.evalComponentSubFrom((*this)[i], i);
                }
                return *this;
            }
    
            template<class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            inline bool operator < (const V & v) const
            {
                if (v.size() != this->size())
                    return false;
                for (typename base_vector<T, ST>::size_type i = 0; i < this->size(); ++i)
                {
                    if ((*this)[i] < v[i])
                        return true;
                    if ((*this)[i] > v[i])
                        return false;
                }
                return false;
            }
    
            template<class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            inline bool operator <= (const V & v) const
            {
                if (v.size() != this->size())
                    return false;
                for (typename base_vector<T, ST>::size_type i = 0; i < this->size(); ++i)
                {
                    if ((*this)[i] < v[i])
                        return true;
                    if ((*this)[i] > v[i])
                        return false;
                }
                return true;
            }
        
            template<class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            inline bool operator > (const V & v) const
            {
                return !(*this <= v);
            }
        
            template<class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            inline bool operator >= (const V & v) const
            {
                return !(*this < v);
            }
        
            inline T normSq() const
            {
                T res;
                linalgold::normSq<T, ST>(res, *this);
                return res;
            }

#ifndef PLLL_INTERNAL_NO_TEMPLATE_FRIENDS
            template<class TT, class STT, class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            friend void mod(const MatrixColumn<TT, STT> & r, const V & a, const TT & b);
            template<class TT, class STT, class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            friend void div(const MatrixColumn<TT, STT> & r, const V & a, const TT & b);
            template<class TT, class STT, class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            friend void mul(const MatrixColumn<TT, STT> & r, const V & a, const TT & b);
            template<class TT, class STT, class V2, class V3> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            friend void add(const MatrixColumn<TT, STT> & r, const V2 & a, const V3 & b);
            template<class TT, class STT, class V2, class V3> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            friend void sub(const MatrixColumn<TT, STT> & r, const V2 & a, const V3 & b);
            template<class TT, class STT, class V2, class TTT> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            friend void addmul(const MatrixColumn<TT, STT> & r, const V2 & a, const TTT & b);
            template<class TT, class STT, class V2, class TTT> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            friend void submul(const MatrixColumn<TT, STT> & r, const V2 & a, const TTT & b);
            template<class TT, class STT, class V> // this can be math_vector<>, base_vector<>, or something else which behaves like a vector
            friend void neg(const MatrixColumn<TT, STT> & r, const V & a);
            template<class TT, class STT, class V>
            friend void linalgold::swap(const MatrixColumn<TT, STT> & v1, V & v2);
            template<class TT, class STT, class V>
            friend void linalgold::swap(const MatrixColumn<TT, STT> & v1, const V & v2);
            template<class TT, class STT>
            friend void linalgold::normSq(TT & r, const MatrixColumn<TT, STT> & a);
            template<class TT, class STT, class V>
            friend void dot(TT & r, const MatrixColumn<TT, STT> & a, const V & b);
#endif
        };

// Constant accessors for rows and columns (constant and non-constant); these behave similar to
// a math_vector<>.
    
        template<class T, class ST>
        class MatrixCRow
        {
            friend class math_matrix<T, ST>;

#ifdef PLLL_INTERNAL_NO_TEMPLATE_FRIENDS
        public:
#else
        private:
#endif
            const math_matrix<T, ST> & d_m;
            typename base_matrix<T, ST>::size_type d_i;
        
            inline MatrixCRow(const math_matrix<T, ST> & m, typename base_matrix<T, ST>::size_type i)
                : d_m(m), d_i(i)
            {
            }

        public:
            enum { has_direct_access = true };
    
            inline const math_matrix<T, ST> & matrix() const
            {
                return d_m;
            }
    
            inline typename ST::pointer_type data() const
            {
                return d_m.data() + d_i * d_m.cols();
            }
    
            inline typename base_matrix<T, ST>::size_type row() const
            {
                return d_i;
            }
        
            inline typename base_matrix<T, ST>::size_type size() const
            {
                return d_m.cols();
            }
        
            inline typename ST::constref_type operator [] (typename base_matrix<T, ST>::size_type j) const
            {
                return d_m(d_i, j);
            }
    
            inline T normSq() const
            {
                T res;
                linalgold::normSq<T, ST>(res, *this);
                return res;
            }

#ifndef PLLL_INTERNAL_NO_TEMPLATE_FRIENDS
            template<class TT, class STT>
            friend void linalgold::normSq(TT & r, const typename math_matrix<TT, STT>::CRow & a);
            template<class TT, class STT, class V>
            friend void dot(TT & r, const typename math_matrix<TT, STT>::CRow & a, const V & b);
#endif
        };
    
        template<class T, class ST>
        class MatrixCColumn
        {
            friend class math_matrix<T, ST>;

#ifdef PLLL_INTERNAL_NO_TEMPLATE_FRIENDS
        public:
#else
        private:
#endif
            const math_matrix<T, ST> & d_m;
            typename base_matrix<T, ST>::size_type d_i;
        
            inline MatrixCColumn(const math_matrix<T, ST> & m, typename base_matrix<T, ST>::size_type i)
                : d_m(m), d_i(i)
            {
            }

        public:
            enum { has_direct_access = false };
    
            inline const math_matrix<T, ST> & matrix() const
            {
                return d_m;
            }
        
            inline typename base_matrix<T, ST>::size_type column() const
            {
                return d_i;
            }
        
            inline typename base_matrix<T, ST>::size_type size() const
            {
                return d_m.rows();
            }
        
            inline typename ST::constref_type operator [] (typename base_matrix<T, ST>::size_type j) const
            {
                return d_m(j, d_i);
            }

            inline T normSq() const
            {
                T res;
                linalgold::normSq<T, ST>(res, *this);
                return res;
            }

#ifndef PLLL_INTERNAL_NO_TEMPLATE_FRIENDS
            template<class TT, class STT>
            friend void linalgold::normSq(TT & r, const MatrixCColumn<TT, STT> & a);
            template<class TT, class STT, class V>
            friend void dot(TT & r, const MatrixCColumn<TT, STT> & a, const V & b);
#endif
        };
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Stream output

#include <iostream> // apparently iosfwd is not enough here...

namespace plll
{
    namespace linalgold
    {
        template<class T, class ST>
        std::ostream & operator << (std::ostream & s, const base_vector<T, ST> & v)
        {
            s << "vec<" << v.size() << ">[";
            for (typename base_vector<T, ST>::size_type i = 0; i < v.size(); ++i)
            {
                if (i > 0)
                    s << ", ";
                s << v[i];
            }
            return s << "]";
        }

        template<class T, class ST>
        std::ostream & operator << (std::ostream & s, const base_matrix<T, ST> & v)
        {
            s << "mat<" << v.rows() << "," << v.cols() << ">[";
            for (typename base_matrix<T, ST>::size_type i = 0; i < v.rows(); ++i)
            {
                if (i > 0)
                    s << ", ";
                s << "[";
                for (typename base_matrix<T, ST>::size_type j = 0; j < v.cols(); ++j)
                {
                    if (j > 0)
                        s << ", ";
                    s << v(i, j);
                }
                s << "]";
            }
            return s << "]";
        }

        template<class T, class ST>
        std::ostream & operator << (std::ostream & s, const MatrixRow<T, ST> & v)
        {
            s << "vec<" << v.size() << ">[";
            for (typename base_vector<T, ST>::size_type i = 0; i < v.size(); ++i)
            {
                if (i > 0)
                    s << ", ";
                s << v[i];
            }
            return s << "]";
        }

        template<class T, class ST>
        std::ostream & operator << (std::ostream & s, const MatrixColumn<T, ST> & v)
        {
            s << "vec<" << v.size() << ">[";
            for (typename base_vector<T, ST>::size_type i = 0; i < v.size(); ++i)
            {
                if (i > 0)
                    s << ", ";
                s << v[i];
            }
            return s << "]";
        }

        template<class T, class ST>
        std::ostream & operator << (std::ostream & s, const MatrixCRow<T, ST> & v)
        {
            s << "vec<" << v.size() << ">[";
            for (typename base_vector<T, ST>::size_type i = 0; i < v.size(); ++i)
            {
                if (i > 0)
                    s << ", ";
                s << v[i];
            }
            return s << "]";
        }

        template<class T, class ST>
        std::ostream & operator << (std::ostream & s, const MatrixCColumn<T, ST> & v)
        {
            s << "vec<" << v.size() << ">[";
            for (typename base_vector<T, ST>::size_type i = 0; i < v.size(); ++i)
            {
                if (i > 0)
                    s << ", ";
                s << v[i];
            }
            return s << "]";
        }

        template<class Op, class Data, class ScT>
        std::ostream & operator << (std::ostream & s, const VecExpression<Op, Data, ScT> & v)
        {
            s << "vec<" << v.size() << ">[";
            for (unsigned int i = 0; i < v.size(); ++i)
            {
                if (i > 0)
                    s << ", ";
                s << v[i];
            }
            return s << "]";
        }
    }
}

#include "vecmat-impl.cpp"
#include "vecmat-ops.hpp"

#endif

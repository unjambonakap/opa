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

#ifndef PLLL_INCLUDE_GUARD__MATRIX_CONVERSION_HPP
#define PLLL_INCLUDE_GUARD__MATRIX_CONVERSION_HPP

namespace plll
{
    template<class IntTypeContext>
    class MatrixConversion
    /* Given an arithmetic::Integer matrix, provides conversions of this matrix to other integer types (specified by
       the given integer context). Conversion is done on-the-fly. In case arithmetic::IntegerContext is given, no
       conversion is done. */
    {
    private:
        IntTypeContext & d_ic;
        linalg::math_matrix<arithmetic::Integer> & d_A;
        linalg::math_matrix<typename IntTypeContext::Integer> d_B;
        bool d_auto_conversion_on_destruction;
    
    public:
        MatrixConversion(IntTypeContext & ic, linalg::math_matrix<arithmetic::Integer> & A, bool auto_conversion_on_destruction = true, bool begin_now = false)
            : d_ic(ic), d_A(A), d_auto_conversion_on_destruction(auto_conversion_on_destruction)
        {
            if (begin_now)
                begin();
        }
    
        ~MatrixConversion()
        {
            if (d_auto_conversion_on_destruction)
                end();
        }
    
        void begin()
        // Makes sure that d_B contains content of d_A (if possible). Until the call of end(), no
        // modifications should be done to d_A.
        {
            d_B.resize(d_A.rows(), d_A.cols(), linalg::Initialize(d_ic));
            for (unsigned i = 0; i < d_A.rows(); ++i)
                for (unsigned j = 0; j < d_A.cols(); ++j)
                    convert(d_B(i, j), d_A(i, j), d_ic);
        }
    
        void end()
        // Makes sure that d_A contains content of d_B. From now on, no calls to matrix() should be
        // made. The content of d_A can be freely modified from now on.
        {
            update_original();
            d_B.resize(0, 0);
        }
    
        inline void update_original()
        // Same as calling end() and then begin(): makes sure that d_A equals d_B.
        {
            d_A.resize(d_B.rows(), d_B.cols());
            for (unsigned i = 0; i < d_A.rows(); ++i)
                for (unsigned j = 0; j < d_A.cols(); ++j)
                    convert(d_A(i, j), d_B(i, j), arithmetic::IntegerContext());
        }
    
        inline linalg::math_matrix<typename IntTypeContext::Integer> & matrix()
        // Valid between begin() and end(). Garuantee: the reference will *not* change between end() and
        // begin().
        {
            return d_B;
        }
    
        inline const linalg::math_matrix<typename IntTypeContext::Integer> & matrix() const
        // Valid between begin() and end(). Garuantee: the reference will *not* change between end() and
        // begin().
        {
            return d_B;
        }
    };

    template<>
    class MatrixConversion<arithmetic::IntegerContext>
    {
    private:
        linalg::math_matrix<arithmetic::Integer> & d_A;
    
    public:
        inline MatrixConversion(arithmetic::IntegerContext &, linalg::math_matrix<arithmetic::Integer> & A, bool = true, bool = false)
            : d_A(A)
        {
        }
    
        inline void begin()
        {
        }
    
        inline void end()
        {
        }
    
        inline void update_original()
        {
        }
    
        inline linalg::math_matrix<arithmetic::Integer> & matrix()
        {
            return d_A;
        }
    
        inline const linalg::math_matrix<arithmetic::Integer> & matrix() const
        {
            return d_A;
        }
    };

    template<class IntTypeContext>
    class VectorConversion
    /* Given an arithmetic::Integer vector, provides conversions of this vector to other integer types (specified by
       the given integer context). Conversion is done on-the-fly. In case arithmetic::IntegerContext is given, no
       conversion is done. */
    {
    private:
        IntTypeContext & d_ic;
        linalg::math_rowvector<arithmetic::Integer> & d_v;
        linalg::math_rowvector<typename IntTypeContext::Integer> d_w;
        bool d_auto_conversion_on_destruction;
    
    public:
        VectorConversion(IntTypeContext & ic, linalg::math_rowvector<arithmetic::Integer> & A, bool auto_conversion_on_destruction = true, bool begin_now = false)
            : d_ic(ic), d_v(A), d_auto_conversion_on_destruction(auto_conversion_on_destruction)
        {
            if (begin_now)
                begin();
        }
    
        ~VectorConversion()
        {
            if (d_auto_conversion_on_destruction)
                end();
        }
    
        void begin()
        // Makes sure that d_w contains content of d_v (if possible). Until the call of end(), no
        // modifications should be done to d_v.
        {
            d_w.resize(d_v.size(), linalg::Initialize(d_ic));
            for (unsigned i = 0; i < d_v.size(); ++i)
                convert(d_w[i], d_v[i], d_ic);
        }
    
        void end()
        // Makes sure that d_v contains content of d_w. From now on, no calls to matrix() should be
        // made. The content of d_v can be freely modified from now on.
        {
            update_original();
            d_w.resize(0);
        }
    
        inline void update_original()
        // Same as calling end() and then begin(): makes sure that d_v equals d_w.
        {
            d_v.resize(d_w.size());
            for (unsigned i = 0; i < d_v.size(); ++i)
                convert(d_v[i], d_w[i], arithmetic::IntegerContext());
        }
    
        inline linalg::math_rowvector<typename IntTypeContext::Integer> & vector()
        // Valid between begin() and end(). Garuantee: the reference will *not* change between end() and
        // begin().
        {
            return d_w;
        }
    
        inline const linalg::math_rowvector<typename IntTypeContext::Integer> & vector() const
        // Valid between begin() and end(). Garuantee: the reference will *not* change between end() and
        // begin().
        {
            return d_w;
        }
    };

    template<>
    class VectorConversion<arithmetic::IntegerContext>
    {
    private:
        linalg::math_rowvector<arithmetic::Integer> & d_v;
    
    public:
        inline VectorConversion(arithmetic::IntegerContext &, linalg::math_rowvector<arithmetic::Integer> & A, bool = true, bool = false)
            : d_v(A)
        {
        }
    
        inline void begin()
        {
        }
    
        inline void end()
        {
        }
    
        inline void update_original()
        {
        }
    
        inline linalg::math_rowvector<arithmetic::Integer> & vector()
        {
            return d_v;
        }
    
        inline const linalg::math_rowvector<arithmetic::Integer> & vector() const
        {
            return d_v;
        }
    };

    template<class IntTypeContext>
    class ValueConversion
    /* Given an arithmetic::Integer value, provides conversions of this value to other integer types (specified by
       the given integer context). Conversion is done on-the-fly. In case arithmetic::IntegerContext is given, no
       conversion is done. */
    {
    private:
        IntTypeContext & d_ic;
        arithmetic::Integer & d_v;
        typename IntTypeContext::Integer d_w;
        bool d_auto_conversion_on_destruction;
    
    public:
        ValueConversion(IntTypeContext & ic, arithmetic::Integer & A, bool auto_conversion_on_destruction = true, bool begin_now = false)
            : d_ic(ic), d_v(A), d_auto_conversion_on_destruction(auto_conversion_on_destruction)
        {
            if (begin_now)
                begin();
        }
    
        ~ValueConversion()
        {
            if (d_auto_conversion_on_destruction)
                end();
        }
    
        void begin()
        // Makes sure that d_w contains content of d_v (if possible). Until the call of end(), no
        // modifications should be done to d_v.
        {
            convert(d_w, d_v, d_ic);
        }
    
        void end()
        // Makes sure that d_v contains content of d_w. From now on, no calls to matrix() should be
        // made. The content of d_v can be freely modified from now on.
        {
            update_original();
        }
    
        inline void update_original()
        // Same as calling end() and then begin(): makes sure that d_v equals d_w.
        {
            convert(d_v, d_w, arithmetic::IntegerContext());
        }
    
        inline typename IntTypeContext::Integer & value()
        // Valid between begin() and end(). Garuantee: the reference will *not* change between end() and
        // begin().
        {
            return d_w;
        }
    
        inline const typename IntTypeContext::Integer & value() const
        // Valid between begin() and end(). Garuantee: the reference will *not* change between end() and
        // begin().
        {
            return d_w;
        }
    };

    template<>
    class ValueConversion<arithmetic::IntegerContext>
    {
    private:
        arithmetic::Integer & d_v;
    
    public:
        inline ValueConversion(arithmetic::IntegerContext &, arithmetic::Integer & A, bool = true, bool = false)
            : d_v(A)
        {
        }
    
        inline void begin()
        {
        }
    
        inline void end()
        {
        }
    
        inline void update_original()
        {
        }
    
        inline arithmetic::Integer & value()
        {
            return d_v;
        }
    
        inline const arithmetic::Integer & value() const
        {
            return d_v;
        }
    };

    template<class IntTypeContext>
    class MatrixConversion2
    /* Given an IntTypeContext matrix, provides conversions of this matrix to arithmetic::Integer. Conversion is done
       on-the-fly. In case arithmetic::IntegerContext is given, no conversion is done. */
    {
    private:
        IntTypeContext & d_ic;
        linalg::math_matrix<typename IntTypeContext::Integer> & d_A;
        linalg::math_matrix<arithmetic::Integer> d_B;
        bool d_auto_conversion_on_destruction;
    
    public:
        MatrixConversion2(IntTypeContext & ic, linalg::math_matrix<typename IntTypeContext::Integer> & A, bool auto_conversion_on_destruction = true, bool begin_now = false)
            : d_ic(ic), d_A(A), d_auto_conversion_on_destruction(auto_conversion_on_destruction)
        {
            if (begin_now)
                begin();
        }
    
        ~MatrixConversion2()
        {
            if (d_auto_conversion_on_destruction)
                end();
        }
    
        void begin()
        // Makes sure that d_B contains content of d_A (if possible). Until the call of end(), no
        // modifications should be done to d_A.
        {
            d_B.resize(d_A.rows(), d_A.cols());
            for (unsigned i = 0; i < d_A.rows(); ++i)
                for (unsigned j = 0; j < d_A.cols(); ++j)
                    convert(d_B(i, j), d_A(i, j), arithmetic::IntegerContext());
        }
    
        void end()
        // Makes sure that d_A contains content of d_B. From now on, no calls to matrix() should be
        // made. The content of d_A can be freely modified from now on.
        {
            update_original();
            d_B.resize(0, 0);
        }
    
        inline void update_original()
        // Same as calling end() and then begin(): makes sure that d_A equals d_B.
        {
            d_A.resize(d_B.rows(), d_B.cols(), linalg::Initialize(d_ic));
            for (unsigned i = 0; i < d_A.rows(); ++i)
                for (unsigned j = 0; j < d_A.cols(); ++j)
                    convert(d_A(i, j), d_B(i, j), d_ic);
        }
    
        inline linalg::math_matrix<arithmetic::Integer> & matrix()
        // Valid between begin() and end(). Garuantee: the reference will *not* change between end() and
        // begin().
        {
            return d_B;
        }
    
        inline const linalg::math_matrix<arithmetic::Integer> & matrix() const
        // Valid between begin() and end(). Garuantee: the reference will *not* change between end() and
        // begin().
        {
            return d_B;
        }
    };

    template<>
    class MatrixConversion2<arithmetic::IntegerContext>
    {
    private:
        linalg::math_matrix<arithmetic::Integer> & d_A;
    
    public:
        inline MatrixConversion2(arithmetic::IntegerContext &, linalg::math_matrix<arithmetic::Integer> & A, bool = true, bool = false)
            : d_A(A)
        {
        }
    
        inline void begin()
        {
        }
    
        inline void end()
        {
        }
    
        inline void update_original()
        {
        }
    
        inline linalg::math_matrix<arithmetic::Integer> & matrix()
        {
            return d_A;
        }
    
        inline const linalg::math_matrix<arithmetic::Integer> & matrix() const
        {
            return d_A;
        }
    };

    template<class IntTypeContext>
    class VectorConversion2
    /* Given an IntTypeContext vector, provides conversions of this vector to arithmetic::Integer. Conversion is done
       on-the-fly. In case arithmetic::IntegerContext is given, no conversion is done. */
    {
    private:
        IntTypeContext & d_ic;
        linalg::math_rowvector<typename IntTypeContext::Integer> & d_v;
        linalg::math_rowvector<arithmetic::Integer> d_w;
        bool d_auto_conversion_on_destruction;
    
    public:
        VectorConversion2(IntTypeContext & ic, linalg::math_rowvector<typename IntTypeContext::Integer> & A, bool auto_conversion_on_destruction = true, bool begin_now = false)
            : d_ic(ic), d_v(A), d_auto_conversion_on_destruction(auto_conversion_on_destruction)
        {
            if (begin_now)
                begin();
        }
    
        ~VectorConversion2()
        {
            if (d_auto_conversion_on_destruction)
                end();
        }
    
        void begin()
        // Makes sure that d_w contains content of d_v (if possible). Until the call of end(), no
        // modifications should be done to d_v.
        {
            d_w.resize(d_v.size());
            for (unsigned i = 0; i < d_v.size(); ++i)
                convert(d_w[i], d_v[i], arithmetic::IntegerContext());
        }
    
        void end()
        // Makes sure that d_v contains content of d_w. From now on, no calls to matrix() should be
        // made. The content of d_v can be freely modified from now on.
        {
            update_original();
            d_w.resize(0);
        }
    
        inline void update_original()
        // Same as calling end() and then begin(): makes sure that d_v equals d_w.
        {
            d_v.resize(d_w.size(), linalg::Initialize(d_ic));
            for (unsigned i = 0; i < d_v.size(); ++i)
                convert(d_v[i], d_w[i], d_ic);
        }
    
        inline linalg::math_rowvector<arithmetic::Integer> & vector()
        // Valid between begin() and end(). Garuantee: the reference will *not* change between end() and
        // begin().
        {
            return d_w;
        }
    
        inline const linalg::math_rowvector<arithmetic::Integer> & vector() const
        // Valid between begin() and end(). Garuantee: the reference will *not* change between end() and
        // begin().
        {
            return d_w;
        }
    };

    template<>
    class VectorConversion2<arithmetic::IntegerContext>
    {
    private:
        linalg::math_rowvector<arithmetic::Integer> & d_v;
    
    public:
        inline VectorConversion2(arithmetic::IntegerContext &, linalg::math_rowvector<arithmetic::Integer> & A, bool = true, bool = false)
            : d_v(A)
        {
        }
    
        inline void begin()
        {
        }
    
        inline void end()
        {
        }
    
        inline void update_original()
        {
        }
    
        inline linalg::math_rowvector<arithmetic::Integer> & vector()
        {
            return d_v;
        }
    
        inline const linalg::math_rowvector<arithmetic::Integer> & vector() const
        {
            return d_v;
        }
    };

    template<class IntTypeContext>
    class ValueConversion2
    /* Given an IntTypeContext value, provides conversions of this value to arithmetic::Integer. Conversion is done
       on-the-fly. In case arithmetic::IntegerContext is given, no conversion is done. */
    {
    private:
        IntTypeContext & d_ic;
        typename IntTypeContext::Integer & d_v;
        arithmetic::Integer d_w;
        bool d_auto_conversion_on_destruction;
    
    public:
        ValueConversion2(IntTypeContext & ic, typename IntTypeContext::Integer & A, bool auto_conversion_on_destruction = true, bool begin_now = false)
            : d_ic(ic), d_v(A), d_auto_conversion_on_destruction(auto_conversion_on_destruction)
        {
            if (begin_now)
                begin();
        }
    
        ~ValueConversion2()
        {
            if (d_auto_conversion_on_destruction)
                end();
        }
    
        void begin()
        // Makes sure that d_w contains content of d_v (if possible). Until the call of end(), no
        // modifications should be done to d_v.
        {
            convert(d_w, d_v, arithmetic::IntegerContext());
        }
    
        void end()
        // Makes sure that d_v contains content of d_w. From now on, no calls to matrix() should be
        // made. The content of d_v can be freely modified from now on.
        {
            update_original();
        }
    
        inline void update_original()
        // Same as calling end() and then begin(): makes sure that d_v equals d_w.
        {
            convert(d_v, d_w, d_ic);
        }
    
        inline arithmetic::Integer & value()
        // Valid between begin() and end(). Garuantee: the reference will *not* change between end() and
        // begin().
        {
            return d_w;
        }
    
        inline const arithmetic::Integer & value() const
        // Valid between begin() and end(). Garuantee: the reference will *not* change between end() and
        // begin().
        {
            return d_w;
        }
    };

    template<>
    class ValueConversion2<arithmetic::IntegerContext>
    {
    private:
        arithmetic::Integer & d_v;
    
    public:
        inline ValueConversion2(arithmetic::IntegerContext &, arithmetic::Integer & A, bool = true, bool = false)
            : d_v(A)
        {
        }
    
        inline void begin()
        {
        }
    
        inline void end()
        {
        }
    
        inline void update_original()
        {
        }
    
        inline arithmetic::Integer & value()
        {
            return d_v;
        }
    
        inline const arithmetic::Integer & value() const
        {
            return d_v;
        }
    };

    template<class IntTypeContext>
    class MatrixConversionConst
    /* Given an arithmetic::Integer matrix, provides conversions of this matrix to other integer types (specified by
       the given integer context). ConversionConst is done on-the-fly. In case arithmetic::IntegerContext is given, no
       conversion is done. */
    {
    private:
        linalg::math_matrix<typename IntTypeContext::Integer> d_B;
    
    public:
        MatrixConversionConst(IntTypeContext & ic, const linalg::math_matrix<arithmetic::Integer> & A)
        {
            d_B.resize(A.rows(), A.cols());
            for (unsigned i = 0; i < A.rows(); ++i)
                for (unsigned j = 0; j < A.cols(); ++j)
                    convert(d_B(i, j), A(i, j), ic);
        }
    
        inline const linalg::math_matrix<typename IntTypeContext::Integer> & matrix() const
        {
            return d_B;
        }
    };

    template<>
    class MatrixConversionConst<arithmetic::IntegerContext>
    {
    private:
        const linalg::math_matrix<arithmetic::Integer> & d_A;
    
    public:
        inline MatrixConversionConst(arithmetic::IntegerContext &, const linalg::math_matrix<arithmetic::Integer> & A)
            : d_A(A)
        {
        }
    
        inline const linalg::math_matrix<arithmetic::Integer> & matrix() const
        {
            return d_A;
        }
    };

    template<class IntTypeContext>
    class VectorConversionConst
    /* Given an arithmetic::Integer vector, provides conversions of this vector to other integer types (specified by
       the given integer context). ConversionConst is done on-the-fly. In case arithmetic::IntegerContext is given, no
       conversion is done. */
    {
    private:
        linalg::math_rowvector<typename IntTypeContext::Integer> d_w;
    
    public:
        VectorConversionConst(IntTypeContext & ic, const linalg::math_rowvector<arithmetic::Integer> & v)
        {
            d_w.resize(v.size());
            for (unsigned i = 0; i < v.size(); ++i)
                convert(d_w[i], v[i], ic);
        }
    
        inline const linalg::math_rowvector<typename IntTypeContext::Integer> & vector() const
        {
            return d_w;
        }
    };

    template<>
    class VectorConversionConst<arithmetic::IntegerContext>
    {
    private:
        const linalg::math_rowvector<arithmetic::Integer> & d_v;
    
    public:
        inline VectorConversionConst(arithmetic::IntegerContext &, const linalg::math_rowvector<arithmetic::Integer> & A)
            : d_v(A)
        {
        }
    
        inline const linalg::math_rowvector<arithmetic::Integer> & vector() const
        {
            return d_v;
        }
    };

    template<class IntTypeContext>
    class ValueConversionConst
    /* Given an arithmetic::Integer value, provides conversions of this value to other integer types (specified by
       the given integer context). ConversionConst is done on-the-fly. In case arithmetic::IntegerContext is given, no
       conversion is done. */
    {
    private:
        typename IntTypeContext::Integer d_w;
    
    public:
        ValueConversionConst(IntTypeContext & ic, const arithmetic::Integer & v)
        {
            convert(d_w, v, ic);
        }
    
        inline const typename IntTypeContext::Integer & value() const
        {
            return d_w;
        }
    };

    template<>
    class ValueConversionConst<arithmetic::IntegerContext>
    {
    private:
        const arithmetic::Integer & d_v;
    
    public:
        inline ValueConversionConst(arithmetic::IntegerContext &, const arithmetic::Integer & A)
            : d_v(A)
        {
        }
    
        inline const arithmetic::Integer & value() const
        {
            return d_v;
        }
    };

    template<class IntTypeContext>
    class MatrixConversionConst2
    /* Given an IntTypeContext matrix, provides conversions of this matrix to arithmetic::Integer. ConversionConst is done
       on-the-fly. In case arithmetic::IntegerContext is given, no conversion is done. */
    {
    private:
        linalg::math_matrix<arithmetic::Integer> d_B;
    
    public:
        MatrixConversionConst2(IntTypeContext & ic, const linalg::math_matrix<typename IntTypeContext::Integer> & A)
        {
            d_B.resize(A.rows(), A.cols());
            for (unsigned i = 0; i < A.rows(); ++i)
                for (unsigned j = 0; j < A.cols(); ++j)
                    convert(d_B(i, j), A(i, j), arithmetic::IntegerContext());
        }
    
        inline const linalg::math_matrix<arithmetic::Integer> & matrix() const
        {
            return d_B;
        }
    };

    template<>
    class MatrixConversionConst2<arithmetic::IntegerContext>
    {
    private:
        const linalg::math_matrix<arithmetic::Integer> & d_A;
    
    public:
        inline MatrixConversionConst2(arithmetic::IntegerContext &, const linalg::math_matrix<arithmetic::Integer> & A)
            : d_A(A)
        {
        }
    
        inline const linalg::math_matrix<arithmetic::Integer> & matrix() const
        {
            return d_A;
        }
    };

    template<class IntTypeContext>
    class VectorConversionConst2
    /* Given an IntTypeContext vector, provides conversions of this vector to arithmetic::Integer. ConversionConst is done
       on-the-fly. In case arithmetic::IntegerContext is given, no conversion is done. */
    {
    private:
        linalg::math_rowvector<arithmetic::Integer> d_w;
    
    public:
        VectorConversionConst2(IntTypeContext & ic, const linalg::math_rowvector<typename IntTypeContext::Integer> & v)
        {
            d_w.resize(v.size());
            for (unsigned i = 0; i < v.size(); ++i)
                convert(d_w[i], v[i], arithmetic::IntegerContext());
        }
    
        inline const linalg::math_rowvector<arithmetic::Integer> & vector() const
        {
            return d_w;
        }
    };

    template<>
    class VectorConversionConst2<arithmetic::IntegerContext>
    {
    private:
        const linalg::math_rowvector<arithmetic::Integer> & d_v;
    
    public:
        inline VectorConversionConst2(arithmetic::IntegerContext &, const linalg::math_rowvector<arithmetic::Integer> & A)
            : d_v(A)
        {
        }
    
        inline const linalg::math_rowvector<arithmetic::Integer> & vector() const
        {
            return d_v;
        }
    };

    template<class IntTypeContext>
    class ValueConversionConst2
    /* Given an IntTypeContext value, provides conversions of this value to arithmetic::Integer. ConversionConst is done
       on-the-fly. In case arithmetic::IntegerContext is given, no conversion is done. */
    {
    private:
        arithmetic::Integer d_w;
    
    public:
        ValueConversionConst2(IntTypeContext & ic, const typename IntTypeContext::Integer & v)
        {
            convert(d_w, v, arithmetic::IntegerContext());
        }
    
        inline const arithmetic::Integer & value() const
        {
            return d_w;
        }
    };

    template<>
    class ValueConversionConst2<arithmetic::IntegerContext>
    {
    private:
        const arithmetic::Integer & d_v;
    
    public:
        inline ValueConversionConst2(arithmetic::IntegerContext &, const arithmetic::Integer & A)
            : d_v(A)
        {
        }
    
        inline const arithmetic::Integer & value() const
        {
            return d_v;
        }
    };
}

#endif

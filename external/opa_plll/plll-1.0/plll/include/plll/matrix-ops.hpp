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

#ifndef PLLL_INCLUDE_GUARD__MATRIX_OPS_HPP
#define PLLL_INCLUDE_GUARD__MATRIX_OPS_HPP

#include "arithmetic.hpp"
#include <memory>

/**
   \file
   \brief Operator definitions for matrices and vectors.
   
   This header contains abstract templates used to implement all operations on vectors and
   matrices. These templates are instantiated in `matrix-ops2.hpp`.
*/
namespace plll
{
    namespace linalg
    {
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        //// expression framework
        ////////////////////////////////////////////////////////////////////////////////////////////////////

        namespace implementation
        {
            namespace expressions
            {
                template<template<typename DataType> class Operator, typename Data>
                class expr : private Operator<Data>
                {
                private:
                    Data d_data;
                    
                public:
                    inline expr(const Data & data)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Data(data)))
                        : d_data(data)
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::expr::expr(" << getAddress(*this) << "; " << getAddress(data) << ")");
                    }
                    
                    inline expr(const Operator<Data> & op, const Data & data)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Operator<Data>(op)) && noexcept(Data(data)))
                        : Operator<Data>(op), d_data(data)
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::expr::expr(" << getAddress(*this) << "; " << getAddress(op) << "; " << getAddress(data) << ")");
                    }
                    
                    inline const expr<Operator, Data> & operator = (const expr<Operator, Data> & a) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(assign(helper::make_type_lvalue<expr<Operator, Data>>(), a)))
                    {
                        assign(*this, a);
                        return *this;
                    }
                    
                    enum { has_direct_access = Operator<Data>::has_direct_access };
                    
                    typedef typename Operator<Data>::CoeffType CoeffType;
                    typedef typename Operator<Data>::CoeffType_Get CoeffType_Get;
                    typedef typename Operator<Data>::GetCoeffSteps_Type GetCoeffSteps_Type;
                    typedef typename Operator<Data>::ValueType ValueType;
                    typedef typename helper::SelectFirstType<Operator<Data>::use_temporary_on_evaluate,
                                                             MatrixTemporaryWrapper<ValueType>,
                                                             expr<Operator, Data> >::result LazyEvalType;
                    typedef typename Operator<Data>::data_pointer_type data_pointer_type;
                    typedef typename Operator<Data>::data_ref_type data_ref_type;
                    
                    template<typename ResultType>
                    inline void get_coeff(ResultType & result, size_type i, size_type j) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(helper::make_type_lvalue<expr<Operator, Data> >().Operator<Data>::get_coeff(result, i, j, d_data)))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::expr::get_coeff(" << getAddress(*this) << "; " << getAddress(result) << " " << i << " " << j << ") const");
                        PLLL_INTERNAL_STATIC_CHECK(MatrixInfo<ValueType>::is_matrix, OnlyDefinedForMatrices);
                        Operator<Data>::get_coeff(result, i, j, d_data);
                    }
                    
                    inline bool get_coeff_alwayszero() const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(helper::make_type_lvalue<expr<Operator, Data> >().Operator<Data>::get_coeff_alwayszero(d_data)))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::expr::get_coeff_alwayszero(" << getAddress(*this) << ")");
                        PLLL_INTERNAL_STATIC_CHECK(MatrixInfo<ValueType>::is_matrix, OnlyDefinedForMatrices);
                        return Operator<Data>::get_coeff_alwayszero(d_data);
                    }
                    
                    inline GetCoeffSteps_Type get_coeff_steps(size_type i, size_type j) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(helper::make_type_lvalue<expr<Operator, Data> >().Operator<Data>::get_coeff_steps(i, j, d_data)))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::expr::get_coeff_steps(" << getAddress(*this) << ")");
                        PLLL_INTERNAL_STATIC_CHECK(MatrixInfo<ValueType>::is_matrix, OnlyDefinedForMatrices);
                        return Operator<Data>::get_coeff_steps(i, j, d_data);
                    }
                    
                    inline CoeffType_Get operator () (size_type i, size_type j) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(helper::make_type_lvalue<expr<Operator, Data> >().Operator<Data>::operator()(i, j, d_data)))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::expr::operator () (" << getAddress(*this) << "; " << i << " " << j << ")");
                        PLLL_INTERNAL_STATIC_CHECK(MatrixInfo<ValueType>::is_matrix, OnlyDefinedForMatrices);
                        return Operator<Data>::operator()(i, j, d_data);
                    }
                    
                    inline CoeffType_Get operator [] (size_type i) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL((MatrixInfo<ValueType>::cols == 1) ?
                                                                  noexcept(helper::make_type_lvalue<expr<Operator, Data> >().Operator<Data>::operator()(i, 0, d_data)) :
                                                                  noexcept(helper::make_type_lvalue<expr<Operator, Data> >().Operator<Data>::operator()(0, i, d_data)))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::expr::operator [] (" << getAddress(*this) << "; " << i << ")");
                        PLLL_INTERNAL_STATIC_CHECK(MatrixInfo<ValueType>::is_matrix, OnlyDefinedForMatrices);
                        PLLL_INTERNAL_STATIC_CHECK((MatrixInfo<ValueType>::cols == 1) ||
                                                   (MatrixInfo<ValueType>::rows == 1), OnlyDefinedForVectors);
                        return (MatrixInfo<ValueType>::cols == 1) ? Operator<Data>::operator()(i, 0, d_data) : Operator<Data>::operator()(0, i, d_data);
                    }
                    
                    inline data_pointer_type data() const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(helper::make_type_lvalue<expr<Operator, Data> >().Operator<Data>::data(d_data)))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::expr::data(" << getAddress(*this) << ") const");
                        return Operator<Data>::data(d_data);
                    }
                    
                    inline data_ref_type data(size_type i) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(helper::make_type_lvalue<expr<Operator, Data> >().Operator<Data>::data(d_data, i)))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::expr::data(" << getAddress(*this) << "; " << i << ")");
                        return Operator<Data>::data(d_data, i);
                    }
                    
                    inline size_type rows() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::expr::rows(" << getAddress(*this) << ")");
                        return Operator<Data>::rows(d_data);
                    }
                    
                    inline size_type cols() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::expr::cols(" << getAddress(*this) << ")");
                        return Operator<Data>::cols(d_data);
                    }
                    
                    inline size_type size() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::expr::size(" << getAddress(*this) << ")");
                        return Operator<Data>::size(d_data);
                    }
                    
                    template<unsigned SRows, unsigned SCols>
                    inline expr<sub<SRows, SCols>::template operation_generic, expr<Operator, Data> >
                    block(size_type r_ofs, size_type c_ofs) const
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::expr::block<" << SRows << " " << SCols << ">(" << getAddress(*this) << "; " << r_ofs << " " << c_ofs << ") const");
                        return expr<sub<SRows, SCols>::template operation_generic, expr<Operator, Data> >
                            (typename sub<SRows, SCols>::template operation_generic<expr<Operator, Data> >(SRows, SCols, r_ofs, c_ofs, *this), *this);
                    }
                    
                    inline expr<sub<Flexible, Flexible>::template operation_generic, expr<Operator, Data> >
                    block(size_type r_ofs, size_type c_ofs, size_type rows, size_type cols) const
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::expr::block(" << getAddress(*this) << "; " << r_ofs << " " << c_ofs << " " << rows << " " << cols << ")");
                        return expr<sub<Flexible, Flexible>::template operation_generic, expr<Operator, Data> >
                            (typename sub<Flexible, Flexible>::template operation_generic<expr<Operator, Data> >(rows, cols, r_ofs, c_ofs, *this), *this);
                    }
                    
                    inline expr<sub_1d<MatrixInfo<typename Operator<Data>::ValueType>::cols>::template operation_row, expr<Operator, Data> >
                    row(size_type row_index) const
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::expr::row(" << getAddress(*this) << "; " << row_index << ") const");
                        return expr<sub_1d<MatrixInfo<typename Operator<Data>::ValueType>::cols>::template operation_row, expr<Operator, Data> >
                            (typename sub_1d<MatrixInfo<typename Operator<Data>::ValueType>::cols>::template operation_row<expr<Operator, Data> >(cols(), row_index, *this), *this);
                    }
                    
                    inline expr<sub_1d<MatrixInfo<typename Operator<Data>::ValueType>::rows>::template operation_col, expr<Operator, Data> >
                    col(size_type col_index) const
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::expr::col(" << getAddress(*this) << "; " << col_index << ") const");
                        return expr<sub_1d<MatrixInfo<typename Operator<Data>::ValueType>::rows>::template operation_col, expr<Operator, Data> >
                            (typename sub_1d<MatrixInfo<typename Operator<Data>::ValueType>::rows>::template operation_col<expr<Operator, Data> >(rows(), col_index, *this), *this);
                    }
                    
                    inline expr<expressions::transpose, expr<Operator, Data> > transpose() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::expr::transpose(" << getAddress(*this) << ") const");
                        return expr<expressions::transpose, expr<Operator, Data> >(*this);
                    }
                    
                    inline ValueType evaluate() const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(ValueType(helper::make_type_lvalue<const expr<Operator, Data> >())))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::expr::evaluate(" << getAddress(*this) << ")");
                        return ValueType(*this);
                    }
                    
                    inline LazyEvalType lazy_evaluate() const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(LazyEvalType(helper::make_type_lvalue<const expr<Operator, Data> >())))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::expr::lazy_evaluate(" << getAddress(*this) << ")");
                        return LazyEvalType(*this);
                    }
                    
                    template<typename MatrixType>
                    inline bool involves_this_matrix(const MatrixType & A) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::expr::involves_this_matrix(" << getAddress(*this) << "; " << getAddress(A) << ")");
                        return Operator<Data>::involves_this_matrix(A, d_data);
                    }
                    
                    template<typename MatrixType>
                    inline bool test_involvement(const MatrixType & A) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::expr::test_involvement(" << getAddress(*this) << "; " << getAddress(A) << ")");
                        return Operator<Data>::test_involvement(A, d_data);
                    }
                    
                    template<typename T, int R, int C, typename ST, bool MO>
                    const expr<Operator, Data> & operator = (const base_matrix<T, R, C, ST, MO> & m) const
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::operator = (" << getAddress(*this) << "; " << getAddress(m) << ")");
                        PLLL_INTERNAL_STATIC_CHECK((MatrixInfo<base_matrix<T, R, C, ST, MO> >::is_matrix), RequiresMatrixType);
                        assign(*this, m);
                        return *this;
                    }
                    
                    template<template<typename> class Op_, typename Data_>
                    const expr<Operator, Data> & operator = (const expr<Op_, Data_> & m) const
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::operator = (" << getAddress(*this) << "; " << getAddress(m) << ")");
                        PLLL_INTERNAL_STATIC_CHECK((MatrixInfo<expr<Op_, Data_> >::is_matrix), RequiresMatrixType);
                        assign(*this, m);
                        return *this;
                    }
                    
#if __cplusplus >= 201103L
                    template<typename T, int R, int C, typename ST, bool MO>
                    const expr<Operator, Data> & operator = (base_matrix<T, R, C, ST, MO> && m) const
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::operator = (" << getAddress(*this) << "; " << getAddress(m) << ") [&&]");
                        PLLL_INTERNAL_STATIC_CHECK((MatrixInfo<base_matrix<T, R, C, ST, MO> >::is_matrix), RequiresMatrixType);
                        assign(*this, std::move(m));
                        return *this;
                    }
                    
                    template<template<typename> class Op_, typename Data_>
                    const expr<Operator, Data> & operator = (expr<Op_, Data_> && m) const
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::operator = (" << getAddress(*this) << "; " << getAddress(m) << ") [&&]");
                        PLLL_INTERNAL_STATIC_CHECK((MatrixInfo<expr<Op_, Data_> >::is_matrix), RequiresMatrixType);
                        assign(*this, std::move(m));
                        return *this;
                    }
#endif
                    
                    template<typename MatrixType>
                    void assign_resize(const MatrixType & mat) const
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::expr::assign_resize(" << getAddress(*this) << "; " << getAddress(mat) << ")");
                        PLLL_INTERNAL_STATIC_CHECK(MatrixInfo<MatrixType>::is_matrix, RequiresMatrixType);
                        PLLL_INTERNAL_STATIC_CHECK(Operator<Data>::can_resize_rows || Operator<Data>::can_resize_cols, RequiresResizeableMatrixType);
                        PLLL_INTERNAL_STATIC_CHECK(Operator<Data>::can_assign_to, RequiresAssignableMatrixType);
                        Operator<Data>::assign_resize(mat, d_data);
                    }
                    
#if __cplusplus >= 201103L
                    template<typename MatrixType>
                    void move_resize(MatrixType && mat) const
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::expr::move_resize(" << getAddress(*this) << "; " << getAddress(mat) << ")");
                        PLLL_INTERNAL_STATIC_CHECK(MatrixInfo<MatrixType>::is_matrix, RequiresMatrixType);
                        PLLL_INTERNAL_STATIC_CHECK(MatrixInfo<MatrixType>::can_move_from, RequiresMoveableMatrixType);
                        PLLL_INTERNAL_STATIC_CHECK(Operator<Data>::can_resize_rows || Operator<Data>::can_resize_cols, RequiresResizeableMatrixType);
                        PLLL_INTERNAL_STATIC_CHECK(Operator<Data>::can_assign_to, RequiresAssignableMatrixType);
                        Operator<Data>::move_resize(std::move(mat), d_data);
                    }
#else
                    template<typename MatrixType>
                    void move_resize(const MatrixType & mat) const
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::expr::move_resize(" << getAddress(*this) << "; " << getAddress(mat) << ")");
                        PLLL_INTERNAL_STATIC_CHECK(MatrixInfo<MatrixType>::is_matrix, RequiresMatrixType);
                        PLLL_INTERNAL_STATIC_CHECK(MatrixInfo<MatrixType>::can_move_from, RequiresMoveableMatrixType);
                        PLLL_INTERNAL_STATIC_CHECK(Operator<Data>::can_resize_rows || Operator<Data>::can_resize_cols, RequiresResizeableMatrixType);
                        PLLL_INTERNAL_STATIC_CHECK(Operator<Data>::can_assign_to, RequiresAssignableMatrixType);
                        Operator<Data>::move_resize(mat, d_data);
                    }
#endif
                    
                    typedef typename Operator<Data>::Enumerator Enumerator;
                    typedef typename Operator<Data>::RowEnumerator RowEnumerator;
                    typedef typename Operator<Data>::ColEnumerator ColEnumerator;
                    typedef typename Operator<Data>::ConstEnumerator ConstEnumerator;
                    typedef typename Operator<Data>::ConstRowEnumerator ConstRowEnumerator;
                    typedef typename Operator<Data>::ConstColEnumerator ConstColEnumerator;
                    typedef typename Operator<Data>::DefaultEnumerator DefaultEnumerator;
                    typedef typename Operator<Data>::DefaultRowEnumerator DefaultRowEnumerator;
                    typedef typename Operator<Data>::DefaultColEnumerator DefaultColEnumerator;
                    
                    inline DefaultEnumerator enumerate() const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(helper::make_type_lvalue<expr<Operator, Data> >().Operator<Data>::enumerate(d_data)))
                    {
                        return Operator<Data>::enumerate(d_data);
                    }
                    
                    inline DefaultRowEnumerator enumerate_row(size_type row) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(helper::make_type_lvalue<expr<Operator, Data> >().Operator<Data>::enumerate_row(row, d_data)))
                    {
                        return Operator<Data>::enumerate_row(row, d_data);
                    }
                    
                    inline DefaultColEnumerator enumerate_col(size_type col) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(helper::make_type_lvalue<expr<Operator, Data> >().Operator<Data>::enumerate_col(col, d_data)))
                    {
                        return Operator<Data>::enumerate_col(col, d_data);
                    }
                    
                    typedef typename Operator<Data>::RowsEnumerator RowsEnumerator;
                    typedef typename Operator<Data>::ColsEnumerator ColsEnumerator;
                    typedef typename Operator<Data>::ConstRowsEnumerator ConstRowsEnumerator;
                    typedef typename Operator<Data>::ConstColsEnumerator ConstColsEnumerator;
                    typedef typename Operator<Data>::DefaultRowsEnumerator DefaultRowsEnumerator;
                    typedef typename Operator<Data>::DefaultColsEnumerator DefaultColsEnumerator;
                    
                    inline DefaultRowsEnumerator enumerate_rows() const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(helper::make_type_lvalue<expr<Operator, Data> >().Operator<Data>::enumerate_rows(d_data)))
                    {
                        return Operator<Data>::enumerate_rows(d_data);
                    }
                    
                    inline DefaultColsEnumerator enumerate_cols() const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(helper::make_type_lvalue<expr<Operator, Data> >().Operator<Data>::enumerate_cols(d_data)))
                    {
                        return Operator<Data>::enumerate_cols(d_data);
                    }
                };
                
                template<typename Data>
                class identity_operation
                {
                public:
                    typedef typename Data::CoeffType CoeffType;
                    typedef typename Data::CoeffType_Get CoeffType_Get;
                    typedef typename Data::GetCoeffSteps_Type GetCoeffSteps_Type;
                    typedef typename Data::ValueType ValueType;
                    typedef typename Data::data_pointer_type data_pointer_type;
                    typedef typename Data::data_ref_type data_ref_type;
                    
                    enum { use_temporary_on_evaluate = MatrixInfo<Data>::use_temporary_on_evaluate,
                           coeffs_are_simple_expressions = MatrixInfo<Data>::coeffs_are_simple_expressions,
                           is_const = MatrixInfo<Data>::is_const,
                           can_assign_to = MatrixInfo<Data>::can_assign_to,
                           can_move_from = MatrixInfo<Data>::can_move_from,
                           can_resize_rows = MatrixInfo<Data>::can_resize_rows,
                           can_resize_cols = MatrixInfo<Data>::can_resize_cols,
                           has_direct_access = Data::has_direct_access };
                    
                    template<typename ResultType>
                    static inline void get_coeff(ResultType & result, size_type i, size_type j, const Data & D)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(D.get_coeff(result, i, j)))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::identity_operation::get_coeff(" << getAddress(*this) << "; " << getAddress(result) << " " << i << " " << j << " " << getAddress(D) << ")");
                        D.get_coeff(result, i, j);
                    }
                    
                    static inline bool get_coeff_alwayszero(const Data & D)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(D.get_coeff_alwayszero()))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::identity_operation::get_coeff_alwayszero(" << getAddress(*this) << "; " << getAddress(D) << ")");
                        return D.get_coeff_alwayszero();
                    }
                    
                    static inline GetCoeffSteps_Type get_coeff_steps(size_type i, size_type j, const Data & D)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(D.get_coeff_steps(i, j)))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::identity_operation::get_coeff_steps(" << getAddress(*this) << "; " << i << " " << j << " " << getAddress(D) << ")");
                        return D.get_coeff_steps(i, j);
                    }
                    
                    inline CoeffType_Get operator () (size_type i, size_type j, const Data & D) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(D(i, j)))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::identity_operation::operator () (" << getAddress(*this) << "; " << i << " " << j << " " << getAddress(D) << ")");
                        return D(i, j);
                    }
                    
                    static inline data_pointer_type data(const Data & D)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(D.data()))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::identity_operation::data(" << getAddress(*this) << ";" << getAddress(D) << ") const");
                        return D.data();
                    }
                    
                    static inline data_ref_type data(const Data & D, size_type i)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(D.data(i)))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::identity_operation::data(" << getAddress(*this) << "; " << getAddress(D) << " " << i << ")");
                        return D.data(i);
                    }
                    
                    static inline size_type rows(const Data & D) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::identity_operation::rows(" << getAddress(*this) << "; " << getAddress(D) << ")");
                        return D.rows();
                    }
                    
                    static inline size_type cols(const Data & D) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::identity_operation::cols(" << getAddress(*this) << "; " << getAddress(D) << ")");
                        return D.cols();
                    }
                    
                    static inline size_type size(const Data & D) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::identity_operation::size(" << getAddress(*this) << "; " << getAddress(D) << ")");
                        return D.size();
                    }
                    
                    template<typename MatrixType>
                    static inline bool involves_this_matrix(const MatrixType & matrix, const Data & D) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::identity_operation::involves_this_matrix(" << getAddress(*this) << "; " << getAddress(matrix) << " " << getAddress(D) << ")");
                        return D.involves_this_matrix(matrix);
                    }
                
                    template<typename MatrixType>
                    static inline bool test_involvement(const MatrixType & matrix, const Data & D) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::identity_operation::test_involvement(" << getAddress(*this) << "; " << getAddress(matrix) << " " << getAddress(D) << ")");
                        return D.test_involvement(matrix);
                    }
                
                    template<typename MatrixType>
                    static void assign_resize(const MatrixType & mat, const Data & D)
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::identity_operation::assign_resize(" << getAddress(*this) << "; " << getAddress(mat) << " " << getAddress(D) << ")");
                        D.assign_resize(mat);
                    }
                
#if __cplusplus >= 201103L
                    template<typename MatrixType>
                    static void move_resize(MatrixType && mat, const Data & D)
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::identity_operation::move_resize(" << getAddress(*this) << "; " << getAddress(mat) << " " << getAddress(D) << ")");
                        D.move_resize(std::move(mat));
                    }
#else
                    template<typename MatrixType>
                    static void move_resize(const MatrixType & mat, const Data & D)
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::identity_operation::move_resize(" << getAddress(*this) << "; " << getAddress(mat) << " " << getAddress(D) << ")");
                        D.move_resize(mat);
                    }
#endif
                    
                    typedef typename Data::Enumerator Enumerator;
                    typedef typename Data::RowEnumerator RowEnumerator;
                    typedef typename Data::ColEnumerator ColEnumerator;
                    typedef typename Data::ConstEnumerator ConstEnumerator;
                    typedef typename Data::ConstRowEnumerator ConstRowEnumerator;
                    typedef typename Data::ConstColEnumerator ConstColEnumerator;
                    typedef typename Data::DefaultEnumerator DefaultEnumerator;
                    typedef typename Data::DefaultRowEnumerator DefaultRowEnumerator;
                    typedef typename Data::DefaultColEnumerator DefaultColEnumerator;
                
                    static inline DefaultEnumerator enumerate(const Data & data)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(data.enumerate()))
                    {
                        return data.enumerate();
                    }
                
                    static inline DefaultRowEnumerator enumerate_row(size_type row, const Data & data)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(data.enumerate_row(row)))
                    {
                        return data.enumerate_row(row);
                    }
                
                    static inline DefaultColEnumerator enumerate_col(size_type col, const Data & data)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(data.enumerate_col(col)))
                    {
                        return data.enumerate_col(col);
                    }
                
                    typedef typename Data::RowsEnumerator RowsEnumerator;
                    typedef typename Data::ColsEnumerator ColsEnumerator;
                    typedef typename Data::ConstRowsEnumerator ConstRowsEnumerator;
                    typedef typename Data::ConstColsEnumerator ConstColsEnumerator;
                    typedef typename Data::DefaultRowsEnumerator DefaultRowsEnumerator;
                    typedef typename Data::DefaultColsEnumerator DefaultColsEnumerator;
                    
                    static inline DefaultRowsEnumerator enumerate_rows(const Data & data)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(data.enumerate_rows()))
                    {
                        return data.enumerate_rows();
                    }
                    
                    static inline DefaultColsEnumerator enumerate_cols(const Data & data)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(data.enumerate_cols()))
                    {
                        return data.enumerate_cols();
                    }
                };
                
                template<typename MatrixType>
                class MatrixWrapper
                {
                private:
                    MatrixType & d_matrix;
                
                public:
                    enum { has_direct_access = MatrixType::has_direct_access };
                
                    typedef typename MatrixInfo<MatrixType>::Type CoeffType;
                    typedef typename MatrixInfo<MatrixType>::Type & CoeffType_Get;
                    typedef typename MatrixType::GetCoeffSteps_Type GetCoeffSteps_Type;
                    typedef MatrixType ValueType;
                    typedef MatrixType & LazyEvalType;
                    typedef typename MatrixInfo<MatrixType>::StorageTraits::pointer_type data_pointer_type;
                    typedef typename MatrixInfo<MatrixType>::StorageTraits::ref_type data_ref_type;
                
                    template<typename ResultType>
                    inline void get_coeff(ResultType & result, size_type i, size_type j) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_matrix.get_coeff(result, i, j)))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::MatrixWrapper::get_coeff(" << getAddress(*this) << "; " << getAddress(result) << " " << i << " " << j << ") const");
                        d_matrix.get_coeff(result, i, j);
                    }
                
                    static inline bool get_coeff_alwayszero() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::MatrixWrapper::get_coeff_alwayszero(" << getAddress(*this) << ")");
                        return false;
                    }
                
                    inline GetCoeffSteps_Type get_coeff_steps(size_type i, size_type j) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_matrix.get_coeff_steps(i, j)))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::MatrixWrapper::get_zero_steps(" << getAddress(*this) << "; " << i << " " << j << ")");
                        return d_matrix.get_coeff_steps(i, j);
                    }
                
                    inline CoeffType_Get operator () (size_type i, size_type j) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_matrix(i, j)))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::MatrixWrapper::operator () (" << getAddress(*this) << "; " << i << " " << j << ")");
                        return d_matrix(i, j);
                    }
                
                    inline CoeffType_Get operator [] (size_type i) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_matrix[i]))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::MatrixWrapper::operator [] (" << getAddress(*this) << "; " << i << ")");
                        return d_matrix[i];
                    }
                
                    inline size_type rows() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::MatrixWrapper::rows(" << getAddress(*this) << ")");
                        return d_matrix.rows();
                    }
                
                    inline size_type cols() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::MatrixWrapper::cols(" << getAddress(*this) << ")");
                        return d_matrix.cols();
                    }
                
                    inline size_type size() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::MatrixWrapper::size(" << getAddress(*this) << ")");
                        return d_matrix.size();
                    }
                
                    inline data_pointer_type data() const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_matrix.data()))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::MatrixWrapper::data(" << getAddress(*this) << ") const");
                        return d_matrix.data();
                    }
                
                    inline data_ref_type data(size_type i) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_matrix.data(i)))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::MatrixWrapper::data(" << getAddress(*this) << "; " << i << ")");
                        return d_matrix.data(i);
                    }
                
                    inline MatrixWrapper(MatrixType & matrix) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                        : d_matrix(matrix)
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::MatrixWrapper::MatrixWrapper(" << getAddress(*this) << "; " << getAddress(matrix) << ")");
                    }
                    
                    inline MatrixType & evaluate() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::MatrixWrapper::evaluate(" << getAddress(*this) << ")");
                        return d_matrix;
                    }
                    
                    inline LazyEvalType lazy_evaluate() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::MatrixWrapper::lazy_evaluate(" << getAddress(*this) << ")");
                        return d_matrix;
                    }
                    
                    template<typename QMatrixType>
                    inline bool involves_this_matrix(const QMatrixType & A) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::MatrixWrapper::involves_this_matrix(" << getAddress(*this) << "; " << getAddress(A) << ")");
                        return d_matrix.involves_this_matrix(A);
                    }
                    
                    template<typename MatrixType_>
                    inline bool test_involvement(const MatrixType_ & A) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::MatrixWrapper::test_involvement(" << getAddress(*this) << "; " << getAddress(A) << ")");
                        return d_matrix.test_involvement(A);
                    }
                    
                    template<typename MatrixType_>
                    void assign_resize(const MatrixType_ & mat) const
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::MatrixWrapper::assign_resize(" << getAddress(*this) << "; " << getAddress(mat) << ")");
                        d_matrix.assign_resize(mat);
                    }
                    
#if __cplusplus >= 201103L
                    template<typename MatrixType_>
                    void move_resize(MatrixType_ && mat) const
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::MatrixWrapper::move_resize(" << getAddress(*this) << "; " << getAddress(mat) << ")");
                        d_matrix.move_resize(std::move(mat));
                    }
#else
                    template<typename MatrixType_>
                    void move_resize(const MatrixType_ & mat) const
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::MatrixWrapper::move_resize(" << getAddress(*this) << "; " << getAddress(mat) << ")");
                        d_matrix.move_resize(mat);
                    }
#endif
                
                    typedef typename MatrixType::Enumerator Enumerator;
                    typedef typename MatrixType::Enumerator DefaultEnumerator;
                    typedef typename MatrixType::ConstEnumerator ConstEnumerator;
                    typedef typename MatrixType::RowEnumerator RowEnumerator;
                    typedef typename MatrixType::RowEnumerator DefaultRowEnumerator;
                    typedef typename MatrixType::ConstRowEnumerator ConstRowEnumerator;
                    typedef typename MatrixType::ColEnumerator ColEnumerator;
                    typedef typename MatrixType::ColEnumerator DefaultColEnumerator;
                    typedef typename MatrixType::ConstColEnumerator ConstColEnumerator;
                
                    inline DefaultEnumerator enumerate() const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_matrix.enumerate()))
                    {
                        return d_matrix.enumerate();
                    }
                
                    inline DefaultRowEnumerator enumerate_row(size_type row) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_matrix.enumerate_row(row)))
                    {
                        return d_matrix.enumerate_row(row);
                    }
                
                    inline DefaultColEnumerator enumerate_col(size_type col) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_matrix.enumerate_col(col)))
                    {
                        return d_matrix.enumerate_col(col);
                    }
                
                    typedef typename MatrixType::RowsEnumerator RowsEnumerator;
                    typedef typename MatrixType::ColsEnumerator ColsEnumerator;
                    typedef typename MatrixType::RowsEnumerator DefaultRowsEnumerator;
                    typedef typename MatrixType::ColsEnumerator DefaultColsEnumerator;
                    typedef typename MatrixType::ConstRowsEnumerator ConstRowsEnumerator;
                    typedef typename MatrixType::ConstColsEnumerator ConstColsEnumerator;
                
                    inline DefaultRowsEnumerator enumerate_rows() const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_matrix.enumerate_rows()))
                    {
                        return d_matrix.enumerate_rows();
                    }
                
                    inline DefaultColsEnumerator enumerate_cols() const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_matrix.enumerate_cols()))
                    {
                        return d_matrix.enumerate_cols();
                    }
                };
            
                template<typename MatrixType>
                class ConstMatrixWrapper
                {
                private:
                    const MatrixType & d_matrix;
                
                public:
                    enum { has_direct_access = MatrixType::has_direct_access };
                
                    typedef typename MatrixInfo<MatrixType>::Type CoeffType;
                    typedef const typename MatrixInfo<MatrixType>::Type & CoeffType_Get;
                    typedef typename MatrixType::GetCoeffSteps_Type GetCoeffSteps_Type;
                    typedef MatrixType ValueType;
                    typedef MatrixType & LazyEvalType;
                    typedef typename MatrixInfo<MatrixType>::StorageTraits::constpointer_type data_pointer_type;
                    typedef typename MatrixInfo<MatrixType>::StorageTraits::constref_type data_ref_type;
                
                    template<typename ResultType>
                    inline void get_coeff(ResultType & result, size_type i, size_type j) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_matrix.get_coeff(result, i, j)))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::ConstMatrixWrapper::get_coeff(" << getAddress(*this) << "; " << getAddress(result) << " " << i << " " << j << ") const");
                        d_matrix.get_coeff(result, i, j);
                    }
                
                    static inline bool get_coeff_alwayszero() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::ConstMatrixWrapper::get_coeff_alwayszero(" << getAddress(*this) << ")");
                        return false;
                    }
                
                    inline GetCoeffSteps_Type get_coeff_steps(size_type i, size_type j) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_matrix.get_coeff_steps(i, j)))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::ConstMatrixWrapper::get_zero_steps(" << getAddress(*this) << "; " << i << " " << j << ")");
                        return d_matrix.get_coeff_steps(i, j);
                    }
                
                    inline CoeffType_Get operator () (size_type i, size_type j) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_matrix(i, j)))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::ConstMatrixWrapper::operator () (" << getAddress(*this) << "; " << i << " " << j << ")");
                        return d_matrix(i, j);
                    }
                
                    inline CoeffType_Get operator [] (size_type i) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_matrix[i]))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::ConstMatrixWrapper::operator [] (" << getAddress(*this) << "; " << i << ")");
                        return d_matrix[i];
                    }
                
                    inline size_type rows() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::ConstMatrixWrapper::rows(" << getAddress(*this) << ")");
                        return d_matrix.rows();
                    }
                
                    inline size_type cols() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::ConstMatrixWrapper::cols(" << getAddress(*this) << ")");
                        return d_matrix.cols();
                    }
                
                    inline size_type size() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::ConstMatrixWrapper::size(" << getAddress(*this) << ")");
                        return d_matrix.size();
                    }
                
                    inline data_pointer_type data() const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_matrix.data()))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::ConstMatrixWrapper::data(" << getAddress(*this) << ") const");
                        return d_matrix.data();
                    }
                
                    inline data_ref_type data(size_type i) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_matrix.data(i)))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::ConstMatrixWrapper::data(" << getAddress(*this) << "; " << i << ") const");
                        return d_matrix.data(i);
                    }
                    
                    inline ConstMatrixWrapper(const MatrixType & matrix) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                        : d_matrix(matrix)
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::ConstMatrixWrapper::ConstMatrixWrapper(" << getAddress(*this) << "; " << getAddress(matrix) << ")");
                    }
                    
                    inline const MatrixType & evaluate() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::ConstMatrixWrapper::evaluate(" << getAddress(*this) << ")");
                        return d_matrix;
                    }
                    
                    inline LazyEvalType lazy_evaluate() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::ConstMatrixWrapper::lazy_evaluate(" << getAddress(*this) << ")");
                        return d_matrix;
                    }
                    
                    template<typename QMatrixType>
                    inline bool involves_this_matrix(const QMatrixType & A) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::ConstMatrixWrapper::involves_this_matrix(" << getAddress(*this) << "; " << getAddress(A) << ")");
                        return d_matrix.involves_this_matrix(A);
                    }
                
                    template<typename MatrixType_>
                    inline bool test_involvement(const MatrixType_ & A) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::ConstMatrixWrapper::test_involvement(" << getAddress(*this) << "; " << getAddress(A) << ")");
                        return d_matrix.test_involvement(A);
                    }
                
                    typedef void Enumerator;
                    typedef void RowEnumerator;
                    typedef void ColEnumerator;
                    typedef typename MatrixType::ConstEnumerator ConstEnumerator;
                    typedef typename MatrixType::ConstRowEnumerator ConstRowEnumerator;
                    typedef typename MatrixType::ConstColEnumerator ConstColEnumerator;
                    typedef typename MatrixType::ConstEnumerator DefaultEnumerator;
                    typedef typename MatrixType::ConstRowEnumerator DefaultRowEnumerator;
                    typedef typename MatrixType::ConstColEnumerator DefaultColEnumerator;
                
                    inline DefaultEnumerator enumerate() const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_matrix.enumerate()))
                    {
                        return d_matrix.enumerate();
                    }
                
                    inline DefaultRowEnumerator enumerate_row(size_type row) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_matrix.enumerate_row(row)))
                    {
                        return d_matrix.enumerate_row(row);
                    }
                
                    inline DefaultColEnumerator enumerate_col(size_type col) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_matrix.enumerate_col(col)))
                    {
                        return d_matrix.enumerate_col(col);
                    }
                
                    typedef void RowsEnumerator;
                    typedef void ColsEnumerator;
                    typedef typename MatrixType::ConstRowsEnumerator ConstRowsEnumerator;
                    typedef typename MatrixType::ConstColsEnumerator ConstColsEnumerator;
                    typedef typename MatrixType::ConstRowsEnumerator DefaultRowsEnumerator;
                    typedef typename MatrixType::ConstColsEnumerator DefaultColsEnumerator;
                
                    inline DefaultRowsEnumerator enumerate_rows() const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_matrix.enumerate_rows()))
                    {
                        return d_matrix.enumerate_rows();
                    }
                
                    inline DefaultColsEnumerator enumerate_cols() const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_matrix.enumerate_cols()))
                    {
                        return d_matrix.enumerate_cols();
                    }
                };
                
                template<typename MatrixType>
                inline MatrixWrapper<MatrixType>
                make_matrix_wrapper(MatrixType & matrix) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    PLLL_DEBUG_OUTPUT_MESSAGE("expressions::make_matrix_wrapper(" << getAddress(matrix) << ") [1]");
                    return MatrixWrapper<MatrixType>(matrix);
                }
                
                template<typename MatrixType>
                inline ConstMatrixWrapper<MatrixType>
                make_matrix_wrapper(const MatrixType & matrix) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    PLLL_DEBUG_OUTPUT_MESSAGE("expressions::make_matrix_wrapper(" << getAddress(matrix) << ") [2]");
                    return ConstMatrixWrapper<MatrixType>(matrix);
                }
                
                template<typename MatrixType>
                inline expr<identity_operation, MatrixWrapper<MatrixType> >
                make_matrix_expression(MatrixType & matrix) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    PLLL_DEBUG_OUTPUT_MESSAGE("expressions::make_matrix_expression(" << getAddress(matrix) << ") [1]");
                    return expr<identity_operation, MatrixWrapper<MatrixType> >(MatrixWrapper<MatrixType>(matrix));
                }
                
                template<typename MatrixType>
                inline expr<identity_operation, ConstMatrixWrapper<MatrixType> >
                make_matrix_expression(const MatrixType & matrix) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    PLLL_DEBUG_OUTPUT_MESSAGE("expressions::make_matrix_expression(" << getAddress(matrix) << ") [2]");
                    return expr<identity_operation, ConstMatrixWrapper<MatrixType> >(ConstMatrixWrapper<MatrixType>(matrix));
                }
                
                template<typename MatrixType>
                class MatrixTemporaryWrapper
                {
                private:
#if __cplusplus >= 201103L
                    mutable std::unique_ptr<MatrixType> d_matrix;
#else
                    mutable std::auto_ptr<MatrixType> d_matrix;
#endif
                
                public:
                    enum { has_direct_access = MatrixType::has_direct_access };
                
                    typedef typename MatrixInfo<MatrixType>::Type CoeffType;
#if __cplusplus >= 201103L
                    typedef typename MatrixInfo<MatrixType>::Type && CoeffType_Get;
#else
                    // make them accessible to allow "move by swap"
                    typedef typename MatrixInfo<MatrixType>::Type & CoeffType_Get;
#endif
                    typedef typename MatrixType::GetCoeffSteps_Type GetCoeffSteps_Type;
                    typedef MatrixType ValueType;
                    typedef const MatrixType & LazyEvalType;
                    typedef typename MatrixInfo<MatrixType>::StorageTraits::pointer_type data_pointer_type;
                    typedef typename MatrixInfo<MatrixType>::StorageTraits::ref_type data_ref_type;
                    
                    template<typename ResultType>
                    inline void get_coeff(ResultType & result, size_type i, size_type j) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_matrix->get_coeff(result, i, j)))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::MatrixTemporaryWrapper::get_coeff(" << getAddress(*this) << "; " << getAddress(result) << " " << i << " " << j << ") const");
                        d_matrix->get_coeff(result, i, j);
                    }
                    
                    static inline bool get_coeff_alwayszero() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::MatrixTemporaryWrapper::get_coeff_alwayszero(" << getAddress(*this) << ")");
                        return false;
                    }
                    
                    inline GetCoeffSteps_Type get_coeff_steps(size_type i, size_type j) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_matrix->get_coeff_steps(i, j)))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::MatrixTemporaryWrapper::get_zero_steps(" << getAddress(*this) << "; " << i << " " << j << ")");
                        return d_matrix->get_coeff_steps(i, j);
                    }
                    
                    inline CoeffType_Get operator () (size_type i, size_type j) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept((*d_matrix)(i, j)))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::MatrixTemporaryWrapper::operator () (" << getAddress(*this) << "; " << i << " " << j << ")");
#if __cplusplus >= 201103L
                        return std::move((*d_matrix)(i, j));
#else
                        return (*d_matrix)(i, j);
#endif
                    }
                    
                    inline CoeffType_Get operator [] (size_type i) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept((*d_matrix)[i]))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::MatrixTemporaryWrapper::operator [] (" << getAddress(*this) << "; " << i << ")");
#if __cplusplus >= 201103L
                        return std::move((*d_matrix)[i]);
#else
                        return (*d_matrix)[i];
#endif
                    }
                    
                    inline size_type rows() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::MatrixTemporaryWrapper::rows(" << getAddress(*this) << ")");
                        return d_matrix->rows();
                    }
                    
                    inline size_type cols() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::MatrixTemporaryWrapper::cols(" << getAddress(*this) << ")");
                        return d_matrix->cols();
                    }
                
                    inline size_type size() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::MatrixTemporaryWrapper::size(" << getAddress(*this) << ")");
                        return d_matrix->size();
                    }
                    
                    inline data_pointer_type data() const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_matrix->data()))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::MatrixTemporaryWrapper::data(" << getAddress(*this) << ") const");
                        return d_matrix->data();
                    }
                    
                    inline data_ref_type data(size_type i) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_matrix->data(i)))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::MatrixTemporaryWrapper::data(" << getAddress(*this) << "; " << i << ") const");
                        return d_matrix->data(i);
                    }
                    
                    template<typename MatrixType_>
                    inline MatrixTemporaryWrapper(const MatrixType_ & matrix)
                        : d_matrix(new MatrixType(matrix))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::MatrixTemporaryWrapper::MatrixTemporaryWrapper(" << getAddress(*this) << "; " << getAddress(matrix) << ") [1]");
                    }
                    
#if __cplusplus >= 201103L
                    inline MatrixTemporaryWrapper(const MatrixTemporaryWrapper<MatrixType> & mtw) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                        : d_matrix(std::move(mtw.d_matrix))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::MatrixTemporaryWrapper::MatrixTemporaryWrapper(" << getAddress(*this) << "; " << getAddress(mtw) << ") [2]");
                    }
                    
                    inline MatrixTemporaryWrapper(MatrixTemporaryWrapper<MatrixType> && mtw) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                        : d_matrix(std::move(mtw.d_matrix))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::MatrixTemporaryWrapper::MatrixTemporaryWrapper(" << getAddress(*this) << "; " << getAddress(mtw) << ") [3&&]");
                    }
                    
                    template<typename MatrixType_>
                    inline MatrixTemporaryWrapper(MatrixType_ && matrix)
                        : d_matrix(new MatrixType(std::move(matrix)))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::MatrixTemporaryWrapper::MatrixTemporaryWrapper(" << getAddress(*this) << "; " << getAddress(matrix) << ") [4&&]");
                    }
#else
                    inline MatrixTemporaryWrapper(const MatrixTemporaryWrapper<MatrixType> & mtw) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                        : d_matrix(mtw.d_matrix)
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::MatrixTemporaryWrapper::MatrixTemporaryWrapper(" << getAddress(*this) << "; " << getAddress(mtw) << ") [2]");
                    }
#endif
                    
#if __cplusplus >= 201103L
                    inline MatrixType && evaluate() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::MatrixTemporaryWrapper::evaluate(" << getAddress(*this) << ")");
                        return std::move(*d_matrix);
                    }
#else
                    inline const MatrixType & evaluate() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::MatrixTemporaryWrapper::evaluate(" << getAddress(*this) << ")");
                        return *d_matrix;
                    }
#endif
                    
                    inline LazyEvalType lazy_evaluate() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::MatrixTemporaryWrapper::lazy_evaluate(" << getAddress(*this) << ")");
                        return *d_matrix;
                    }
                    
                    template<typename QMatrixType>
                    static inline bool involves_this_matrix(const QMatrixType & A) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::MatrixTemporaryWrapper::involves_this_matrix(" << getAddress(*this) << "; " << getAddress(A) << ")");
                        return false;
                    }
                    
                    template<typename MatrixType_>
                    static inline bool test_involvement(const MatrixType_ & A) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::MatrixTemporaryWrapper::test_involvement(" << getAddress(*this) << "; " << getAddress(A) << ")");
                        return false;
                    }
                    
                    typedef void Enumerator;
                    typedef void RowEnumerator;
                    typedef void ColEnumerator;
                    typedef typename MatrixType::ConstEnumerator ConstEnumerator;
                    typedef typename MatrixType::ConstRowEnumerator ConstRowEnumerator;
                    typedef typename MatrixType::ConstColEnumerator ConstColEnumerator;
                    typedef typename MatrixType::ConstEnumerator DefaultEnumerator;
                    typedef typename MatrixType::ConstRowEnumerator DefaultRowEnumerator;
                    typedef typename MatrixType::ConstColEnumerator DefaultColEnumerator;
                
                    inline DefaultEnumerator enumerate() const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_matrix->enumerate()))
                    {
                        return d_matrix->enumerate();
                    }
                
                    inline DefaultRowEnumerator enumerate_row(size_type row) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_matrix->enumerate_row(row)))
                    {
                        return d_matrix->enumerate_row(row);
                    }
                
                    inline DefaultColEnumerator enumerate_col(size_type col) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_matrix->enumerate_col(col)))
                    {
                        return d_matrix->enumerate_col(col);
                    }
                
                    typedef void RowsEnumerator;
                    typedef void ColsEnumerator;
                    typedef typename MatrixType::ConstRowsEnumerator ConstRowsEnumerator;
                    typedef typename MatrixType::ConstColsEnumerator ConstColsEnumerator;
                    typedef typename MatrixType::ConstRowsEnumerator DefaultRowsEnumerator;
                    typedef typename MatrixType::ConstColsEnumerator DefaultColsEnumerator;
                
                    inline DefaultRowsEnumerator enumerate_rows() const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_matrix->enumerate_rows()))
                    {
                        return d_matrix->enumerate_rows();
                    }
                
                    inline DefaultColsEnumerator enumerate_cols() const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_matrix->enumerate_cols()))
                    {
                        return d_matrix->enumerate_cols();
                    }
                };
                
                template<template<typename> class Op, typename Data>
                inline expr<identity_operation, MatrixTemporaryWrapper<typename expr<Op, Data>::ValueType> >
                make_matrix_temporary_expression(const expr<Op, Data> & matrix)
                {
                    PLLL_DEBUG_OUTPUT_MESSAGE("expressions::make_matrix_temporary_expression(" << getAddress(matrix) << ") [1]");
                    return expr<identity_operation, MatrixTemporaryWrapper<typename expr<Op, Data>::ValueType> >
                        (MatrixTemporaryWrapper<typename expr<Op, Data>::ValueType>(matrix));
                }
                
                template<typename MatrixType>
                inline expr<identity_operation, MatrixTemporaryWrapper<MatrixType> >
                make_matrix_temporary_expression(MatrixType & matrix)
                {
                    PLLL_DEBUG_OUTPUT_MESSAGE("expressions::make_matrix_temporary_expression(" << getAddress(matrix) << ") [2]");
                    return expr<identity_operation, MatrixTemporaryWrapper<MatrixType> >(MatrixTemporaryWrapper<MatrixType>(matrix));
                }
                
                template<typename ScalarType>
                class ScalarWrapper
                {
                private:
                    const ScalarType & d_scalar;
                
                public:
                    typedef ScalarType ValueType;
                    typedef const ScalarType & LazyEvalType;
                
                    inline ScalarWrapper(const ScalarType & scalar) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                        : d_scalar(scalar)
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::ScalarWrapper::ScalarWrapper(" << getAddress(*this) << "; " << getAddress(scalar) << ")");
                    }
                
                    inline const ScalarType & evaluate() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::ScalarWrapper::evaluate(" << getAddress(*this) << ")");
                        return d_scalar;
                    }
                
                    inline LazyEvalType lazy_evaluate() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::ScalarWrapper::lazy_evaluate(" << getAddress(*this) << ")");
                        return d_scalar;
                    }
                
                    template<typename MatrixType>
                    static inline bool involves_this_matrix(const MatrixType & A) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::ScalarWrapper::involves_this_matrix(" << getAddress(*this) << "; " << getAddress(A) << ")");
                        return false;
                    }
                
                    template<typename MatrixType>
                    static inline bool test_involvement(const MatrixType & A) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::ScalarWrapper::test_involvement(" << getAddress(*this) << "; " << getAddress(A) << ")");
                        return false;
                    }
                };
                
                template<typename ScalarType>
                inline ScalarWrapper<ScalarType> make_scalar_wrapper(const ScalarType & scalar) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    PLLL_DEBUG_OUTPUT_MESSAGE("expressions::make_scalar_wrapper(" << getAddress(scalar) << ")");
                    return ScalarWrapper<ScalarType>(scalar);
                }
                
                template<typename PairMatrices>
                class matrix_matrix_multiplication
                {
                private:
                    typedef typename PairMatrices::first_type MTA;
                    typedef typename PairMatrices::second_type MTB;
                
                public:
                    typedef typename arithmetic::binary_operation<typename arithmetic::binary_operation<typename MTA::CoeffType_Get, typename MTB::CoeffType_Get,
                                                                                                        arithmetic::op::multiplication>::ResultType,
                                                                  typename arithmetic::binary_operation<typename MTA::CoeffType_Get, typename MTB::CoeffType_Get,
                                                                                                        arithmetic::op::multiplication>::ResultType,
                                                                  arithmetic::op::addition>::ResultType CoeffType;
                    typedef CoeffType CoeffType_Get;
                    typedef base_matrix<CoeffType,
                                        MatrixInfo<MTA>::rows,
                                        MatrixInfo<MTB>::cols,
                                        typename MatrixInfo<MTA>::StorageTraits,
                                        true> ValueType;
                    typedef void data_pointer_type;
                    typedef void data_ref_type;
                
                    enum { use_temporary_on_evaluate = true/*MatrixInfo<MTA>::use_temporary_on_evaluate | MatrixInfo<MTB>::use_temporary_on_evaluate*/,
                           coeffs_are_simple_expressions = false,
                           is_const = true,
                           can_assign_to = false,
                           can_move_from = false,
                           can_resize_rows = false,
                           can_resize_cols = false,
                           has_direct_access = false };
                    
                    matrix_matrix_multiplication() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_INTERNAL_STATIC_CHECK((MatrixInfo<MTA>::cols < 0) || (MatrixInfo<MTB>::rows < 0) ||
                                                   (static_cast<size_type>(MatrixInfo<MTA>::cols) == static_cast<size_type>(MatrixInfo<MTB>::rows)), FormatsDoNotMatch);
                        PLLL_INTERNAL_STATIC_CHECK(MatrixInfo<MTA>::is_matrix && MatrixInfo<MTA>::is_math_object, RequiresMathMatrix);
                        PLLL_INTERNAL_STATIC_CHECK(MatrixInfo<MTB>::is_matrix && MatrixInfo<MTB>::is_math_object, RequiresMathMatrix);
                    }
                    
                    template<typename ResultType>
                    static inline void get_coeff(ResultType & result, size_type i, size_type j, const std::pair<MTA, MTB> & AB)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(setZero(result)) &&
                                                                  noexcept(CoeffType()) &&
                                                                  (((MatrixInfo<MTA>::cols != 0) && (MatrixInfo<MTB>::rows != 0)) ?
                                                                   (noexcept(AB.first.enumerate_row(i)) && noexcept(AB.second.enumerate_col(j)) &&
                                                                    noexcept(helper::make_type_lvalue<typename MTA::ConstRowEnumerator>().next()) &&
                                                                    noexcept(helper::make_type_lvalue<typename MTB::ConstColEnumerator>().next()) &&
                                                                    noexcept(result = helper::make_type_lvalue<typename MTA::ConstRowEnumerator>().current() * helper::make_type_lvalue<typename MTB::ConstColEnumerator>().current()) &&
                                                                    noexcept(result += helper::make_type_lvalue<typename MTA::ConstRowEnumerator>().current() * helper::make_type_lvalue<typename MTB::ConstColEnumerator>().current())) : true))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::matrix_matrix_multiplication::get_coeff(" << getAddress(*this) << "; " << getAddress(result) << " " << i << " " << j << " " << getAddress(AB) << ") const");
                        assert(i < rows(AB));
                        assert(j < cols(AB));
                        size_type tt = MatrixInfo<MTA>::cols < 0 ? AB.second.rows() : AB.first.cols();
                        if (tt == 0)
                            setZero(result);
                        else
                        {
                            typename MTA::ConstRowEnumerator eA = AB.first.enumerate_row(i);
                            typename MTB::ConstColEnumerator eB = AB.second.enumerate_col(j);
                            result = eA.current() * eB.current();
                            for (eA.next(), eB.next(); eA.has_current(); eA.next(), eB.next())
                                result += eA.current() * eB.current();
                        }
                    }
                    
                    static inline bool get_coeff_alwayszero(const std::pair<MTA, MTB> & AB) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::matrix_matrix_multiplication::get_coeff_alwayszero(" << getAddress(*this) << "; " << getAddress(AB) << ")");
                        return (MatrixInfo<MTA>::cols < 0 ? AB.second.rows() : AB.first.cols()) == 0;
                    }
                    
                    class GetCoeffSteps_Type
                    {
                    protected:
                        mutable typename MTA::ConstRowEnumerator d_eA;
                        mutable typename MTB::ConstColEnumerator d_eB;
                        
                    public:
                        inline GetCoeffSteps_Type(typename MTA::ConstRowEnumerator eA, typename MTB::ConstColEnumerator eB)
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(typename MTA::ConstRowEnumerator(eA)) &&
                                                                      noexcept(typename MTB::ConstColEnumerator(eB)))
                            : d_eA(eA), d_eB(eB)
                        {
                        }
                        
                        typedef typename arithmetic::binary_operation<typename MTA::CoeffType_Get, typename MTB::CoeffType_Get,
                                                                      arithmetic::op::multiplication>::IntermediateType Step1_Type;
                        
                        inline Step1_Type step1() const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_eA.current() * d_eB.current()))
                        {
                            return d_eA.current() * d_eB.current();
                        }
                        
                        template<typename ResultType>
                        inline void step2(ResultType & result) const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(result += d_eA.current() * d_eB.current()))
                        {
                            for (d_eA.next(), d_eB.next(); d_eA.has_current(); d_eA.next(), d_eB.next())
                                result += d_eA.current() * d_eB.current();
                        }
                    };
                    
                    static inline GetCoeffSteps_Type get_coeff_steps(size_type i, size_type j, const std::pair<MTA, MTB> & AB)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(GetCoeffSteps_Type(AB.first.enumerate_row(i), AB.second.enumerate_col(j))))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::matrix_matrix_multiplication::get_coeff_steps(" << getAddress(*this) << "; " << i << " " << j << " " << getAddress(AB) << ")");
                        assert(i < rows(AB));
                        assert(j < cols(AB));
                        return GetCoeffSteps_Type(AB.first.enumerate_row(i), AB.second.enumerate_col(j));
                    }
                    
                    inline CoeffType_Get operator () (size_type i, size_type j, const std::pair<MTA, MTB> & AB) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(setZero(helper::make_type_lvalue<CoeffType>())) &&
                                                                  noexcept(CoeffType()) &&
                                                                  (((MatrixInfo<MTA>::cols != 0) && (MatrixInfo<MTB>::rows != 0)) ?
                                                                   (noexcept(AB.first.enumerate_row(i)) && noexcept(AB.second.enumerate_col(j)) &&
                                                                    noexcept(helper::make_type_lvalue<typename MTA::ConstRowEnumerator>().next()) &&
                                                                    noexcept(helper::make_type_lvalue<typename MTB::ConstColEnumerator>().next()) &&
                                                                    noexcept(helper::make_type_lvalue<CoeffType>() = helper::make_type_lvalue<typename MTA::ConstRowEnumerator>().current() * helper::make_type_lvalue<typename MTB::ConstColEnumerator>().current()) &&
                                                                    noexcept(helper::make_type_lvalue<CoeffType>() += helper::make_type_lvalue<typename MTA::ConstRowEnumerator>().current() * helper::make_type_lvalue<typename MTB::ConstColEnumerator>().current())) : true))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::matrix_matrix_multiplication::operator () (" << getAddress(*this) << "; " << i << " " << j << " " << getAddress(AB) << ")");
                        assert(i < rows(AB));
                        assert(j < cols(AB));
                        size_type tt = MatrixInfo<MTA>::cols < 0 ? AB.second.rows() : AB.first.cols();
                        if (tt == 0)
                        {
                            CoeffType res;
                            setZero(res);
                            return res;
                        }
                        else
                        {
                            typename MTA::ConstRowEnumerator eA = AB.first.enumerate_row(i);
                            typename MTB::ConstColEnumerator eB = AB.second.enumerate_col(j);
                            CoeffType res = eA.current() * eB.current();
                            for (eA.next(), eB.next(); eA.has_current(); eA.next(), eB.next())
                                res += eA.current() * eB.current();
                            return res;
                        }
                    }
                    
                    static inline size_type rows(const std::pair<MTA, MTB> & AB) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::matrix_matrix_multiplication::rows(" << getAddress(*this) << "; " << getAddress(AB) << ")");
                        return AB.first.rows();
                    }
                    
                    static inline size_type cols(const std::pair<MTA, MTB> & AB) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::matrix_matrix_multiplication::cols(" << getAddress(*this) << "; " << getAddress(AB) << ")");
                        return AB.second.cols();
                    }
                    
                    static inline size_type size(const std::pair<MTA, MTB> & AB) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::matrix_matrix_multiplication::size(" << getAddress(*this) << "; " << getAddress(AB) << ")");
                        return AB.first.rows() * AB.second.cols();
                    }
                    
                    template<typename MatrixType>
                    static inline bool involves_this_matrix(const MatrixType & matrix, const std::pair<MTA, MTB> & AB) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::matrix_matrix_multiplication::involves_this_matrix(" << getAddress(*this) << "; " << getAddress(matrix) << " " << getAddress(AB) << ")");
                        return AB.first.involves_this_matrix(matrix) || AB.second.involves_this_matrix(matrix);
                    }
                    
                    template<typename MatrixType>
                    static inline bool test_involvement(const MatrixType & matrix, const std::pair<MTA, MTB> & AB) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::matrix_matrix_multiplication::test_involvement(" << getAddress(*this) << "; " << getAddress(matrix) << " " << getAddress(AB) << ")");
                        return AB.first.test_involvement(matrix) || AB.second.test_involvement(matrix);
                    }
                    
                    class ConstRowEnumerator
                    {
                    private:
                        typename MTA::ConstRowEnumerator d_row;
                        typename MTB::ConstColsEnumerator d_cols;
                    
                    public:
                        typedef typename matrix_matrix_multiplication<std::pair<MTA, MTB> >::CoeffType_Get Type;
                        typedef typename matrix_matrix_multiplication<std::pair<MTA, MTB> >::GetCoeffSteps_Type GetCoeffSteps_Type;
                        
                        inline ConstRowEnumerator(typename MTA::ConstRowEnumerator row, typename MTB::ConstColsEnumerator cols)
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(typename MTA::ConstRowEnumerator(row)) &&
                                                                      noexcept(typename MTB::ConstColsEnumerator(cols)))
                            : d_row(row), d_cols(cols)
                        {
                        }
                        
                        inline Type current() const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(setZero(helper::make_type_lvalue<Type>())) &&
                                                                      noexcept(helper::make_type_lvalue<Type>() = d_row.current() * d_cols.current().current()) &&
                                                                      noexcept(helper::make_type_lvalue<Type>() += d_row.current() * d_cols.current().current()))
                        {
                            typename MTA::ConstRowEnumerator er = d_row;
                            if (er.has_current())
                            {
                                typename MTB::ConstColsEnumerator::Type ec = d_cols.current();
                                Type result = er.current() * ec.current();
                                for (er.next(), ec.next(); er.has_current(); er.next(), ec.next())
                                    result += er.current() * ec.current();
                                return result;
                            }
                            else
                            {
                                Type result;
                                setZero(result);
                                return result;
                            }
                        }
                        
                        template<typename Result>
                        inline void get_current(Result & result) const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(setZero(result)) &&
                                                                      noexcept(result = d_row.current() * d_cols.current().current()) &&
                                                                      noexcept(result += d_row.current() * d_cols.current().current()))
                        {
                            typename MTA::ConstRowEnumerator er = d_row;
                            if (er.has_current())
                            {
                                typename MTB::ConstColsEnumerator::Type ec = d_cols.current();
                                result = er.current() * ec.current();
                                for (er.next(), ec.next(); er.has_current(); er.next(), ec.next())
                                    result += er.current() * ec.current();
                            }
                            else
                                setZero(result);
                        }
                        
                        inline GetCoeffSteps_Type get_current_steps() const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(GetCoeffSteps_Type(d_row, d_cols.current())))
                        {
                            return GetCoeffSteps_Type(d_row, d_cols.current());
                        }
                        
                        inline bool has_current() const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_cols.has_current()))
                        {
                            return d_cols.has_current();
                        }
                        
                        inline void next()
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_cols.next()))
                        {
                            d_cols.next();
                        }
                    };
                    
                    class ConstColEnumerator
                    {
                    private:
                        typename MTA::ConstRowsEnumerator d_rows;
                        typename MTB::ConstColEnumerator d_col;
                        
                    public:
                        typedef typename matrix_matrix_multiplication<std::pair<MTA, MTB> >::CoeffType_Get Type;
                        typedef typename matrix_matrix_multiplication<std::pair<MTA, MTB> >::GetCoeffSteps_Type GetCoeffSteps_Type;
                        
                        inline ConstColEnumerator(typename MTA::ConstRowsEnumerator rows, typename MTB::ConstColEnumerator col)
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(typename MTA::ConstRowsEnumerator(rows)) &&
                                                                      noexcept(typename MTB::ConstColEnumerator(col)))
                            : d_rows(rows), d_col(col)
                        {
                        }
                        
                        inline Type current() const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(setZero(helper::make_type_lvalue<Type>())) &&
                                                                      noexcept(helper::make_type_lvalue<Type>() = d_rows.current().current() * d_col.current()) &&
                                                                      noexcept(helper::make_type_lvalue<Type>() += d_rows.current().current() * d_col.current()))
                        {
                            typename MTA::ConstRowsEnumerator::Type er = d_rows.current();
                            if (er.has_current())
                            {
                                typename MTB::ConstColEnumerator ec = d_col;
                                Type result = er.current() * ec.current();
                                for (er.next(), ec.next(); er.has_current(); er.next(), ec.next())
                                    result += er.current() * ec.current();
                                return result;
                            }
                            else
                            {
                                Type result;
                                setZero(result);
                                return result;
                            }
                        }
                        
                        template<typename Result>
                        inline void get_current(Result & result) const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(setZero(result)) &&
                                                                      noexcept(result = d_rows.current().current() * d_col.current()) &&
                                                                      noexcept(result += d_rows.current().current() * d_col.current()))
                        {
                            typename MTA::ConstRowsEnumerator::Type er = d_rows.current();
                            if (er.has_current())
                            {
                                typename MTB::ConstColEnumerator ec = d_col;
                                result = er.current() * ec.current();
                                for (er.next(), ec.next(); er.has_current(); er.next(), ec.next())
                                    result += er.current() * ec.current();
                            }
                            else
                                setZero(result);
                        }
                        
                        inline GetCoeffSteps_Type get_current_steps() const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(GetCoeffSteps_Type(d_rows.current(), d_col)))
                        {
                            return GetCoeffSteps_Type(d_rows.current(), d_col);
                        }
                        
                        inline bool has_current() const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_rows.has_current()))
                        {
                            return d_rows.has_current();
                        }
                        
                        inline void next()
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_rows.next()))
                        {
                            d_rows.next();
                        }
                    };
                    
                    class ConstEnumerator_Generic
                    {
                    private:
                        typename MTA::ConstRowsEnumerator d_rows;
                        typename MTB::ConstColsEnumerator d_cols, d_cols_init;
                        
                    public:
                        typedef typename matrix_matrix_multiplication<std::pair<MTA, MTB> >::CoeffType_Get Type;
                        typedef typename matrix_matrix_multiplication<std::pair<MTA, MTB> >::GetCoeffSteps_Type GetCoeffSteps_Type;
                        
                        inline ConstEnumerator_Generic(typename MTA::ConstRowsEnumerator rows, typename MTB::ConstColsEnumerator cols)
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(typename MTA::ConstRowsEnumerator(rows)) &&
                                                                      noexcept(typename MTB::ConstColsEnumerator(cols)))
                            : d_rows(rows), d_cols(cols), d_cols_init(d_cols)
                        {
                        }
                        
                        inline Type current() const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(setZero(helper::make_type_lvalue<Type>())) &&
                                                                      noexcept(helper::make_type_lvalue<Type>() = d_rows.current().current() * d_cols.current().current()) &&
                                                                      noexcept(helper::make_type_lvalue<Type>() += d_rows.current().current() * d_cols.current().current()))
                        {
                            typename MTA::ConstRowsEnumerator::Type er = d_rows.current();
                            if (er.has_current())
                            {
                                typename MTB::ConstColsEnumerator::Type ec = d_cols.current();
                                Type result = er.current() * ec.current();
                                for (er.next(), ec.next(); er.has_current(); er.next(), ec.next())
                                    result += er.current() * ec.current();
                                return result;
                            }
                            else
                            {
                                Type result;
                                setZero(result);
                                return result;
                            }
                        }
                        
                        template<typename Result>
                        inline void get_current(Result & result) const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(setZero(result)) &&
                                                                      noexcept(result = d_rows.current().current() * d_cols.current().current()) &&
                                                                      noexcept(result += d_rows.current().current() * d_cols.current().current()))
                        {
                            typename MTA::ConstRowsEnumerator::Type er = d_rows.current();
                            if (er.has_current())
                            {
                                typename MTB::ConstColsEnumerator::Type ec = d_cols.current();
                                result = er.current() * ec.current();
                                for (er.next(), ec.next(); er.has_current(); er.next(), ec.next())
                                    result += er.current() * ec.current();
                            }
                            else
                                setZero(result);
                        }
                        
                        inline GetCoeffSteps_Type get_current_steps() const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(GetCoeffSteps_Type(d_rows.current(), d_cols.current())))
                        {
                            return GetCoeffSteps_Type(d_rows.current(), d_cols.current());
                        }
                        
                        inline bool has_current() const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_cols.has_current()))
                        {
                            return d_cols.has_current();
                        }
                        
                        inline void next()
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_cols.next()) && noexcept(d_cols.has_current()) &&
                                                                      noexcept(d_rows.next()) && noexcept(d_rows.has_current()) &&
                                                                      noexcept(d_cols = d_cols_init))
                        {
                            d_cols.next();
                            if (!d_cols.has_current())
                            {
                                d_rows.next();
                                if (d_rows.has_current())
                                    d_cols = d_cols_init;
                            }
                        }
                    };
                    
                    class ConstEnumerator_OneRow : public ConstRowEnumerator
                    {
                    public:
                        inline ConstEnumerator_OneRow(typename MTA::ConstRowsEnumerator rows, typename MTB::ConstColsEnumerator cols)
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(rows.current()))
                            : ConstRowEnumerator(rows.current(), cols)
                        {
                        }
                    };
                    
                    class ConstEnumerator_OneCol : public ConstColEnumerator
                    {
                    public:
                        inline ConstEnumerator_OneCol(typename MTA::ConstRowsEnumerator rows, typename MTB::ConstColsEnumerator cols)
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(cols.current()))
                            : ConstColEnumerator(rows, cols.current())
                        {
                        }
                    };
                    
                    typedef typename helper::SelectFirstType<MatrixInfo<MTB>::cols == 1,
                                                             ConstEnumerator_OneCol,
                                                             typename helper::SelectFirstType<MatrixInfo<MTA>::rows == 1,
                                                                                              ConstEnumerator_OneRow,
                                                                                              ConstEnumerator_Generic>::result>::result ConstEnumerator;
                    
                    typedef void Enumerator;
                    typedef void RowEnumerator;
                    typedef void ColEnumerator;
                    
                    typedef ConstEnumerator DefaultEnumerator;
                    typedef ConstRowEnumerator DefaultRowEnumerator;
                    typedef ConstColEnumerator DefaultColEnumerator;
                    
                    static inline DefaultEnumerator enumerate(const std::pair<MTA, MTB> & AB)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(AB.first.enumerate_rows()) && noexcept(AB.second.enumerate_cols()))
                    {
                        return DefaultEnumerator(AB.first.enumerate_rows(), AB.second.enumerate_cols());
                    }
                    
                    static inline DefaultRowEnumerator enumerate_row(size_type row, const std::pair<MTA, MTB> & AB)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(AB.second.enumerate_row(row)) && noexcept(AB.second.enumerate_cols()))
                    {
                        return DefaultRowEnumerator(AB.first.enumerate_row(row), AB.second.enumerate_cols());
                    }
                    
                    static inline DefaultColEnumerator enumerate_col(size_type col, const std::pair<MTA, MTB> & AB)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(AB.first.enumerate_rows()) && noexcept(AB.second.enumerate_col(col)))
                    {
                        return DefaultColEnumerator(AB.first.enumerate_rows(), AB.second.enumerate_col(col));
                    }
                    
                    class ConstRowsEnumerator
                    {
                    private:
                        typename MTA::ConstRowsEnumerator d_rows;
                        typename MTB::ConstColsEnumerator d_cols;
                        
                    public:
                        typedef ConstRowEnumerator Type;
                        
                        ConstRowsEnumerator(const std::pair<MTA, MTB> & AB)
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(AB.first.enumerate_rows()) && noexcept(AB.second.enumerate_cols()))
                            : d_rows(AB.first.enumerate_rows()), d_cols(AB.second.enumerate_cols())
                        {
                        }
                        
                        Type current() const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(ConstRowEnumerator(d_rows.current(), d_cols)))
                        {
                            return ConstRowEnumerator(d_rows.current(), d_cols);
                        }
                        
                        inline bool has_current() const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_rows.has_current()))
                        {
                            return d_rows.has_current();
                        }
                        
                        inline void next() PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_rows.next()))
                        {
                            d_rows.next();
                        }
                    };
                    
                    class ConstColsEnumerator
                    {
                    private:
                        typename MTA::ConstRowsEnumerator d_rows;
                        typename MTB::ConstColsEnumerator d_cols;
                        
                    public:
                        typedef ConstColEnumerator Type;
                        
                        ConstColsEnumerator(const std::pair<MTA, MTB> & AB)
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(AB.first.enumerate_rows()) && noexcept(AB.second.enumerate_cols()))
                            : d_rows(AB.first.enumerate_rows()), d_cols(AB.second.enumerate_cols())
                        {
                        }
                        
                        Type current() const PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(ConstColEnumerator(d_rows, d_cols.current())))
                        {
                            return ConstColEnumerator(d_rows, d_cols.current());
                        }
                        
                        inline bool has_current() const PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_cols.has_current()))
                        {
                            return d_cols.has_current();
                        }
                        
                        inline void next() PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_cols.next()))
                        {
                            d_cols.next();
                        }
                    };
                    
                    typedef void RowsEnumerator;
                    typedef void ColsEnumerator;
                    typedef ConstRowsEnumerator DefaultRowsEnumerator;
                    typedef ConstColsEnumerator DefaultColsEnumerator;
                    
                    inline DefaultRowsEnumerator enumerate_rows(const std::pair<MTA, MTB> & AB) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(DefaultRowsEnumerator(AB)))
                    {
                        return DefaultRowsEnumerator(AB);
                    }
                    
                    inline DefaultColsEnumerator enumerate_cols(const std::pair<MTA, MTB> & AB) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(DefaultColsEnumerator(AB)))
                    {
                        return DefaultColsEnumerator(AB);
                    }
                };
                
                template<typename OpA, typename OpB>
                class MSMul
                {
                public:
                    typedef typename arithmetic::binary_operation<OpA, OpB, arithmetic::op::multiplication>::IntermediateType ReturnType;
                    typedef typename arithmetic::binary_operation<OpA, OpB, arithmetic::op::multiplication>::ResultType ResultType;
                    
                    static inline ReturnType op(const OpA & a, const OpB & b)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(ReturnType(a * b)))
                    {
                        return a * b;
                    }
                    
                    template<typename ResultType_>
                    static inline void op(ResultType_ & result, const OpA & a, const OpB & b)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(result = a * b))
                    {
                        result = a * b;
                    }
                };
                
                template<typename OpA, typename OpB>
                class SMMul
                {
                public:
                    typedef typename arithmetic::binary_operation<OpB, OpA, arithmetic::op::multiplication>::IntermediateType ReturnType;
                    typedef typename arithmetic::binary_operation<OpB, OpA, arithmetic::op::multiplication>::ResultType ResultType;
                    
                    static inline ReturnType op(const OpA & a, const OpB & b)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(ReturnType(b * a)))
                    {
                        return b * a;
                    }
                    
                    template<typename ResultType_>
                    static inline void op(ResultType_ & result, const OpA & a, const OpB & b)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(result = b * a))
                    {
                        result = b * a;
                    }
                };
                
                template<typename OpA, typename OpB>
                class MSDiv
                {
                public:
                    typedef typename arithmetic::binary_operation<OpA, OpB, arithmetic::op::division>::IntermediateType ReturnType;
                    typedef typename arithmetic::binary_operation<OpA, OpB, arithmetic::op::division>::ResultType ResultType;
                    
                    static inline ReturnType op(const OpA & a, const OpB & b)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(ReturnType(a / b)))
                    {
                        return a / b;
                    }
                    
                    template<typename ResultType_>
                    static inline void op(ResultType_ & result, const OpA & a, const OpB & b)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(result = a / b))
                    {
                        result = a / b;
                    }
                };
                
                template<typename OpA, typename OpB>
                class MSMod
                {
                public:
                    typedef typename arithmetic::binary_operation<OpA, OpB, arithmetic::op::modulo>::IntermediateType ReturnType;
                    typedef typename arithmetic::binary_operation<OpA, OpB, arithmetic::op::modulo>::ResultType ResultType;
                    
                    static inline ReturnType op(const OpA & a, const OpB & b)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(ReturnType(a % b)))
                    {
                        return a % b;
                    }
                    
                    template<typename ResultType_>
                    static inline void op(ResultType_ & result, const OpA & a, const OpB & b)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(result = a % b))
                    {
                        result = a % b;
                    }
                };
                
                template<template<typename OpA, typename OpB> class Operation>
                class matrix_scalar_operation
                {
                public:
                    template<typename MatrixScalarPair>
                    class operation
                    {
                    private:
                        typedef typename MatrixScalarPair::first_type MT;
                        typedef typename MatrixScalarPair::second_type ST;
                        typedef Operation<typename plll::helper::remove_decorations<typename MT::CoeffType_Get>::Result,
                                          typename plll::helper::remove_decorations<typename ST::ValueType>::Result> Op;
                        
                    public:
                        typedef typename Op::ResultType CoeffType;
                        typedef typename Op::ReturnType CoeffType_Get;
                        typedef base_matrix<CoeffType,
                                            MatrixInfo<MT>::rows,
                                            MatrixInfo<MT>::cols,
                                            typename MatrixInfo<MT>::StorageTraits,
                                            true> ValueType;
                        typedef void data_pointer_type;
                        typedef void data_ref_type;
                        
                        enum { use_temporary_on_evaluate = MatrixInfo<MT>::use_temporary_on_evaluate,
                               coeffs_are_simple_expressions = true,
                               is_const = true,
                               can_assign_to = false,
                               can_move_from = false,
                               can_resize_rows = false,
                               can_resize_cols = false,
                               has_direct_access = false };
                        
                        operation() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                        {
                            PLLL_INTERNAL_STATIC_CHECK(MatrixInfo<MT>::is_matrix && MatrixInfo<MT>::is_math_object, RequiresMathMatrix);
                        }
                        
                        template<typename ResultType>
                        static inline void get_coeff(ResultType & result, size_type i, size_type j, const std::pair<MT, ST> & AB)
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Op::op(result, AB.first(i, j), AB.second.evaluate())))
                        {
                            PLLL_DEBUG_OUTPUT_MESSAGE("expressions::matrix_scalar_operation::operation::get_coeff(" << getAddress(*this) << "; " << getAddress(result) << " " << i << " " << j << " " << getAddress(AB) << ") const");
                            assert(i < rows(AB));
                            assert(j < cols(AB));
                            Op::op(result, AB.first(i, j), AB.second.evaluate());
                        }
                        
                        static inline bool get_coeff_alwayszero(const std::pair<MT, ST> & AB) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                        {
                            PLLL_DEBUG_OUTPUT_MESSAGE("expressions::matrix_scalar_operation::operation::get_coeff_alwayszero(" << getAddress(*this) << "; " << getAddress(AB) << ")");
                            return false;
                        }
                        
                        class GetCoeffSteps_Type
                        {
                        protected:
                            typename Op::ReturnType d_value;
                            
                        public:
                            inline GetCoeffSteps_Type(const typename Op::ReturnType & value)
                                PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(typename Op::ReturnType(value)))
                                : d_value(value)
                            {
                            }
                            
#if __cplusplus >= 201103L
                            inline GetCoeffSteps_Type(typename Op::ReturnType && value)
                                PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(typename Op::ReturnType(std::move(value))))
                                : d_value(std::move(value))
                            {
                            }
#endif
                            
                            typedef const typename Op::ReturnType & Step1_Type;
                            
                            inline Step1_Type step1() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                            {
                                return d_value;
                            }
                            
                            template<typename ResultType>
                            inline void step2(ResultType & result) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                            {
                            }
                        };
                        
                        static inline GetCoeffSteps_Type get_coeff_steps(size_type i, size_type j, const std::pair<MT, ST> & AB)
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(GetCoeffSteps_Type(Op::op(AB.first(i, j), AB.second.evaluate()))))
                        {
                            PLLL_DEBUG_OUTPUT_MESSAGE("expressions::matrix_scalar_operation::operation::get_coeff_steps(" << getAddress(*this) << "; " << i << " " << j << " " << getAddress(AB) << ")");
                            assert(i < rows(AB));
                            assert(j < cols(AB));
                            return GetCoeffSteps_Type(Op::op(AB.first(i, j), AB.second.evaluate()));
                        }
                        
                        inline CoeffType_Get operator () (size_type i, size_type j, const std::pair<MT, ST> & AB) const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Op::op(AB.first(i, j), AB.second.evaluate())))
                        {
                            PLLL_DEBUG_OUTPUT_MESSAGE("expressions::matrix_scalar_operation::operation::operator () (" << getAddress(*this) << "; " << i << " " << j << " " << getAddress(AB) << ")");
                            assert(i < rows(AB));
                            assert(j < cols(AB));
                            return Op::op(AB.first(i, j), AB.second.evaluate());
                        }
                        
                        static inline size_type rows(const std::pair<MT, ST> & AB) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                        {
                            PLLL_DEBUG_OUTPUT_MESSAGE("expressions::matrix_scalar_operation::operation::rows(" << getAddress(*this) << "; " << getAddress(AB) << ")");
                            return AB.first.rows();
                        }
                        
                        static inline size_type cols(const std::pair<MT, ST> & AB) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                        {
                            PLLL_DEBUG_OUTPUT_MESSAGE("expressions::matrix_scalar_operation::operation::cols(" << getAddress(*this) << "; " << getAddress(AB) << ")");
                            return AB.first.cols();
                        }
                        
                        static inline size_type size(const std::pair<MT, ST> & AB) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                        {
                            PLLL_DEBUG_OUTPUT_MESSAGE("expressions::matrix_scalar_operation::operation::size(" << getAddress(*this) << "; " << getAddress(AB) << ")");
                            return AB.first.size();
                        }
                        
                        template<typename MatrixType>
                        static inline bool involves_this_matrix(const MatrixType & matrix, const std::pair<MT, ST> & AB) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                        {
                            PLLL_DEBUG_OUTPUT_MESSAGE("expressions::matrix_scalar_operation::operation::involves_this_matrix(" << getAddress(*this) << "; " << getAddress(matrix) << " " << getAddress(AB) << ")");
                            return AB.first.involves_this_matrix(matrix);
                        }
                        
                        template<typename MatrixType>
                        static inline bool test_involvement(const MatrixType & matrix, const std::pair<MT, ST> & AB) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                        {
                            PLLL_DEBUG_OUTPUT_MESSAGE("expressions::matrix_scalar_operation::operation::test_involvement(" << getAddress(*this) << "; " << getAddress(matrix) << " " << getAddress(AB) << ")");
                            return AB.first.test_involvement(matrix);
                        }
                        
                        template<typename Enum>
                        class MSEnumerator
                        {
                        private:
                            Enum d_enum;
                            const ST & d_B;
                            
                        public:
                            typedef CoeffType_Get Type;
                            typedef operation<MatrixScalarPair>::GetCoeffSteps_Type GetCoeffSteps_Type;
                            
                            inline MSEnumerator(const Enum & e, const ST & B)
                                PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Enum(e)))
                                : d_enum(e), d_B(B)
                            {
                            }
                            
#if __cplusplus >= 201103L
                            inline MSEnumerator(Enum && e, const ST & B)
                                PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Enum(std::move(e))))
                                : d_enum(std::move(e)), d_B(B)
                            {
                            }
#endif
                            
                            inline Type current() const
                                PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Op::op(d_enum.current(), d_B.evaluate())))
                            {
                                return Op::op(d_enum.current(), d_B.evaluate());
                            }
                            
                            template<typename Result>
                            inline void get_current(Result & result) const
                                PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Op::op(result, d_enum.current(), d_B.evaluate())))
                            {
                                Op::op(result, d_enum.current(), d_B.evaluate());
                            }
                            
                            inline GetCoeffSteps_Type get_current_steps() const
                                PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(GetCoeffSteps_Type(Op::op(d_enum.current(), d_B.evaluate()))))
                            {
                                return GetCoeffSteps_Type(Op::op(d_enum.current(), d_B.evaluate()));
                            }
                            
                            inline bool has_current() const
                                PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_enum.has_current()))
                            {
                                return d_enum.has_current();
                            }
                            
                            inline void next()
                                PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_enum.next()))
                            {
                                d_enum.next();
                            }
                        };
                        
                        typedef void Enumerator;
                        typedef MSEnumerator<typename MT::ConstEnumerator> ConstEnumerator;
                        typedef MSEnumerator<typename MT::ConstEnumerator> DefaultEnumerator;
                        typedef void RowEnumerator;
                        typedef MSEnumerator<typename MT::ConstRowEnumerator> ConstRowEnumerator;
                        typedef MSEnumerator<typename MT::ConstRowEnumerator> DefaultRowEnumerator;
                        typedef void ColEnumerator;
                        typedef MSEnumerator<typename MT::ConstColEnumerator> ConstColEnumerator;
                        typedef MSEnumerator<typename MT::ConstColEnumerator> DefaultColEnumerator;
                        
                        static inline DefaultEnumerator enumerate(const std::pair<MT, ST> & AB)
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(DefaultEnumerator(AB.first.enumerate(), AB.second)))
                        {
                            return DefaultEnumerator(AB.first.enumerate(), AB.second);
                        }
                        
                        static inline DefaultRowEnumerator enumerate_row(size_type row, const std::pair<MT, ST> & AB)
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(DefaultRowEnumerator(AB.first.enumerate_row(row), AB.second)))
                        {
                            return DefaultRowEnumerator(AB.first.enumerate_row(row), AB.second);
                        }
                        
                        static inline DefaultColEnumerator enumerate_col(size_type col, const std::pair<MT, ST> & AB)
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(DefaultColEnumerator(AB.first.enumerate_col(col), AB.second)))
                        {
                            return DefaultColEnumerator(AB.first.enumerate_col(col), AB.second);
                        }
                        
                        template<typename Enum_meta, typename Enum>
                        class MSEnumerator_meta
                        {
                        private:
                            Enum_meta d_enum;
                            const ST & d_B;
                            
                        public:
                            typedef MSEnumerator<Enum> Type;
                            
                            inline MSEnumerator_meta(const Enum_meta & e, const ST & B)
                                PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Enum_meta(e)))
                                : d_enum(e), d_B(B)
                            {
                            }
                            
#if __cplusplus >= 201103L
                            inline MSEnumerator_meta(Enum_meta && e, const ST & B)
                                PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Enum_meta(std::move(e))))
                                : d_enum(std::move(e)), d_B(B)
                            {
                            }
#endif
                            
                            inline Type current() const
                                PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Type(d_enum.current(), d_B)))
                            {
                                return Type(d_enum.current(), d_B);
                            }
                            
                            inline bool has_current() const
                                PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_enum.has_current()))
                            {
                                return d_enum.has_current();
                            }
                            
                            inline void next()
                                PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_enum.next()))
                            {
                                d_enum.next();
                            }
                        };
                        
                        typedef void RowsEnumerator;
                        typedef void ColsEnumerator;
                        typedef MSEnumerator_meta<typename MT::ConstRowsEnumerator, typename MT::ConstRowEnumerator> ConstRowsEnumerator;
                        typedef MSEnumerator_meta<typename MT::ConstColsEnumerator, typename MT::ConstColEnumerator> ConstColsEnumerator;
                        typedef ConstRowsEnumerator DefaultRowsEnumerator;
                        typedef ConstColsEnumerator DefaultColsEnumerator;
                        
                        inline DefaultRowsEnumerator enumerate_rows(const std::pair<MT, ST> & AB) const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(DefaultRowsEnumerator(AB.first.enumerate_rows(), AB.second)))
                        {
                            return DefaultRowsEnumerator(AB.first.enumerate_rows(), AB.second);
                        }
                        
                        inline DefaultColsEnumerator enumerate_cols(const std::pair<MT, ST> & AB) const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(DefaultColsEnumerator(AB.first.enumerate_cols(), AB.second)))
                        {
                            return DefaultColsEnumerator(AB.first.enumerate_cols(), AB.second);
                        }
                    };
                };
                
                template<typename OpA, typename OpB>
                class CWAdd
                {
                public:
                    typedef typename arithmetic::binary_operation<OpA, OpB, arithmetic::op::addition>::IntermediateType ReturnType;
                    typedef typename arithmetic::binary_operation<OpA, OpB, arithmetic::op::addition>::ResultType ResultType;
                    
                    static inline ReturnType op(const OpA & a, const OpB & b)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(ReturnType(a + b)))
                    {
                        return a + b;
                    }
                    
                    template<typename ResultType_>
                    static inline void op(ResultType_ & result, const OpA & a, const OpB & b)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(result = a + b))
                    {
                        result = a + b;
                    }
                };
                
                template<typename OpA, typename OpB>
                class CWSub
                {
                public:
                    typedef typename arithmetic::binary_operation<OpA, OpB, arithmetic::op::subtraction>::IntermediateType ReturnType;
                    typedef typename arithmetic::binary_operation<OpA, OpB, arithmetic::op::subtraction>::ResultType ResultType;
                    
                    static inline ReturnType op(const OpA & a, const OpB & b)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(ReturnType(a - b)))
                    {
                        return a - b;
                    }
                    
                    template<typename ResultType_>
                    static inline void op(ResultType_ & result, const OpA & a, const OpB & b)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(result = a - b))
                    {
                        result = a - b;
                    }
                };
                
                template<typename OpA, typename OpB>
                class CWMul
                {
                public:
                    typedef typename arithmetic::binary_operation<OpA, OpB, arithmetic::op::multiplication>::IntermediateType ReturnType;
                    typedef typename arithmetic::binary_operation<OpA, OpB, arithmetic::op::multiplication>::ResultType ResultType;
                    
                    static inline ReturnType op(const OpA & a, const OpB & b)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(ReturnType(a * b)))
                    {
                        return a * b;
                    }
                    
                    template<typename ResultType_>
                    static inline void op(ResultType_ & result, const OpA & a, const OpB & b)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(result = a * b))
                    {
                        result = a * b;
                    }
                };
                
                template<typename OpA, typename OpB>
                class CWDiv
                {
                public:
                    typedef typename arithmetic::binary_operation<OpA, OpB, arithmetic::op::division>::IntermediateType ReturnType;
                    typedef typename arithmetic::binary_operation<OpA, OpB, arithmetic::op::division>::ResultType ResultType;
                    
                    static inline ReturnType op(const OpA & a, const OpB & b)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(ReturnType(a / b)))
                    {
                        return a / b;
                    }
                    
                    template<typename ResultType_>
                    static inline void op(ResultType_ & result, const OpA & a, const OpB & b)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(result = a / b))
                    {
                        result = a / b;
                    }
                };
                
                template<typename OpA, typename OpB>
                class CWMod
                {
                public:
                    typedef typename arithmetic::binary_operation<OpA, OpB, arithmetic::op::modulo>::IntermediateType ReturnType;
                    typedef typename arithmetic::binary_operation<OpA, OpB, arithmetic::op::modulo>::ResultType ResultType;
                    
                    static inline ReturnType op(const OpA & a, const OpB & b)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(ReturnType(a % b)))
                    {
                        return a % b;
                    }
                    
                    template<typename ResultType_>
                    static inline void op(ResultType_ & result, const OpA & a, const OpB & b)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(result = a % b))
                    {
                        result = a % b;
                    }
                };
                
                template<template<typename OpA, typename OpB> class Operation>
                class componentwise_operation
                {
                public:
                    template<typename PairMatrices>
                    class operation
                    {
                    private:
                        typedef typename PairMatrices::first_type MTA;
                        typedef typename PairMatrices::second_type MTB;
                        typedef Operation<typename plll::helper::remove_decorations<typename MTA::CoeffType_Get>::Result,
                                          typename plll::helper::remove_decorations<typename MTB::CoeffType_Get>::Result> Op;
                        
                    public:
                        typedef typename Op::ResultType CoeffType;
                        typedef typename Op::ReturnType CoeffType_Get;
                        typedef base_matrix<CoeffType,
                                            MatrixInfo<MTA>::rows < 0 ? static_cast<int>(MatrixInfo<MTB>::rows) : static_cast<int>(MatrixInfo<MTA>::rows),
                                                                    MatrixInfo<MTA>::cols < 0 ? static_cast<int>(MatrixInfo<MTB>::cols) : static_cast<int>(MatrixInfo<MTA>::cols),
                                                                                            typename MatrixInfo<MTA>::StorageTraits,
                                                                                            true> ValueType;
                        typedef void data_pointer_type;
                        typedef void data_ref_type;
                        
                        enum { use_temporary_on_evaluate = MatrixInfo<MTA>::use_temporary_on_evaluate | MatrixInfo<MTB>::use_temporary_on_evaluate,
                               coeffs_are_simple_expressions = true,
                               is_const = true,
                               can_assign_to = false,
                               can_move_from = false,
                               can_resize_rows = false,
                               can_resize_cols = false,
                               has_direct_access = false };
                        
                        operation() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                        {
                            PLLL_INTERNAL_STATIC_CHECK((MatrixInfo<MTA>::cols < 0) || (MatrixInfo<MTB>::cols < 0) ||
                                                       (static_cast<size_type>(MatrixInfo<MTA>::cols) == static_cast<size_type>(MatrixInfo<MTB>::cols)), FormatsDoNotMatch);
                            PLLL_INTERNAL_STATIC_CHECK((MatrixInfo<MTA>::rows < 0) || (MatrixInfo<MTB>::rows < 0) ||
                                                       (static_cast<size_type>(MatrixInfo<MTA>::rows) == static_cast<size_type>(MatrixInfo<MTB>::rows)), FormatsDoNotMatch);
                            PLLL_INTERNAL_STATIC_CHECK(MatrixInfo<MTA>::is_matrix && MatrixInfo<MTA>::is_math_object, RequiresMathMatrix);
                            PLLL_INTERNAL_STATIC_CHECK(MatrixInfo<MTB>::is_matrix && MatrixInfo<MTB>::is_math_object, RequiresMathMatrix);
                        }
                        
                        template<typename ResultType>
                        static inline void get_coeff(ResultType & result, size_type i, size_type j, const std::pair<MTA, MTB> & AB)
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Op::op(result, AB.first(i, j), AB.second(i, j))))
                        {
                            PLLL_DEBUG_OUTPUT_MESSAGE("expressions::componentwise_operation::operation::get_coeff(" << getAddress(*this) << "; " << getAddress(result) << " " << i << " " << j << " " << getAddress(AB) << ")");
                            assert(i < rows(AB));
                            assert(j < cols(AB));
                            Op::op(result, AB.first(i, j), AB.second(i, j));
                        }
                        
                        static inline bool get_coeff_alwayszero(const std::pair<MTA, MTB> & AB) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                        {
                            PLLL_DEBUG_OUTPUT_MESSAGE("expressions::componentwise_operation::operation::get_coeff_alwayszero(" << getAddress(*this) << "; " << getAddress(AB) << ")");
                            return false;
                        }
                        
                        class GetCoeffSteps_Type
                        {
                        protected:
                            typename Op::ReturnType d_value;
                            
                        public:
                            inline GetCoeffSteps_Type(const typename Op::ReturnType & value)
                                PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(typename Op::ReturnType(value)))
                                : d_value(value)
                            {
                            }
                            
#if __cplusplus >= 201103L
                            inline GetCoeffSteps_Type(typename Op::ReturnType && value)
                                PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(typename Op::ReturnType(std::move(value))))
                                : d_value(std::move(value))
                            {
                            }
#endif
                            
                            typedef const typename Op::ReturnType & Step1_Type;
                            
                            inline Step1_Type step1() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                            {
                                return d_value;
                            }
                            
                            template<typename ResultType>
                            inline void step2(ResultType & result) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                            {
                            }
                        };
                        
                        static inline GetCoeffSteps_Type get_coeff_steps(size_type i, size_type j, const std::pair<MTA, MTB> & AB)
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(GetCoeffSteps_Type(Op::op(AB.first(i, j), AB.second(i, j)))))
                        {
                            PLLL_DEBUG_OUTPUT_MESSAGE("expressions::componentwise_operation::operation::get_coeff_steps(" << getAddress(*this) << "; " << i << " " << j << " " << getAddress(AB) << ")");
                            assert(i < rows(AB));
                            assert(j < cols(AB));
                            return GetCoeffSteps_Type(Op::op(AB.first(i, j), AB.second(i, j)));
                        }
                        
                        inline CoeffType_Get operator () (size_type i, size_type j, const std::pair<MTA, MTB> & AB) const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Op::op(AB.first(i, j), AB.second(i, j))))
                        {
                            PLLL_DEBUG_OUTPUT_MESSAGE("expressions::componentwise_operation::operation::operator () (" << getAddress(*this) << "; " << i << " " << j << " " << getAddress(AB) << ")");
                            assert(i < rows(AB));
                            assert(j < cols(AB));
                            return Op::op(AB.first(i, j), AB.second(i, j));
                        }
                        
                        static inline size_type rows(const std::pair<MTA, MTB> & AB) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                        {
                            PLLL_DEBUG_OUTPUT_MESSAGE("expressions::componentwise_operation::operation::rows(" << getAddress(*this) << "; " << getAddress(AB) << ")");
                            return AB.first.rows();
                        }
                        
                        static inline size_type cols(const std::pair<MTA, MTB> & AB) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                        {
                            PLLL_DEBUG_OUTPUT_MESSAGE("expressions::componentwise_operation::operation::cols(" << getAddress(*this) << "; " << getAddress(AB) << ")");
                            return AB.first.cols();
                        }
                        
                        static inline size_type size(const std::pair<MTA, MTB> & AB) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                        {
                            PLLL_DEBUG_OUTPUT_MESSAGE("expressions::componentwise_operation::operation::size(" << getAddress(*this) << "; " << getAddress(AB) << ")");
                            return AB.first.rows() * AB.first.cols();
                        }
                        
                        template<typename MatrixType>
                        static inline bool involves_this_matrix(const MatrixType & matrix, const std::pair<MTA, MTB> & AB) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                        {
                            PLLL_DEBUG_OUTPUT_MESSAGE("expressions::componentwise_operation::operation::involves_this_matrix(" << getAddress(*this) << "; " << getAddress(matrix) << " " << getAddress(AB) << ")");
                            return AB.first.involves_this_matrix(matrix) || AB.second.involves_this_matrix(matrix);
                        }
                        
                        template<typename MatrixType>
                        static inline bool test_involvement(const MatrixType & matrix, const std::pair<MTA, MTB> & AB) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                        {
                            PLLL_DEBUG_OUTPUT_MESSAGE("expressions::componentwise_operation::operation::test_involvement(" << getAddress(*this) << "; " << getAddress(matrix) << " " << getAddress(AB) << ")");
                            return AB.first.test_involvement(matrix) || AB.second.test_involvement(matrix);
                        }
                        
                        template<typename EnumA, typename EnumB>
                        class CWEnumerator
                        {
                        private:
                            EnumA d_enum_A;
                            EnumB d_enum_B;
                            
                        public:
                            typedef CoeffType_Get Type;
                            typedef operation<PairMatrices>::GetCoeffSteps_Type GetCoeffSteps_Type;
                            
                            inline CWEnumerator(const EnumA & eA, const EnumB & eB)
                                PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(EnumA(eA)) && noexcept(EnumB(eB)))
                                : d_enum_A(eA), d_enum_B(eB)
                            {
                            }
                            
#if __cplusplus >= 201103L
                            inline CWEnumerator(EnumA && eA, EnumB && eB)
                                PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(EnumA(std::move(eA))) && noexcept(EnumB(std::move(eB))))
                                : d_enum_A(std::move(eA)), d_enum_B(std::move(eB))
                            {
                            }
#endif
                            
                            inline Type current() const
                                PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Op::op(d_enum_A.current(), d_enum_B.current())))
                            {
                                return Op::op(d_enum_A.current(), d_enum_B.current());
                            }
                            
                            template<typename Result>
                            inline void get_current(Result & result) const
                                PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Op::op(result, d_enum_A.current(), d_enum_B.current())))
                            {
                                Op::op(result, d_enum_A.current(), d_enum_B.current());
                            }
                            
                            inline GetCoeffSteps_Type get_current_steps() const
                                PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(GetCoeffSteps_Type(Op::op(d_enum_A.current(), d_enum_B.current()))))
                            {
                                return GetCoeffSteps_Type(Op::op(d_enum_A.current(), d_enum_B.current()));
                            }
                            
                            inline bool has_current() const
                                PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_enum_A.has_current()))
                            {
                                return d_enum_A.has_current();
                            }
                            
                            inline void next()
                                PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_enum_A.next()) && noexcept(d_enum_B.next()))
                            {
                                d_enum_A.next();
                                d_enum_B.next();
                            }
                        };
                        
                        typedef void Enumerator;
                        typedef CWEnumerator<typename MTA::ConstEnumerator, typename MTB::ConstEnumerator> ConstEnumerator;
                        typedef CWEnumerator<typename MTA::ConstEnumerator, typename MTB::ConstEnumerator> DefaultEnumerator;
                        typedef void RowEnumerator;
                        typedef CWEnumerator<typename MTA::ConstRowEnumerator, typename MTB::ConstRowEnumerator> ConstRowEnumerator;
                        typedef CWEnumerator<typename MTA::ConstRowEnumerator, typename MTB::ConstRowEnumerator> DefaultRowEnumerator;
                        typedef void ColEnumerator;
                        typedef CWEnumerator<typename MTA::ConstColEnumerator, typename MTB::ConstColEnumerator> ConstColEnumerator;
                        typedef CWEnumerator<typename MTA::ConstColEnumerator, typename MTB::ConstColEnumerator> DefaultColEnumerator;
                        
                        static inline DefaultEnumerator enumerate(const std::pair<MTA, MTB> & AB)
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(DefaultEnumerator(AB.first.enumerate(), AB.second.enumerate())))
                        {
                            return DefaultEnumerator(AB.first.enumerate(), AB.second.enumerate());
                        }
                        
                        static inline DefaultRowEnumerator enumerate_row(size_type row, const std::pair<MTA, MTB> & AB)
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(DefaultRowEnumerator(AB.first.enumerate_row(row), AB.second.enumerate_row(row))))
                        {
                            return DefaultRowEnumerator(AB.first.enumerate_row(row), AB.second.enumerate_row(row));
                        }
                        
                        static inline DefaultColEnumerator enumerate_col(size_type col, const std::pair<MTA, MTB> & AB)
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(DefaultColEnumerator(AB.first.enumerate_col(col), AB.second.enumerate_col(col))))
                        {
                            return DefaultColEnumerator(AB.first.enumerate_col(col), AB.second.enumerate_col(col));
                        }
                        
                        template<typename Enum_meta_A, typename Enum_meta_B, typename Enum_A, typename Enum_B>
                        class CWEnumerator_meta
                        {
                        private:
                            Enum_meta_A d_enum_A;
                            Enum_meta_B d_enum_B;
                            
                        public:
                            typedef CWEnumerator<Enum_A, Enum_B> Type;
                            
                            inline CWEnumerator_meta(const Enum_meta_A & eA, const Enum_meta_B & eB)
                                PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Enum_meta_A(eA)) && noexcept(Enum_mbeta_B(eB)))
                                : d_enum_A(eA), d_enum_B(eB)
                            {
                            }
                            
#if __cplusplus >= 201103L
                            inline CWEnumerator_meta(Enum_meta_A && eA, Enum_meta_B && eB)
                                PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Enum_meta_A(std::move(eA))) && noexcept(Enum_mbeta_B(std::move(eB))))
                                : d_enum_A(std::move(eA)), d_enum_B(std::move(eB))
                            {
                            }
#endif
                            
                            inline Type current() const
                                PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Type(d_enum_A.current(), d_enum_B.current())))
                            {
                                return Type(d_enum_A.current(), d_enum_B.current());
                            }
                            
                            inline bool has_current() const
                                PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_enum_A.has_current()))
                            {
                                return d_enum_A.has_current();
                            }
                            
                            inline void next()
                                PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_enum_A.next()) && noexcept(d_enum_B.next()))
                            {
                                d_enum_A.next();
                                d_enum_B.next();
                            }
                        };
                        
                        typedef void RowsEnumerator;
                        typedef void ColsEnumerator;
                        typedef CWEnumerator_meta<typename MTA::ConstRowsEnumerator, typename MTB::ConstRowsEnumerator,
                                                  typename MTA::ConstRowEnumerator,  typename MTB::ConstRowEnumerator> ConstRowsEnumerator;
                        typedef CWEnumerator_meta<typename MTA::ConstColsEnumerator, typename MTB::ConstColsEnumerator,
                                                  typename MTA::ConstColEnumerator,  typename MTB::ConstColEnumerator> ConstColsEnumerator;
                        typedef ConstRowsEnumerator DefaultRowsEnumerator;
                        typedef ConstColsEnumerator DefaultColsEnumerator;
                        
                        inline DefaultRowsEnumerator enumerate_rows(const std::pair<MTA, MTB> & AB) const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(DefaultRowsEnumerator(AB.first.enumerate_rows(), AB.second.enumerate_rows())))
                        {
                            return DefaultRowsEnumerator(AB.first.enumerate_rows(), AB.second.enumerate_rows());
                        }
                        
                        inline DefaultColsEnumerator enumerate_cols(const std::pair<MTA, MTB> & AB) const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(DefaultColsEnumerator(AB.first.enumerate_cols(), AB.second.enumerate_cols())))
                        {
                            return DefaultColsEnumerator(AB.first.enumerate_cols(), AB.second.enumerate_cols());
                        }
                    };
                };
                
                template<typename MTIn>
                class matrix_negation
                {
                public:
                    typedef typename arithmetic::unary_operation<typename MTIn::CoeffType_Get, arithmetic::op::negation>::ResultType CoeffType;
                    typedef typename arithmetic::unary_operation<typename MTIn::CoeffType_Get, arithmetic::op::negation>::IntermediateType CoeffType_Get;
                    typedef base_matrix<CoeffType,
                                        MatrixInfo<MTIn>::rows,
                                        MatrixInfo<MTIn>::cols,
                                        typename MatrixInfo<MTIn>::StorageTraits,
                                        true> ValueType;
                    typedef void data_pointer_type;
                    typedef void data_ref_type;
                    
                    enum { use_temporary_on_evaluate = MatrixInfo<MTIn>::use_temporary_on_evaluate,
                           coeffs_are_simple_expressions = true,
                           is_const = true,
                           can_assign_to = false,
                           can_move_from = false,
                           can_resize_rows = false,
                           can_resize_cols = false,
                           has_direct_access = false };
                    
                    matrix_negation() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_INTERNAL_STATIC_CHECK(MatrixInfo<MTIn>::is_matrix && MatrixInfo<MTIn>::is_math_object, RequiresMathMatrix);
                    }
                
                    template<typename ResultType>
                    static inline void get_coeff(ResultType & result, size_type i, size_type j, const MTIn & A)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(result = -A(i, j)))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::matrix_negation::get_coeff(" << getAddress(*this) << "; " << getAddress(result) << " " << i << " " << j << " " << getAddress(A) << ") const");
                        assert(i < rows(A));
                        assert(j < cols(A));
                        result = -A(i, j);
                    }
                
                    static inline bool get_coeff_alwayszero(const MTIn & A) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::matrix_negation::get_coeff_alwayszero(" << getAddress(*this) << "; " << getAddress(A) << ")");
                        return false;
                    }
                    
                    class GetCoeffSteps_Type
                    {
                    protected:
                        typedef typename arithmetic::unary_operation<typename MTIn::CoeffType_Get, arithmetic::op::negation>::IntermediateType Type;
                        Type d_value;
                        
                    public:
                        inline GetCoeffSteps_Type(const Type & value)
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Type(value)))
                            : d_value(value)
                        {
                        }
                        
#if __cplusplus >= 201103L
                        inline GetCoeffSteps_Type(Type && value)
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Type(std::move(value))))
                            : d_value(std::move(value))
                        {
                        }
#endif
                        
                        typedef const Type & Step1_Type;
                        
                        inline Step1_Type step1() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                        {
                            return d_value;
                        }
                        
                        template<typename ResultType>
                        inline void step2(ResultType & result) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                        {
                        }
                    };
                    
                    static inline GetCoeffSteps_Type get_coeff_steps(size_type i, size_type j, const MTIn & A)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(GetCoeffSteps_Type(-A(i, j))))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::matrix_negation::get_coeff_steps(" << getAddress(*this) << "; " << i << " " << j << " " << getAddress(A) << ")");
                        assert(i < rows(A));
                        assert(j < cols(A));
                        return GetCoeffSteps_Type(-A(i, j));
                    }
                
                    inline CoeffType_Get operator () (size_type i, size_type j, const MTIn & A) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(-A(i, j)))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::matrix_negation::operator () (" << getAddress(*this) << "; " << i << " " << j << " " << getAddress(A) << ")");
                        assert(i < rows(A));
                        assert(j < cols(A));
                        return -A(i, j);
                    }
                
                    static inline size_type rows(const MTIn & A) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::matrix_negation::rows(" << getAddress(*this) << "; " << getAddress(A) << ")");
                        return A.rows();
                    }
                
                    static inline size_type cols(const MTIn & A) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::matrix_negation::cols(" << getAddress(*this) << "; " << getAddress(A) << ")");
                        return A.cols();
                    }
                    
                    static inline size_type size(const MTIn & A) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::matrix_negation::size(" << getAddress(*this) << "; " << getAddress(A) << ")");
                        return A.rows() * A.cols();
                    }
                    
                    template<typename MatrixType>
                    static inline bool involves_this_matrix(const MatrixType & matrix, const MTIn & A) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::matrix_negation::involves_this_matrix(" << getAddress(*this) << "; " << getAddress(matrix) << " " << getAddress(A) << ")");
                        return A.involves_this_matrix(matrix);
                    }
                    
                    template<typename MatrixType>
                    static inline bool test_involvement(const MatrixType & matrix, const MTIn & A) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::matrix_negation::test_involvement(" << getAddress(*this) << "; " << getAddress(matrix) << " " << getAddress(A) << ")");
                        return A.test_involvement(matrix);
                    }
                    
                    template<typename Enum>
                    class NegEnumerator
                    {
                    private:
                        Enum d_enum;
                        
                    public:
                        typedef CoeffType_Get Type;
                        typedef matrix_negation<MTIn>::GetCoeffSteps_Type GetCoeffSteps_Type;
                        
                        inline NegEnumerator(const Enum & e)
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Enum(e)))
                            : d_enum(e)
                        {
                        }
                    
#if __cplusplus >= 201103L
                        NegEnumerator(Enum && e)
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Enum(std::move(e))))
                            : d_enum(std::move(e))
                        {
                        }
#endif
                        
                        inline Type current() const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Type(-d_enum.current())))
                        {
                            return -d_enum.current();
                        }
                        
                        template<typename Result>
                        inline void get_current(Result & result) const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(result = -d_enum.current()))
                        {
                            result = -d_enum.current();
                        }
                        
                        inline GetCoeffSteps_Type get_current_steps() const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(GetCoeffSteps_Type(-d_enum.current())))
                        {
                            return GetCoeffSteps_Type(-d_enum.current());
                        }
                        
                        inline bool has_current() const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_enum.has_current()))
                        {
                            return d_enum.has_current();
                        }
                    
                        inline void next()
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_enum.next()))
                        {
                            d_enum.next();
                        }
                    };
                    
                    typedef void Enumerator;
                    typedef NegEnumerator<typename MTIn::ConstEnumerator> ConstEnumerator;
                    typedef NegEnumerator<typename MTIn::ConstEnumerator> DefaultEnumerator;
                    typedef void RowEnumerator;
                    typedef NegEnumerator<typename MTIn::ConstRowEnumerator> ConstRowEnumerator;
                    typedef NegEnumerator<typename MTIn::ConstRowEnumerator> DefaultRowEnumerator;
                    typedef void ColEnumerator;
                    typedef NegEnumerator<typename MTIn::ConstColEnumerator> ConstColEnumerator;
                    typedef NegEnumerator<typename MTIn::ConstColEnumerator> DefaultColEnumerator;
                    
                    static inline DefaultEnumerator enumerate(const MTIn & A)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(DefaultEnumerator(A.enumerate())))
                    {
                        return DefaultEnumerator(A.enumerate());
                    }
                    
                    static inline DefaultRowEnumerator enumerate_row(size_type row, const MTIn & A)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(DefaultRowEnumerator(A.enumerate_row(row))))
                    {
                        return DefaultRowEnumerator(A.enumerate_row(row));
                    }
                    
                    static inline DefaultColEnumerator enumerate_col(size_type col, const MTIn & A)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(DefaultColEnumerator(A.enumerate_col(col))))
                    {
                        return DefaultColEnumerator(A.enumerate_col(col));
                    }
                    
                    template<typename Enum_meta, typename Enum>
                    class NegEnumerator_meta
                    {
                    private:
                        Enum_meta d_enum;
                        
                    public:
                        typedef NegEnumerator<Enum> Type;
                        
                        inline NegEnumerator_meta(const Enum_meta & e)
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Enum_meta(e)))
                            : d_enum(e)
                        {
                        }
                        
#if __cplusplus >= 201103L
                        inline NegEnumerator_meta(Enum_meta && e)
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Enum_meta(std::move(e))))
                            : d_enum(std::move(e))
                        {
                        }
#endif
                        
                        inline Type current() const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Type(d_enum.current())))
                        {
                            return Type(d_enum.current());
                        }
                        
                        inline bool has_current() const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_enum.has_current()))
                        {
                            return d_enum.has_current();
                        }
                        
                        inline void next()
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_enum.next()))
                        {
                            d_enum.next();
                        }
                    };
                    
                    typedef void RowsEnumerator;
                    typedef void ColsEnumerator;
                    typedef NegEnumerator_meta<typename MTIn::ConstRowsEnumerator, typename MTIn::ConstRowEnumerator> ConstRowsEnumerator;
                    typedef NegEnumerator_meta<typename MTIn::ConstColsEnumerator, typename MTIn::ConstColEnumerator> ConstColsEnumerator;
                    typedef ConstRowsEnumerator DefaultRowsEnumerator;
                    typedef ConstColsEnumerator DefaultColsEnumerator;
                    
                    inline DefaultRowsEnumerator enumerate_rows(const MTIn & A) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(DefaultRowsEnumerator(A.enumerate_rows())))
                    {
                        return DefaultRowsEnumerator(A.enumerate_rows());
                    }
                    
                    inline DefaultColsEnumerator enumerate_cols(const MTIn & A) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(DefaultColsEnumerator(A.enumerate_cols())))
                    {
                        return DefaultColsEnumerator(A.enumerate_cols());
                    }
                };
                
                template<int Rows, int Cols>
                template<typename MT>
                class sub<Rows, Cols>::operation_generic : private row_count_storage<Rows>, private col_count_storage<Cols>
                {
                private:
                    size_type d_r_ofs, d_c_ofs;
                    
                public:
                    typedef typename MT::CoeffType CoeffType;
                    typedef typename MT::CoeffType_Get CoeffType_Get;
                    typedef typename MT::GetCoeffSteps_Type GetCoeffSteps_Type;
                    typedef base_matrix<typename MatrixInfo<MT>::Type,
                                        Rows, Cols,
                                        typename MatrixInfo<MT>::StorageTraits,
                                        MatrixInfo<MT>::is_math_object> ValueType;
                    typedef void data_pointer_type;
                    typedef void data_ref_type;
                    
                    enum { use_temporary_on_evaluate = MatrixInfo<MT>::use_temporary_on_evaluate,
                           coeffs_are_simple_expressions = true,
                           is_const = MatrixInfo<MT>::is_const,
                           can_assign_to = MatrixInfo<MT>::can_assign_to,
                           can_move_from = MatrixInfo<MT>::can_move_from,
                           can_resize_rows = false,
                           can_resize_cols = false,
                           has_direct_access = false };
                    
                    template<typename MT_>
                    operation_generic(size_type rows, size_type cols, size_type r_ofs, size_type c_ofs, const MT_ & data) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                        : row_count_storage<Rows>(rows), col_count_storage<Cols>(cols), d_r_ofs(r_ofs), d_c_ofs(c_ofs)
                          // data is given only for verification reasons!
                    {
                        PLLL_INTERNAL_STATIC_CHECK((Rows < 0) || (MatrixInfo<MT>::rows < 0) ||
                                                   (static_cast<size_type>(Rows) == static_cast<size_type>(MatrixInfo<MT>::rows)), TooFewRows);
                        PLLL_INTERNAL_STATIC_CHECK((Cols < 0) || (MatrixInfo<MT>::cols < 0) ||
                                                   (static_cast<size_type>(Cols) == static_cast<size_type>(MatrixInfo<MT>::cols)), TooFewCols);
                        PLLL_INTERNAL_STATIC_CHECK(MatrixInfo<MT>::is_matrix, RequiresMatrix);
                        if (Rows > 0) assert(row_count_storage<Rows>::rows() == static_cast<size_type>(Rows));
                        if (Cols > 0) assert(col_count_storage<Cols>::cols() == static_cast<size_type>(Cols));
                        assert(row_count_storage<Rows>::rows() + d_r_ofs <= data.rows());
                        assert(col_count_storage<Cols>::cols() + d_c_ofs <= data.cols());
                    }
                    
                    template<typename ResultType>
                    inline void get_coeff(ResultType & result, size_type i, size_type j, const MT & A)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(result = A(d_r_ofs + i, d_c_ofs + j)))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::sub::operation_generic::get_coeff(" << getAddress(*this) << "; " << getAddress(result) << " " << i << " " << j << " " << getAddress(A) << ") const");
                        assert(i < rows(A));
                        assert(j < cols(A));
                        result = A(d_r_ofs + i, d_c_ofs + j);
                    }
                    
                    inline bool get_coeff_alwayszero(const MT & A) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(A.get_coeff_alwayszero()))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::sub::operation_generic::get_coeff_alwayszero(" << getAddress(*this) << "; " << getAddress(A) << ")");
                        return A.get_coeff_alwayszero();
                    }
                    
                    inline GetCoeffSteps_Type get_coeff_steps(size_type i, size_type j, const MT & A) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(A.get_coeff_steps(d_r_ofs + i, d_c_ofs + j)))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::sub::operation_generic::get_coeff_steps(" << getAddress(*this) << "; " << i << " " << j << " " << getAddress(A) << ")");
                        assert(i < rows(A));
                        assert(j < cols(A));
                        return A.get_coeff_steps(d_r_ofs + i, d_c_ofs + j);
                    }
                    
                    inline CoeffType_Get operator () (size_type i, size_type j, const MT & A) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(A(d_r_ofs + i, d_c_ofs + j)))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::sub::operation_generic::operator () (" << getAddress(*this) << "; " << i << " " << j << " " << getAddress(A) << ")");
                        assert(i < rows(A));
                        assert(j < cols(A));
                        return A(d_r_ofs + i, d_c_ofs + j);
                    }
                    
                    inline size_type rows(const MT & A) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::sub::operation_generic::rows(" << getAddress(*this) << "; " << getAddress(A) << ")");
                        return row_count_storage<Rows>::rows();
                    }
                    
                    inline size_type cols(const MT & A) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::sub::operation_generic::cols(" << getAddress(*this) << "; " << getAddress(A) << ")");
                        return col_count_storage<Cols>::cols();
                    }
                    
                    inline size_type size(const MT & A) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::sub::operation_generic::size(" << getAddress(*this) << "; " << getAddress(A) << ")");
                        return row_count_storage<Rows>::rows() * col_count_storage<Cols>::cols();
                    }
                    
                    template<typename MatrixType>
                    inline bool involves_this_matrix(const MatrixType & matrix, const MT & A) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::sub::operation_generic::involves_this_matrix(" << getAddress(*this) << "; " << getAddress(matrix) << " " << getAddress(A) << ")");
                        return A.involves_this_matrix(matrix);
                    }
                    
                    template<typename MatrixType>
                    inline bool test_involvement(const MatrixType & matrix, const MT & A) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::sub::operation_generic::test_involvement(" << getAddress(*this) << "; " << getAddress(matrix) << " " << getAddress(A) << ")");
                        return A.test_involvement(matrix);
                    }
                    
                    class Enumerator : private row_count_storage<Rows>, private col_count_storage<Cols>
                    {
                    private:
                        size_type d_r_ofs, d_c_ofs;
                        size_type d_row, d_col;
                        MT d_A;
                        
                    public:
                        typedef CoeffType_Get Type;
                        typedef typename MT::GetCoeffSteps_Type GetCoeffSteps_Type;
                        
                        template<typename SType>
                        inline Enumerator(const SType & s, const MT & A)
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(MT(A)))
                            : row_count_storage<Rows>(s.rows(A)), col_count_storage<Cols>(s.cols(A)),
                              d_r_ofs(s.d_r_ofs), d_c_ofs(s.d_c_ofs), d_row(0), d_col(0), d_A(A)
                        {
                        }
                        
                        inline Type current() const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_A(d_r_ofs + d_row, d_c_ofs + d_col)))
                        {
                            return d_A(d_r_ofs + d_row, d_c_ofs + d_col);
                        }
                        
                        template<typename Result>
                        inline void get_current(Result & result) const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_A.get_coeffs(result, d_r_ofs + d_row, d_c_ofs + d_col)))
                        {
                            d_A.get_coeffs(result, d_r_ofs + d_row, d_c_ofs + d_col);
                        }
                        
                        inline GetCoeffSteps_Type get_current_steps() const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_A.get_coeff_steps(d_r_ofs + d_row, d_c_ofs + d_col)))
                        {
                            return d_A.get_coeff_steps(d_r_ofs + d_row, d_c_ofs + d_col);
                        }
                        
                        inline bool has_current() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                        {
                            return d_row < row_count_storage<Rows>::rows();
                        }
                        
                        inline void next() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                        {
                            if (++d_col >= col_count_storage<Cols>::cols())
                            {
                                d_col = 0;
                                ++d_row;
                            }
                        }
                    };
                    
                    class RowEnumerator
                    {
                    private:
                        MT d_A;
                        size_type d_row, d_col;
                        size_type d_left;
                        
                    public:
                        typedef CoeffType_Get Type;
                        typedef typename MT::GetCoeffSteps_Type GetCoeffSteps_Type;
                        
                        template<typename SType>
                        inline RowEnumerator(const SType & s, size_type row, const MT & A)
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(MT(A)))
                            : d_A(A), d_row(s.d_r_ofs + row), d_col(s.d_c_ofs), d_left(s.cols(A))
                        {
                        }
                        
                        inline Type current() const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_A(d_row, d_col)))
                        {
                            return d_A(d_row, d_col);
                        }
                        
                        template<typename Result>
                        inline void get_current(Result & result) const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_A.get_coeff(result, d_row, d_col)))
                        {
                            d_A.get_coeff(result, d_row, d_col);
                        }
                        
                        inline GetCoeffSteps_Type get_current_steps() const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_A.get_coeff_steps(d_row, d_col)))
                        {
                            return d_A.get_coeff_steps(d_row, d_col);
                        }
                        
                        inline bool has_current() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                        {
                            return d_left > 0;
                        }
                    
                        inline void next() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                        {
                            --d_left;
                            ++d_col;
                        }
                    };
                    
                    class ColEnumerator
                    {
                    private:
                        MT d_A;
                        size_type d_row, d_col;
                        size_type d_left;
                        
                    public:
                        typedef CoeffType_Get Type;
                        typedef typename MT::GetCoeffSteps_Type GetCoeffSteps_Type;
                        
                        template<typename SType>
                        inline ColEnumerator(const SType & s, size_type col, const MT & A)
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(MT(A)))
                            : d_A(A), d_row(s.d_r_ofs), d_col(s.d_c_ofs + col), d_left(s.rows(A))
                        {
                        }
                        
                        inline Type current() const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_A(d_row, d_col)))
                        {
                            return d_A(d_row, d_col);
                        }
                        
                        template<typename Result>
                        inline void get_current(Result & result) const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_A.get_coeff(result, d_row, d_col)))
                        {
                            d_A.get_coeff(result, d_row, d_col);
                        }
                        
                        inline GetCoeffSteps_Type get_current_steps() const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_A.get_coeff_steps(d_row, d_col)))
                        {
                            return d_A.get_coeff_steps(d_row, d_col);
                        }
                        
                        inline bool has_current() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                        {
                            return d_left > 0;
                        }
                    
                        inline void next() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                        {
                            --d_left;
                            ++d_row;
                        }
                    };
                    
                    friend class Enumerator;
                    friend class RowEnumerator;
                    friend class ColEnumerator;
                    
                    typedef Enumerator ConstEnumerator;
                    typedef RowEnumerator ConstRowEnumerator;
                    typedef ColEnumerator ConstColEnumerator;
                    typedef Enumerator DefaultEnumerator;
                    typedef RowEnumerator DefaultRowEnumerator;
                    typedef ColEnumerator DefaultColEnumerator;
                    
                    inline DefaultEnumerator enumerate(const MT & A) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(DefaultEnumerator(helper::make_type_lvalue<const operation_generic<MT> >(), A)))
                    {
                        return DefaultEnumerator(*this, A);
                    }
                    
                    inline DefaultRowEnumerator enumerate_row(size_type row, const MT & A)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(DefaultRowEnumerator(helper::make_type_lvalue<const operation_generic<MT> >(), row, A)))
                    {
                        return DefaultRowEnumerator(*this, row, A);
                    }
                    
                    inline DefaultColEnumerator enumerate_col(size_type col, const MT & A) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(DefaultColEnumerator(helper::make_type_lvalue<const operation_generic<MT> >(), col, A)))
                    {
                        return DefaultColEnumerator(*this, col, A);
                    }
                    
                    class RowsEnumerator
                    {
                    private:
                        typename sub<Rows, Cols>::template operation_generic<MT> * d_data;
                        const MT & d_A;
                        size_type d_row;
                        
                    public:
                        typedef RowEnumerator Type;
                        
                        RowsEnumerator(typename sub<Rows, Cols>::template operation_generic<MT> * data, const MT & A) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                            : d_data(data), d_A(A), d_row(0)
                        {
                        }
                        
                        inline Type current() const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Type(*d_data, d_data->d_r_ofs + d_row, d_A)))
                        {
                            return Type(*d_data, d_data->d_r_ofs + d_row, d_A);
                        }
                        
                        inline bool has_current() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                        {
                            return d_row < d_data->rows(d_A);
                        }
                        
                        inline void next() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                        {
                            ++d_row;
                        }
                    };
                    
                    class ColsEnumerator
                    {
                    private:
                        typename sub<Rows, Cols>::template operation_generic<MT> * d_data;
                        const MT & d_A;
                        size_type d_col;
                        
                    public:
                        typedef ColEnumerator Type;
                        
                        ColsEnumerator(typename sub<Rows, Cols>::template operation_generic<MT> * data, const MT & A) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                            : d_data(data), d_A(A), d_col(0)
                        {
                        }
                        
                        inline Type current() const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Type(*d_data, d_data->d_c_ofs + d_col, d_A)))
                        {
                            return Type(*d_data, d_data->d_c_ofs + d_col, d_A);
                        }
                        
                        inline bool has_current() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                        {
                            return d_col < d_data->cols(d_A);
                        }
                    
                        inline void next()
                            PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                        {
                            ++d_col;
                        }
                    };
                    
                    friend class RowsEnumerator;
                    friend class ColsEnumerator;
                    
                    typedef RowsEnumerator ConstRowsEnumerator;
                    typedef ColsEnumerator ConstColsEnumerator;
                    typedef RowsEnumerator DefaultRowsEnumerator;
                    typedef ColsEnumerator DefaultColsEnumerator;
                    
                    inline DefaultRowsEnumerator enumerate_row(const MT & A) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(DefaultRowsEnumerator(helper::make_type_lvalue<const operation_generic<MT> *>(), A)))
                    {
                        return DefaultRowsEnumerator(this, A);
                    }
                    
                    inline DefaultColsEnumerator enumerate_col(const MT & A) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(DefaultColsEnumerator(helper::make_type_lvalue<const operation_generic<MT> *>(), A)))
                    {
                        return DefaultColsEnumerator(this, A);
                    }
                };
                
                template<int Cols>
                template<typename MT>
                class sub_1d<Cols>::operation_row : private col_count_storage<Cols>
                {
                private:
                    size_type d_r_ofs;
                    
                public:
                    typedef typename MT::CoeffType CoeffType;
                    typedef typename MT::CoeffType_Get CoeffType_Get;
                    typedef typename MT::GetCoeffSteps_Type GetCoeffSteps_Type;
                    typedef base_matrix<typename MatrixInfo<MT>::Type,
                                        1, Cols,
                                        typename MatrixInfo<MT>::StorageTraits,
                                        MatrixInfo<MT>::is_math_object> ValueType;
                    typedef typename MT::data_pointer_type data_pointer_type;
                    typedef typename MT::data_ref_type data_ref_type;
                    
                    enum { use_temporary_on_evaluate = MatrixInfo<MT>::use_temporary_on_evaluate,
                           coeffs_are_simple_expressions = true,
                           is_const = MatrixInfo<MT>::is_const,
                           can_assign_to = MatrixInfo<MT>::can_assign_to,
                           can_move_from = MatrixInfo<MT>::can_move_from,
                           can_resize_rows = false,
                           can_resize_cols = false,
                           has_direct_access = MT::has_direct_access };
                    
                    template<typename MT_>
                    operation_row(size_type cols, size_type r_ofs, const MT_ & data)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                        : col_count_storage<Cols>(cols), d_r_ofs(r_ofs)
                          // data is given only for verification reasons!
                    {
                        PLLL_INTERNAL_STATIC_CHECK((Cols < 0) || (MatrixInfo<MT>::cols < 0) ||
                                                   (static_cast<size_type>(Cols) == static_cast<size_type>(MatrixInfo<MT>::cols)), TooFewCols);
                        PLLL_INTERNAL_STATIC_CHECK(MatrixInfo<MT>::is_matrix, RequiresMatrix);
                        if (Cols > 0) assert(col_count_storage<Cols>::cols() == static_cast<size_type>(Cols));
                        assert(d_r_ofs < data.rows());
                    }
                    
                    template<typename ResultType>
                    inline void get_coeff(ResultType & result, size_type i, size_type j, const MT & A) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(result = A(d_r_ofs, j)))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::sub_1d::operation_row::get_coeff(" << getAddress(*this) << "; " << getAddress(result) << " " << i << " " << j << " " << getAddress(A) << ") const");
                        assert(i == 0);
                        assert(j < cols(A));
                        result = A(d_r_ofs, j);
                    }
                    
                    inline bool get_coeff_alwayszero(const MT & A) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(A.get_coeff_alwayszero()))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::sub_1d::operation_row::get_coeff_alwayszero(" << getAddress(*this) << "; " << getAddress(A) << ")");
                        return A.get_coeff_alwayszero();
                    }
                    
                    inline GetCoeffSteps_Type get_coeff_steps(size_type i, size_type j, const MT & A) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(A.get_coeff_steps(d_r_ofs, j)))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::sub_1d::operation_row::get_coeff_steps(" << getAddress(*this) << "; " << i << " " << j << " " << getAddress(A) << ")");
                        assert(i == 0);
                        assert(j < cols(A));
                        return A.get_coeff_steps(d_r_ofs, j);
                    }
                    
                    inline CoeffType_Get operator () (size_type i, size_type j, const MT & A) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(A(d_r_ofs, j)))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::sub_1d::operation_row::operator () (" << getAddress(*this) << "; " << i << " " << j << " " << getAddress(A) << ")");
                        assert(i == 0);
                        assert(j < cols(A));
                        return A(d_r_ofs, j);
                    }
                    
                    inline size_type rows(const MT & A) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::sub_1d::operation_row::rows(" << getAddress(*this) << "; " << getAddress(A) << ")");
                        return 1;
                    }
                    
                    inline size_type cols(const MT & A) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::sub_1d::operation_row::cols(" << getAddress(*this) << "; " << getAddress(A) << ")");
                        return col_count_storage<Cols>::cols();
                    }
                    
                    inline size_type size(const MT & A) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::sub_1d::operation_row::size(" << getAddress(*this) << "; " << getAddress(A) << ")");
                        return col_count_storage<Cols>::cols();
                    }
                    
                    inline data_pointer_type data(const MT & A) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(A.data() + helper::make_type_lvalue<const operation_row>().d_r_ofs *
                                                                           helper::make_type_lvalue<const operation_row>().col_count_storage<Cols>::cols()))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::sub_1d::operation_row::data(" << getAddress(*this) << "; " << getAddress(A) << ") const");
                        return A.data() + d_r_ofs * col_count_storage<Cols>::cols();
                    }
                    
                    inline data_ref_type data(size_type i, const MT & A) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(A.data(i + helper::make_type_lvalue<const operation_row>().d_r_ofs *
                                                                                  helper::make_type_lvalue<const operation_row>().col_count_storage<Cols>::cols())))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::sub_1d::operation_row::data(" << getAddress(*this) << "; " << i << " " << getAddress(A) << ") const");
                        return A.data(i + d_r_ofs * col_count_storage<Cols>::cols());
                    }
                    
                    template<typename MatrixType>
                    inline bool involves_this_matrix(const MatrixType & matrix, const MT & A) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::sub_1d::operation_row::involves_this_matrix(" << getAddress(*this) << "; " << getAddress(matrix) << " " << getAddress(A) << ")");
                        return A.involves_this_matrix(matrix);
                    }
                
                    template<typename MatrixType>
                    inline bool test_involvement(const MatrixType & matrix, const MT & A) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::sub_1d::operation_row::test_involvement(" << getAddress(*this) << "; " << getAddress(matrix) << " " << getAddress(A) << ")");
                        return A.test_involvement(matrix);
                    }
                    
                    typedef typename MT::RowEnumerator Enumerator;
                    typedef typename MT::ConstRowEnumerator ConstEnumerator;
                    typedef typename MT::DefaultRowEnumerator DefaultEnumerator;
                    typedef typename MT::RowEnumerator RowEnumerator;
                    typedef typename MT::ConstRowEnumerator ConstRowEnumerator;
                    typedef typename MT::DefaultRowEnumerator DefaultRowEnumerator;
                    
                    class ColEnumerator
                    {
                    private:
                        MT d_A;
                        size_type d_row, d_col;
                        bool d_first;
                        
                    public:
                        typedef CoeffType_Get Type;
                        typedef typename operation_row<MT>::GetCoeffSteps_Type GetCoeffSteps_Type;
                        
                        template<typename SType>
                        inline ColEnumerator(const SType & s, size_type col, const MT & A) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                            : d_A(A), d_row(s.d_r_ofs), d_col(col), d_first(true)
                        {
                        }
                        
                        inline Type current() const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_A(d_row, d_col)))
                        {
                            return d_A(d_row, d_col);
                        }
                        
                        template<typename Result>
                        inline void get_current(Result & result) const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_A.get_coeff(result, d_row, d_col)))
                        {
                            d_A.get_coeff(result, d_row, d_col);
                        }
                        
                        inline GetCoeffSteps_Type get_current_steps() const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_A.get_coeff_steps(d_row, d_col)))
                        {
                            return d_A.get_coeff_steps(d_row, d_col);
                        }
                        
                        inline bool has_current() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                        {
                            return d_first;
                        }
                        
                        inline void next() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                        {
                            d_first = false;
                        }
                    };
                    
                    friend class ColEnumerator;
                    
                    typedef ColEnumerator ConstColEnumerator;
                    typedef ColEnumerator DefaultColEnumerator;
                    
                    inline DefaultEnumerator enumerate(const MT & A) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(A.enumerate_row(d_r_ofs)))
                    {
                        return A.enumerate_row(d_r_ofs);
                    }
                    
                    inline DefaultRowEnumerator enumerate_row(size_type row, const MT & A) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(A.enumerate_row(d_r_ofs)))
                    {
                        assert(row == 0);
                        return A.enumerate_row(d_r_ofs);
                    }
                    
                    inline DefaultColEnumerator enumerate_col(size_type col, const MT & A) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(DefaultColEnumerator(helper::make_type_lvalue<operation_row<MT> >(), col, A)))
                    {
                        return DefaultColEnumerator(*this, col, A);
                    }
                    
                    class RowsEnumerator
                    {
                    private:
                        DefaultRowEnumerator d_row;
                        bool d_first;
                        
                    public:
                        typedef DefaultRowEnumerator Type;
                        
                        RowsEnumerator(const DefaultRowEnumerator & e)
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(DefaultRowEnumerator(e)))
                            : d_row(e), d_first(true)
                        {
                        }
                        
                        inline Type current() const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Type(d_row)))
                        {
                            return d_row;
                        }
                        
                        inline bool has_current() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                        {
                            return d_first;
                        }
                        
                        inline void next() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                        {
                            d_first = false;
                        }
                    };
                    
                    class ColsEnumerator
                    {
                    private:
                        DefaultRowEnumerator d_row;
                        
                    public:
                        class CustomColEnumerator
                        {
                        private:
                            const DefaultRowEnumerator * d_row;
                            bool d_first;
                            
                        public:
                            typedef CoeffType_Get Type;
                            typedef typename DefaultRowEnumerator::GetCoeffSteps_Type GetCoeffSteps_Type;
                            
                            CustomColEnumerator(const DefaultRowEnumerator & row) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                                : d_row(&row), d_first(true)
                            {
                            }
                            
                            inline Type current() const
                                PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_row->current()))
                            {
                                return d_row->current();
                            }
                            
                            template<typename Result>
                            inline void get_current(Result & result) const
                                PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_row->get_current(result)))
                            {
                                d_row->get_current(result);
                            }
                            
                            inline GetCoeffSteps_Type get_current_steps() const
                                PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_row->get_current_steps()))
                            {
                                return d_row->get_current_steps();
                            }
                            
                            inline bool has_current() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                            {
                                return d_first;
                            }
                        
                            inline void next() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                            {
                                d_first = false;
                            }
                        };
                        
                        typedef CustomColEnumerator Type;
                        
                        ColsEnumerator(const RowEnumerator & row)
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(RowEnumerator(row)))
                            : d_row(row)
                        {
                        }
                        
                        inline Type current() const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Type(d_row)))
                        {
                            return Type(d_row);
                        }
                        
                        inline bool has_current() const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_row.has_current()))
                        {
                            return d_row.has_current();
                        }
                        
                        inline void next()
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_row.next()))
                        {
                            d_row.next();
                        }
                    };
                    
                    friend class RowsEnumerator;
                    friend class ColsEnumerator;
                    
                    typedef RowsEnumerator ConstRowsEnumerator;
                    typedef ColsEnumerator ConstColsEnumerator;
                    typedef RowsEnumerator DefaultRowsEnumerator;
                    typedef ColsEnumerator DefaultColsEnumerator;
                    
                    inline DefaultRowsEnumerator enumerate_rows(const MT & A) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(DefaultRowsEnumerator(helper::make_type_lvalue<const operation_row<MT> >().enumerate_row(0, A))))
                    {
                        return DefaultRowsEnumerator(enumerate_row(0, A));
                    }
                    
                    inline DefaultColsEnumerator enumerate_cols(const MT & A) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(DefaultColsEnumerator(helper::make_type_lvalue<const operation_row<MT> >().enumerate_row(0, A))))
                    {
                        return DefaultColsEnumerator(enumerate_row(0, A));
                    }
                };
                
                template<int Rows>
                template<typename MT>
                class sub_1d<Rows>::operation_col : private row_count_storage<Rows>
                {
                private:
                    size_type d_c_ofs;
                    
                public:
                    typedef typename MT::CoeffType CoeffType;
                    typedef typename MT::CoeffType_Get CoeffType_Get;
                    typedef typename MT::GetCoeffSteps_Type GetCoeffSteps_Type;
                    typedef base_matrix<typename MatrixInfo<MT>::Type,
                                        Rows, 1,
                                        typename MatrixInfo<MT>::StorageTraits,
                                        MatrixInfo<MT>::is_math_object> ValueType;
                    typedef void data_pointer_type;
                    typedef void data_ref_type;
                    
                    enum { use_temporary_on_evaluate = MatrixInfo<MT>::use_temporary_on_evaluate,
                           coeffs_are_simple_expressions = true,
                           is_const = MatrixInfo<MT>::is_const,
                           can_assign_to = MatrixInfo<MT>::can_assign_to,
                           can_move_from = MatrixInfo<MT>::can_move_from,
                           can_resize_rows = false,
                           can_resize_cols = false,
                           has_direct_access = false };
                    
                    template<typename MT_>
                    operation_col(size_type rows, size_type c_ofs, const MT_ & data) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                        : row_count_storage<Rows>(rows), d_c_ofs(c_ofs)
                          // data is given only for verification reasons!
                    {
                        PLLL_INTERNAL_STATIC_CHECK((Rows < 0) || (MatrixInfo<MT>::rows < 0) ||
                                                   (static_cast<size_type>(Rows) == static_cast<size_type>(MatrixInfo<MT>::rows)), TooFewRows);
                        PLLL_INTERNAL_STATIC_CHECK(MatrixInfo<MT>::is_matrix, RequiresMatrix);
                        if (Rows > 0) assert(row_count_storage<Rows>::rows() == static_cast<size_type>(Rows));
                        assert(d_c_ofs < data.cols());
                    }
                    
                    template<typename ResultType>
                    inline void get_coeff(ResultType & result, size_type i, size_type j, const MT & A) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(result = A(i, d_c_ofs)))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::sub_1d::operation_col::get_coeff(" << getAddress(*this) << "; " << getAddress(result) << " " << i << " " << j << " " << getAddress(A) << ") const");
                        assert(i == 0);
                        assert(j < cols(A));
                        result = A(i, d_c_ofs);
                    }
                    
                    inline bool get_coeff_alwayszero(const MT & A) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(A.get_coeff_alwayszero()))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::sub_1d::operation_col::get_coeff_alwayszero(" << getAddress(*this) << "; " << getAddress(A) << ")");
                        return A.get_coeff_alwayszero();
                    }
                    
                    inline GetCoeffSteps_Type get_coeff_steps(size_type i, size_type j, const MT & A) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(A.get_coeff_steps(i, d_c_ofs)))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::sub_1d::operation_col::get_coeff_steps(" << getAddress(*this) << "; " << i << " " << j << " " << getAddress(A) << ")");
                        assert(i == 0);
                        assert(j < cols(A));
                        return A.get_coeff_steps(i, d_c_ofs);
                    }
                    
                    inline CoeffType_Get operator () (size_type i, size_type j, const MT & A) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(A(i, d_c_ofs)))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::sub_1d::operation_col::operator () (" << getAddress(*this) << "; " << i << " " << j << " " << getAddress(A) << ")");
                        assert(i == 0);
                        assert(j < cols(A));
                        return A(i, d_c_ofs);
                    }
                    
                    inline size_type rows(const MT & A) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::sub_1d::operation_col::rows(" << getAddress(*this) << "; " << getAddress(A) << ")");
                        return row_count_storage<Rows>::rows();
                    }
                    
                    inline size_type cols(const MT & A) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::sub_1d::operation_col::cols(" << getAddress(*this) << "; " << getAddress(A) << ")");
                        return 1;
                    }
                    
                    inline size_type size(const MT & A) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::sub_1d::operation_col::size(" << getAddress(*this) << "; " << getAddress(A) << ")");
                        return row_count_storage<Rows>::rows();
                    }
                    
                    template<typename MatrixType>
                    inline bool involves_this_matrix(const MatrixType & matrix, const MT & A) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::sub_1d::operation_col::involves_this_matrix(" << getAddress(*this) << "; " << getAddress(matrix) << " " << getAddress(A) << ")");
                        return A.involves_this_matrix(matrix);
                    }
                    
                    template<typename MatrixType>
                    inline bool test_involvement(const MatrixType & matrix, const MT & A) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::sub_1d::operation_col::test_involvement(" << getAddress(*this) << "; " << getAddress(matrix) << " " << getAddress(A) << ")");
                        return A.test_involvement(matrix);
                    }
                    
                    typedef typename MT::ColEnumerator Enumerator;
                    typedef typename MT::ConstColEnumerator ConstEnumerator;
                    typedef typename MT::DefaultColEnumerator DefaultEnumerator;
                    typedef typename MT::ColEnumerator ColEnumerator;
                    typedef typename MT::ConstColEnumerator ConstColEnumerator;
                    typedef typename MT::DefaultColEnumerator DefaultColEnumerator;
                    
                    class RowEnumerator
                    {
                    private:
                        MT d_A;
                        size_type d_row, d_col;
                        bool d_first;
                        
                    public:
                        typedef CoeffType_Get Type;
                        typedef typename operation_col<MT>::GetCoeffSteps_Type GetCoeffSteps_Type;
                        
                        template<typename SType>
                        inline RowEnumerator(const SType & s, size_type row, const MT & A)
                             PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                            : d_A(A), d_row(row), d_col(s.d_c_ofs), d_first(true)
                        {
                        }
                        
                        inline Type current() const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_A(d_row, d_col)))
                        {
                            return d_A(d_row, d_col);
                        }
                        
                        template<typename Result>
                        inline void get_current(Result & result) const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_A.get_coeff(result, d_row, d_col)))
                        {
                            d_A.get_coeff(result, d_row, d_col);
                        }
                        
                        inline GetCoeffSteps_Type get_current_steps() const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_A.get_coeff_steps(d_row, d_col)))
                        {
                            return d_A.get_coeff_steps(d_row, d_col);
                        }
                        
                        inline bool has_current() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                        {
                            return d_first;
                        }
                        
                        inline void next() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                        {
                            d_first = false;
                        }
                    };
                    
                    friend class RowEnumerator;
                    
                    typedef RowEnumerator ConstRowEnumerator;
                    typedef RowEnumerator DefaultRowEnumerator;
                    
                    inline DefaultEnumerator enumerate(const MT & A) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(A.enumerate_col(d_c_ofs)))
                    {
                        return A.enumerate_col(d_c_ofs);
                    }
                    
                    inline DefaultRowEnumerator enumerate_row(size_type row, const MT & A) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(DefaultRowEnumerator(helper::make_type_lvalue<operation_col<MT> >(), row, A)))
                    {
                        return DefaultRowEnumerator(*this, row, A);
                    }
                    
                    inline DefaultColEnumerator enumerate_col(size_type col, const MT & A) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(A.enumerate_col(d_c_ofs)))
                    {
                        assert(col == 0);
                        return A.enumerate_col(d_c_ofs);
                    }
                    
                    class RowsEnumerator
                    {
                    private:
                        DefaultColEnumerator d_col;
                        
                    public:
                        class CustomRowEnumerator
                        {
                        private:
                            const DefaultColEnumerator & d_col;
                            bool d_first;
                            
                        public:
                            typedef CoeffType_Get Type;
                            typedef typename DefaultColEnumerator::GetCoeffSteps_Type GetCoeffSteps_Type;
                            
                            CustomRowEnumerator(const DefaultColEnumerator & col) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                                : d_col(col), d_first(true)
                            {
                            }
                            
                            inline Type current() const
                                PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_col.current()))
                            {
                                return d_col.current();
                            }
                            
                            template<typename Result>
                            inline void get_current(Result & result) const
                                PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_col.get_current(result)))
                            {
                                d_col.get_current(result);
                            }
                            
                            inline GetCoeffSteps_Type get_current_steps() const
                                PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_col.get_current_steps()))
                            {
                                return d_col.get_current_steps();
                            }
                            
                            inline bool has_current() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                            {
                                return d_first;
                            }
                            
                            inline void next() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                            {
                                d_first = false;
                            }
                        };
                        
                        typedef CustomRowEnumerator Type;
                        
                        RowsEnumerator(const ColEnumerator & col)
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(DefaultColEnumerator(col)))
                            : d_col(col)
                        {
                        }
                        
                        inline Type current() const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Type(d_col)))
                        {
                            return Type(d_col);
                        }
                        
                        inline bool has_current() const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_col.has_current()))
                        {
                            return d_col.has_current();
                        }
                        
                        inline void next()
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_col.next()))
                        {
                            d_col.next();
                        }
                    };
                    
                    class ColsEnumerator
                    {
                    private:
                        DefaultColEnumerator d_col;
                        bool d_first;
                        
                    public:
                        typedef DefaultColEnumerator Type;
                        
                        ColsEnumerator(const DefaultColEnumerator & e)
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(DefaultColEnumerator(e)))
                            : d_col(e), d_first(true)
                        {
                        }
                        
                        inline Type current() const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Type(d_col)))
                        {
                            return d_col;
                        }
                        
                        inline bool has_current() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                        {
                            return d_first;
                        }
                        
                        inline void next() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                        {
                            d_first = false;
                        }
                    };
                    
                    friend class RowsEnumerator;
                    friend class ColsEnumerator;
                    
                    typedef RowsEnumerator ConstRowsEnumerator;
                    typedef ColsEnumerator ConstColsEnumerator;
                    typedef RowsEnumerator DefaultRowsEnumerator;
                    typedef ColsEnumerator DefaultColsEnumerator;
                    
                    inline DefaultRowsEnumerator enumerate_rows(const MT & A) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(DefaultRowsEnumerator(helper::make_type_lvalue<const operation_col<MT> >().enumerate_col(0, A))))
                    {
                        return DefaultRowsEnumerator(enumerate_col(0, A));
                    }
                    
                    inline DefaultColsEnumerator enumerate_cols(const MT & A) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(DefaultColsEnumerator(helper::make_type_lvalue<const operation_col<MT> >().enumerate_col(0, A))))
                    {
                        return DefaultColsEnumerator(enumerate_col(0, A));
                    }
                };
                
                template<typename MT>
                class transpose
                {
                public:
                    typedef typename MT::CoeffType CoeffType;
                    typedef typename MT::CoeffType_Get CoeffType_Get;
                    typedef typename MT::GetCoeffSteps_Type GetCoeffSteps_Type;
                    typedef base_matrix<CoeffType,
                                        MatrixInfo<MT>::cols, MatrixInfo<MT>::rows,
                                        typename MatrixInfo<MT>::StorageTraits,
                                        MatrixInfo<MT>::is_math_object> ValueType;
                    typedef void data_pointer_type;
                    typedef void data_ref_type;
                    
                    enum { use_temporary_on_evaluate = MatrixInfo<MT>::use_temporary_on_evaluate,
                           coeffs_are_simple_expressions = true,
                           is_const = MatrixInfo<MT>::is_const,
                           can_assign_to = MatrixInfo<MT>::can_assign_to,
                           can_move_from = MatrixInfo<MT>::can_move_from,
                           can_resize_rows = false,
                           can_resize_cols = false,
                           has_direct_access = false };
                    
                    transpose() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_INTERNAL_STATIC_CHECK(MatrixInfo<MT>::is_matrix, RequiresMatrix);
                    }
                
                    template<typename ResultType>
                    static inline void get_coeff(ResultType & result, size_type i, size_type j, const MT & A)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(result = A(j, i)))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::transpose::get_coeff(" << getAddress(*this) << "; " << getAddress(result) << " " << i << " " << j << " " << getAddress(A) << ") const");
                        assert(i < rows(A));
                        assert(j < cols(A));
                        result = A(j, i);
                    }
                
                    static inline bool get_coeff_alwayszero(const MT & A)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(A.get_coeff_alwayszero()))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::transpose::get_coeff_alwayszero(" << getAddress(*this) << "; " << getAddress(A) << ")");
                        return A.get_coeff_alwayszero();
                    }
                
                    static inline GetCoeffSteps_Type get_coeff_steps(size_type i, size_type j, const MT & A)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(A.get_coeff_steps(j, i)))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::transpose::get_coeff_steps(" << getAddress(*this) << "; " << i << " " << j << " " << getAddress(A) << ")");
                        assert(i < rows(A));
                        assert(j < cols(A));
                        return A.get_coeff_steps(j, i);
                    }
                
                    inline CoeffType_Get operator () (size_type i, size_type j, const MT & A) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(CoeffType_Get(A(j, i))))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::transpose::operator () (" << getAddress(*this) << "; " << i << " " << j << " " << getAddress(A) << ")");
                        assert(i < rows(A));
                        assert(j < cols(A));
                        return A(j, i);
                    }
                
                    static inline size_type rows(const MT & A) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::transpose::rows(" << getAddress(*this) << "; " << getAddress(A) << ")");
                        return A.cols();
                    }
                    
                    static inline size_type cols(const MT & A) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::transpose::cols(" << getAddress(*this) << "; " << getAddress(A) << ")");
                        return A.rows();
                    }
                    
                    static inline size_type size(const MT & A) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::transpose::size(" << getAddress(*this) << "; " << getAddress(A) << ")");
                        return A.size();
                    }
                    
                    template<typename MatrixType>
                    static inline bool involves_this_matrix(const MatrixType & matrix, const MT & A) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::transpose::involves_this_matrix(" << getAddress(*this) << "; " << getAddress(matrix) << " " << getAddress(A) << ")");
                        return A.involves_this_matrix(matrix);
                    }
                    
                    template<typename MatrixType>
                    static inline bool test_involvement(const MatrixType & matrix, const MT & A) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("expressions::transpose::test_involvement(" << getAddress(*this) << "; " << getAddress(matrix) << " " << getAddress(A) << ")");
                        return A.test_involvement(matrix);
                    }
                    
                    class Enumerator : private row_count_storage<MatrixInfo<MT>::rows>, private col_count_storage<MatrixInfo<MT>::cols>
                    {
                    private:
                        typename MT::DefaultColsEnumerator d_cols;
                        typename MT::DefaultColsEnumerator::Type d_col;
                        
                    public:
                        typedef typename MT::DefaultColsEnumerator::Type::Type Type;
                        typedef typename MT::DefaultColsEnumerator::Type::GetCoeffSteps_Type GetCoeffSteps_Type;
                        
                        inline Enumerator(const MT & matrix)
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(typename MT::DefaultColsEnumerator(matrix.enumerate_cols())) &&
                                                                      noexcept(typename MT::DefaultColsEnumerator::Type(d_cols.current())))
                            : d_cols(matrix.enumerate_cols()), d_col(d_cols.current())
                        {
                        }
                        
                        inline Type current() const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Type(d_col.current())))
                        {
                            return d_col.current();
                        }
                        
                        template<typename Result>
                        inline void get_current(Result & result) const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_col.current(result)))
                        {
                            d_col.current(result);
                        }
                        
                        inline GetCoeffSteps_Type get_current_steps() const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_col.get_current_steps()))
                        {
                            return d_col.get_current_steps();
                        }
                        
                        inline bool has_current() const
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_col.has_current()))
                        {
                            return d_col.has_current();
                        }
                        
                        inline void next()
                            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(d_col.next()) &&
                                                                      noexcept(d_col.has_current()) &&
                                                                      noexcept(d_cols.next()) &&
                                                                      noexcept(d_cols.has_current()) &&
                                                                      noexcept(d_col = d_cols.current()))
                        {
                            d_col.next();
                            if (!d_col.has_current())
                            {
                                d_cols.next();
                                if (d_cols.has_current())
                                    d_col = d_cols.current();
                            }
                        }
                    };
                    
                    typedef Enumerator ConstEnumerator;
                    typedef Enumerator DefaultEnumerator;
                    typedef typename MT::ColEnumerator RowEnumerator;
                    typedef typename MT::RowEnumerator ColEnumerator;
                    typedef typename MT::ConstColEnumerator ConstRowEnumerator;
                    typedef typename MT::ConstRowEnumerator ConstColEnumerator;
                    typedef typename MT::DefaultColEnumerator DefaultRowEnumerator;
                    typedef typename MT::DefaultRowEnumerator DefaultColEnumerator;
                    
                    static inline DefaultEnumerator enumerate(const MT & A)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(DefaultEnumerator(A)))
                    {
                        return DefaultEnumerator(A);
                    }
                    
                    static inline DefaultRowEnumerator enumerate_row(size_type row, const MT & A)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(A.enumerate_col(row)))
                    {
                        return A.enumerate_col(row);
                    }
                    
                    static inline DefaultColEnumerator enumerate_col(size_type col, const MT & A)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(A.enumerate_row(col)))
                    {
                        return A.enumerate_row(col);
                    }
                    
                    typedef typename MT::ColsEnumerator RowsEnumerator;
                    typedef typename MT::RowsEnumerator ColsEnumerator;
                    typedef typename MT::ConstColsEnumerator ConstRowsEnumerator;
                    typedef typename MT::ConstRowsEnumerator ConstColsEnumerator;
                    typedef typename MT::DefaultColsEnumerator DefaultRowsEnumerator;
                    typedef typename MT::DefaultRowsEnumerator DefaultColsEnumerator;
                    
                    inline DefaultRowsEnumerator enumerate_rows(const MT & A) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(DefaultRowsEnumerator(A.enumerate_cols())))
                    {
                        return DefaultRowsEnumerator(A.enumerate_cols());
                    }
                    
                    inline DefaultColsEnumerator enumerate_cols(const MT & A) const
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(DefaultColsEnumerator(A.enumerate_rows())))
                    {
                        return DefaultColsEnumerator(A.enumerate_rows());
                    }
                };
            }
            
            template<typename MatrixType>
            struct MatrixInfo<expressions::ConstMatrixWrapper<MatrixType> >
            {
                enum { is_matrix = MatrixInfo<const MatrixType>::is_matrix,
                       is_const = true,
                       is_math_object = MatrixInfo<const MatrixType>::is_math_object,
                       is_expression = true,
                       only_defined_for_matrices = 1,
                       can_move_from = false,
                       can_assign_to = false,
                       can_resize_rows = false,
                       can_resize_cols = false,
                       use_temporary_on_evaluate = MatrixInfo<const MatrixType>::use_temporary_on_evaluate,
                       coeffs_are_simple_expressions = MatrixInfo<const MatrixType>::coeffs_are_simple_expressions };
                enum { rows = MatrixInfo<MatrixType>::rows, cols = MatrixInfo<MatrixType>::cols };
                typedef typename MatrixInfo<MatrixType>::Type Type;
                typedef typename MatrixInfo<MatrixType>::StorageTraits StorageTraits;
            };
            
            template<typename MatrixType>
            struct MatrixInfo<const expressions::ConstMatrixWrapper<MatrixType> >
            {
                enum { is_matrix = MatrixInfo<const MatrixType>::is_matrix,
                       is_const = true,
                       is_math_object = MatrixInfo<const MatrixType>::is_math_object,
                       is_expression = true,
                       only_defined_for_matrices = 1,
                       can_move_from = false,
                       can_assign_to = false,
                       can_resize_rows = false,
                       can_resize_cols = false,
                       use_temporary_on_evaluate = MatrixInfo<const MatrixType>::use_temporary_on_evaluate,
                       coeffs_are_simple_expressions = MatrixInfo<const MatrixType>::coeffs_are_simple_expressions };
                enum { rows = MatrixInfo<MatrixType>::rows, cols = MatrixInfo<MatrixType>::cols };
                typedef typename MatrixInfo<MatrixType>::Type Type;
                typedef typename MatrixInfo<MatrixType>::StorageTraits StorageTraits;
            };
            
            template<typename MatrixType>
            struct MatrixInfo<expressions::MatrixWrapper<MatrixType> >
            {
                enum { is_matrix = MatrixInfo<MatrixType>::is_matrix,
                       is_const = MatrixInfo<MatrixType>::is_const,
                       is_math_object = MatrixInfo<MatrixType>::is_math_object,
                       is_expression = true,
                       only_defined_for_matrices = 1,
                       can_move_from = MatrixInfo<MatrixType>::can_move_from,
                       can_assign_to = MatrixInfo<MatrixType>::can_assign_to,
                       can_resize_rows = MatrixInfo<MatrixType>::can_resize_rows,
                       can_resize_cols = MatrixInfo<MatrixType>::can_resize_cols,
                       use_temporary_on_evaluate = MatrixInfo<MatrixType>::use_temporary_on_evaluate,
                       coeffs_are_simple_expressions = MatrixInfo<MatrixType>::coeffs_are_simple_expressions };
                enum { rows = MatrixInfo<MatrixType>::rows, cols = MatrixInfo<MatrixType>::cols };
                typedef typename MatrixInfo<MatrixType>::Type Type;
                typedef typename MatrixInfo<MatrixType>::StorageTraits StorageTraits;
            };
            
            template<typename MatrixType>
            struct MatrixInfo<const expressions::MatrixWrapper<MatrixType> >
            {
                enum { is_matrix = MatrixInfo<MatrixType>::is_matrix,
                       is_const = MatrixInfo<MatrixType>::is_const,
                       is_math_object = MatrixInfo<MatrixType>::is_math_object,
                       is_expression = true,
                       only_defined_for_matrices = 1,
                       can_move_from = MatrixInfo<MatrixType>::can_move_from,
                       can_assign_to = MatrixInfo<MatrixType>::can_assign_to,
                       can_resize_rows = MatrixInfo<MatrixType>::can_resize_rows,
                       can_resize_cols = MatrixInfo<MatrixType>::can_resize_cols,
                       use_temporary_on_evaluate = MatrixInfo<MatrixType>::use_temporary_on_evaluate,
                       coeffs_are_simple_expressions = MatrixInfo<MatrixType>::coeffs_are_simple_expressions };
                enum { rows = MatrixInfo<MatrixType>::rows, cols = MatrixInfo<MatrixType>::cols };
                typedef typename MatrixInfo<MatrixType>::Type Type;
                typedef typename MatrixInfo<MatrixType>::StorageTraits StorageTraits;
            };
            
            template<typename MatrixType>
            struct MatrixInfo<expressions::MatrixTemporaryWrapper<MatrixType> >
            {
                enum { is_matrix = MatrixInfo<const MatrixType>::is_matrix,
                       is_const = true,
                       is_math_object = MatrixInfo<const MatrixType>::is_math_object,
                       is_expression = true,
                       only_defined_for_matrices = 1,
                       can_move_from = true,
                       can_assign_to = false,
                       can_resize_rows = false,
                       can_resize_cols = false,
                       use_temporary_on_evaluate = false,
                       coeffs_are_simple_expressions = true };
                enum { rows = MatrixInfo<MatrixType>::rows, cols = MatrixInfo<MatrixType>::cols };
                typedef typename MatrixInfo<MatrixType>::Type Type;
                typedef typename MatrixInfo<MatrixType>::StorageTraits StorageTraits;
            };
            
            template<typename MatrixType>
            struct MatrixInfo<const expressions::MatrixTemporaryWrapper<MatrixType> >
            {
                enum { is_matrix = MatrixInfo<const MatrixType>::is_matrix,
                       is_const = true,
                       is_math_object = MatrixInfo<const MatrixType>::is_math_object,
                       is_expression = true,
                       only_defined_for_matrices = 1,
                       can_move_from = true,
                       can_assign_to = false,
                       can_resize_rows = false,
                       can_resize_cols = false,
                       use_temporary_on_evaluate = false,
                       coeffs_are_simple_expressions = true };
                enum { rows = MatrixInfo<MatrixType>::rows, cols = MatrixInfo<MatrixType>::cols };
                typedef typename MatrixInfo<MatrixType>::Type Type;
                typedef typename MatrixInfo<MatrixType>::StorageTraits StorageTraits;
            };
            
            template<typename ScalarType>
            struct MatrixInfo<expressions::ScalarWrapper<ScalarType> >
            {
                enum { is_matrix = false,
                       is_const = true,
                       is_math_object = true,
                       is_expression = true,
                       can_move_from = false,
                       can_assign_to = false,
                       can_resize_rows = false,
                       can_resize_cols = false,
                       use_temporary_on_evaluate = false,
                       coeffs_are_simple_expressions = true };
                typedef ScalarType Type;
            };
            
            template<typename ScalarType>
            struct MatrixInfo<const expressions::ScalarWrapper<ScalarType> >
            {
                enum { is_matrix = false,
                       is_const = true,
                       is_math_object = true,
                       is_expression = true,
                       can_move_from = false,
                       can_assign_to = false,
                       can_resize_rows = false,
                       can_resize_cols = false,
                       use_temporary_on_evaluate = false,
                       coeffs_are_simple_expressions = true };
                typedef ScalarType Type;
            };
            
            // Operations: binary matrix operations
            
            template<typename A, bool A_is_matrix, bool A_needs_wrapper, typename B, bool B_is_matrix, bool B_needs_wrapper>
            class BinaryMatrixOperationImpl_Impl;
        
            template<typename A, typename B>
            class BinaryMatrixOperationImpl_Impl<A, true, false, B, true, false>
            {
            public:
                typedef expressions::expr<expressions::matrix_matrix_multiplication, std::pair<typename A::LazyEvalType, typename B::LazyEvalType> > Mul_ReturnType;
                
                inline static Mul_ReturnType multiply(const A & a, const B & b)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Mul_ReturnType(std::make_pair(a.lazy_evaluate(), b.lazy_evaluate()))))
                {
                    assert(a.cols() == b.rows());
                    return Mul_ReturnType(std::make_pair(a.lazy_evaluate(), b.lazy_evaluate()));
                }
                
                typedef void Div_ReturnType;
                typedef void Mod_ReturnType;
                
                typedef expressions::expr<expressions::componentwise_operation<expressions::CWAdd>::operation, std::pair<A, B> > Add_ReturnType;
                
                inline static Add_ReturnType add(const A & a, const B & b)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Add_ReturnType(std::make_pair(a, b))))
                {
                    assert(a.rows() == b.rows());
                    assert(a.cols() == b.cols());
                    return Add_ReturnType(std::make_pair(a, b));
                }
                
                typedef expressions::expr<expressions::componentwise_operation<expressions::CWSub>::operation, std::pair<A, B> > Sub_ReturnType;
                
                inline static Sub_ReturnType sub(const A & a, const B & b)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Sub_ReturnType(std::make_pair(a, b))))
                {
                    assert(a.rows() == b.rows());
                    assert(a.cols() == b.cols());
                    return Sub_ReturnType(std::make_pair(a, b));
                }
                
                typedef expressions::expr<expressions::componentwise_operation<expressions::CWMul>::operation, std::pair<A, B> > CwMul_ReturnType;
                
                inline static CwMul_ReturnType componentwise_mul(const A & a, const B & b)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(CwMul_ReturnType(std::make_pair(a, b))))
                {
                    assert(a.rows() == b.rows());
                    assert(a.cols() == b.cols());
                    return CwMul_ReturnType(std::make_pair(a, b));
                }
                
                typedef expressions::expr<expressions::componentwise_operation<expressions::CWDiv>::operation, std::pair<A, B> > CwDiv_ReturnType;
                
                inline static CwDiv_ReturnType componentwise_div(const A & a, const B & b)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(CwDiv_ReturnType(std::make_pair(a, b))))
                {
                    assert(a.rows() == b.rows());
                    assert(a.cols() == b.cols());
                    return CwDiv_ReturnType(std::make_pair(a, b));
                }
                
                typedef expressions::expr<expressions::componentwise_operation<expressions::CWMod>::operation, std::pair<A, B> > CwMod_ReturnType;
                
                inline static CwMod_ReturnType componentwise_mod(const A & a, const B & b)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(CwMod_ReturnType(std::make_pair(a, b))))
                {
                    assert(a.rows() == b.rows());
                    assert(a.cols() == b.cols());
                    return CwMod_ReturnType(std::make_pair(a, b));
                }
            };
            
            template<typename A, typename B>
            class BinaryMatrixOperationImpl_Impl<A, true, true, B, true, false>
            {
            public:
                typedef expressions::expr<expressions::matrix_matrix_multiplication,
                                          std::pair<expressions::ConstMatrixWrapper<A>, typename B::LazyEvalType> > Mul_ReturnType;
                
                inline static Mul_ReturnType multiply(const A & a, const B & b)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Mul_ReturnType(std::make_pair(expressions::make_matrix_wrapper(a), b.lazy_evaluate()))))
                {
                    assert(a.cols() == b.rows());
                    return Mul_ReturnType(std::make_pair(expressions::make_matrix_wrapper(a), b.lazy_evaluate()));
                }
                
                typedef void Div_ReturnType;
                typedef void Mod_ReturnType;
                
                typedef expressions::expr<expressions::componentwise_operation<expressions::CWAdd>::operation,
                                          std::pair<expressions::ConstMatrixWrapper<A>, B> > Add_ReturnType;
                
                inline static Add_ReturnType add(const A & a, const B & b)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Add_ReturnType(std::make_pair(expressions::make_matrix_wrapper(a), b))))
                {
                    assert(a.rows() == b.rows());
                    assert(a.cols() == b.cols());
                    return Add_ReturnType(std::make_pair(expressions::make_matrix_wrapper(a), b));
                }
                
                typedef expressions::expr<expressions::componentwise_operation<expressions::CWSub>::operation,
                                          std::pair<expressions::ConstMatrixWrapper<A>, B> > Sub_ReturnType;
                
                inline static Sub_ReturnType sub(const A & a, const B & b)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Sub_ReturnType(std::make_pair(expressions::make_matrix_wrapper(a), b))))
                {
                    assert(a.rows() == b.rows());
                    assert(a.cols() == b.cols());
                    return Sub_ReturnType(std::make_pair(expressions::make_matrix_wrapper(a), b));
                }
                
                typedef expressions::expr<expressions::componentwise_operation<expressions::CWMul>::operation,
                                          std::pair<expressions::ConstMatrixWrapper<A>, B> > CwMul_ReturnType;
                
                inline static CwMul_ReturnType componentwise_mul(const A & a, const B & b)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(CwMul_ReturnType(std::make_pair(expressions::make_matrix_wrapper(a), b))))
                {
                    assert(a.rows() == b.rows());
                    assert(a.cols() == b.cols());
                    return CwMul_ReturnType(std::make_pair(expressions::make_matrix_wrapper(a), b));
                }
                
                typedef expressions::expr<expressions::componentwise_operation<expressions::CWDiv>::operation,
                                          std::pair<expressions::ConstMatrixWrapper<A>, B> > CwDiv_ReturnType;
                
                inline static CwDiv_ReturnType componentwise_div(const A & a, const B & b)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(CwDiv_ReturnType(std::make_pair(expressions::make_matrix_wrapper(a), b))))
                {
                    assert(a.rows() == b.rows());
                    assert(a.cols() == b.cols());
                    return CwDiv_ReturnType(std::make_pair(expressions::make_matrix_wrapper(a), b));
                }
                
                typedef expressions::expr<expressions::componentwise_operation<expressions::CWMod>::operation,
                                          std::pair<expressions::ConstMatrixWrapper<A>, B> > CwMod_ReturnType;
                
                inline static CwMod_ReturnType componentwise_mod(const A & a, const B & b)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(CwMod_ReturnType(std::make_pair(expressions::make_matrix_wrapper(a), b))))
                {
                    assert(a.rows() == b.rows());
                    assert(a.cols() == b.cols());
                    return CwMod_ReturnType(std::make_pair(expressions::make_matrix_wrapper(a), b));
                }
            };
            
            template<typename A, typename B>
            class BinaryMatrixOperationImpl_Impl<A, true, false, B, true, true>
            {
            public:
                typedef expressions::expr<expressions::matrix_matrix_multiplication,
                                          std::pair<typename A::LazyEvalType, expressions::ConstMatrixWrapper<B> > > Mul_ReturnType;
                
                inline static Mul_ReturnType multiply(const A & a, const B & b)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Mul_ReturnType(std::make_pair(a.lazy_evaluate(), expressions::make_matrix_wrapper(b)))))
                {
                    assert(a.cols() == b.rows());
                    return Mul_ReturnType(std::make_pair(a.lazy_evaluate(), expressions::make_matrix_wrapper(b)));
                }
                
                typedef void Div_ReturnType;
                typedef void Mod_ReturnType;
                
                typedef expressions::expr<expressions::componentwise_operation<expressions::CWAdd>::operation,
                                          std::pair<A, expressions::ConstMatrixWrapper<B> > > Add_ReturnType;
                
                inline static Add_ReturnType add(const A & a, const B & b)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Add_ReturnType(std::make_pair(a, expressions::make_matrix_wrapper(b)))))
                {
                    assert(a.rows() == b.rows());
                    assert(a.cols() == b.cols());
                    return Add_ReturnType(std::make_pair(a, expressions::make_matrix_wrapper(b)));
                }
                
                typedef expressions::expr<expressions::componentwise_operation<expressions::CWSub>::operation,
                                          std::pair<A, expressions::ConstMatrixWrapper<B> > > Sub_ReturnType;
                
                inline static Sub_ReturnType sub(const A & a, const B & b)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Sub_ReturnType(std::make_pair(a, expressions::make_matrix_wrapper(b)))))
                {
                    assert(a.rows() == b.rows());
                    assert(a.cols() == b.cols());
                    return Sub_ReturnType(std::make_pair(a, expressions::make_matrix_wrapper(b)));
                }
                
                typedef expressions::expr<expressions::componentwise_operation<expressions::CWMul>::operation,
                                          std::pair<A, expressions::ConstMatrixWrapper<B> > > CwMul_ReturnType;
                
                inline static CwMul_ReturnType componentwise_mul(const A & a, const B & b)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(CwMul_ReturnType(std::make_pair(a, expressions::make_matrix_wrapper(b)))))
                {
                    assert(a.rows() == b.rows());
                    assert(a.cols() == b.cols());
                    return CwMul_ReturnType(std::make_pair(a, expressions::make_matrix_wrapper(b)));
                }
                
                typedef expressions::expr<expressions::componentwise_operation<expressions::CWDiv>::operation,
                                          std::pair<A, expressions::ConstMatrixWrapper<B> > > CwDiv_ReturnType;
                
                inline static CwDiv_ReturnType componentwise_div(const A & a, const B & b)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(CwDiv_ReturnType(std::make_pair(a, expressions::make_matrix_wrapper(b)))))
                {
                    assert(a.rows() == b.rows());
                    assert(a.cols() == b.cols());
                    return CwDiv_ReturnType(std::make_pair(a, expressions::make_matrix_wrapper(b)));
                }
                
                typedef expressions::expr<expressions::componentwise_operation<expressions::CWMod>::operation,
                                          std::pair<A, expressions::ConstMatrixWrapper<B> > > CwMod_ReturnType;
                
                inline static CwMod_ReturnType componentwise_mod(const A & a, const B & b)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(CwMod_ReturnType(std::make_pair(a, expressions::make_matrix_wrapper(b)))))
                {
                    assert(a.rows() == b.rows());
                    assert(a.cols() == b.cols());
                    return CwMod_ReturnType(std::make_pair(a, expressions::make_matrix_wrapper(b)));
                }
            };
            
            template<typename A, typename B>
            class BinaryMatrixOperationImpl_Impl<A, true, true, B, true, true>
            {
            public:
                typedef expressions::expr<expressions::matrix_matrix_multiplication,
                                          std::pair<expressions::ConstMatrixWrapper<A>,
                                                    expressions::ConstMatrixWrapper<B> > > Mul_ReturnType;
                
                inline static Mul_ReturnType multiply(const A & a, const B & b)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Mul_ReturnType(std::make_pair(expressions::make_matrix_wrapper(a),
                                                                                                     expressions::make_matrix_wrapper(b)))))
                {
                    assert(a.cols() == b.rows());
                    return Mul_ReturnType(std::make_pair(expressions::make_matrix_wrapper(a),
                                                         expressions::make_matrix_wrapper(b)));
                }
                
                typedef void Div_ReturnType;
                typedef void Mod_ReturnType;
                
                typedef expressions::expr<expressions::componentwise_operation<expressions::CWAdd>::operation,
                                          std::pair<expressions::ConstMatrixWrapper<A>,
                                                    expressions::ConstMatrixWrapper<B> > > Add_ReturnType;
                
                inline static Add_ReturnType add(const A & a, const B & b)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Add_ReturnType(std::make_pair(expressions::make_matrix_wrapper(a),
                                                                                                     expressions::make_matrix_wrapper(b)))))
                {
                    assert(a.rows() == b.rows());
                    assert(a.cols() == b.cols());
                    return Add_ReturnType(std::make_pair(expressions::make_matrix_wrapper(a),
                                                         expressions::make_matrix_wrapper(b)));
                }
                
                typedef expressions::expr<expressions::componentwise_operation<expressions::CWSub>::operation,
                                          std::pair<expressions::ConstMatrixWrapper<A>,
                                                    expressions::ConstMatrixWrapper<B> > > Sub_ReturnType;
                
                inline static Sub_ReturnType sub(const A & a, const B & b)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Sub_ReturnType(std::make_pair(expressions::make_matrix_wrapper(a),
                                                                                                     expressions::make_matrix_wrapper(b)))))
                {
                    assert(a.rows() == b.rows());
                    assert(a.cols() == b.cols());
                    return Sub_ReturnType(std::make_pair(expressions::make_matrix_wrapper(a),
                                                         expressions::make_matrix_wrapper(b)));
                }
                
                typedef expressions::expr<expressions::componentwise_operation<expressions::CWMul>::operation,
                                          std::pair<expressions::ConstMatrixWrapper<A>,
                                                    expressions::ConstMatrixWrapper<B> > > CwMul_ReturnType;
                
                inline static CwMul_ReturnType componentwise_mul(const A & a, const B & b)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(CwMul_ReturnType(std::make_pair(expressions::make_matrix_wrapper(a),
                                                                                                       expressions::make_matrix_wrapper(b)))))
                {
                    assert(a.rows() == b.rows());
                    assert(a.cols() == b.cols());
                    return CwMul_ReturnType(std::make_pair(expressions::make_matrix_wrapper(a),
                                                           expressions::make_matrix_wrapper(b)));
                }
                
                typedef expressions::expr<expressions::componentwise_operation<expressions::CWDiv>::operation,
                                          std::pair<expressions::ConstMatrixWrapper<A>,
                                                    expressions::ConstMatrixWrapper<B> > > CwDiv_ReturnType;
                
                inline static CwDiv_ReturnType componentwise_div(const A & a, const B & b)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(CwDiv_ReturnType(std::make_pair(expressions::make_matrix_wrapper(a),
                                                                                                       expressions::make_matrix_wrapper(b)))))
                {
                    assert(a.rows() == b.rows());
                    assert(a.cols() == b.cols());
                    return CwDiv_ReturnType(std::make_pair(expressions::make_matrix_wrapper(a),
                                                           expressions::make_matrix_wrapper(b)));
                }
                
                typedef expressions::expr<expressions::componentwise_operation<expressions::CWMod>::operation,
                                          std::pair<expressions::ConstMatrixWrapper<A>,
                                                    expressions::ConstMatrixWrapper<B> > > CwMod_ReturnType;
                
                inline static CwMod_ReturnType componentwise_mod(const A & a, const B & b)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(CwMod_ReturnType(std::make_pair(expressions::make_matrix_wrapper(a),
                                                                                                       expressions::make_matrix_wrapper(b)))))
                {
                    assert(a.rows() == b.rows());
                    assert(a.cols() == b.cols());
                    return CwMod_ReturnType(std::make_pair(expressions::make_matrix_wrapper(a),
                                                           expressions::make_matrix_wrapper(b)));
                }
            };
            
            template<typename A, typename B>
            class BinaryMatrixOperationImpl_Impl<A, false, true, B, true, false>
            {
            public:
                typedef expressions::expr<expressions::matrix_scalar_operation<expressions::SMMul>::operation,
                                          std::pair<B, expressions::ScalarWrapper<A> > > Mul_ReturnType;
            
                inline static Mul_ReturnType multiply(const A & a, const B & b)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Mul_ReturnType(std::make_pair(b, expressions::make_scalar_wrapper(a)))))
                {
                    return Mul_ReturnType(std::make_pair(b, expressions::make_scalar_wrapper(a)));
                }
            
                typedef void Div_ReturnType;
                typedef void Mod_ReturnType;
                typedef void Add_ReturnType;
                typedef void Sub_ReturnType;
                typedef void CwMul_ReturnType;
                typedef void CwDiv_ReturnType;
                typedef void CwMod_ReturnType;
            };
            
            template<typename A, typename B>
            class BinaryMatrixOperationImpl_Impl<A, false, true, B, true, true>
            {
            public:
                typedef expressions::expr<expressions::matrix_scalar_operation<expressions::SMMul>::operation,
                                          std::pair<expressions::ConstMatrixWrapper<B>,
                                                    expressions::ScalarWrapper<A> > > Mul_ReturnType;
                
                inline static Mul_ReturnType multiply(const A & a, const B & b)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Mul_ReturnType(std::make_pair(expressions::make_matrix_wrapper(b),
                                                                                                     expressions::make_scalar_wrapper(a)))))
                {
                    return Mul_ReturnType(std::make_pair(expressions::make_matrix_wrapper(b),
                                                         expressions::make_scalar_wrapper(a)));
                }
                
                typedef void Div_ReturnType;
                typedef void Mod_ReturnType;
                typedef void Add_ReturnType;
                typedef void Sub_ReturnType;
                typedef void CwMul_ReturnType;
                typedef void CwDiv_ReturnType;
                typedef void CwMod_ReturnType;
            };
            
            template<typename A, typename B>
            class BinaryMatrixOperationImpl_Impl<A, true, false, B, false, true>
            {
            public:
                typedef expressions::expr<expressions::matrix_scalar_operation<expressions::MSMul>::operation,
                                          std::pair<A, expressions::ScalarWrapper<B> > > Mul_ReturnType;
                
                inline static Mul_ReturnType multiply(const A & a, const B & b)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Mul_ReturnType(std::make_pair(a, expressions::make_scalar_wrapper(b)))))
                {
                    return Mul_ReturnType(std::make_pair(a, expressions::make_scalar_wrapper(b)));
                }
                
                typedef expressions::expr<expressions::matrix_scalar_operation<expressions::MSDiv>::operation,
                                          std::pair<A, expressions::ScalarWrapper<B> > > Div_ReturnType;
                
                inline static Div_ReturnType divide(const A & a, const B & b)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Div_ReturnType(std::make_pair(a, expressions::make_scalar_wrapper(b)))))
                {
                    return Div_ReturnType(std::make_pair(a, expressions::make_scalar_wrapper(b)));
                }
                
                typedef expressions::expr<expressions::matrix_scalar_operation<expressions::MSMod>::operation,
                                          std::pair<A, expressions::ScalarWrapper<B> > > Mod_ReturnType;
                
                inline static Mod_ReturnType modulo(const A & a, const B & b)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Mod_ReturnType(std::make_pair(a, expressions::make_scalar_wrapper(b)))))
                {
                    return Mod_ReturnType(std::make_pair(a, expressions::make_scalar_wrapper(b)));
                }
                
                typedef void Add_ReturnType;
                typedef void Sub_ReturnType;
                typedef void CwMul_ReturnType;
                typedef void CwDiv_ReturnType;
                typedef void CwMod_ReturnType;
            };
            
            template<typename A, typename B>
            class BinaryMatrixOperationImpl_Impl<A, true, true, B, false, true>
            {
            public:
                typedef expressions::expr<expressions::matrix_scalar_operation<expressions::MSMul>::operation,
                                          std::pair<expressions::ConstMatrixWrapper<A>,
                                                    expressions::ScalarWrapper<B> > > Mul_ReturnType;
                
                inline static Mul_ReturnType multiply(const A & a, const B & b)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Mul_ReturnType(std::make_pair(expressions::make_matrix_wrapper(a),
                                                                                                     expressions::make_scalar_wrapper(b)))))
                {
                    return Mul_ReturnType(std::make_pair(expressions::make_matrix_wrapper(a),
                                                         expressions::make_scalar_wrapper(b)));
                }
                
                typedef expressions::expr<expressions::matrix_scalar_operation<expressions::MSDiv>::operation,
                                          std::pair<expressions::ConstMatrixWrapper<A>,
                                                    expressions::ScalarWrapper<B> > > Div_ReturnType;
                
                inline static Div_ReturnType divide(const A & a, const B & b)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Div_ReturnType(std::make_pair(expressions::make_matrix_wrapper(a),
                                                                                                     expressions::make_scalar_wrapper(b)))))
                {
                    return Div_ReturnType(std::make_pair(expressions::make_matrix_wrapper(a),
                                                         expressions::make_scalar_wrapper(b)));
                }
                
                typedef expressions::expr<expressions::matrix_scalar_operation<expressions::MSMod>::operation,
                                          std::pair<expressions::ConstMatrixWrapper<A>,
                                                    expressions::ScalarWrapper<B> > > Mod_ReturnType;
                
                inline static Mod_ReturnType modulo(const A & a, const B & b)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Mod_ReturnType(std::make_pair(expressions::make_matrix_wrapper(a),
                                                                                                     expressions::make_scalar_wrapper(b)))))
                {
                    return Mod_ReturnType(std::make_pair(expressions::make_matrix_wrapper(a),
                                                         expressions::make_scalar_wrapper(b)));
                }
                
                typedef void Add_ReturnType;
                typedef void Sub_ReturnType;
                typedef void CwMul_ReturnType;
                typedef void CwDiv_ReturnType;
                typedef void CwMod_ReturnType;
            };
            
            template<typename A, typename B>
            class BinaryMatrixOperationImpl
            {
            public:
                typedef BinaryMatrixOperationImpl_Impl<A, MatrixInfo<A>::is_matrix, !MatrixInfo<A>::is_expression,
                                                       B, MatrixInfo<B>::is_matrix, !MatrixInfo<B>::is_expression> BMOII;
                typedef typename BMOII::Mul_ReturnType Mul_ReturnType;
                typedef typename BMOII::Div_ReturnType Div_ReturnType;
                typedef typename BMOII::Mod_ReturnType Mod_ReturnType;
                typedef typename BMOII::Add_ReturnType Add_ReturnType;
                typedef typename BMOII::Sub_ReturnType Sub_ReturnType;
                typedef typename BMOII::CwMul_ReturnType CwMul_ReturnType;
                typedef typename BMOII::CwDiv_ReturnType CwDiv_ReturnType;
                typedef typename BMOII::CwMod_ReturnType CwMod_ReturnType;
                
                inline static Mul_ReturnType multiply(const A & a, const B & b)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(BMOII::multiply(a, b)))
                {
                    return BMOII::multiply(a, b);
                }
                
                inline static Div_ReturnType divide(const A & a, const B & b)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(BMOII::divide(a, b)))
                {
                    return BMOII::divide(a, b);
                }
                
                inline static Mod_ReturnType modulo(const A & a, const B & b)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(BMOII::modulo(a, b)))
                {
                    return BMOII::modulo(a, b);
                }
                
                inline static Add_ReturnType add(const A & a, const B & b)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(BMOII::add(a, b)))
                {
                    return BMOII::add(a, b);
                }
                
                inline static Sub_ReturnType sub(const A & a, const B & b)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(BMOII::sub(a, b)))
                {
                    return BMOII::sub(a, b);
                }
                
                inline static CwMul_ReturnType componentwise_mul(const A & a, const B & b)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(BMOII::componentwise_mul(a, b)))
                {
                    return BMOII::componentwise_mul(a, b);
                }
                
                inline static CwDiv_ReturnType componentwise_div(const A & a, const B & b)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(BMOII::componentwise_div(a, b)))
                {
                    return BMOII::componentwise_div(a, b);
                }
                
                inline static CwMod_ReturnType componentwise_mod(const A & a, const B & b)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(BMOII::componentwise_mod(a, b)))
                {
                    return BMOII::componentwise_mod(a, b);
                }
            };
            
            // operations: unary operations
            
            template<typename A, bool A_needs_wrapper>
            class UnaryMatrixOperationImpl_Impl;
            
            template<typename A>
            class UnaryMatrixOperationImpl_Impl<A, false>
            {
            public:
                typedef expressions::expr<expressions::matrix_negation, A> Neg_ReturnType;
                
                inline static Neg_ReturnType negate(const A & a)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Neg_ReturnType(a)))
                {
                    return Neg_ReturnType(a);
                }
            };
            
            template<typename A>
            class UnaryMatrixOperationImpl_Impl<A, true>
            {
            public:
                typedef expressions::expr<expressions::matrix_negation, expressions::ConstMatrixWrapper<A> > Neg_ReturnType;
                
                inline static Neg_ReturnType negate(const A & a)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(Neg_ReturnType(expressions::make_matrix_wrapper(a))))
                {
                    return Neg_ReturnType(expressions::make_matrix_wrapper(a));
                }
            };
            
            template<typename A>
            class UnaryMatrixOperationImpl
            {
            public:
                typedef UnaryMatrixOperationImpl_Impl<A, !MatrixInfo<A>::is_expression> UMOII;
                typedef typename UMOII::Neg_ReturnType Neg_ReturnType;
                
                inline static Neg_ReturnType negate(const A & a)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(UMOII::negate(a)))
                {
                    return UMOII::negate(a);
                }
            };
            
            // operations: assign-operations
            
            template<typename MT1, typename MT2>
            class AssignmentMatrixOperationImpl
            {
            private:
                // A += B
                
                template<typename MT2_>
                static inline void add_assign_impl(const MT1 & A, const MT2_ & B)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(A.enumerate()) && noexcept(B.enumerate()) &&
                                                              noexcept(helper::make_type_lvalue<typename MT1::Enumerator>().has_current()) &&
                                                              noexcept(helper::make_type_lvalue<typename MT2_::ConstEnumerator>().has_current()) &&
                                                              noexcept(helper::make_type_lvalue<typename MT1::Enumerator>().next()) &&
                                                              noexcept(helper::make_type_lvalue<typename MT2_::ConstEnumerator>().next()) &&
                                                              noexcept(helper::make_type_lvalue<typename MT1::Enumerator>().current() += helper::make_type_lvalue<typename MT2_::ConstEnumerator>().current()))
                {
                    typename MT1::Enumerator eA = A.enumerate();
                    typename MT2_::ConstEnumerator eB = B.enumerate();
                    for (; eA.has_current(); eA.next(), eB.next())
                        eA.current() += eB.current();
                }
                
                // A -= B
                
                template<typename MT2_>
                static inline void sub_assign_impl(const MT1 & A, const MT2_ & B)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(A.enumerate()) && noexcept(B.enumerate()) &&
                                                              noexcept(helper::make_type_lvalue<typename MT1::Enumerator>().has_current()) &&
                                                              noexcept(helper::make_type_lvalue<typename MT2_::ConstEnumerator>().has_current()) &&
                                                              noexcept(helper::make_type_lvalue<typename MT1::Enumerator>().next()) &&
                                                              noexcept(helper::make_type_lvalue<typename MT2_::ConstEnumerator>().next()) &&
                                                              noexcept(helper::make_type_lvalue<typename MT1::Enumerator>().current() -= helper::make_type_lvalue<typename MT2_::ConstEnumerator>().current()))
                {
                    typename MT1::Enumerator eA = A.enumerate();
                    typename MT2_::ConstEnumerator eB = B.enumerate();
                    for (; eA.has_current(); eA.next(), eB.next())
                        eA.current() -= eB.current();
                }
                
                // A *= B
                
                template<typename MT2_>
                static inline void mul_assign(helper::BoolToType<true>, const MT1 & A, const MT2_ & B) // A *= B with B matrix
#if __cplusplus >= 201103L
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(math_matrix<typename MatrixInfo<MT1>::Type, MatrixInfo<MT1>::rows, MatrixInfo<MT2_>::cols, typename MatrixInfo<MT1>::StorageTraits>(A * B)) &&
                                                              noexcept(A = std::move(math_matrix<typename MatrixInfo<MT1>::Type, MatrixInfo<MT1>::rows, MatrixInfo<MT2_>::cols, typename MatrixInfo<MT1>::StorageTraits>(A * B))))
#endif
                {
                    assert(A.cols() == B.rows());
                    // matrix-matrix multiplication
                    PLLL_INTERNAL_STATIC_CHECK((MatrixInfo<MT1>::cols < 0) || (MatrixInfo<MT2_>::rows < 0) ||
                                               (static_cast<size_type>(MatrixInfo<MT1>::cols) == static_cast<size_type>(MatrixInfo<MT2_>::rows)), FormatsDoNotMatch);
                    math_matrix<typename MatrixInfo<MT1>::Type, MatrixInfo<MT1>::rows, MatrixInfo<MT2_>::cols, typename MatrixInfo<MT1>::StorageTraits>
                        temporary = A * B;
#if __cplusplus >= 201103L
                    A = std::move(temporary);
#else
                    swap(A, temporary);
#endif
                }
                
                template<typename MT2_>
                static inline void mul_assign_cw_impl(const MT1 & A, const MT2_ & B)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(A.enumerate()) &&
                                                              noexcept(helper::make_type_lvalue<typename MT1::Enumerator>().has_current()) &&
                                                              noexcept(helper::make_type_lvalue<typename MT1::Enumerator>().next()) &&
                                                              noexcept(helper::make_type_lvalue<typename MT1::Enumerator>().current() /= B))
                {
                    typename MT1::Enumerator eA = A.enumerate();
                    for (; eA.has_current(); eA.next())
                        eA.current() *= B;
                }
                
                template<typename MT2_>
                static inline void mul_assign(helper::BoolToType<false>, const MT1 & A, const MT2_ & B) // A *= B with B scalar
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(mul_assign_cw_impl(A, B)))
                {
                    // matrix-scalar multiplication
                    mul_assign_cw_impl(A, B);
                }
                
                // A /= B
                
                template<typename MT2_>
                static inline void div_assign_impl(const MT1 & A, const MT2_ & B)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(A.enumerate()) &&
                                                              noexcept(helper::make_type_lvalue<typename MT1::Enumerator>().has_current()) &&
                                                              noexcept(helper::make_type_lvalue<typename MT1::Enumerator>().next()) &&
                                                              noexcept(helper::make_type_lvalue<typename MT1::Enumerator>().current() /= B))
                {
                    typename MT1::Enumerator eA = A.enumerate();
                    for (; eA.has_current(); eA.next())
                        eA.current() /= B;
                }
                
                // A %= B
                
                template<typename MT2_>
                static inline void mod_assign_impl(const MT1 & A, const MT2_ & B)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(A.enumerate()) &&
                                                              noexcept(helper::make_type_lvalue<typename MT1::Enumerator>().has_current()) &&
                                                              noexcept(helper::make_type_lvalue<typename MT1::Enumerator>().next()) &&
                                                              noexcept(helper::make_type_lvalue<typename MT1::Enumerator>().current() %= B))
                {
                    typename MT1::Enumerator eA = A.enumerate();
                    for (; eA.has_current(); eA.next())
                        eA.current() /= B;
                }
                
            public:
                static inline const MT1 & add_assign(const MT1 & A, const MT2 & B) // A += B
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(MatrixInfo<MT2>::use_temporary_on_evaluate ?
                                                              (noexcept(math_matrix<typename MatrixInfo<MT1>::Type, MatrixInfo<MT1>::rows, MatrixInfo<MT1>::cols, typename MatrixInfo<MT1>::StorageTraits>(A * B)) &&
                                                               noexcept(add_assign_impl(A, helper::make_type_lvalue<math_matrix<typename MatrixInfo<MT1>::Type, MatrixInfo<MT1>::rows, MatrixInfo<MT1>::cols, typename MatrixInfo<MT1>::StorageTraits> >()))) :
                                                              (noexcept(add_assign_impl(A, B)) &&
                                                               noexcept(math_matrix<typename MatrixInfo<MT1>::Type, MatrixInfo<MT1>::rows, MatrixInfo<MT1>::cols, typename MatrixInfo<MT1>::StorageTraits>(A * B)) &&
                                                               noexcept(add_assign_impl(A, helper::make_type_lvalue<math_matrix<typename MatrixInfo<MT1>::Type, MatrixInfo<MT1>::rows, MatrixInfo<MT1>::cols, typename MatrixInfo<MT1>::StorageTraits> >()))))
                {
                    PLLL_INTERNAL_STATIC_CHECK(!MatrixInfo<MT1>::is_const && MatrixInfo<MT1>::can_assign_to, RequiresAssignNonconstMatrix);
                    PLLL_INTERNAL_STATIC_CHECK(MatrixInfo<MT1>::is_matrix && MatrixInfo<MT1>::is_math_object, RequiresMathMatrix);
                    PLLL_INTERNAL_STATIC_CHECK(MatrixInfo<MT2>::is_matrix && MatrixInfo<MT2>::is_math_object, RequiresMathMatrix);
                    PLLL_INTERNAL_STATIC_CHECK((MatrixInfo<MT1>::rows < 0) || (MatrixInfo<MT2>::rows < 0) ||
                                               (static_cast<size_type>(MatrixInfo<MT1>::rows) == static_cast<size_type>(MatrixInfo<MT2>::rows)), FormatsDoNotMatch);
                    PLLL_INTERNAL_STATIC_CHECK((MatrixInfo<MT1>::cols < 0) || (MatrixInfo<MT2>::cols < 0) ||
                                               (static_cast<size_type>(MatrixInfo<MT1>::cols) == static_cast<size_type>(MatrixInfo<MT2>::cols)), FormatsDoNotMatch);
                    assert(A.rows() == B.rows());
                    assert(A.cols() == B.cols());
                    if (MatrixInfo<MT2>::use_temporary_on_evaluate && A.test_involvement(B))
                    {
                        math_matrix<typename MatrixInfo<MT1>::Type, MatrixInfo<MT1>::rows, MatrixInfo<MT1>::cols, typename MatrixInfo<MT1>::StorageTraits>
                            temporary = B;
                        add_assign_impl(A, temporary);
                    }
                    else
                        add_assign_impl(A, B);
                    return A;
                }
                
                static inline const MT1 & sub_assign(const MT1 & A, const MT2 & B) // A -= B
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(MatrixInfo<MT2>::use_temporary_on_evaluate ?
                                                              (noexcept(math_matrix<typename MatrixInfo<MT1>::Type, MatrixInfo<MT1>::rows, MatrixInfo<MT1>::cols, typename MatrixInfo<MT1>::StorageTraits>(A * B)) &&
                                                               noexcept(sub_assign_impl(A, helper::make_type_lvalue<math_matrix<typename MatrixInfo<MT1>::Type, MatrixInfo<MT1>::rows, MatrixInfo<MT1>::cols, typename MatrixInfo<MT1>::StorageTraits> >()))) :
                                                              (noexcept(sub_assign_impl(A, B)) &&
                                                               noexcept(math_matrix<typename MatrixInfo<MT1>::Type, MatrixInfo<MT1>::rows, MatrixInfo<MT1>::cols, typename MatrixInfo<MT1>::StorageTraits>(A * B)) &&
                                                               noexcept(sub_assign_impl(A, helper::make_type_lvalue<math_matrix<typename MatrixInfo<MT1>::Type, MatrixInfo<MT1>::rows, MatrixInfo<MT1>::cols, typename MatrixInfo<MT1>::StorageTraits> >()))))
                {
                    PLLL_INTERNAL_STATIC_CHECK(!MatrixInfo<MT1>::is_const && MatrixInfo<MT1>::can_assign_to, RequiresAssignNonconstMatrix);
                    PLLL_INTERNAL_STATIC_CHECK(MatrixInfo<MT1>::is_matrix && MatrixInfo<MT1>::is_math_object, RequiresMathMatrix);
                    PLLL_INTERNAL_STATIC_CHECK(MatrixInfo<MT2>::is_matrix && MatrixInfo<MT2>::is_math_object, RequiresMathMatrix);
                    PLLL_INTERNAL_STATIC_CHECK((MatrixInfo<MT1>::rows < 0) || (MatrixInfo<MT2>::rows < 0) ||
                                               (static_cast<size_type>(MatrixInfo<MT1>::rows) == static_cast<size_type>(MatrixInfo<MT2>::rows)), FormatsDoNotMatch);
                    PLLL_INTERNAL_STATIC_CHECK((MatrixInfo<MT1>::cols < 0) || (MatrixInfo<MT2>::cols < 0) ||
                                               (static_cast<size_type>(MatrixInfo<MT1>::cols) == static_cast<size_type>(MatrixInfo<MT2>::cols)), FormatsDoNotMatch);
                    assert(A.rows() == B.rows());
                    assert(A.cols() == B.cols());
                    if (MatrixInfo<MT2>::use_temporary_on_evaluate && A.test_involvement(B))
                    {
                        math_matrix<typename MatrixInfo<MT1>::Type, MatrixInfo<MT1>::rows, MatrixInfo<MT1>::cols, typename MatrixInfo<MT1>::StorageTraits>
                            temporary = B;
                        sub_assign_impl(A, temporary);
                    }
                    else
                        sub_assign_impl(A, B);
                    return A;
                }
                
                static inline const MT1 & mul_assign(const MT1 & A, const MT2 & B) // A *= B
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(mul_assign(helper::BoolToType<MatrixInfo<MT2>::is_matrix>(), A, B)))
                {
                    PLLL_INTERNAL_STATIC_CHECK(!MatrixInfo<MT1>::is_const && MatrixInfo<MT1>::can_assign_to, RequiresAssignNonconstMatrix);
                    PLLL_INTERNAL_STATIC_CHECK(MatrixInfo<MT1>::is_matrix && MatrixInfo<MT1>::is_math_object, RequiresMathMatrix);
                    PLLL_INTERNAL_STATIC_CHECK(!MatrixInfo<MT2>::is_matrix || MatrixInfo<MT2>::is_math_object, RequiresMathObject);
                    mul_assign(helper::BoolToType<MatrixInfo<MT2>::is_matrix>(), A, B);
                    return A;
                }
                
                static inline const MT1 & div_assign(const MT1 & A, const MT2 & B) // A /= B
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(div_assign_impl(A, B)))
                {
                    PLLL_INTERNAL_STATIC_CHECK(!MatrixInfo<MT1>::is_const && MatrixInfo<MT1>::can_assign_to, RequiresAssignNonconstMatrix);
                    PLLL_INTERNAL_STATIC_CHECK(MatrixInfo<MT1>::is_matrix && MatrixInfo<MT1>::is_math_object, RequiresMathMatrix);
                    PLLL_INTERNAL_STATIC_CHECK(!MatrixInfo<MT2>::is_matrix, RequiresScalar);
                    div_assign_impl(A, B);
                    return A;
                }
                
                static inline const MT1 & mod_assign(const MT1 & A, const MT2 & B) // A %= B
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(mod_assign_impl(A, B)))
                {
                    PLLL_INTERNAL_STATIC_CHECK(!MatrixInfo<MT1>::is_const && MatrixInfo<MT1>::can_assign_to, RequiresAssignNonconstMatrix);
                    PLLL_INTERNAL_STATIC_CHECK(MatrixInfo<MT1>::is_matrix && MatrixInfo<MT1>::is_math_object, RequiresMathMatrix);
                    PLLL_INTERNAL_STATIC_CHECK(!MatrixInfo<MT2>::is_matrix, RequiresScalar);
                    mod_assign_impl(A, B);
                    return A;
                }
            };
        }
    }
}

#endif

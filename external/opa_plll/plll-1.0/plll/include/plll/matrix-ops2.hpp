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

#ifndef PLLL_INCLUDE_GUARD__MATRIX_OPS2_HPP
#define PLLL_INCLUDE_GUARD__MATRIX_OPS2_HPP

/**
   \file
   \brief Operator instantiations for matrices and vectors.
   
   This header contains instantiations of the abstract templates in `matrix-ops.hpp`. They provide
   implementations to essentially all operations on matrices.
*/
namespace plll
{
    namespace implemenation
    {
    }
    
    namespace linalg
    {
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        //// swap
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        /**@{
           \name Swap functions.
        */
        
        template<typename T, int R, int C, typename ST, bool MO>
        void swap(base_matrix<T, R, C, ST, MO> & A, base_matrix<T, R, C, ST, MO> & B) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("swap(" << getAddress(A) << " " << getAddress(B) << ") [0]");
            std::swap(A.d_data, B.d_data);
            A.implementation::template row_count_storage<R>::swap(B);
            A.implementation::template col_count_storage<C>::swap(B);
        }
        
        template<typename T, int R1, int C1, int R2, int C2, typename ST, bool MO1, bool MO2>
        void swap(base_matrix<T, R1, C1, ST, MO1> & A, base_matrix<T, R2, C2, ST, MO2> & B) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("swap(" << getAddress(A) << " " << getAddress(B) << ") [1]");
            PLLL_INTERNAL_STATIC_CHECK((R1 < 0) || (R2 < 0) || (R1 == R2), TypeAShouldNotBeConst);
            PLLL_INTERNAL_STATIC_CHECK((C1 < 0) || (C2 < 0) || (C1 == C2), TypeAShouldNotBeConst);
            if (((R1 >= 0) && (R2 < 0)) || ((R2 >= 0) && (R1 < 0))) assert(A.rows() == B.rows());
            if (((C1 >= 0) && (C2 < 0)) || ((C2 >= 0) && (C1 < 0))) assert(A.cols() == B.cols());
            std::swap(A.d_data, B.d_data);
            A.implementation::template row_count_storage<R1>::swap(B);
            A.implementation::template col_count_storage<C1>::swap(B);
        }
        
        namespace implementation
        {
            using std::swap;
            
            template<template<typename AData> class AOp, typename AData, template<typename BData> class BOp, typename BData>
            void do_swap(const implementation::expressions::expr<AOp, AData> & A, const implementation::expressions::expr<BOp, BData> & B)
                PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(A.enumerate()) && noexcept(B.enumerate()) &&
                                                          noexcept(helper::make_type_lvalue<typename implementation::expressions::expr<AOp, AData>::Enumerator>().has_current()) &&
                                                          noexcept(helper::make_type_lvalue<typename implementation::expressions::expr<AOp, AData>::Enumerator>().next()) &&
                                                          noexcept(helper::make_type_lvalue<typename implementation::expressions::expr<BOp, BData>::Enumerator>().next()) &&
                                                          noexcept(swap(helper::make_type_lvalue<typename implementation::expressions::expr<AOp, AData>::Enumerator>().current(),
                                                                        helper::make_type_lvalue<typename implementation::expressions::expr<BOp, BData>::Enumerator>().current())))
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("swap(" << getAddress(A) << " " << getAddress(B) << ") [2]");
                PLLL_INTERNAL_STATIC_CHECK(!(implementation::MatrixInfo<implementation::expressions::expr<AOp, AData> >::is_const), TypeAShouldNotBeConst);
                PLLL_INTERNAL_STATIC_CHECK(!(implementation::MatrixInfo<implementation::expressions::expr<BOp, BData> >::is_const), TypeBShouldNotBeConst);
                PLLL_INTERNAL_STATIC_CHECK((implementation::MatrixInfo<implementation::expressions::expr<AOp, AData> >::rows < 0) || (implementation::MatrixInfo<implementation::expressions::expr<BOp, BData> >::rows < 0) ||
                                           (static_cast<size_type>(implementation::MatrixInfo<implementation::expressions::expr<AOp, AData> >::rows) == static_cast<size_type>(implementation::MatrixInfo<implementation::expressions::expr<BOp, BData> >::rows)), NeedSameNumberOfRows);
                PLLL_INTERNAL_STATIC_CHECK((implementation::MatrixInfo<implementation::expressions::expr<AOp, AData> >::cols < 0) || (implementation::MatrixInfo<implementation::expressions::expr<BOp, BData> >::cols < 0) ||
                                           (static_cast<size_type>(implementation::MatrixInfo<implementation::expressions::expr<AOp, AData> >::cols) == static_cast<size_type>(implementation::MatrixInfo<implementation::expressions::expr<BOp, BData> >::cols)), NeedSameNumberOfCols);
                assert(A.rows() == B.rows());
                assert(A.cols() == B.cols());
                typename implementation::expressions::expr<AOp, AData>::Enumerator eA = A.enumerate();
                typename implementation::expressions::expr<BOp, BData>::Enumerator eB = B.enumerate();
                for (; eA.has_current(); eA.next(), eB.next())
                {
                    using std::swap;
                    swap(eA.current(), eB.current());
                }
            }
        }
        
        template<template<typename AData> class AOp, typename AData, template<typename BData> class BOp, typename BData>
        void swap(const implementation::expressions::expr<AOp, AData> & A, const implementation::expressions::expr<BOp, BData> & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::do_swap(A, B)))
        {
            implementation::do_swap(A, B);
        }
        
        template<template<typename AData> class AOp, typename AData,
                 typename BT, int BRows, int BCols, typename BST, bool BMO>
        void swap(const implementation::expressions::expr<AOp, AData> & A, base_matrix<BT, BRows, BCols, BST, BMO> & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(linalg::swap(A, implementation::expressions::make_matrix_expression(B))))
        {
            linalg::swap(A, implementation::expressions::make_matrix_expression(B));
        }
        
        template<typename AT, int ARows, int ACols, typename AST, bool AMO,
                 template<typename BData> class BOp, typename BData>
        void swap(base_matrix<AT, ARows, ACols, AST, AMO> & A, const implementation::expressions::expr<BOp, BData> & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(linalg::swap(implementation::expressions::make_matrix_expression(A), B)))
        {
            linalg::swap(implementation::expressions::make_matrix_expression(A), B);
        }
        
        template<typename T,
                 int ARows, int ACols, typename AST, bool AMO,
                 int BRows, int BCols, typename BST, bool BMO>
        void swap(base_matrix<T, ARows, ACols, AST, AMO> & A, base_matrix<T, BRows, BCols, BST, BMO> & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(linalg::swap(implementation::expressions::make_matrix_expression(A), implementation::expressions::make_matrix_expression(B))))
        {
            linalg::swap(implementation::expressions::make_matrix_expression(A), implementation::expressions::make_matrix_expression(B));
        }
        
        ///@}
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        //// assignment
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        /**@{
           \name Assignment functions.
        */
        
        namespace implementation
        {
            template<template<typename SourceData> class SourceOp, typename SourceData, template<typename DestData> class DestOp, typename DestData>
            inline void assign_impl(const expressions::expr<DestOp, DestData> & destination,
                                    const expressions::expr<SourceOp, SourceData> & source, helper::BoolToType<false>)
                PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(destination.enumerate()) && noexcept(source.enumerate()) &&
                                                          noexcept(helper::make_type_lvalue<typename expressions::expr<DestOp, DestData>::Enumerator>().has_current()) &&
                                                          noexcept(helper::make_type_lvalue<typename expressions::expr<DestOp, DestData>::Enumerator>().next()) &&
                                                          noexcept(helper::make_type_lvalue<typename expressions::expr<SourceOp, SourceData>::ConstEnumerator>().next()) &&
                                                          noexcept(helper::make_type_lvalue<typename expressions::expr<DestOp, DestData>::Enumerator>().current() =
                                                                   helper::make_type_lvalue<typename expressions::expr<SourceOp, SourceData>::ConstEnumerator>().current()))
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("assign_impl(" << getAddress(destination) << " " << getAddress(source) << ") [1-1]");
                typename expressions::expr<DestOp, DestData>::Enumerator eDest = destination.enumerate();
                typename expressions::expr<SourceOp, SourceData>::ConstEnumerator eSrc = source.enumerate();
                for (; eDest.has_current(); eDest.next(), eSrc.next())
                    eDest.current() = eSrc.current();
            }
            
#if __cplusplus >= 201103L
            template<template<typename SourceData> class SourceOp, typename SourceData, template<typename DestData> class DestOp, typename DestData>
            inline void assign_impl(const expressions::expr<DestOp, DestData> & destination,
                                    const expressions::expr<SourceOp, SourceData> & source, helper::BoolToType<true> move)
                PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(destination.enumerate()) && noexcept(source.enumerate()) &&
                                                          noexcept(helper::make_type_lvalue<typename expressions::expr<DestOp, DestData>::Enumerator>().has_current()) &&
                                                          noexcept(helper::make_type_lvalue<typename expressions::expr<DestOp, DestData>::Enumerator>().next()) &&
                                                          noexcept(helper::make_type_lvalue<typename expressions::expr<SourceOp, SourceData>::ConstEnumerator>().next()) &&
                                                          (MatrixInfo<expressions::expr<SourceOp, SourceData> >::can_move_from ? 
                                                           noexcept(helper::make_type_lvalue<typename expressions::expr<DestOp, DestData>::Enumerator>().current() =
                                                                    std::move(helper::make_type_lvalue<typename expressions::expr<SourceOp, SourceData>::ConstEnumerator>().current())) :
                                                           noexcept(helper::make_type_lvalue<typename expressions::expr<DestOp, DestData>::Enumerator>().current() =
                                                                    helper::make_type_lvalue<typename expressions::expr<SourceOp, SourceData>::ConstEnumerator>().current())))
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("assign_impl(" << getAddress(destination) << " " << getAddress(source) << ") [1-1]");
                typename expressions::expr<DestOp, DestData>::Enumerator eDest = destination.enumerate();
                typename expressions::expr<SourceOp, SourceData>::ConstEnumerator eSrc = source.enumerate();
                for (; eDest.has_current(); eDest.next(), eSrc.next())
                    if (MatrixInfo<expressions::expr<SourceOp, SourceData> >::can_move_from)
                        eDest.current() = std::move(eSrc.current());
                    else
                        eDest.current() = eSrc.current();
            }
#else
            template<template<typename SourceData> class SourceOp, typename SourceData, template<typename DestData> class DestOp, typename DestData>
            inline void assign_impl(const expressions::expr<DestOp, DestData> & destination,
                                    const expressions::expr<SourceOp, SourceData> & source, helper::BoolToType<true> move)
                PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(destination.enumerate()) && noexcept(source.enumerate()) &&
                                                          noexcept(helper::make_type_lvalue<typename expressions::expr<DestOp, DestData>::Enumerator>().has_current()) &&
                                                          noexcept(helper::make_type_lvalue<typename expressions::expr<DestOp, DestData>::Enumerator>().next()) &&
                                                          noexcept(helper::make_type_lvalue<typename expressions::expr<SourceOp, SourceData>::ConstEnumerator>().next()) &&
                                                          noexcept(helper::make_type_lvalue<typename expressions::expr<DestOp, DestData>::Enumerator>().current() =
                                                                   helper::make_type_lvalue<typename expressions::expr<SourceOp, SourceData>::ConstEnumerator>().current()))
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("assign_impl(" << getAddress(destination) << " " << getAddress(source) << ") [1-1]");
                typename expressions::expr<DestOp, DestData>::Enumerator eDest = destination.enumerate();
                typename expressions::expr<SourceOp, SourceData>::ConstEnumerator eSrc = source.enumerate();
                for (; eDest.has_current(); eDest.next(), eSrc.next())
                {
                    eDest.current() = eSrc.current();
                    // TODO: add "compatibility move" !!! ??? ...
                }
            }
#endif
            
            template<template<typename SourceData> class SourceOp, typename SourceData, template<typename DestData> class DestOp, typename DestData>
            void assign_impl_resize(const expressions::expr<DestOp, DestData> & destination,
                                    const expressions::expr<SourceOp, SourceData> & source, helper::BoolToType<true> move, helper::BoolToType<true> canresize)
            {
#if __cplusplus >= 201103L
                destination.move_resize(std::move(source));
#else
                destination.move_resize(source);
#endif
            }
            
            template<template<typename SourceData> class SourceOp, typename SourceData, template<typename DestData> class DestOp, typename DestData>
            void assign_impl_resize(const expressions::expr<DestOp, DestData> & destination,
                                    const expressions::expr<SourceOp, SourceData> & source, helper::BoolToType<false> move, helper::BoolToType<true> canresize)
            {
                destination.assign_resize(source);
            }
            
            template<template<typename SourceData> class SourceOp, typename SourceData, template<typename DestData> class DestOp, typename DestData, bool move>
            void assign_impl_resize(const expressions::expr<DestOp, DestData> & destination,
                                    const expressions::expr<SourceOp, SourceData> & source, helper::BoolToType<move>, helper::BoolToType<false> canresize)
            {
                assert(!"Trying to resize a matrix object which is not resizeable!");
            }
            
            template<template<typename SourceData> class SourceOp, typename SourceData, template<typename DestData> class DestOp, typename DestData, bool move>
            void assign_with_temporary(const expressions::expr<DestOp, DestData> & destination,
                                       const expressions::expr<SourceOp, SourceData> & source, helper::BoolToType<move> ittm, helper::BoolToType<true> with_temporary)
            {
                if (destination.test_involvement(source))
                    assign_impl(destination, make_matrix_temporary_expression(source), helper::BoolToType<true>());
                else
                    assign_impl(destination, source, ittm);
            }
            
            template<template<typename SourceData> class SourceOp, typename SourceData, template<typename DestData> class DestOp, typename DestData, bool move>
            void assign_with_temporary(const expressions::expr<DestOp, DestData> & destination,
                                       const expressions::expr<SourceOp, SourceData> & source, helper::BoolToType<move> ittm, helper::BoolToType<false> with_temporary)
                PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(assign_impl(destination, source, ittm)))
            {
                assign_impl(destination, source, ittm);
            }
        }
        
        template<template<typename SourceData> class SourceOp, typename SourceData, template<typename DestData> class DestOp, typename DestData, bool move>
        void assign(const implementation::expressions::expr<DestOp, DestData> & destination,
                    const implementation::expressions::expr<SourceOp, SourceData> & source, helper::BoolToType<move> ittm)
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("assign(" << getAddress(destination) << " " << getAddress(source) << ") [1]");
            PLLL_INTERNAL_STATIC_CHECK((!implementation::MatrixInfo<implementation::expressions::expr<DestOp, DestData> >::is_const), DestinationShouldNotBeConst);
            PLLL_INTERNAL_STATIC_CHECK((implementation::MatrixInfo<implementation::expressions::expr<DestOp, DestData> >::rows < 0) || (implementation::MatrixInfo<implementation::expressions::expr<SourceOp, SourceData> >::rows < 0) ||
                                       (static_cast<int>(implementation::MatrixInfo<implementation::expressions::expr<DestOp, DestData> >::rows) == static_cast<int>(implementation::MatrixInfo<implementation::expressions::expr<SourceOp, SourceData> >::rows)) ||
                                       ((static_cast<int>(implementation::MatrixInfo<implementation::expressions::expr<DestOp, DestData> >::rows) != static_cast<int>(implementation::MatrixInfo<implementation::expressions::expr<SourceOp, SourceData> >::rows))
                                          && implementation::MatrixInfo<implementation::expressions::expr<DestOp, DestData> >::can_resize_rows),
                                       NeedToResizeRowsOfDestination);
            PLLL_INTERNAL_STATIC_CHECK((implementation::MatrixInfo<implementation::expressions::expr<DestOp, DestData> >::cols < 0) || (implementation::MatrixInfo<implementation::expressions::expr<SourceOp, SourceData> >::cols < 0) ||
                                       (static_cast<int>(implementation::MatrixInfo<implementation::expressions::expr<DestOp, DestData> >::cols) == static_cast<int>(implementation::MatrixInfo<implementation::expressions::expr<SourceOp, SourceData> >::cols)) ||
                                       ((static_cast<int>(implementation::MatrixInfo<implementation::expressions::expr<DestOp, DestData> >::cols) != static_cast<int>(implementation::MatrixInfo<implementation::expressions::expr<SourceOp, SourceData> >::cols))
                                          && implementation::MatrixInfo<implementation::expressions::expr<DestOp, DestData> >::can_resize_cols),
                                       NeedToResizeColsOfDestination);
            assert((destination.rows() == source.rows()) || (implementation::MatrixInfo<implementation::expressions::expr<DestOp, DestData> >::can_resize_rows));
            assert((destination.cols() == source.cols()) || (implementation::MatrixInfo<implementation::expressions::expr<DestOp, DestData> >::can_resize_cols));
            if ((destination.rows() != source.rows()) || (destination.cols() != source.cols()))
                implementation::assign_impl_resize(destination, source,
                                                   helper::BoolToType<move && implementation::MatrixInfo<implementation::expressions::expr<SourceOp, SourceData> >::can_move_from>(),
                                                   helper::BoolToType<implementation::MatrixInfo<implementation::expressions::expr<DestOp, DestData> >::can_resize_rows ||
                                                   implementation::MatrixInfo<implementation::expressions::expr<DestOp, DestData> >::can_resize_cols>());
            else
                implementation::assign_with_temporary(destination, source, ittm, helper::BoolToType<implementation::MatrixInfo<implementation::expressions::expr<SourceOp, SourceData> >::use_temporary_on_evaluate>());
        }
        
        template<template<typename SourceData> class SourceOp, typename SourceData, template<typename DestData> class DestOp, typename DestData>
        inline void assign(const implementation::expressions::expr<DestOp, DestData> & destination,
                           const implementation::expressions::expr<SourceOp, SourceData> & source)
        {
            assign(destination, source, helper::BoolToType<false>());
        }
        
        template<typename SourceT, int SourceRows, int SourceCols, typename SourceST, bool SourceMO, template<typename DestData> class DestOp, typename DestData>
        inline void assign(const implementation::expressions::expr<DestOp, DestData> & destination,
                           const base_matrix<SourceT, SourceRows, SourceCols, SourceST, SourceMO> & source)
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("assign(" << getAddress(destination) << " " << getAddress(source) << ") [2]");
            assign(destination, implementation::expressions::make_matrix_expression(source));
        }
        
        template<template<typename SourceData> class SourceOp, typename SourceData, typename DestT, int DestRows, int DestCols, typename DestST, bool DestMO>
        inline void assign(base_matrix<DestT, DestRows, DestCols, DestST, DestMO> & destination,
                           const implementation::expressions::expr<SourceOp, SourceData> & source)
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("assign(" << getAddress(destination) << " " << getAddress(source) << ") [3]");
            assign(implementation::expressions::make_matrix_expression(destination), source);
        }
        
        template<typename SourceT, int SourceRows, int SourceCols, typename SourceST, bool SourceMO,
                 typename DestT, int DestRows, int DestCols, typename DestST, bool DestMO>
        inline void assign(base_matrix<DestT, DestRows, DestCols, DestST, DestMO> & destination,
                           const base_matrix<SourceT, SourceRows, SourceCols, SourceST, SourceMO> & source)
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("assign(" << getAddress(destination) << " " << getAddress(source) << ") [4]");
            assign(implementation::expressions::make_matrix_expression(destination), implementation::expressions::make_matrix_expression(source));
        }
        
#if __cplusplus >= 201103L
        template<typename SourceT, int SourceRows, int SourceCols, typename SourceST, bool SourceMO, template<typename DestData> class DestOp, typename DestData>
        inline void assign(const implementation::expressions::expr<DestOp, DestData> & destination,
                           base_matrix<SourceT, SourceRows, SourceCols, SourceST, SourceMO> && source)
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("assign(" << getAddress(destination) << " " << getAddress(source) << ") [2]");
            assign(destination, implementation::expressions::make_matrix_expression(source), helper::BoolToType<true>());
        }
        
        template<typename SourceT, int SourceRows, int SourceCols, typename SourceST, bool SourceMO,
                 typename DestT, int DestRows, int DestCols, typename DestST, bool DestMO>
        inline void assign(base_matrix<DestT, DestRows, DestCols, DestST, DestMO> & destination,
                           base_matrix<SourceT, SourceRows, SourceCols, SourceST, SourceMO> && source)
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("assign(" << getAddress(destination) << " " << getAddress(source) << ") [4]");
            assign(implementation::expressions::make_matrix_expression(destination), implementation::expressions::make_matrix_expression(source), helper::BoolToType<true>());
        }
#endif
        
        ///@}
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        //// transpose (functional)
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        /**@{
           \name Transpose functions.
        */
        
        template<template<typename SourceData> class SourceOp, typename SourceData, template<typename DestData> class DestOp, typename DestData, bool move>
        inline void transpose(const implementation::expressions::expr<DestOp, DestData> & destination,
                              const implementation::expressions::expr<SourceOp, SourceData> & source, helper::BoolToType<move> ittm)
        {
            if (implementation::expressions::expr<SourceOp, SourceData>::has_direct_access)
                assign(implementation::expressions::expr<implementation::expressions::transpose, implementation::expressions::expr<DestOp, DestData> >(destination), source, ittm);
            else
                assign(destination, implementation::expressions::expr<implementation::expressions::transpose, implementation::expressions::expr<SourceOp, SourceData> >(source), ittm);
        }
        
        template<template<typename SourceData> class SourceOp, typename SourceData, template<typename DestData> class DestOp, typename DestData>
        inline void transpose(const implementation::expressions::expr<DestOp, DestData> & destination,
                              const implementation::expressions::expr<SourceOp, SourceData> & source)
        {
            transpose(destination, source, helper::BoolToType<false>());
        }
        
        template<typename SourceT, int SourceRows, int SourceCols, typename SourceST, bool SourceMO, template<typename DestData> class DestOp, typename DestData>
        inline void transpose(const implementation::expressions::expr<DestOp, DestData> & destination,
                              const base_matrix<SourceT, SourceRows, SourceCols, SourceST, SourceMO> & source)
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("transpose(" << getAddress(destination) << " " << getAddress(source) << ") [2]");
            transpose(destination, implementation::expressions::make_matrix_expression(source));
        }
        
        template<template<typename SourceData> class SourceOp, typename SourceData, typename DestT, int DestRows, int DestCols, typename DestST, bool DestMO>
        inline void transpose(base_matrix<DestT, DestRows, DestCols, DestST, DestMO> & destination,
                              const implementation::expressions::expr<SourceOp, SourceData> & source)
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("transpose(" << getAddress(destination) << " " << getAddress(source) << ") [3]");
            transpose(implementation::expressions::make_matrix_expression(destination), source);
        }
        
        template<typename SourceT, int SourceRows, int SourceCols, typename SourceST, bool SourceMO,
                 typename DestT, int DestRows, int DestCols, typename DestST, bool DestMO>
        inline void transpose(base_matrix<DestT, DestRows, DestCols, DestST, DestMO> & destination,
                              const base_matrix<SourceT, SourceRows, SourceCols, SourceST, SourceMO> & source)
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("transpose(" << getAddress(destination) << " " << getAddress(source) << ") [4]");
            transpose(implementation::expressions::make_matrix_expression(destination), implementation::expressions::make_matrix_expression(source));
        }
        
#if __cplusplus >= 201103L
        template<typename SourceT, int SourceRows, int SourceCols, typename SourceST, bool SourceMO, template<typename DestData> class DestOp, typename DestData>
        inline void transpose(const implementation::expressions::expr<DestOp, DestData> & destination,
                              base_matrix<SourceT, SourceRows, SourceCols, SourceST, SourceMO> && source)
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("transpose(" << getAddress(destination) << " " << getAddress(source) << ") [2]");
            transpose(destination, implementation::expressions::make_matrix_expression(source), helper::BoolToType<true>());
        }
        
        template<typename SourceT, int SourceRows, int SourceCols, typename SourceST, bool SourceMO,
                 typename DestT, int DestRows, int DestCols, typename DestST, bool DestMO>
        inline void transpose(base_matrix<DestT, DestRows, DestCols, DestST, DestMO> & destination,
                              base_matrix<SourceT, SourceRows, SourceCols, SourceST, SourceMO> && source)
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("transpose(" << getAddress(destination) << " " << getAddress(source) << ") [4]");
            transpose(implementation::expressions::make_matrix_expression(destination), implementation::expressions::make_matrix_expression(source), helper::BoolToType<true>());
        }
#endif
        
        ///@}
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        //// (usual) operators
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        /**@{
           \name Operators.
        */
        
        // Operations: base_matrix [op] other

        /**
           \brief Computes the product of `A` with `B`.
           
           \param A A matrix.
           \param B A matrix or scalar.
           \return The product of `A` and `B`.
         */
        template<typename T, int Rows, int Cols, typename ST, typename MT>
        inline typename implementation::BinaryMatrixOperationImpl<base_matrix<T, Rows, Cols, ST, true>, MT>::Mul_ReturnType
        operator * (const base_matrix<T, Rows, Cols, ST, true> & A, const MT & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::BinaryMatrixOperationImpl<base_matrix<T, Rows, Cols, ST, true>, MT>::multiply(A, B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator * (" << getAddress(A) << " " << getAddress(B) << ")");
            return implementation::BinaryMatrixOperationImpl<base_matrix<T, Rows, Cols, ST, true>, MT>::multiply(A, B);
        }
        
        /**
           \brief Computes the product of `A` with `B`.
           
           \param A A scalar.
           \param B A matrix.
           \return The product of `A` and `B`.
         */
        template<typename T, int Rows, int Cols, typename ST>
        inline typename implementation::BinaryMatrixOperationImpl<T, base_matrix<T, Rows, Cols, ST, true> >::Mul_ReturnType
        operator * (const T & A, const base_matrix<T, Rows, Cols, ST, true> & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::BinaryMatrixOperationImpl<T, base_matrix<T, Rows, Cols, ST, true> >::multiply(A, B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator * (" << getAddress(A) << " " << getAddress(B) << ")");
            return implementation::BinaryMatrixOperationImpl<T, base_matrix<T, Rows, Cols, ST, true> >::multiply(A, B);
        }
        
        /**
           \brief Computes the componentwise division of `A` with the scalar `B`.
           
           \param A A matrix.
           \param B A scalar.
           \return The componentwise division of `A` by `B`.
         */
        template<typename T, int Rows, int Cols, typename ST, typename MT>
        inline typename implementation::BinaryMatrixOperationImpl<base_matrix<T, Rows, Cols, ST, true>, MT>::Div_ReturnType
        operator / (const base_matrix<T, Rows, Cols, ST, true> & A, const MT & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::BinaryMatrixOperationImpl<base_matrix<T, Rows, Cols, ST, true>, MT>::divide(A, B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator / (" << getAddress(A) << " " << getAddress(B) << ")");
            return implementation::BinaryMatrixOperationImpl<base_matrix<T, Rows, Cols, ST, true>, MT>::divide(A, B);
        }
        
        /**
           \brief Computes the componentwise modulo of `A` with the scalar `B`.
           
           \param A A matrix.
           \param B A scalar.
           \return The componentwise modulo of `A` by `B`.
         */
        template<typename T, int Rows, int Cols, typename ST, typename MT>
        inline typename implementation::BinaryMatrixOperationImpl<base_matrix<T, Rows, Cols, ST, true>, MT>::Mod_ReturnType
        operator % (const base_matrix<T, Rows, Cols, ST, true> & A, const MT & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::BinaryMatrixOperationImpl<base_matrix<T, Rows, Cols, ST, true>, MT>::modulo(A, B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator % (" << getAddress(A) << " " << getAddress(B) << ")");
            return implementation::BinaryMatrixOperationImpl<base_matrix<T, Rows, Cols, ST, true>, MT>::modulo(A, B);
        }
        
        /**
           \brief Computes the componentwise multiplication of `A` with `B`.
           
           \param A A matrix.
           \param B A matrix.
           \return The componentwise multiplication of `A` and `B`.
         */
        template<typename T, int Rows, int Cols, typename ST, typename MT>
        inline typename implementation::BinaryMatrixOperationImpl<base_matrix<T, Rows, Cols, ST, true>, MT>::CwMul_ReturnType
        componentwise_mul(const base_matrix<T, Rows, Cols, ST, true> & A, const MT & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::BinaryMatrixOperationImpl<base_matrix<T, Rows, Cols, ST, true>, MT>::componentwise_mul(A, B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("componentwise_mul(" << getAddress(A) << " " << getAddress(B) << ")");
            return implementation::BinaryMatrixOperationImpl<base_matrix<T, Rows, Cols, ST, true>, MT>::componentwise_mul(A, B);
        }
        
        /**
           \brief Computes the componentwise division of `A` with `B`.
           
           \param A A matrix.
           \param B A matrix.
           \return The componentwise division of `A` by `B`.
         */
        template<typename T, int Rows, int Cols, typename ST, typename MT>
        inline typename implementation::BinaryMatrixOperationImpl<base_matrix<T, Rows, Cols, ST, true>, MT>::CwDiv_ReturnType
        componentwise_div(const base_matrix<T, Rows, Cols, ST, true> & A, const MT & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::BinaryMatrixOperationImpl<base_matrix<T, Rows, Cols, ST, true>, MT>::componentwise_div(A, B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("componentwise_div(" << getAddress(A) << " " << getAddress(B) << ")");
            return implementation::BinaryMatrixOperationImpl<base_matrix<T, Rows, Cols, ST, true>, MT>::componentwise_div(A, B);
        }
        
        /**
           \brief Computes the componentwise modulo of `A` with `B`.
           
           \param A A matrix.
           \param B A matrix.
           \return The componentwise modulo of `A` by `B`.
         */
        template<typename T, int Rows, int Cols, typename ST, typename MT>
        inline typename implementation::BinaryMatrixOperationImpl<base_matrix<T, Rows, Cols, ST, true>, MT>::CwMod_ReturnType
        componentwise_mod(const base_matrix<T, Rows, Cols, ST, true> & A, const MT & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::BinaryMatrixOperationImpl<base_matrix<T, Rows, Cols, ST, true>, MT>::componentwise_mod(A, B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("componentwise_mod(" << getAddress(A) << " " << getAddress(B) << ")");
            return implementation::BinaryMatrixOperationImpl<base_matrix<T, Rows, Cols, ST, true>, MT>::componentwise_mod(A, B);
        }
        
        /**
           \brief Computes the sum of `A` and `B`.
           
           \param A A matrix.
           \param B A matrix.
           \return The sum of `A` and `B`.
         */
        template<typename T, int Rows, int Cols, typename ST, typename MT>
        inline typename implementation::BinaryMatrixOperationImpl<base_matrix<T, Rows, Cols, ST, true>, MT>::Add_ReturnType
        operator + (const base_matrix<T, Rows, Cols, ST, true> & A, const MT & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::BinaryMatrixOperationImpl<base_matrix<T, Rows, Cols, ST, true>, MT>::add(A, B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator + (" << getAddress(A) << " " << getAddress(B) << ")");
            return implementation::BinaryMatrixOperationImpl<base_matrix<T, Rows, Cols, ST, true>, MT>::add(A, B);
        }
        
        /**
           \brief Computes the difference of `A` and `B`.
           
           \param A A matrix.
           \param B A matrix.
           \return The difference of `A` and `B`.
         */
        template<typename T, int Rows, int Cols, typename ST, typename MT>
        inline typename implementation::BinaryMatrixOperationImpl<base_matrix<T, Rows, Cols, ST, true>, MT>::Sub_ReturnType
        operator - (const base_matrix<T, Rows, Cols, ST, true> & A, const MT & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::BinaryMatrixOperationImpl<base_matrix<T, Rows, Cols, ST, true>, MT>::sub(A, B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator - (" << getAddress(A) << " " << getAddress(B) << ")");
            return implementation::BinaryMatrixOperationImpl<base_matrix<T, Rows, Cols, ST, true>, MT>::sub(A, B);
        }
        
        /**
           \brief Computes the negation of `A`.
           
           \param A A matrix.
           \return The componentwise negation of `A`.
         */
        template<typename T, int Rows, int Cols, typename ST>
        inline typename implementation::UnaryMatrixOperationImpl<base_matrix<T, Rows, Cols, ST, true> >::Neg_ReturnType
        operator - (const base_matrix<T, Rows, Cols, ST, true> & A)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::UnaryMatrixOperationImpl<base_matrix<T, Rows, Cols, ST, true> >::negate(A)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator - (" << getAddress(A) << ")");
            return implementation::UnaryMatrixOperationImpl<base_matrix<T, Rows, Cols, ST, true> >::negate(A);
        }
        
        // Operations: implementation::expressions::expr [op] other
        
        template<template<typename DataType> class Operator, typename Data, typename MT>
        inline typename implementation::BinaryMatrixOperationImpl<implementation::expressions::expr<Operator, Data>, MT>::Mul_ReturnType
        operator * (const implementation::expressions::expr<Operator, Data> & A, const MT & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::BinaryMatrixOperationImpl<implementation::expressions::expr<Operator, Data>, MT>::multiply(A, B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator * (" << getAddress(A) << " " << getAddress(B) << ")");
            return implementation::BinaryMatrixOperationImpl<implementation::expressions::expr<Operator, Data>, MT>::multiply(A, B);
        }
        
        template<template<typename DataType> class Operator, typename Data>
        inline typename implementation::BinaryMatrixOperationImpl<typename implementation::MatrixInfo<implementation::expressions::expr<Operator, Data> >::Type,
                                                 implementation::expressions::expr<Operator, Data> >::Mul_ReturnType
        operator * (const typename implementation::MatrixInfo<implementation::expressions::expr<Operator, Data> >::Type & A,
                    const implementation::expressions::expr<Operator, Data> & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::BinaryMatrixOperationImpl<
                                                               typename implementation::MatrixInfo<implementation::expressions::expr<Operator, Data> >::Type,
                                                               implementation::expressions::expr<Operator, Data> >::multiply(A, B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator * (" << getAddress(A) << " " << getAddress(B) << ")");
            return implementation::BinaryMatrixOperationImpl<typename implementation::MatrixInfo<implementation::expressions::expr<Operator, Data> >::Type,
                                                             implementation::expressions::expr<Operator, Data> >::multiply(A, B);
        }
        
        template<template<typename DataType> class Operator, typename Data, typename MT>
        inline typename implementation::BinaryMatrixOperationImpl<implementation::expressions::expr<Operator, Data>, MT>::Div_ReturnType
        operator / (const implementation::expressions::expr<Operator, Data> & A, const MT & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::BinaryMatrixOperationImpl<implementation::expressions::expr<Operator, Data>, MT>::divide(A, B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator / (" << getAddress(A) << " " << getAddress(B) << ")");
            return implementation::BinaryMatrixOperationImpl<implementation::expressions::expr<Operator, Data>, MT>::divide(A, B);
        }
        
        template<template<typename DataType> class Operator, typename Data, typename MT>
        inline typename implementation::BinaryMatrixOperationImpl<implementation::expressions::expr<Operator, Data>, MT>::Mod_ReturnType
        operator % (const implementation::expressions::expr<Operator, Data> & A, const MT & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::BinaryMatrixOperationImpl<implementation::expressions::expr<Operator, Data>, MT>::modulo(A, B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator % (" << getAddress(A) << " " << getAddress(B) << ")");
            return implementation::BinaryMatrixOperationImpl<implementation::expressions::expr<Operator, Data>, MT>::modulo(A, B);
        }
        
        template<template<typename DataType> class Operator, typename Data, typename MT>
        inline typename implementation::BinaryMatrixOperationImpl<implementation::expressions::expr<Operator, Data>, MT>::CwMul_ReturnType
        componentwise_mul(const implementation::expressions::expr<Operator, Data> & A, const MT & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::BinaryMatrixOperationImpl<implementation::expressions::expr<Operator, Data>, MT>::componentwise_mul(A, B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("componentwise_mul(" << getAddress(A) << " " << getAddress(B) << ")");
            return implementation::BinaryMatrixOperationImpl<implementation::expressions::expr<Operator, Data>, MT>::componentwise_mul(A, B);
        }
        
        template<template<typename DataType> class Operator, typename Data, typename MT>
        inline typename implementation::BinaryMatrixOperationImpl<implementation::expressions::expr<Operator, Data>, MT>::CwDiv_ReturnType
        componentwise_div(const implementation::expressions::expr<Operator, Data> & A, const MT & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::BinaryMatrixOperationImpl<implementation::expressions::expr<Operator, Data>, MT>::componentwise_div(A, B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("componentwise_div(" << getAddress(A) << " " << getAddress(B) << ")");
            return implementation::BinaryMatrixOperationImpl<implementation::expressions::expr<Operator, Data>, MT>::componentwise_div(A, B);
        }
        
        template<template<typename DataType> class Operator, typename Data, typename MT>
        inline typename implementation::BinaryMatrixOperationImpl<implementation::expressions::expr<Operator, Data>, MT>::CwMod_ReturnType
        componentwise_mod(const implementation::expressions::expr<Operator, Data> & A, const MT & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::BinaryMatrixOperationImpl<implementation::expressions::expr<Operator, Data>, MT>::componentwise_mod(A, B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("componentwise_mod(" << getAddress(A) << " " << getAddress(B) << ")");
            return implementation::BinaryMatrixOperationImpl<implementation::expressions::expr<Operator, Data>, MT>::componentwise_mod(A, B);
        }
        
        template<template<typename DataType> class Operator, typename Data, typename MT>
        inline typename implementation::BinaryMatrixOperationImpl<implementation::expressions::expr<Operator, Data>, MT>::Add_ReturnType
        operator + (const implementation::expressions::expr<Operator, Data> & A, const MT & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::BinaryMatrixOperationImpl<implementation::expressions::expr<Operator, Data>, MT>::add(A, B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator + (" << getAddress(A) << " " << getAddress(B) << ")");
            return implementation::BinaryMatrixOperationImpl<implementation::expressions::expr<Operator, Data>, MT>::add(A, B);
        }
        
        template<template<typename DataType> class Operator, typename Data, typename MT>
        inline typename implementation::BinaryMatrixOperationImpl<implementation::expressions::expr<Operator, Data>, MT>::Mul_ReturnType
        operator - (const implementation::expressions::expr<Operator, Data> & A, const MT & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::BinaryMatrixOperationImpl<implementation::expressions::expr<Operator, Data>, MT>::sub(A, B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator - (" << getAddress(A) << " " << getAddress(B) << ")");
            return implementation::BinaryMatrixOperationImpl<implementation::expressions::expr<Operator, Data>, MT>::sub(A, B);
        }
        
        template<template<typename DataType> class Operator, typename Data>
        inline typename implementation::UnaryMatrixOperationImpl<implementation::expressions::expr<Operator, Data> >::Neg_ReturnType
        operator - (const implementation::expressions::expr<Operator, Data> & A)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::UnaryMatrixOperationImpl<implementation::expressions::expr<Operator, Data> >::negate(A)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator - (" << getAddress(A) << ")");
            return implementation::UnaryMatrixOperationImpl<implementation::expressions::expr<Operator, Data> >::negate(A);
        }
        
        ///@}
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        //// functional versions of operators, i.e. something along "void op(result, opA, ...)"
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        /**@{
           \name Functional versions of operators.
        */
        
        // Operations: base_matrix [op] other

        /**
           \brief Computes the product of `A` and `B` and stores the result in `result`.
           
           \param result The variable to store the result in.
           \param A The first operand.
           \param B The second operand.
         */
        template<typename T, int Rows, int Cols, typename ST, typename MT1, typename MT2>
        void mul(base_matrix<T, Rows, Cols, ST, true> & result, const MT1 & A, const MT2 & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(assign(result, implementation::BinaryMatrixOperationImpl<MT1, MT2>::multiply(A, B))))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("mul(" << getAddress(result) << " " << getAddress(A) << " " << getAddress(B) << ")");
            assign(result, implementation::BinaryMatrixOperationImpl<MT1, MT2>::multiply(A, B));
        }
        
        /**
           \brief Computes the componentwise division of `A` by `B` and stores the result in
                  `result`.
           
           \param result The variable to store the result in.
           \param A The first operand. Must be a matrix.
           \param B The second operand. Must be a scalar.
         */
        template<typename T, int Rows, int Cols, typename ST, typename MT1, typename MT2>
        void div(base_matrix<T, Rows, Cols, ST, true> & result, const MT1 & A, const MT2 & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(assign(result, implementation::BinaryMatrixOperationImpl<MT1, MT2>::divide(A, B))))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("div(" << getAddress(result) << " " << getAddress(A) << " " << getAddress(B) << ")");
            assign(result, implementation::BinaryMatrixOperationImpl<MT1, MT2>::divide(A, B));
        }
        
        /**
           \brief Computes the componentwise modulo of `A` by `B` and stores the result in `result`.
           
           \param result The variable to store the result in.
           \param A The first operand. Must be a matrix.
           \param B The second operand. Must be a scalar.
         */
        template<typename T, int Rows, int Cols, typename ST, typename MT1, typename MT2>
        void mod(base_matrix<T, Rows, Cols, ST, true> & result, const MT1 & A, const MT2 & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(assign(result, implementation::BinaryMatrixOperationImpl<MT1, MT2>::modulo(A, B))))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("mod(" << getAddress(result) << " " << getAddress(A) << " " << getAddress(B) << ")");
            assign(result, implementation::BinaryMatrixOperationImpl<MT1, MT2>::modulo(A, B));
        }
        
        /**
           \brief Computes the componentwise product of `A` and `B` and stores the result in
                  `result`.
           
           \param result The variable to store the result in.
           \param A The first operand. Must be a matrix.
           \param B The second operand. Must be a matrix.
         */
        template<typename T, int Rows, int Cols, typename ST, typename MT1, typename MT2>
        void componentwise_mul(base_matrix<T, Rows, Cols, ST, true> & result, const MT1 & A, const MT2 & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(assign(result, implementation::BinaryMatrixOperationImpl<MT1, MT2>::componentwise_mul(A, B))))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("componentwise_mul(" << getAddress(result) << " " << getAddress(A) << " " << getAddress(B) << ")");
            assign(result, implementation::BinaryMatrixOperationImpl<MT1, MT2>::componentwise_mul(A, B));
        }
        
        /**
           \brief Computes the componentwise division of `A` by `B` and stores the result in
                  `result`.
           
           \param result The variable to store the result in.
           \param A The first operand. Must be a matrix.
           \param B The second operand. Must be a matrix.
         */
        template<typename T, int Rows, int Cols, typename ST, typename MT1, typename MT2>
        void componentwise_div(base_matrix<T, Rows, Cols, ST, true> & result, const MT1 & A, const MT2 & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(assign(result, implementation::BinaryMatrixOperationImpl<MT1, MT2>::componentwise_div(A, B))))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("componentwise_div(" << getAddress(result) << " " << getAddress(A) << " " << getAddress(B) << ")");
            assign(result, implementation::BinaryMatrixOperationImpl<MT1, MT2>::componentwise_div(A, B));
        }
        
        /**
           \brief Computes the componentwise modulo of `A` by `B` and stores the result in `result`.
           
           \param result The variable to store the result in.
           \param A The first operand. Must be a matrix.
           \param B The second operand. Must be a matrix.
         */
        template<typename T, int Rows, int Cols, typename ST, typename MT1, typename MT2>
        void componentwise_mod(base_matrix<T, Rows, Cols, ST, true> & result, const MT1 & A, const MT2 & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(assign(result, implementation::BinaryMatrixOperationImpl<MT1, MT2>::componentwise_mod(A, B))))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("componentwise_mod(" << getAddress(result) << " " << getAddress(A) << " " << getAddress(B) << ")");
            assign(result, implementation::BinaryMatrixOperationImpl<MT1, MT2>::componentwise_mod(A, B));
        }
        
        /**
           \brief Computes the sum of `A` and `B` and stores the result in `result`.
           
           \param result The variable to store the result in.
           \param A The first operand. Must be a matrix.
           \param B The second operand. Must be a matrix.
         */
        template<typename T, int Rows, int Cols, typename ST, typename MT1, typename MT2>
        void add(base_matrix<T, Rows, Cols, ST, true> & result, const MT1 & A, const MT2 & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(assign(result, implementation::BinaryMatrixOperationImpl<MT1, MT2>::add(A, B))))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("add(" << getAddress(result) << " " << getAddress(A) << " " << getAddress(B) << ")");
            assign(result, implementation::BinaryMatrixOperationImpl<MT1, MT2>::add(A, B));
        }
        
        /**
           \brief Computes the difference of `A` and `B` and stores the result in `result`.
           
           \param result The variable to store the result in.
           \param A The first operand. Must be a matrix.
           \param B The second operand. Must be a matrix.
         */
        template<typename T, int Rows, int Cols, typename ST, typename MT1, typename MT2>
        void sub(base_matrix<T, Rows, Cols, ST, true> & result, const MT1 & A, const MT2 & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(assign(result, implementation::BinaryMatrixOperationImpl<MT1, MT2>::sub(A, B))))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("sub(" << getAddress(result) << " " << getAddress(A) << " " << getAddress(B) << ")");
            assign(result, implementation::BinaryMatrixOperationImpl<MT1, MT2>::sub(A, B));
        }
        
        /**
           \brief Computes the negation of `A` and stores the result in `result`.
           
           \param result The variable to store the result in.
           \param A The operand. Must be a matrix.
         */
        template<typename T, int Rows, int Cols, typename ST, typename MT>
        void neg(base_matrix<T, Rows, Cols, ST, true> & result, const MT & A)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(assign(result, implementation::UnaryMatrixOperationImpl<MT>::negate(A))))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("neg(" << getAddress(result) << " " << getAddress(A) << ")");
            assign(result, implementation::UnaryMatrixOperationImpl<MT>::negate(A));
        }
        
        // Operations: implementation::expressions::expr [op] other
        
        template<template<typename DataType> class Operator, typename Data, typename MT1, typename MT2>
        void mul(const implementation::expressions::expr<Operator, Data> & result, const MT1 & A, const MT2 & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(assign(result, implementation::BinaryMatrixOperationImpl<MT1, MT2>::multiply(A, B))))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("mul(" << getAddress(result) << " " << getAddress(A) << " " << getAddress(B) << ")");
            assign(result, implementation::BinaryMatrixOperationImpl<MT1, MT2>::multiply(A, B));
        }
        
        template<template<typename DataType> class Operator, typename Data, typename MT1, typename MT2>
        void div(const implementation::expressions::expr<Operator, Data> & result, const MT1 & A, const MT2 & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(assign(result, implementation::BinaryMatrixOperationImpl<MT1, MT2>::divide(A, B))))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("div(" << getAddress(result) << " " << getAddress(A) << " " << getAddress(B) << ")");
            assign(result, implementation::BinaryMatrixOperationImpl<MT1, MT2>::divide(A, B));
        }
        
        template<template<typename DataType> class Operator, typename Data, typename MT1, typename MT2>
        void mod(const implementation::expressions::expr<Operator, Data> & result, const MT1 & A, const MT2 & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(assign(result, implementation::BinaryMatrixOperationImpl<MT1, MT2>::modulo(A, B))))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("mod(" << getAddress(result) << " " << getAddress(A) << " " << getAddress(B) << ")");
            assign(result, implementation::BinaryMatrixOperationImpl<MT1, MT2>::modulo(A, B));
        }
        
        template<template<typename DataType> class Operator, typename Data, typename MT1, typename MT2>
        void componentwise_mul(const implementation::expressions::expr<Operator, Data> & result, const MT1 & A, const MT2 & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(assign(result, implementation::BinaryMatrixOperationImpl<MT1, MT2>::componentwise_mul(A, B))))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("componentwise_mul(" << getAddress(result) << " " << getAddress(A) << " " << getAddress(B) << ")");
            assign(result, implementation::BinaryMatrixOperationImpl<MT1, MT2>::componentwise_mul(A, B));
        }
        
        template<template<typename DataType> class Operator, typename Data, typename MT1, typename MT2>
        void componentwise_div(const implementation::expressions::expr<Operator, Data> & result, const MT1 & A, const MT2 & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(assign(result, implementation::BinaryMatrixOperationImpl<MT1, MT2>::componentwise_div(A, B))))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("componentwise_div(" << getAddress(result) << " " << getAddress(A) << " " << getAddress(B) << ")");
            assign(result, implementation::BinaryMatrixOperationImpl<MT1, MT2>::componentwise_div(A, B));
        }
        
        template<template<typename DataType> class Operator, typename Data, typename MT1, typename MT2>
        void componentwise_mod(const implementation::expressions::expr<Operator, Data> & result, const MT1 & A, const MT2 & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(assign(result, implementation::BinaryMatrixOperationImpl<MT1, MT2>::componentwise_mod(A, B))))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("componentwise_mod(" << getAddress(result) << " " << getAddress(A) << " " << getAddress(B) << ")");
            assign(result, implementation::BinaryMatrixOperationImpl<MT1, MT2>::componentwise_mod(A, B));
        }
        
        template<template<typename DataType> class Operator, typename Data, typename MT1, typename MT2>
        void add(const implementation::expressions::expr<Operator, Data> & result, const MT1 & A, const MT2 & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(assign(result, implementation::BinaryMatrixOperationImpl<MT1, MT2>::add(A, B))))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("add(" << getAddress(result) << " " << getAddress(A) << " " << getAddress(B) << ")");
            assign(result, implementation::BinaryMatrixOperationImpl<MT1, MT2>::add(A, B));
        }
        
        template<template<typename DataType> class Operator, typename Data, typename MT1, typename MT2>
        void sub(const implementation::expressions::expr<Operator, Data> & result, const MT1 & A, const MT2 & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(assign(result, implementation::BinaryMatrixOperationImpl<MT1, MT2>::sub(A, B))))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("sub(" << getAddress(result) << " " << getAddress(A) << " " << getAddress(B) << ")");
            assign(result, implementation::BinaryMatrixOperationImpl<MT1, MT2>::sub(A, B));
        }
        
        template<template<typename DataType> class Operator, typename Data, typename MT>
        void neg(const implementation::expressions::expr<Operator, Data> & result, const MT & A)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(assign(result, implementation::UnaryMatrixOperationImpl<MT>::negate(A))))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("neg(" << getAddress(result) << " " << getAddress(A) << ")");
            assign(result, implementation::UnaryMatrixOperationImpl<MT>::negate(A));
        }
        
        ///@}
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        //// (usual) assignment-operation operators
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        /**@{
           \name Assignment-operation operators.
        */
        
        // Operations: base_matrix [op]= other
        
        /**
           \brief Adds `B` to the current matrix.
           
           \param A The current matrix.
           \param B The matrix to add to the current matrix.
           \return A reference to the current matrix.
         */
        template<typename T, int Rows, int Cols, typename ST, typename MT>
        inline base_matrix<T, Rows, Cols, ST, true> & operator += (base_matrix<T, Rows, Cols, ST, true> & A, const MT & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::expressions::make_matrix_expression(A) += B))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator += (" << getAddress(A) << " " << getAddress(B) << ")");
            implementation::expressions::make_matrix_expression(A) += B;
            return A;
        }
        
        /**
           \brief Subtracts `B` from the current matrix.
           
           \param A The current matrix.
           \param B The matrix to subtract from the current matrix.
           \return A reference to the current matrix.
         */
        template<typename T, int Rows, int Cols, typename ST, typename MT>
        inline base_matrix<T, Rows, Cols, ST, true> & operator -= (base_matrix<T, Rows, Cols, ST, true> & A, const MT & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::expressions::make_matrix_expression(A) -= B))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator -= (" << getAddress(A) << " " << getAddress(B) << ")");
            implementation::expressions::make_matrix_expression(A) -= B;
            return A;
        }
        
        /**
           \brief Multiplies the current matrix with `B` and stores the result in the current matrix.
           
           \param A The current matrix.
           \param B The matrix or scalar to multiply with.
           \return A reference to the current matrix.
         */
        template<typename T, int Rows, int Cols, typename ST, typename MT>
        inline base_matrix<T, Rows, Cols, ST, true> & operator *= (base_matrix<T, Rows, Cols, ST, true> & A, const MT & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::expressions::make_matrix_expression(A) *= B))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator *= (" << getAddress(A) << " " << getAddress(B) << ")");
            implementation::expressions::make_matrix_expression(A) *= B;
            return A;
        }
        
        /**
           \brief Divides the current matrix by the scalar `B` and stores the result in the current matrix.
           
           \param A The current matrix.
           \param B The scalar to divide by.
           \return A reference to the current matrix.
         */
        template<typename T, int Rows, int Cols, typename ST, typename MT>
        inline base_matrix<T, Rows, Cols, ST, true> & operator /= (base_matrix<T, Rows, Cols, ST, true> & A, const MT & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::expressions::make_matrix_expression(A) /= B))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator /= (" << getAddress(A) << " " << getAddress(B) << ")");
            implementation::expressions::make_matrix_expression(A) /= B;
            return A;
        }
        
        /**
           \brief Takes the current matrix modulo the scalar `B` and stores the result in the
                  current matrix.
           
           \param A The current matrix.
           \param B The scalar to mod by.
           \return A reference to the current matrix.
         */
        template<typename T, int Rows, int Cols, typename ST, typename MT>
        inline base_matrix<T, Rows, Cols, ST, true> & operator %= (base_matrix<T, Rows, Cols, ST, true> & A, const MT & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::expressions::make_matrix_expression(A) %= B))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator %= (" << getAddress(A) << " " << getAddress(B) << ")");
            implementation::expressions::make_matrix_expression(A) %= B;
            return A;
        }
        
        // Operations: expr [op]= other
        
        template<template<typename DataType> class Op, typename Data, typename MT>
        inline const implementation::expressions::expr<Op, Data> & operator += (const implementation::expressions::expr<Op, Data> & A, const MT & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::AssignmentMatrixOperationImpl<implementation::expressions::expr<Op, Data>, MT>::add_assign(A, B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator += (" << getAddress(A) << " " << getAddress(B) << ")");
            return implementation::AssignmentMatrixOperationImpl<implementation::expressions::expr<Op, Data>, MT>::add_assign(A, B);
        }
        
        template<template<typename DataType> class Op, typename Data, typename MT>
        inline const implementation::expressions::expr<Op, Data> & operator -= (const implementation::expressions::expr<Op, Data> & A, const MT & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::AssignmentMatrixOperationImpl<implementation::expressions::expr<Op, Data>, MT>::sub_assign(A, B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator -= (" << getAddress(A) << " " << getAddress(B) << ")");
            return implementation::AssignmentMatrixOperationImpl<implementation::expressions::expr<Op, Data>, MT>::sub_assign(A, B);
        }
        
        template<template<typename DataType> class Op, typename Data, typename MT>
        inline const implementation::expressions::expr<Op, Data> & operator *= (const implementation::expressions::expr<Op, Data> & A, const MT & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::AssignmentMatrixOperationImpl<implementation::expressions::expr<Op, Data>, MT>::mul_assign(A, B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator *= (" << getAddress(A) << " " << getAddress(B) << ")");
            return implementation::AssignmentMatrixOperationImpl<implementation::expressions::expr<Op, Data>, MT>::mul_assign(A, B);
        }
        
        template<template<typename DataType> class Op, typename Data, typename MT>
        inline const implementation::expressions::expr<Op, Data> & operator /= (const implementation::expressions::expr<Op, Data> & A, const MT & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::AssignmentMatrixOperationImpl<implementation::expressions::expr<Op, Data>, MT>::div_assign(A, B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator /= (" << getAddress(A) << " " << getAddress(B) << ")");
            return implementation::AssignmentMatrixOperationImpl<implementation::expressions::expr<Op, Data>, MT>::div_assign(A, B);
        }
        
        template<template<typename DataType> class Op, typename Data, typename MT>
        inline const implementation::expressions::expr<Op, Data> & operator %= (const implementation::expressions::expr<Op, Data> & A, const MT & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::AssignmentMatrixOperationImpl<implementation::expressions::expr<Op, Data>, MT>::mod_assign(A, B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator %= (" << getAddress(A) << " " << getAddress(B) << ")");
            return implementation::AssignmentMatrixOperationImpl<implementation::expressions::expr<Op, Data>, MT>::mod_assign(A, B);
        }
        
        ///@}
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        //// other operations
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        /**@{
           \name Other operations.
        */
        
        // sum: sums all elements in the matrix
        
        namespace implementation
        {
            template<typename MatrixType>
            class SumImpl
            {
            public:
                typedef typename MatrixInfo<MatrixType>::Type RetValue;
            
                template<typename ReturnType>
                inline static void sum(ReturnType & result, const MatrixType & A)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(setZero(result)) &&
                                                              noexcept(A.enumerate()) &&
                                                              noexcept(helper::make_type_lvalue<typename MatrixType::ConstEnumerator>().has_current()) &&
                                                              noexcept(helper::make_type_lvalue<typename MatrixType::ConstEnumerator>().next()) &&
                                                              noexcept(helper::make_type_lvalue<typename MatrixType::ConstEnumerator>().get_current(result)) &&
                                                              noexcept(result += helper::make_type_lvalue<typename MatrixType::ConstEnumerator>().current()))
                {
                    PLLL_DEBUG_OUTPUT_MESSAGE("SumImpl::sum(" << getAddress(result) << " " << getAddress(A) << ")");
                    PLLL_INTERNAL_STATIC_CHECK(MatrixInfo<MatrixType>::is_matrix && MatrixInfo<MatrixType>::is_math_object, RequiresMathMatrix);
                    if ((A.rows() == 0) || (A.cols() == 0))
                        setZero(result);
                    else
                    {
                        typename MatrixType::ConstEnumerator e = A.enumerate();
                        e.get_current(result);
                        for (e.next(); e.has_current(); e.next())
                            result += e.current();
                    }
                }
            
                inline static RetValue sum(const MatrixType & A)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(setZero(helper::make_type_lvalue<RetValue>())) &&
                                                              noexcept(A.enumerate()) && noexcept(RetValue()) &&
                                                              noexcept(helper::make_type_lvalue<typename MatrixType::ConstEnumerator>().has_current()) &&
                                                              noexcept(helper::make_type_lvalue<typename MatrixType::ConstEnumerator>().next()) &&
                                                              noexcept(helper::make_type_lvalue<RetValue>() = helper::make_type_lvalue<typename MatrixType::ConstEnumerator>().current()) &&
                                                              noexcept(helper::make_type_lvalue<RetValue>() += helper::make_type_lvalue<typename MatrixType::ConstEnumerator>().current()))
                {
                    PLLL_DEBUG_OUTPUT_MESSAGE("SumImpl::sum(" << getAddress(A) << ")");
                    PLLL_INTERNAL_STATIC_CHECK(MatrixInfo<MatrixType>::is_matrix && MatrixInfo<MatrixType>::is_math_object, RequiresMathMatrix);
                    size_type rr = A.rows(), cc = A.cols();
                    if ((rr == 0) || (cc == 0))
                    {
                        RetValue result;
                        setZero(result);
                        return result;
                    }
                    else
                    {
                        typename MatrixType::ConstEnumerator e = A.enumerate();
                        RetValue result = e.current();
                        for (e.next(); e.has_current(); e.next())
                            result += e.current();
                        return result;
                    }
                }
            };
        }
        
        /**
           \brief Computes the sum of all entries of `A` and stores the result in `result`.
           
           \param result A scalar where to compute the result in.
           \param A The matrix whose entries to sum up.
         */
        template<typename RetType, typename T, int Rows, int Cols, typename ST>
        void sum(RetType & result, const base_matrix<T, Rows, Cols, ST, true> & A)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::SumImpl<base_matrix<T, Rows, Cols, ST, true> >::sum(result, A)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("sum(" << getAddress(result) << " " << getAddress(A) << ")");
            implementation::SumImpl<base_matrix<T, Rows, Cols, ST, true> >::sum(result, A);
        }
        
        /**
           \brief Computes the sum of all entries of `A` and returns the result.
           
           \return A scalar whose value is the sum of all entries of `A`.
           \param A The matrix whose entries to sum up.
         */
        template<typename T, int Rows, int Cols, typename ST>
        typename implementation::SumImpl<base_matrix<T, Rows, Cols, ST, true> >::RetValue sum(const base_matrix<T, Rows, Cols, ST, true> & A)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::SumImpl<base_matrix<T, Rows, Cols, ST, true> >::sum(A)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("sum(" << getAddress(A) << ")");
            return implementation::SumImpl<base_matrix<T, Rows, Cols, ST, true> >::sum(A);
        }
        
        template<typename RetType, template<typename DataType> class Operator, typename Data>
        void sum(RetType & result, const implementation::expressions::expr<Operator, Data> & A)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::SumImpl<implementation::expressions::expr<Operator, Data> >::sum(result, A)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("sum(" << getAddress(result) << " " << getAddress(A) << ")");
            implementation::SumImpl<implementation::expressions::expr<Operator, Data> >::sum(result, A);
        }
        
        template<template<typename DataType> class Operator, typename Data>
        typename implementation::SumImpl<implementation::expressions::expr<Operator, Data> >::RetValue sum(const implementation::expressions::expr<Operator, Data> & A)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::SumImpl<implementation::expressions::expr<Operator, Data> >::sum(A)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("sum(" << getAddress(A) << ")");
            return implementation::SumImpl<implementation::expressions::expr<Operator, Data> >::sum(A);
        }
        
        // dot: computes the componentwise product of two matrices and sums over the result

        namespace implementation
        {
            template<typename MatrixType, typename MatrixType2>
            class DotImpl
            {
            public:
                typedef typename MatrixInfo<MatrixType>::Type RetValue;
            
                template<typename ReturnType>
                inline static void dot(ReturnType & result, const MatrixType & A, const MatrixType2 & B)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(noexcept(setZero(result))) &&
                                                              noexcept(A.enumerate()) && noexcept(B.enumerate()) &&
                                                              noexcept(helper::make_type_lvalue<typename MatrixType::ConstEnumerator>().has_current()) &&
                                                              noexcept(helper::make_type_lvalue<typename MatrixType::ConstEnumerator>().next()) &&
                                                              noexcept(helper::make_type_lvalue<typename MatrixType2::ConstEnumerator>().next()) &&
                                                              noexcept(result = helper::make_type_lvalue<typename MatrixType::ConstEnumerator>().current() *
                                                                       helper::make_type_lvalue<typename MatrixType2::ConstEnumerator>().current()) &&
                                                              noexcept(result += helper::make_type_lvalue<typename MatrixType::ConstEnumerator>().current() *
                                                                       helper::make_type_lvalue<typename MatrixType2::ConstEnumerator>().current()))
                {
                    PLLL_DEBUG_OUTPUT_MESSAGE("DotImpl::dot(" << getAddress(result) << " " << getAddress(A) << " " << getAddress(B) << ")");
                    PLLL_INTERNAL_STATIC_CHECK(MatrixInfo<MatrixType>::is_matrix && MatrixInfo<MatrixType>::is_math_object, RequiresMathMatrix);
                    PLLL_INTERNAL_STATIC_CHECK(MatrixInfo<MatrixType2>::is_matrix && MatrixInfo<MatrixType2>::is_math_object, RequiresMathMatrix);
                    PLLL_INTERNAL_STATIC_CHECK((MatrixInfo<MatrixType>::rows < 0) || (MatrixInfo<MatrixType2>::rows < 0) ||
                                               (static_cast<int>(MatrixInfo<MatrixType>::rows) == static_cast<int>(MatrixInfo<MatrixType2>::rows)), FormatsDoNotMatch);
                    PLLL_INTERNAL_STATIC_CHECK((MatrixInfo<MatrixType>::cols < 0) || (MatrixInfo<MatrixType2>::cols < 0) ||
                                               (static_cast<int>(MatrixInfo<MatrixType>::cols) == static_cast<int>(MatrixInfo<MatrixType2>::cols)), FormatsDoNotMatch);
                    size_type rr = A.rows(), cc = A.cols();
                    if ((rr == 0) || (cc == 0))
                        setZero(result);
                    else
                    {
                        typename MatrixType::ConstEnumerator e = A.enumerate();
                        typename MatrixType2::ConstEnumerator e2 = B.enumerate();
                        result = e.current() * e2.current();
                        for (e.next(), e2.next(); e.has_current(); e.next(), e2.next())
                            result += e.current() * e2.current();
                    }
                }
            
                inline static RetValue dot(const MatrixType & A, const MatrixType2 & B)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(setZero(helper::make_type_lvalue<RetValue>())) &&
                                                              noexcept(A.enumerate()) && noexcept(B.enumerate()) && noexcept(RetValue()) &&
                                                              noexcept(helper::make_type_lvalue<typename MatrixType::ConstEnumerator>().has_current()) &&
                                                              noexcept(helper::make_type_lvalue<typename MatrixType::ConstEnumerator>().next()) &&
                                                              noexcept(helper::make_type_lvalue<typename MatrixType2::ConstEnumerator>().next()) &&
                                                              noexcept(helper::make_type_lvalue<RetValue>() =
                                                                       helper::make_type_lvalue<typename MatrixType::ConstEnumerator>().current() *
                                                                       helper::make_type_lvalue<typename MatrixType2::ConstEnumerator>().current()) &&
                                                              noexcept(helper::make_type_lvalue<RetValue>() +=
                                                                       helper::make_type_lvalue<typename MatrixType::ConstEnumerator>().current() *
                                                                       helper::make_type_lvalue<typename MatrixType2::ConstEnumerator>().current()))
                {
                    PLLL_DEBUG_OUTPUT_MESSAGE("DotImpl::dot(" << getAddress(A) << " " << getAddress(B) << ")");
                    PLLL_INTERNAL_STATIC_CHECK(MatrixInfo<MatrixType>::is_matrix && MatrixInfo<MatrixType>::is_math_object, RequiresMathMatrix);
                    PLLL_INTERNAL_STATIC_CHECK(MatrixInfo<MatrixType2>::is_matrix && MatrixInfo<MatrixType2>::is_math_object, RequiresMathMatrix);
                    PLLL_INTERNAL_STATIC_CHECK((MatrixInfo<MatrixType>::rows < 0) || (MatrixInfo<MatrixType2>::rows < 0) ||
                                               (static_cast<int>(MatrixInfo<MatrixType>::rows) == static_cast<int>(MatrixInfo<MatrixType2>::rows)), FormatsDoNotMatch);
                    PLLL_INTERNAL_STATIC_CHECK((MatrixInfo<MatrixType>::cols < 0) || (MatrixInfo<MatrixType2>::cols < 0) ||
                                               (static_cast<int>(MatrixInfo<MatrixType>::cols) == static_cast<int>(MatrixInfo<MatrixType2>::cols)), FormatsDoNotMatch);
                    size_type rr = A.rows(), cc = A.cols();
                    if ((rr == 0) || (cc == 0))
                    {
                        RetValue result;
                        setZero(result);
                        return result;
                    }
                    else
                    {
                        typename MatrixType::ConstEnumerator e = A.enumerate();
                        typename MatrixType2::ConstEnumerator e2 = B.enumerate();
                        RetValue result = e.current() * e2.current();
                        for (e.next(), e2.next(); e.has_current(); e.next(), e2.next())
                            result += e.current() * e2.current();
                        return result;
                    }
                }
            };
        }
        
        /**
           \brief Computes the dot product of `A` and `B` and stores the result in `result`.
           
           What this does is componentwise multiplying `A` and `B` and then calling `sum()` of the
           resulting matrix.
           
           \param result A scalar where to compute the result in.
           \param A The first operand.
           \param B The second operand.
         */
        template<typename RetType, typename T, int Rows, int Cols, typename ST, typename MT>
        void dot(RetType & result, const base_matrix<T, Rows, Cols, ST, true> & A, const MT & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::DotImpl<base_matrix<T, Rows, Cols, ST, true>, MT>::dot(result, A, B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("dot(" << getAddress(result) << " " << getAddress(A) << " " << getAddress(B) << ")");
            implementation::DotImpl<base_matrix<T, Rows, Cols, ST, true>, MT>::dot(result, A, B);
        }
        
        /**
           \brief Computes the dot product of `A` and `B` and returns the result.
           
           What this does is componentwise multiplying `A` and `B` and then calling `sum()` of the
           resulting matrix.
           
           \return A scalar whose value is the dot product of `A` and `B`.
           \param A The first operand.
           \param B The second operand.
         */
        template<typename T, int Rows, int Cols, typename ST, typename MT>
        typename implementation::DotImpl<base_matrix<T, Rows, Cols, ST, true>, MT>::RetValue dot(const base_matrix<T, Rows, Cols, ST, true> & A, const MT & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::DotImpl<base_matrix<T, Rows, Cols, ST, true>, MT>::dot(A, B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("dot(" << getAddress(A) << " " << getAddress(B) << ")");
            return implementation::DotImpl<base_matrix<T, Rows, Cols, ST, true>, MT>::dot(A, B);
        }
        
        template<typename RetType, template<typename DataType> class Operator, typename Data, typename MT>
        void dot(RetType & result, const implementation::expressions::expr<Operator, Data> & A, const MT & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::DotImpl<implementation::expressions::expr<Operator, Data>, MT>::dot(result, A, B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("dot(" << getAddress(result) << " " << getAddress(A) << " " << getAddress(B) << ")");
            implementation::DotImpl<implementation::expressions::expr<Operator, Data>, MT>::dot(result, A, B);
        }
        
        template<template<typename DataType> class Operator, typename Data, typename MT>
        typename implementation::DotImpl<implementation::expressions::expr<Operator, Data>, MT>::RetValue dot(const implementation::expressions::expr<Operator, Data> & A, const MT & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::DotImpl<implementation::expressions::expr<Operator, Data>, MT>::dot(A, B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("dot(" << getAddress(A) << " " << getAddress(B) << ")");
            return implementation::DotImpl<implementation::expressions::expr<Operator, Data>, MT>::dot(A, B);
        }
        
        // normSq: computes the squared l^2 norm of the matrix, i.e. sums up the squares of all entries
        
        namespace implementation
        {
            template<typename MatrixType>
            class NormSqImpl
            {
            public:
                typedef typename MatrixInfo<MatrixType>::Type RetValue;
            
                template<typename ReturnType>
                inline static void normSq(ReturnType & result, const MatrixType & A)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(setZero(result)) &&
                                                              noexcept(A.enumerate()) &&
                                                              noexcept(helper::make_type_lvalue<typename MatrixType::ConstEnumerator>().has_current()) &&
                                                              noexcept(helper::make_type_lvalue<typename MatrixType::ConstEnumerator>().next()) &&
                                                              noexcept(result = square(helper::make_type_lvalue<typename MatrixType::ConstEnumerator>().current())) &&
                                                              noexcept(result += square(helper::make_type_lvalue<typename MatrixType::ConstEnumerator>().current())))
                {
                    PLLL_DEBUG_OUTPUT_MESSAGE("NormSqImpl::normSq(" << getAddress(result) << " " << getAddress(A) << ")");
                    PLLL_INTERNAL_STATIC_CHECK(MatrixInfo<MatrixType>::is_matrix && MatrixInfo<MatrixType>::is_math_object, RequiresMathMatrix);
                    if ((0 == A.rows()) || (0 == A.cols()))
                        setZero(result);
                    else
                    {
                        typename MatrixType::ConstEnumerator e = A.enumerate();
                        result = square(e.current());
                        for (e.next(); e.has_current(); e.next())
                            result += square(e.current());
                    }
                }
            
                inline static RetValue normSq(const MatrixType & A)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(setZero(helper::make_type_lvalue<RetValue>())) &&
                                                              noexcept(A.enumerate()) && noexcept(RetValue()) &&
                                                              noexcept(helper::make_type_lvalue<typename MatrixType::ConstEnumerator>().has_current()) &&
                                                              noexcept(helper::make_type_lvalue<typename MatrixType::ConstEnumerator>().next()) &&
                                                              noexcept(RetValue(square(helper::make_type_lvalue<typename MatrixType::ConstEnumerator>().current()))) &&
                                                              noexcept(helper::make_type_lvalue<RetValue>() += square(helper::make_type_lvalue<typename MatrixType::ConstEnumerator>().current())))
                {
                    PLLL_DEBUG_OUTPUT_MESSAGE("NormSqImpl::normSq(" << getAddress(A) << ")");
                    PLLL_INTERNAL_STATIC_CHECK(MatrixInfo<MatrixType>::is_matrix && MatrixInfo<MatrixType>::is_math_object, RequiresMathMatrix);
                    if ((A.rows() == 0) || (A.cols() == 0))
                    {
                        RetValue result;
                        setZero(result);
                        return result;
                    }
                    else
                    {
                        typename MatrixType::ConstEnumerator e = A.enumerate();
                        RetValue result = square(e.current());
                        for (e.next(); e.has_current(); e.next())
                            result += square(e.current());
                        return result;
                    }
                }
            };
        }
        
        /**
           \brief Computes the squared norm of `A`, i.e. the sum of the squares of all entries of
                  `A`, and stores the result in `result`.
           
           \param result A scalar where to compute the result in.
           \param A The matrix whose square norm to compute.
         */
        template<typename RetType, typename T, int Rows, int Cols, typename ST>
        void normSq(RetType & result, const base_matrix<T, Rows, Cols, ST, true> & A)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::NormSqImpl<base_matrix<T, Rows, Cols, ST, true> >::normSq(result, A)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("normSq(" << getAddress(result) << " " << getAddress(A) << ")");
            implementation::NormSqImpl<base_matrix<T, Rows, Cols, ST, true> >::normSq(result, A);
        }
        
        /**
           \brief Computes the squared norm of `A`, i.e. the sum of the squares of all entries of
                  `A`, and returns the result.
           
           \return A scalar whose value is the squared norm of `A`.
           \param A The matrix whose square norm to compute.
         */
        template<typename T, int Rows, int Cols, typename ST>
        typename implementation::NormSqImpl<base_matrix<T, Rows, Cols, ST, true> >::RetValue normSq(const base_matrix<T, Rows, Cols, ST, true> & A)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::NormSqImpl<base_matrix<T, Rows, Cols, ST, true> >::normSq(A)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("normSq(" << getAddress(A) << ")");
            return implementation::NormSqImpl<base_matrix<T, Rows, Cols, ST, true> >::normSq(A);
        }
        
        template<typename RetType, template<typename DataType> class Operator, typename Data>
        void normSq(RetType & result, const implementation::expressions::expr<Operator, Data> & A)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::NormSqImpl<implementation::expressions::expr<Operator, Data> >::normSq(result, A)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("normSq(" << getAddress(result) << " " << getAddress(A) << ")");
            implementation::NormSqImpl<implementation::expressions::expr<Operator, Data> >::normSq(result, A);
        }
        
        template<template<typename DataType> class Operator, typename Data>
        typename implementation::NormSqImpl<implementation::expressions::expr<Operator, Data> >::RetValue normSq(const implementation::expressions::expr<Operator, Data> & A)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::NormSqImpl<implementation::expressions::expr<Operator, Data> >::normSq(A)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("normSq(" << getAddress(A) << ")");
            return implementation::NormSqImpl<implementation::expressions::expr<Operator, Data> >::normSq(A);
        }
        
        // addmul() and submul()
        
        /**
           \brief Computes `A += B * C`.
           
           \param A The variable to add the product to.
           \param B The first operand of the product.
           \param C The second operand of the product.
         */
        template<typename T, int Rows, int Cols, typename ST, typename T1, typename T2>
        void addmul(base_matrix<T, Rows, Cols, ST, true> & A, const T1 & B, const T2 & C)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(A += B * C))
        {
            A += B * C;
        }
        
        /**
           \brief Computes `A -= B * C`.
           
           \param A The variable to subtract the product from.
           \param B The first operand of the product.
           \param C The second operand of the product.
         */
        template<typename T, int Rows, int Cols, typename ST, typename T1, typename T2>
        void submul(base_matrix<T, Rows, Cols, ST, true> & A, const T1 & B, const T2 & C)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(A -= B * C))
        {
            A -= B * C;
        }
        
        template<template<typename DataType> class Op, typename Data, typename T1, typename T2>
        void addmul(const implementation::expressions::expr<Op, Data> & A, const T1 & B, const T2 & C)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(A += B * C))
        {
            A += B * C;
        }
        
        template<template<typename DataType> class Op, typename Data, typename T1, typename T2>
        void submul(const implementation::expressions::expr<Op, Data> & A, const T1 & B, const T2 & C)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(A -= B * C))
        {
            A -= B * C;
        }
        
        ///@}
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        //// comparison
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        /**@{
           \name Comparisons.
        */
        
        namespace implementation
        {
            template<typename MT1, typename MT2>
            class comparison_impl
            {
            private:
                class compare_impl
                {
                public:
                    typedef int IMResultType;
                    typedef int ResultType;
                    enum { default_value = 0 };
                    
                    static inline ResultType transform(IMResultType r)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        return r;
                    }
                    
                    static inline IMResultType cmp_dimension(size_type a, size_type b)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        if (a == b)
                            return 0;
                        else
                            return a < b ? -1 : 1;
                    }
                    
                    template<typename A, typename B>
                    static inline IMResultType cmp(const A & a, const B & b)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(a == b) && noexcept(a < b))
                    {
                        if (a == b)
                            return 0;
                        else
                            return a < b ? -1 : 1;
                    }
                };
                
                class eq_impl
                {
                public:
                    typedef bool IMResultType;
                    typedef bool ResultType;
                    enum { default_value = true };
                    
                    static inline ResultType transform(IMResultType r)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        return !r;
                    }
                    
                    static inline IMResultType cmp_dimension(size_type a, size_type b)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        return a != b;
                    }
                    
                    template<typename A, typename B>
                    static inline IMResultType cmp(const A & a, const B & b)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(a == b))
                    {
                        return !(a == b);
                    }
                };
                
                class less_impl
                {
                public:
                    typedef int IMResultType;
                    typedef bool ResultType;
                    enum { default_value = false };
                    
                    static inline ResultType transform(IMResultType r)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        return r < 0;
                    }
                    
                    static inline IMResultType cmp_dimension(size_type a, size_type b)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        if (a == b)
                            return 0;
                        else
                            return a < b ? -1 : 1;
                    }
                    
                    template<typename A, typename B>
                    static inline IMResultType cmp(const A & a, const B & b)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(a == b) && noexcept(a < b))
                    {
                        if (a == b)
                            return 0;
                        else
                            return a < b ? -1 : 1;
                    }
                };
                
                class less_eq_impl
                {
                public:
                    typedef int IMResultType;
                    typedef bool ResultType;
                    enum { default_value = true };
                    
                    static inline ResultType transform(IMResultType r)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        return r <= 0;
                    }
                    
                    static inline IMResultType cmp_dimension(size_type a, size_type b)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        if (a == b)
                            return 0;
                        else
                            return a < b ? -1 : 1;
                    }
                    
                    template<typename A, typename B>
                    static inline IMResultType cmp(const A & a, const B & b)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(a == b) && noexcept(a < b))
                    {
                        if (a == b)
                            return 0;
                        else
                            return a < b ? -1 : 1;
                    }
                };
                
                class greater_impl
                {
                public:
                    typedef int IMResultType;
                    typedef bool ResultType;
                    enum { default_value = false };
                    
                    static inline ResultType transform(IMResultType r)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        return r > 0;
                    }
                    
                    static inline IMResultType cmp_dimension(size_type a, size_type b)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        if (a == b)
                            return 0;
                        else
                            return a < b ? -1 : 1;
                    }
                    
                    template<typename A, typename B>
                    static inline IMResultType cmp(const A & a, const B & b)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(a == b) && noexcept(a < b))
                    {
                        if (a == b)
                            return 0;
                        else
                            return a < b ? -1 : 1;
                    }
                };
                
                class greater_eq_impl
                {
                public:
                    typedef int IMResultType;
                    typedef bool ResultType;
                    enum { default_value = true };
                    
                    static inline ResultType transform(IMResultType r)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        return r >= 0;
                    }
                    
                    static inline IMResultType cmp_dimension(size_type a, size_type b)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    {
                        if (a == b)
                            return 0;
                        else
                            return a < b ? -1 : 1;
                    }
                    
                    template<typename A, typename B>
                    static inline IMResultType cmp(const A & a, const B & b)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(a == b) && noexcept(a < b))
                    {
                        if (a == b)
                            return 0;
                        else
                            return a < b ? -1 : 1;
                    }
                };
                
                template<typename cmp_impl>
                class comparator
                {
                public:
                    static inline typename cmp_impl::ResultType cmp(const MT1 & A, const MT2 & B)
                        PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(A.enumerate()) && noexcept(B.enumerate()) &&
                                                                  noexcept(helper::make_type_lvalue<typename MT1::ConstEnumerator>().next()) &&
                                                                  noexcept(helper::make_type_lvalue<typename MT2::ConstEnumerator>().next()) &&
                                                                  noexcept(helper::make_type_lvalue<typename MT1::ConstEnumerator>().has_current()) &&
                                                                  noexcept(cmp_impl::cmp(helper::make_type_lvalue<typename MT1::ConstEnumerator>().current(),
                                                                                         helper::make_type_lvalue<typename MT2::ConstEnumerator>().current())))
                    {
                        PLLL_DEBUG_OUTPUT_MESSAGE("comparison_impl::comparator::cmp(" << getAddress(A) << " " << getAddress(B) << ") [0]");
                        if (typename cmp_impl::IMResultType res = cmp_impl::cmp_dimension(A.rows(), B.rows()))
                            return cmp_impl::transform(res);
                        if (typename cmp_impl::IMResultType res = cmp_impl::cmp_dimension(A.cols(), B.cols()))
                            return cmp_impl::transform(res);
                        typename MT1::ConstEnumerator eA = A.enumerate();
                        typename MT2::ConstEnumerator eB = B.enumerate();
                        for (; eA.has_current(); eA.next(), eB.next())
                            if (typename cmp_impl::IMResultType res = cmp_impl::cmp(eA.current(), eB.current()))
                                return cmp_impl::transform(res);
                        return cmp_impl::default_value;
                    }
                };
                
            public:
                static inline int compare(const MT1 & A, const MT2 & B)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(comparator<compare_impl>::cmp(A, B)))
                {
                    PLLL_DEBUG_OUTPUT_MESSAGE("comparison_impl::compare(" << getAddress(A) << " " << getAddress(B) << ")");
                    return comparator<compare_impl>::cmp(A, B);
                }
                
                static inline bool less(const MT1 & A, const MT2 & B)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(comparator<less_impl>::cmp(A, B)))
                {
                    PLLL_DEBUG_OUTPUT_MESSAGE("comparison_impl::less(" << getAddress(A) << " " << getAddress(B) << ")");
                    return comparator<less_impl>::cmp(A, B);
                }
                
                static inline bool less_eq(const MT1 & A, const MT2 & B)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(comparator<less_eq_impl>::cmp(A, B)))
                {
                    PLLL_DEBUG_OUTPUT_MESSAGE("comparison_impl::less_eq(" << getAddress(A) << " " << getAddress(B) << ")");
                    return comparator<less_eq_impl>::cmp(A, B);
                }
                
                static inline bool greater(const MT1 & A, const MT2 & B)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(comparator<greater_impl>::cmp(A, B)))
                {
                    PLLL_DEBUG_OUTPUT_MESSAGE("comparison_impl::greater(" << getAddress(A) << " " << getAddress(B) << ")");
                    return comparator<greater_impl>::cmp(A, B);
                }
                
                static inline bool greater_eq(const MT1 & A, const MT2 & B)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(comparator<greater_eq_impl>::cmp(A, B)))
                {
                    PLLL_DEBUG_OUTPUT_MESSAGE("comparison_impl::greater_eq(" << getAddress(A) << " " << getAddress(B) << ")");
                    return comparator<greater_eq_impl>::cmp(A, B);
                }
                
                static inline bool eq(const MT1 & A, const MT2 & B)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(comparator<eq_impl>::cmp(A, B)))
                {
                    PLLL_DEBUG_OUTPUT_MESSAGE("comparison_impl::eq(" << getAddress(A) << " " << getAddress(B) << ")");
                    return comparator<eq_impl>::cmp(A, B);
                }
                
                static inline bool neq(const MT1 & A, const MT2 & B)
                    PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(!eq(A, B)))
                {
                    PLLL_DEBUG_OUTPUT_MESSAGE("comparison_impl::neq(" << getAddress(A) << " " << getAddress(B) << ")");
                    return !eq(A, B);
                }
            };
        }
        
        template<template<typename AData> class AOp, typename AData, template<typename BData> class BOp, typename BData>
        inline int compare(const implementation::expressions::expr<AOp, AData> & A, const implementation::expressions::expr<BOp, BData> & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::comparison_impl<implementation::expressions::expr<AOp, AData>, implementation::expressions::expr<BOp, BData> >::compare(A, B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("compare(" << getAddress(A) << " " << getAddress(B) << ") [1]");
            return implementation::comparison_impl<implementation::expressions::expr<AOp, AData>, implementation::expressions::expr<BOp, BData> >::compare(A, B);
        }
        
        template<template<typename AData> class AOp, typename AData,
                 typename BT, int BRows, int BCols, typename BST, bool BMO>
        inline int compare(const implementation::expressions::expr<AOp, AData> & A, const base_matrix<BT, BRows, BCols, BST, BMO> & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(compare(A, implementation::expressions::make_matrix_expression(B))))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("compare(" << getAddress(A) << " " << getAddress(B) << ") [2]");
            return compare(A, implementation::expressions::make_matrix_expression(B));
        }
        
        template<typename AT, int ARows, int ACols, typename AST, bool AMO,
                 template<typename BData> class BOp, typename BData>
        inline int compare(const base_matrix<AT, ARows, ACols, AST, AMO> & A, const implementation::expressions::expr<BOp, BData> & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(compare(implementation::expressions::make_matrix_expression(A), B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("compare(" << getAddress(A) << " " << getAddress(B) << ") [3]");
            return compare(implementation::expressions::make_matrix_expression(A), B);
        }
        
        /**
           \brief Lexicographically compares the matrices `A` and `B`.

           First compares rows and columns, and then all entries. The exact procedure is described
           in \ref matrixvector_ops_general.
           
           \param A The first operand.
           \param B The second operand.
           \return Returns a negative integer if `A < B`, zero if `A == B` and a positive integer if `A > B`.
         */
        template<typename AT, int ARows, int ACols, typename AST, bool AMO,
                 typename BT, int BRows, int BCols, typename BST, bool BMO>
        inline int compare(const base_matrix<AT, ARows, ACols, AST, AMO> & A, const base_matrix<BT, BRows, BCols, BST, BMO> & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(compare(implementation::expressions::make_matrix_expression(A), implementation::expressions::make_matrix_expression(B))))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("compare(" << getAddress(A) << " " << getAddress(B) << ") [4]");
            return compare(implementation::expressions::make_matrix_expression(A), implementation::expressions::make_matrix_expression(B));
        }
        
        template<template<typename AData> class AOp, typename AData, template<typename BData> class BOp, typename BData>
        inline bool operator == (const implementation::expressions::expr<AOp, AData> & A, const implementation::expressions::expr<BOp, BData> & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::comparison_impl<implementation::expressions::expr<AOp, AData>, implementation::expressions::expr<BOp, BData> >::eq(A, B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator == (" << getAddress(A) << " " << getAddress(B) << ") [1]");
            return implementation::comparison_impl<implementation::expressions::expr<AOp, AData>, implementation::expressions::expr<BOp, BData> >::eq(A, B);
        }
        
        template<template<typename AData> class AOp, typename AData,
                 typename BT, int BRows, int BCols, typename BST, bool BMO>
        inline bool operator == (const implementation::expressions::expr<AOp, AData> & A, const base_matrix<BT, BRows, BCols, BST, BMO> & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(A == implementation::expressions::make_matrix_expression(B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator == (" << getAddress(A) << " " << getAddress(B) << ") [2]");
            return A == implementation::expressions::make_matrix_expression(B);
        }
        
        template<typename AT, int ARows, int ACols, typename AST, bool AMO,
                 template<typename BData> class BOp, typename BData>
        inline bool operator == (const base_matrix<AT, ARows, ACols, AST, AMO> & A, const implementation::expressions::expr<BOp, BData> & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::expressions::make_matrix_expression(A) == B))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator == (" << getAddress(A) << " " << getAddress(B) << ") [3]");
            return implementation::expressions::make_matrix_expression(A) == B;
        }
        
        /**
           \brief Tests `A` and `B` for equality.
           
           \param A The first operand.
           \param B The second operand.
           \return Returns `true` if and only if `A` and `B` have the same format and the same entries.
         */
        template<typename AT, int ARows, int ACols, typename AST, bool AMO,
                 typename BT, int BRows, int BCols, typename BST, bool BMO>
        inline bool operator == (const base_matrix<AT, ARows, ACols, AST, AMO> & A, const base_matrix<BT, BRows, BCols, BST, BMO> & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::expressions::make_matrix_expression(A) == implementation::expressions::make_matrix_expression(B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator == (" << getAddress(A) << " " << getAddress(B) << ") [4]");
            return implementation::expressions::make_matrix_expression(A) == implementation::expressions::make_matrix_expression(B);
        }
        
        template<template<typename AData> class AOp, typename AData, template<typename BData> class BOp, typename BData>
        inline bool operator != (const implementation::expressions::expr<AOp, AData> & A, const implementation::expressions::expr<BOp, BData> & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::comparison_impl<implementation::expressions::expr<AOp, AData>, implementation::expressions::expr<BOp, BData> >::neq(A, B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator != (" << getAddress(A) << " " << getAddress(B) << ") [1]");
            return implementation::comparison_impl<implementation::expressions::expr<AOp, AData>, implementation::expressions::expr<BOp, BData> >::neq(A, B);
        }
        
        template<template<typename AData> class AOp, typename AData,
                 typename BT, int BRows, int BCols, typename BST, bool BMO>
        inline bool operator != (const implementation::expressions::expr<AOp, AData> & A, const base_matrix<BT, BRows, BCols, BST, BMO> & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(A != implementation::expressions::make_matrix_expression(B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator != (" << getAddress(A) << " " << getAddress(B) << ") [2]");
            return A != implementation::expressions::make_matrix_expression(B);
        }
        
        template<typename AT, int ARows, int ACols, typename AST, bool AMO,
                 template<typename BData> class BOp, typename BData>
        inline bool operator != (const base_matrix<AT, ARows, ACols, AST, AMO> & A, const implementation::expressions::expr<BOp, BData> & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::expressions::make_matrix_expression(A) != B))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator != (" << getAddress(A) << " " << getAddress(B) << ") [3]");
            return implementation::expressions::make_matrix_expression(A) != B;
        }
        
        /**
           \brief Tests `A` and `B` for inequality.
           
           \param A The first operand.
           \param B The second operand.
           \return Returns `false` if and only if `A` and `B` have the same format and the same entries.
         */
        template<typename AT, int ARows, int ACols, typename AST, bool AMO,
                 typename BT, int BRows, int BCols, typename BST, bool BMO>
        inline bool operator != (const base_matrix<AT, ARows, ACols, AST, AMO> & A, const base_matrix<BT, BRows, BCols, BST, BMO> & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::expressions::make_matrix_expression(A) != implementation::expressions::make_matrix_expression(B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator != (" << getAddress(A) << " " << getAddress(B) << ") [4]");
            return implementation::expressions::make_matrix_expression(A) != implementation::expressions::make_matrix_expression(B);
        }
        
        template<template<typename AData> class AOp, typename AData, template<typename BData> class BOp, typename BData>
        inline bool operator < (const implementation::expressions::expr<AOp, AData> & A, const implementation::expressions::expr<BOp, BData> & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::comparison_impl<implementation::expressions::expr<AOp, AData>, implementation::expressions::expr<BOp, BData> >::less(A, B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator < (" << getAddress(A) << " " << getAddress(B) << ") [1]");
            return implementation::comparison_impl<implementation::expressions::expr<AOp, AData>, implementation::expressions::expr<BOp, BData> >::less(A, B);
        }
        
        template<template<typename AData> class AOp, typename AData,
                 typename BT, int BRows, int BCols, typename BST, bool BMO>
        inline bool operator < (const implementation::expressions::expr<AOp, AData> & A, const base_matrix<BT, BRows, BCols, BST, BMO> & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(A < implementation::expressions::make_matrix_expression(B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator < (" << getAddress(A) << " " << getAddress(B) << ") [2]");
            return A < implementation::expressions::make_matrix_expression(B);
        }
        
        template<typename AT, int ARows, int ACols, typename AST, bool AMO,
                 template<typename BData> class BOp, typename BData>
        inline bool operator < (const base_matrix<AT, ARows, ACols, AST, AMO> & A, const implementation::expressions::expr<BOp, BData> & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::expressions::make_matrix_expression(A) < B))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator < (" << getAddress(A) << " " << getAddress(B) << ") [3]");
            return implementation::expressions::make_matrix_expression(A) < B;
        }
        
        /**
           \brief Tests `A` and `B` lexicographically for `A` being less than `B`.

           The exact procedure is described in \ref matrixvector_ops_general.
           
           \param A The first operand.
           \param B The second operand.
           \return Returns `true` if and only if \f$A < B\f$ lexicographically.
         */
        template<typename AT, int ARows, int ACols, typename AST, bool AMO,
                 typename BT, int BRows, int BCols, typename BST, bool BMO>
        inline bool operator < (const base_matrix<AT, ARows, ACols, AST, AMO> & A, const base_matrix<BT, BRows, BCols, BST, BMO> & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::expressions::make_matrix_expression(A) < implementation::expressions::make_matrix_expression(B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator < (" << getAddress(A) << " " << getAddress(B) << ") [4]");
            return implementation::expressions::make_matrix_expression(A) < implementation::expressions::make_matrix_expression(B);
        }
        
        template<template<typename AData> class AOp, typename AData, template<typename BData> class BOp, typename BData>
        inline bool operator <= (const implementation::expressions::expr<AOp, AData> & A, const implementation::expressions::expr<BOp, BData> & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::comparison_impl<implementation::expressions::expr<AOp, AData>, implementation::expressions::expr<BOp, BData> >::less_eq(A, B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator <= (" << getAddress(A) << " " << getAddress(B) << ") [1]");
            return implementation::comparison_impl<implementation::expressions::expr<AOp, AData>, implementation::expressions::expr<BOp, BData> >::less_eq(A, B);
        }
        
        template<template<typename AData> class AOp, typename AData,
                 typename BT, int BRows, int BCols, typename BST, bool BMO>
        inline bool operator <= (const implementation::expressions::expr<AOp, AData> & A, const base_matrix<BT, BRows, BCols, BST, BMO> & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(A <= implementation::expressions::make_matrix_expression(B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator <= (" << getAddress(A) << " " << getAddress(B) << ") [2]");
            return A <= implementation::expressions::make_matrix_expression(B);
        }
        
        template<typename AT, int ARows, int ACols, typename AST, bool AMO,
                 template<typename BData> class BOp, typename BData>
        inline bool operator <= (const base_matrix<AT, ARows, ACols, AST, AMO> & A, const implementation::expressions::expr<BOp, BData> & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::expressions::make_matrix_expression(A) <= B))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator <= (" << getAddress(A) << " " << getAddress(B) << ") [3]");
            return implementation::expressions::make_matrix_expression(A) <= B;
        }
        
        /**
           \brief Tests `A` and `B` lexicographically for `A` being less than or equal to `B`.

           The exact procedure is described in \ref matrixvector_ops_general.
           
           \param A The first operand.
           \param B The second operand.
           \return Returns `true` if and only if \f$A \le B\f$ lexicographically.
         */
        template<typename AT, int ARows, int ACols, typename AST, bool AMO,
                 typename BT, int BRows, int BCols, typename BST, bool BMO>
        inline bool operator <= (const base_matrix<AT, ARows, ACols, AST, AMO> & A, const base_matrix<BT, BRows, BCols, BST, BMO> & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::expressions::make_matrix_expression(A) <= implementation::expressions::make_matrix_expression(B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator <= (" << getAddress(A) << " " << getAddress(B) << ") [4]");
            return implementation::expressions::make_matrix_expression(A) <= implementation::expressions::make_matrix_expression(B);
        }
        
        template<template<typename AData> class AOp, typename AData, template<typename BData> class BOp, typename BData>
        inline bool operator > (const implementation::expressions::expr<AOp, AData> & A, const implementation::expressions::expr<BOp, BData> & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::comparison_impl<implementation::expressions::expr<AOp, AData>, implementation::expressions::expr<BOp, BData> >::greater(A, B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator > (" << getAddress(A) << " " << getAddress(B) << ") [1]");
            return implementation::comparison_impl<implementation::expressions::expr<AOp, AData>, implementation::expressions::expr<BOp, BData> >::greater(A, B);
        }
        
        template<template<typename AData> class AOp, typename AData,
                 typename BT, int BRows, int BCols, typename BST, bool BMO>
        inline bool operator > (const implementation::expressions::expr<AOp, AData> & A, const base_matrix<BT, BRows, BCols, BST, BMO> & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(A > implementation::expressions::make_matrix_expression(B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator > (" << getAddress(A) << " " << getAddress(B) << ") [2]");
            return A > implementation::expressions::make_matrix_expression(B);
        }
        
        template<typename AT, int ARows, int ACols, typename AST, bool AMO,
                 template<typename BData> class BOp, typename BData>
        inline bool operator > (const base_matrix<AT, ARows, ACols, AST, AMO> & A, const implementation::expressions::expr<BOp, BData> & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::expressions::make_matrix_expression(A) > B))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator > (" << getAddress(A) << " " << getAddress(B) << ") [3]");
            return implementation::expressions::make_matrix_expression(A) > B;
        }
        
        /**
           \brief Tests `A` and `B` lexicographically for `A` being greater than `B`.

           The exact procedure is described in \ref matrixvector_ops_general.
           
           \param A The first operand.
           \param B The second operand.
           \return Returns `true` if and only if \f$A > B\f$ lexicographically.
         */
        template<typename AT, int ARows, int ACols, typename AST, bool AMO,
                 typename BT, int BRows, int BCols, typename BST, bool BMO>
        inline bool operator > (const base_matrix<AT, ARows, ACols, AST, AMO> & A, const base_matrix<BT, BRows, BCols, BST, BMO> & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::expressions::make_matrix_expression(A) > implementation::expressions::make_matrix_expression(B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator > (" << getAddress(A) << " " << getAddress(B) << ") [4]");
            return implementation::expressions::make_matrix_expression(A) > implementation::expressions::make_matrix_expression(B);
        }
        
        template<template<typename AData> class AOp, typename AData, template<typename BData> class BOp, typename BData>
        inline bool operator >= (const implementation::expressions::expr<AOp, AData> & A, const implementation::expressions::expr<BOp, BData> & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::comparison_impl<implementation::expressions::expr<AOp, AData>, implementation::expressions::expr<BOp, BData> >::greater_eq(A, B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator >= (" << getAddress(A) << " " << getAddress(B) << ") [1]");
            return implementation::comparison_impl<implementation::expressions::expr<AOp, AData>, implementation::expressions::expr<BOp, BData> >::greater_eq(A, B);
        }
        
        template<template<typename AData> class AOp, typename AData,
                 typename BT, int BRows, int BCols, typename BST, bool BMO>
        inline bool operator >= (const implementation::expressions::expr<AOp, AData> & A, const base_matrix<BT, BRows, BCols, BST, BMO> & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(A >= implementation::expressions::make_matrix_expression(B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator >= (" << getAddress(A) << " " << getAddress(B) << ") [2]");
            return A >= implementation::expressions::make_matrix_expression(B);
        }
        
        template<typename AT, int ARows, int ACols, typename AST, bool AMO,
                 template<typename BData> class BOp, typename BData>
        inline bool operator >= (const base_matrix<AT, ARows, ACols, AST, AMO> & A, const implementation::expressions::expr<BOp, BData> & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::expressions::make_matrix_expression(A) >= B))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator >= (" << getAddress(A) << " " << getAddress(B) << ") [3]");
            return implementation::expressions::make_matrix_expression(A) >= B;
        }
        
        /**
           \brief Tests `A` and `B` lexicographically for `A` being greater than or equal to `B`.
           
           The exact procedure is described in \ref matrixvector_ops_general.
           
           \param A The first operand.
           \param B The second operand.
           \return Returns `true` if and only if \f$A \ge B\f$ lexicographically.
         */
        template<typename AT, int ARows, int ACols, typename AST, bool AMO,
                 typename BT, int BRows, int BCols, typename BST, bool BMO>
        inline bool operator >= (const base_matrix<AT, ARows, ACols, AST, AMO> & A, const base_matrix<BT, BRows, BCols, BST, BMO> & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::expressions::make_matrix_expression(A) >= implementation::expressions::make_matrix_expression(B)))
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator >= (" << getAddress(A) << " " << getAddress(B) << ") [4]");
            return implementation::expressions::make_matrix_expression(A) >= implementation::expressions::make_matrix_expression(B);
        }
        ///@}
    }
}

#endif

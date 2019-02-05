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

#ifndef PLLL_INCLUDE_GUARD__MATRIX_HPP
#define PLLL_INCLUDE_GUARD__MATRIX_HPP

// #define PLLL_DEBUG_MATRIX_DEBUG

#include "helper.hpp"
#include "matrix-mem.hpp"
#include <list>
#include <cctype>
#include <iostream>

#if __cplusplus >= 201103L
    #include <utility>
#endif

/**
   \file
   \brief The matrix and vector template library.
   
   This header provides the matrix and vector template library for `plll`. It declares the main
   templates and includes the other header files related to the library. See \ref matrixvector for
   documentation of the matrix and vector library.
*/

namespace plll
{
    #ifndef PLLL_DEBUG_MATRIX_DEBUG
        #define PLLL_DEBUG_OUTPUT_MESSAGE(msg)
    #else
        #define PLLL_DEBUG_OUTPUT_MESSAGE(msg) { std::cout << msg << "\n"; }
        
        template<typename T>
        inline void * getAddress(const T & t) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return &const_cast<T &>(t);
        }
        
        template<typename T>
        inline void * getAddress(const T * t); // should yield linker error
    #endif
    
    /**
       \brief Contains the matrix and vector library as long as other linear algebra functionality.
       
       The `plll::linalg` namespace contains the matrix and vector library as long as other linear
       algebra functionality of `plll`.
    */
    namespace linalg
    {
        enum { Flexible = -1 /**< Identifier to make number of rows and/or columns flexible instead
                                  of determined at compile time. */ };
        
        namespace implementation
        {
            template<typename S, bool def>
            class Initialize_Impl;
            
            template<typename S>
            class Initialize_Impl<S, false>
            {
            private:
                const S & d_ref;
                
            public:
                Initialize_Impl(const S & ref) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_ref(ref)
                {
                }
                
                const S & ref() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    return d_ref;
                }
            };
            
            template<typename S>
            class Initialize_Impl<S, true>
            {
            private:
                S d_obj;
                
            public:
                const S & ref() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    return d_obj;
                }
            };
        }
        
        /**@{
           \name Initializers.
        */
        template<typename S>
        /**
           \brief Creates an initializer object to initialize with the given reference.
           
           \tparam S The type of the object to initialize with.
           \param ref The object to initialize with.
           \return An initializer object.
         */
        inline implementation::Initialize_Impl<S, false> Initialize(const S & ref) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            return implementation::Initialize_Impl<S, false>(ref);
        }
        
        template<typename S>
        /**
           \brief Creates an initializer object to default-initialize with type `S`.
           
           \tparam S The type of the object to default-initialize with.
           \return An initializer object.
         */
        inline implementation::Initialize_Impl<S, true> Initialize() PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(S()))
        {
            return implementation::Initialize_Impl<S, true>();
        }
        ///@}
        
        typedef unsigned int size_type; /**< Type used for indices and sizes. */
        
        namespace implementation
        {
            template<typename MatrixType>
            struct MatrixInfo
            {
                enum { is_matrix = false,
                       is_const = false,
                       is_math_object = false,
                       is_expression = false,
                       can_move_from = false,
                       can_assign_to = false,
                       can_resize_rows = false,
                       can_resize_cols = false,
                       use_temporary_on_evaluate = false,
                       coeffs_are_simple_expressions = false };
                enum { rows = 0, cols = 0 };
                typedef void Type;
                typedef void StorageTraits;
            };
        }

        template<typename T, int Rows, int Cols, typename StorageTraits, bool MathObject>
        class base_matrix;
        
        namespace implementation
        {
            template<typename T, int Rows, int Cols, typename ST, bool MO>
            struct MatrixInfo<base_matrix<T, Rows, Cols, ST, MO> >
            {
                enum { is_matrix = true,
                       is_const = false,
                       is_math_object = MO,
                       is_expression = false,
                       only_defined_for_matrices = 1,
                       can_move_from = true,
                       can_assign_to = true,
                       can_resize_rows = Rows < 0,
                       can_resize_cols = Cols < 0,
                       use_temporary_on_evaluate = false,
                       coeffs_are_simple_expressions = true };
                enum { rows = Rows, cols = Cols };
                typedef T Type;
                typedef ST StorageTraits;
            };
            
            template<typename T, int Rows, int Cols, typename ST, bool MO>
            struct MatrixInfo<const base_matrix<T, Rows, Cols, ST, MO> >
            {
                enum { is_matrix = true,
                       is_const = true,
                       is_math_object = MO,
                       is_expression = false,
                       only_defined_for_matrices = 1,
                       can_move_from = false,
                       can_assign_to = true,
                       can_resize_rows = false, // const matrices cannot be resized
                       can_resize_cols = false,
                       use_temporary_on_evaluate = false,
                       coeffs_are_simple_expressions = true };
                enum { rows = Rows, cols = Cols };
                typedef T Type;
                typedef ST StorageTraits;
            };
        }
        
#if __cplusplus < 201103L
        template<typename T, int Rows, int Cols, typename StorageTraits>
        class math_matrix;
        
        namespace implementation
        {
            template<typename T, int Rows, int Cols, typename ST>
            struct MatrixInfo<math_matrix<T, Rows, Cols, ST> >
            {
                enum { is_matrix = true,
                       is_const = false,
                       is_math_object = true,
                       is_expression = false,
                       only_defined_for_matrices = 1,
                       can_move_from = true,
                       can_assign_to = true,
                       can_resize_rows = Rows < 0,
                       can_resize_cols = Cols < 0,
                       use_temporary_on_evaluate = false,
                       coeffs_are_simple_expressions = true };
                enum { rows = Rows, cols = Cols };
                typedef T Type;
                typedef ST StorageTraits;
            };
            
            template<typename T, int Rows, int Cols, typename ST>
            struct MatrixInfo<const math_matrix<T, Rows, Cols, ST> >
            {
                enum { is_matrix = true,
                       is_const = true,
                       is_math_object = true,
                       is_expression = false,
                       only_defined_for_matrices = 1,
                       can_move_from = false,
                       can_assign_to = true,
                       can_resize_rows = false, // const matrices cannot be resized
                       can_resize_cols = false,
                       use_temporary_on_evaluate = false,
                       coeffs_are_simple_expressions = true };
                enum { rows = Rows, cols = Cols };
                typedef T Type;
                typedef ST StorageTraits;
            };
        }
        
        template<typename T, int Cols, typename StorageTraits, bool MO>
        class base_rowvector;
        
        namespace implementation
        {
            template<typename T, int Cols, typename ST, bool MO>
            struct MatrixInfo<base_rowvector<T, Cols, ST, MO> >
            {
                enum { is_matrix = true,
                       is_const = false,
                       is_math_object = MO,
                       is_expression = false,
                       only_defined_for_matrices = 1,
                       can_move_from = true,
                       can_assign_to = true,
                       can_resize_rows = false,
                       can_resize_cols = Cols < 0,
                       use_temporary_on_evaluate = false,
                       coeffs_are_simple_expressions = true };
                enum { rows = 1, cols = Cols };
                typedef T Type;
                typedef ST StorageTraits;
            };
            
            template<typename T, int Cols, typename ST, bool MO>
            struct MatrixInfo<const base_rowvector<T, Cols, ST, MO> >
            {
                enum { is_matrix = true,
                       is_const = true,
                       is_math_object = MO,
                       is_expression = false,
                       only_defined_for_matrices = 1,
                       can_move_from = false,
                       can_assign_to = true,
                       can_resize_rows = false, // const matrices cannot be resized
                       can_resize_cols = false,
                       use_temporary_on_evaluate = false,
                       coeffs_are_simple_expressions = true };
                enum { rows = 1, cols = Cols };
                typedef T Type;
                typedef ST StorageTraits;
            };
        }
            
        template<typename T, int Rows, typename StorageTraits, bool MO>
        class base_colvector;
        
        namespace implementation
        {
            template<typename T, int Rows, typename ST, bool MO>
            struct MatrixInfo<base_colvector<T, Rows, ST, MO> >
            {
                enum { is_matrix = true,
                       is_const = false,
                       is_math_object = MO,
                       is_expression = false,
                       only_defined_for_matrices = 1,
                       can_move_from = true,
                       can_assign_to = true,
                       can_resize_rows = Rows < 0,
                       can_resize_cols = false,
                       use_temporary_on_evaluate = false,
                       coeffs_are_simple_expressions = true };
                enum { rows = Rows, cols = 1 };
                typedef T Type;
                typedef ST StorageTraits;
            };
            
            template<typename T, int Rows, typename ST, bool MO>
            struct MatrixInfo<const base_colvector<T, Rows, ST, MO> >
            {
                enum { is_matrix = true,
                       is_const = true,
                       is_math_object = MO,
                       is_expression = false,
                       only_defined_for_matrices = 1,
                       can_move_from = false,
                       can_assign_to = true,
                       can_resize_rows = false, // const matrices cannot be resized
                       can_resize_cols = false,
                       use_temporary_on_evaluate = false,
                       coeffs_are_simple_expressions = true };
                enum { rows = Rows, cols = 1 };
                typedef T Type;
                typedef ST StorageTraits;
            };
        }
        
        template<typename T, int Cols, typename StorageTraits>
        class math_rowvector;
        
        namespace implementation
        {
            template<typename T, int Cols, typename ST>
            struct MatrixInfo<math_rowvector<T, Cols, ST> >
            {
                enum { is_matrix = true,
                       is_const = false,
                       is_math_object = true,
                       is_expression = false,
                       only_defined_for_matrices = 1,
                       can_move_from = true,
                       can_assign_to = true,
                       can_resize_rows = false,
                       can_resize_cols = Cols < 0,
                       use_temporary_on_evaluate = false,
                       coeffs_are_simple_expressions = true };
                enum { rows = 1, cols = Cols };
                typedef T Type;
                typedef ST StorageTraits;
            };
            
            template<typename T, int Cols, typename ST>
            struct MatrixInfo<const math_rowvector<T, Cols, ST> >
            {
                enum { is_matrix = true,
                       is_const = true,
                       is_math_object = true,
                       is_expression = false,
                       only_defined_for_matrices = 1,
                       can_move_from = false,
                       can_assign_to = true,
                       can_resize_rows = false, // const matrices cannot be resized
                       can_resize_cols = false,
                       use_temporary_on_evaluate = false,
                       coeffs_are_simple_expressions = true };
                enum { rows = 1, cols = Cols };
                typedef T Type;
                typedef ST StorageTraits;
            };
        }
        
        template<typename T, int Rows, typename StorageTraits>
        class math_colvector;
        
        namespace implementation
        {
            template<typename T, int Rows, typename ST>
            struct MatrixInfo<math_colvector<T, Rows, ST> >
            {
                enum { is_matrix = true,
                       is_const = false,
                       is_math_object = true,
                       is_expression = false,
                       only_defined_for_matrices = 1,
                       can_move_from = true,
                       can_assign_to = true,
                       can_resize_rows = Rows < 0,
                       can_resize_cols = false,
                       use_temporary_on_evaluate = false,
                       coeffs_are_simple_expressions = true };
                enum { rows = Rows, cols = 1 };
                typedef T Type;
                typedef ST StorageTraits;
            };
            
            template<typename T, int Rows, typename ST>
            struct MatrixInfo<const math_colvector<T, Rows, ST> >
            {
                enum { is_matrix = true,
                       is_const = true,
                       is_math_object = true,
                       is_expression = false,
                       only_defined_for_matrices = 1,
                       can_move_from = false,
                       can_assign_to = true,
                       can_resize_rows = false, // const matrices cannot be resized
                       can_resize_cols = false,
                       use_temporary_on_evaluate = false,
                       coeffs_are_simple_expressions = true };
                enum { rows = Rows, cols = 1 };
                typedef T Type;
                typedef ST StorageTraits;
            };
        }
#endif
            
        namespace implementation
        {
            namespace expressions
            {
                template<template<typename DataType> class Operator, typename Data>
                class expr;
                
                template<typename Data>
                class identity_operation;
                
                template<typename MatrixType>
                class MatrixWrapper;
                
                template<typename MatrixType>
                class ConstMatrixWrapper;
                
                template<typename MatrixType>
                inline MatrixWrapper<MatrixType> make_matrix_wrapper(MatrixType &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
                
                template<typename MatrixType>
                inline ConstMatrixWrapper<MatrixType> make_matrix_wrapper(const MatrixType &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
                
                template<typename MatrixType>
                inline expr<identity_operation, MatrixWrapper<MatrixType> > make_matrix_expression(MatrixType &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
                
                template<typename MatrixType>
                inline expr<identity_operation, ConstMatrixWrapper<MatrixType> > make_matrix_expression(const MatrixType &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
                
                template<typename MatrixType>
                class MatrixTemporaryWrapper;
                
                template<typename ScalarType>
                class ScalarWrapper;
                
                template<typename ScalarType>
                inline ScalarWrapper<ScalarType> make_scalar_wrapper(const ScalarType &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
                
                template<typename PairMatrices>
                class matrix_matrix_multiplication;
                
                template<typename OpA, typename OpB>
                class MSMul;
                
                template<typename OpA, typename OpB>
                class SMMul;
                
                template<typename OpA, typename OpB>
                class MSDiv;
                
                template<typename OpA, typename OpB>
                class MSMod;
                
                template<template<typename OpA, typename OpB> class Operation>
                class matrix_scalar_operation;
                
                template<typename OpA, typename OpB>
                class CWAdd;
                
                template<typename OpA, typename OpB>
                class CWSub;
                
                template<typename OpA, typename OpB>
                class CWMul;
                
                template<typename OpA, typename OpB>
                class CWDiv;
                
                template<typename OpA, typename OpB>
                class CWMod;
                
                template<template<typename OpA, typename OpB> class Operation>
                class componentwise_operation;
                
                template<typename MTIn>
                class matrix_negation;
                
                template<int Rows, int Cols>
                class sub
                {
                public:
                    template<typename MT>
                    class operation_generic;
                };
                
                template<int Dim>
                class sub_1d
                {
                public:
                    template<typename MT>
                    class operation_row;
                    
                    template<typename MT>
                    class operation_col;
                };
                
                template<typename MT>
                class transpose;
            }
            
            template<template<typename DataType> class Operator, typename Data>
            struct MatrixInfo<expressions::expr<Operator, Data> >
            {
                enum { is_matrix = MatrixInfo<typename Operator<Data>::ValueType>::is_matrix,
                       is_const = Operator<Data>::is_const,
                       is_math_object = MatrixInfo<typename Operator<Data>::ValueType>::is_math_object,
                       is_expression = true,
                       only_defined_for_matrices = 1,
                       can_move_from = Operator<Data>::can_move_from,
                       can_assign_to = Operator<Data>::can_assign_to,
                       can_resize_rows = Operator<Data>::can_resize_rows,
                       can_resize_cols = Operator<Data>::can_resize_cols,
                       use_temporary_on_evaluate = Operator<Data>::use_temporary_on_evaluate,
                       coeffs_are_simple_expressions = Operator<Data>::coeffs_are_simple_expressions };
                enum { rows = MatrixInfo<typename Operator<Data>::ValueType>::rows, cols = MatrixInfo<typename Operator<Data>::ValueType>::cols };
                typedef typename MatrixInfo<typename Operator<Data>::ValueType>::Type Type;
                typedef typename MatrixInfo<typename Operator<Data>::ValueType>::StorageTraits StorageTraits;
            };
            
            template<template<typename DataType> class Operator, typename Data>
            struct MatrixInfo<const expressions::expr<Operator, Data> >
            {
                enum { is_matrix = MatrixInfo<typename Operator<Data>::ValueType>::is_matrix,
                       is_const = Operator<Data>::is_const,
                       is_math_object = MatrixInfo<typename Operator<Data>::ValueType>::is_math_object,
                       is_expression = true,
                       only_defined_for_matrices = 1,
                       can_move_from = Operator<Data>::can_move_from,
                       can_assign_to = Operator<Data>::can_assign_to,
                       can_resize_rows = Operator<Data>::can_resize_rows,
                       can_resize_cols = Operator<Data>::can_resize_cols,
                       use_temporary_on_evaluate = Operator<Data>::use_temporary_on_evaluate,
                       coeffs_are_simple_expressions = Operator<Data>::coeffs_are_simple_expressions };
                enum { rows = MatrixInfo<typename Operator<Data>::ValueType>::rows, cols = MatrixInfo<typename Operator<Data>::ValueType>::cols };
                typedef typename MatrixInfo<typename Operator<Data>::ValueType>::Type Type;
                typedef typename MatrixInfo<typename Operator<Data>::ValueType>::StorageTraits StorageTraits;
            };
        }
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        //// swap
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        
        template<typename T, int R, int C, typename ST, bool MO>
        /**
           \brief Swaps matrices `A` and `B` efficiently.
           
           \param A First matrix.
           \param B Second matrix.
         */
        void swap(base_matrix<T, R, C, ST, MO> & A, base_matrix<T, R, C, ST, MO> & B) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        // We need this overload so that std::swap() won't be used in some obscure cases...
        
        template<typename T, int R1, int C1, int R2, int C2, typename ST, bool MO1, bool MO2>
        /**
           \brief Swaps matrices `A` and `B` efficiently.
           
           \param A First matrix.
           \param B Second matrix.
         */
        void swap(base_matrix<T, R1, C1, ST, MO1> & A, base_matrix<T, R2, C2, ST, MO2> & B) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
        
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
                                                                        helper::make_type_lvalue<typename implementation::expressions::expr<BOp, BData>::Enumerator>().current())));
        }
        
        template<template<typename AData> class AOp, typename AData, template<typename BData> class BOp, typename BData>
        /**
           \brief Swaps matrices `A` and `B` componentwise.
           
           \param A First matrix.
           \param B Second matrix.
         */
        void swap(const implementation::expressions::expr<AOp, AData> & A, const implementation::expressions::expr<BOp, BData> & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(implementation::do_swap(A, B)));
        
        template<template<typename AData> class AOp, typename AData,
                 typename BT, int BRows, int BCols, typename BST, bool BMO>
        /**
           \brief Swaps matrices `A` and `B` componentwise.
           
           \param A First matrix.
           \param B Second matrix.
         */
        void swap(const implementation::expressions::expr<AOp, AData> & A, base_matrix<BT, BRows, BCols, BST, BMO> & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(linalg::swap(A, implementation::expressions::make_matrix_expression(B))));
        
        template<typename AT, int ARows, int ACols, typename AST, bool AMO,
                 template<typename BData> class BOp, typename BData>
        /**
           \brief Swaps matrices `A` and `B` componentwise.
           
           \param A First matrix.
           \param B Second matrix.
         */
        void swap(base_matrix<AT, ARows, ACols, AST, AMO> & A, const implementation::expressions::expr<BOp, BData> & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(linalg::swap(implementation::expressions::make_matrix_expression(A), B)));
        
        template<typename T,
                 int ARows, int ACols, typename AST, bool AMO,
                 int BRows, int BCols, typename BST, bool BMO>
        /**
           \brief Swaps matrices `A` and `B` componentwise.
           
           \param A First matrix.
           \param B Second matrix.
         */
        void swap(base_matrix<T, ARows, ACols, AST, AMO> & A, base_matrix<T, BRows, BCols, BST, BMO> & B)
            PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(linalg::swap(implementation::expressions::make_matrix_expression(A), implementation::expressions::make_matrix_expression(B))));
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        //// assignment
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        
        template<template<typename SourceData> class SourceOp, typename SourceData, template<typename DestData> class DestOp, typename DestData, bool move>
        /**
           \brief Assigns the content of matrix `source` to matrix `destination`. Allows to force
                  move.
           
           \param source Source matrix, whose contents are copied.
           \param destination Destination matrix, where the contents of the source matrix are stored.
           \param ittm Indicates whether a copy or move (if possible) should be executed.
         */
        void assign(const implementation::expressions::expr<DestOp, DestData> & destination,
                    const implementation::expressions::expr<SourceOp, SourceData> & source, helper::BoolToType<move> ittm);
        
        template<template<typename SourceData> class SourceOp, typename SourceData, template<typename DestData> class DestOp, typename DestData>
        /**
           \brief Assigns the content of matrix `source` to matrix `destination`.
           
           \param source Source matrix, whose contents are copied.
           \param destination Destination matrix, where the contents of the source matrix are stored.
         */
        inline void assign(const implementation::expressions::expr<DestOp, DestData> & destination,
                           const implementation::expressions::expr<SourceOp, SourceData> & source);
        
        template<typename SourceT, int SourceRows, int SourceCols, typename SourceST, bool SourceMO, template<typename DestData> class DestOp, typename DestData>
        /**
           \brief Assigns the content of matrix `source` to matrix `destination`.
           
           \param source Source matrix, whose contents are copied.
           \param destination Destination matrix, where the contents of the source matrix are stored.
         */
        inline void assign(const implementation::expressions::expr<DestOp, DestData> & destination,
                           const base_matrix<SourceT, SourceRows, SourceCols, SourceST, SourceMO> & source);
        
        template<template<typename SourceData> class SourceOp, typename SourceData, typename DestT, int DestRows, int DestCols, typename DestST, bool DestMO>
        /**
           \brief Assigns the content of matrix `source` to matrix `destination`.
           
           \param source Source matrix, whose contents are copied.
           \param destination Destination matrix, where the contents of the source matrix are stored.
         */
        inline void assign(base_matrix<DestT, DestRows, DestCols, DestST, DestMO> & destination,
                           const implementation::expressions::expr<SourceOp, SourceData> & source);
        
        template<typename SourceT, int SourceRows, int SourceCols, typename SourceST, bool SourceMO,
                 typename DestT, int DestRows, int DestCols, typename DestST, bool DestMO>
        /**
           \brief Assigns the content of matrix `source` to matrix `destination`.
           
           \param source Source matrix, whose contents are copied.
           \param destination Destination matrix, where the contents of the source matrix are stored.
         */
        inline void assign(base_matrix<DestT, DestRows, DestCols, DestST, DestMO> & destination,
                           const base_matrix<SourceT, SourceRows, SourceCols, SourceST, SourceMO> & source);
        
#if __cplusplus >= 201103L
        template<typename SourceT, int SourceRows, int SourceCols, typename SourceST, bool SourceMO, template<typename DestData> class DestOp, typename DestData>
        /**
           \brief Assigns the content of matrix `source` to matrix `destination`.
           
           \param source Source matrix, whose contents are copied.
           \param destination Destination matrix, where the contents of the source matrix are stored.
         */
        inline void assign(const implementation::expressions::expr<DestOp, DestData> & destination,
                           base_matrix<SourceT, SourceRows, SourceCols, SourceST, SourceMO> && source);
        
        template<typename SourceT, int SourceRows, int SourceCols, typename SourceST, bool SourceMO,
                 typename DestT, int DestRows, int DestCols, typename DestST, bool DestMO>
        /**
           \brief Assigns the content of matrix `source` to matrix `destination`.
           
           \param source Source matrix, whose contents are copied.
           \param destination Destination matrix, where the contents of the source matrix are stored.
         */
        inline void assign(base_matrix<DestT, DestRows, DestCols, DestST, DestMO> & destination,
                           base_matrix<SourceT, SourceRows, SourceCols, SourceST, SourceMO> && source);
#endif
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        //// transpose (functional)
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        
        template<template<typename SourceData> class SourceOp, typename SourceData, template<typename DestData> class DestOp, typename DestData, bool move>
        /**
           \brief Transposes the content of matrix `source` and stores the result in the matrix
                  `destination`. Allows to force move.
           
           \param source Source matrix, whose contents are transposed copied.
           \param destination Destination matrix, where the contents of the transposed source matrix
                              are stored.
           \param ittm Indicates whether a copy or move (if possible) should be executed.
         */
        inline void transpose(const implementation::expressions::expr<DestOp, DestData> & destination,
                              const implementation::expressions::expr<SourceOp, SourceData> & source, helper::BoolToType<move> ittm);
        
        template<template<typename SourceData> class SourceOp, typename SourceData, template<typename DestData> class DestOp, typename DestData>
        /**
           \brief Transposes the content of matrix `source` and stores the result in the matrix
                  `destination`.
           
           \param source Source matrix, whose contents are transposed copied.
           \param destination Destination matrix, where the contents of the transposed source matrix
                              are stored.
         */
        inline void transpose(const implementation::expressions::expr<DestOp, DestData> & destination,
                              const implementation::expressions::expr<SourceOp, SourceData> & source);
        
        template<typename SourceT, int SourceRows, int SourceCols, typename SourceST, bool SourceMO, template<typename DestData> class DestOp, typename DestData>
        /**
           \brief Transposes the content of matrix `source` and stores the result in the matrix
                  `destination`.
           
           \param source Source matrix, whose contents are transposed copied.
           \param destination Destination matrix, where the contents of the transposed source matrix
                              are stored.
         */
        inline void transpose(const implementation::expressions::expr<DestOp, DestData> & destination,
                              const base_matrix<SourceT, SourceRows, SourceCols, SourceST, SourceMO> & source);
        
        template<template<typename SourceData> class SourceOp, typename SourceData, typename DestT, int DestRows, int DestCols, typename DestST, bool DestMO>
        /**
           \brief Transposes the content of matrix `source` and stores the result in the matrix
                  `destination`.
           
           \param source Source matrix, whose contents are transposed copied.
           \param destination Destination matrix, where the contents of the transposed source matrix
                              are stored.
         */
        inline void transpose(base_matrix<DestT, DestRows, DestCols, DestST, DestMO> & destination,
                              const implementation::expressions::expr<SourceOp, SourceData> & source);
        
        template<typename SourceT, int SourceRows, int SourceCols, typename SourceST, bool SourceMO,
                 typename DestT, int DestRows, int DestCols, typename DestST, bool DestMO>
        /**
           \brief Transposes the content of matrix `source` and stores the result in the matrix
                  `destination`.
           
           \param source Source matrix, whose contents are transposed copied.
           \param destination Destination matrix, where the contents of the transposed source matrix
                              are stored.
         */
        inline void transpose(base_matrix<DestT, DestRows, DestCols, DestST, DestMO> & destination,
                              const base_matrix<SourceT, SourceRows, SourceCols, SourceST, SourceMO> & source);
        
#if __cplusplus >= 201103L
        template<typename SourceT, int SourceRows, int SourceCols, typename SourceST, bool SourceMO, template<typename DestData> class DestOp, typename DestData>
        /**
           \brief Transposes the content of matrix `source` and stores the result in the matrix
                  `destination`.
           
           \param source Source matrix, whose contents are transposed copied.
           \param destination Destination matrix, where the contents of the transposed source matrix
                              are stored.
         */
        inline void transpose(const implementation::expressions::expr<DestOp, DestData> & destination,
                              base_matrix<SourceT, SourceRows, SourceCols, SourceST, SourceMO> && source);
        
        template<typename SourceT, int SourceRows, int SourceCols, typename SourceST, bool SourceMO,
                 typename DestT, int DestRows, int DestCols, typename DestST, bool DestMO>
        /**
           \brief Transposes the content of matrix `source` and stores the result in the matrix
                  `destination`.
           
           \param source Source matrix, whose contents are transposed copied.
           \param destination Destination matrix, where the contents of the transposed source matrix
                              are stored.
         */
        inline void transpose(base_matrix<DestT, DestRows, DestCols, DestST, DestMO> & destination,
                              base_matrix<SourceT, SourceRows, SourceCols, SourceST, SourceMO> && source);
#endif
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        //// row_count_storage and col_count_storage
        ////////////////////////////////////////////////////////////////////////////////////////////////////

        namespace implementation
        {
            template<int Rows>
            class row_count_storage
            {
            public:
                inline row_count_storage() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                }
                
                inline row_count_storage(size_type) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                }
                
                template<int R>
                inline row_count_storage(const row_count_storage<R> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    assert(r.rows() == Rows);
                }
                
                static inline void set_rows(size_type r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    assert(r == Rows);
                }
                
                static inline size_type rows() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    return Rows;
                }
                
                template<int R>
                static void swap(row_count_storage<R> & other) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    other.set_rows(Rows);
                }
            };
            
            template<>
            class row_count_storage<Flexible>
            {
            protected:
                size_type d_rows;
                
            public:
                inline row_count_storage() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    : d_rows(0)
                {
                }
                
                inline row_count_storage(size_type rows) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    : d_rows(rows)
                {
                }
                
                template<int R>
                inline row_count_storage(const row_count_storage<R> & r) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    : d_rows(r.rows())
                {
                }
                
                inline void set_rows(size_type rows) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    d_rows = rows;
                }
                
                inline size_type rows() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    return d_rows;
                }
                
                template<int R>
                void swap(row_count_storage<R> & other) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    unsigned r = other.rows();
                    other.set_rows(d_rows);
                    d_rows = r;
                }
            };
            
            template<int Cols>
            class col_count_storage
            {
            public:
                inline col_count_storage() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                }
                
                inline col_count_storage(size_type) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                }
                
                template<int C>
                inline col_count_storage(const col_count_storage<C> & c) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    assert(c.cols() == Cols);
                }
                
                static inline void set_cols(size_type c) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    assert(c == Cols);
                }
                
                static inline size_type cols() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    return Cols;
                }
                
                template<int C>
                static void swap(col_count_storage<C> & other) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    other.set_cols(Cols);
                }
            };
            
            template<>
            class col_count_storage<Flexible>
            {
            protected:
                size_type d_cols;
                
            public:
                inline col_count_storage() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    : d_cols(0)
                {
                }
                
                inline col_count_storage(size_type cols) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    : d_cols(cols)
                {
                }
                
                template<int C>
                inline col_count_storage(const col_count_storage<C> & c) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    : d_cols(c.cols())
                {
                }
                
                inline void set_cols(size_type cols) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    d_cols = cols;
                }
                
                inline size_type cols() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    return d_cols;
                }
                
                template<int C>
                void swap(col_count_storage<C> & other) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    unsigned r = other.cols();
                    other.set_cols(d_cols);
                    d_cols = r;
                }
            };
        }
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        //// base_matrix
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        
        template<typename T, int Rows = Flexible, int Cols = Flexible, typename StorageTraits = storage_traits<T>, bool MathObject = false>
        /**
           \brief Represents a matrix with coefficients in `T`.
           
           The main template. Stores a matrix with coefficients in `T` with possible compile-time
           row and column numbers enforcements (controlled with `Rows` and `Cols`) and storage
           management (`StorageTraits`). Can store both math objects and non-math objects (via
           parameter `MathObject`).
           
           \tparam T              The coefficient type.
           \tparam Rows           If set to `plll::linalg::Flexible`, allows to flexibly choose
                                  and change the number of rows. Otherwise, enforces fixed number
                                  of rows during compile time. Default value is
                                  `plll::linalg::Flexible`.
           \tparam Cols           If set to `plll::linalg::Flexible`, allows to flexibly choose
                                  and change the number of columns. Otherwise, enforces fixed number
                                  of columns during compile time. Default value is
                                  `plll::linalg::Flexible`.
           \tparam StorageTraits  Allows to determine how the matrix's contents are stored. The
                                  default value is `plll::linalg::storage_traits<T>`, which should
                                  usually not be changed.
           \tparam MathObject     Allows to determine whether the matrix is a math object (`true`)
                                  and thus operations can be performed on it, or whether it is a
                                  non-math object (`false`) which allows no operations except
                                  reading and storing coefficients and resizing.

           There exists more specialized versions of `base_matrix<>` for row and colum vectors
           (`plll::linalg::base_rowvector<>` and `plll::linalg::base_colvector<>`), as well as for
           math objects (`plll::linalg::math_matrix<>`, `plll::linalg::math_rowvector<>` and
           `plll::linalg::math_colvector<>`).
         */
        class base_matrix : private implementation::row_count_storage<Rows>, private implementation::col_count_storage<Cols>
        {
        public:
            enum { has_direct_access = true };
            
#ifdef PLLL_INTERNAL_NO_TEMPLATE_FRIENDS
        public:
#else
            template<typename TT, int R, int C, typename ST, bool MO>
            friend void swap(base_matrix<TT, R, C, ST, MO> &, base_matrix<TT, R, C, ST, MO> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            
            template<typename TT, int R1, int C1, int R2, int C2, typename ST, bool MO1, bool MO2>
            friend void swap(base_matrix<TT, R1, C1, ST, MO1> &, base_matrix<TT, R2, C2, ST, MO2> &) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE;
            
            template<typename TT, int RR, int CC, typename SSTT, bool MMOO>
            friend class base_matrix;
            
        private:
            inline static void set_zero(typename StorageTraits::ref_type e, helper::BoolToType<true>)
                PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(setZero(e)))
            {
                setZero(e);
            }
            
            inline static void set_zero(typename StorageTraits::ref_type, helper::BoolToType<false>)
                PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
            }
            
            inline static void set_zero(typename StorageTraits::ref_type e)
                PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(set_zero(e, helper::BoolToType<MathObject>())))
            {
                set_zero(e, helper::BoolToType<MathObject>());
            }
            
#endif
            typename StorageTraits::pointer_type d_data;
            
            template<typename MatrixType>
            inline void create_from(const MatrixType & mat, helper::BoolToType<true>)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::create_from(" << getAddress(*this) << "; " << getAddress(mat) << ") [1]");
                d_data = StorageTraits::clone(mat.data(), implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
            }
            
            template<typename MatrixType>
            inline void create_from(const MatrixType & mat, helper::BoolToType<false>)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::create_from(" << getAddress(*this) << "; " << getAddress(mat) << ") [2]");
                d_data = StorageTraits::alloc_dontconstruct(implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                if (mat.get_coeff_alwayszero())
                {
                    size_type s = implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols();
                    for (size_type index = 0; index < s; ++index)
                    {
                        StorageTraits::construct(d_data[index]);
                        set_zero(d_data[index]);
                    }
                }
                else
                {
                    typename MatrixType::ConstEnumerator e = mat.enumerate();
                    for (size_type index = 0; e.has_current(); e.next(), ++index)
                        if (implementation::MatrixInfo<MatrixType>::coeffs_are_simple_expressions)
                            StorageTraits::copy_construct(d_data[index], e.current());
                        else
                        {
                            typename MatrixType::ConstEnumerator::GetCoeffSteps_Type coeff = e.get_current_steps();
                            StorageTraits::copy_construct(d_data[index], coeff.step1());
                            coeff.step2(d_data[index]);
                        }
                }
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
            }
            
#if __cplusplus >= 201103L
            template<typename MatrixType>
            inline void move_from(MatrixType && mat, helper::BoolToType<true>)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::create_from(" << getAddress(*this) << "; " << getAddress(mat) << ") [1&&]");
                d_data = StorageTraits::clone_move(mat.data(), implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
            }
            
            template<typename MatrixType>
            inline void move_from(MatrixType && mat, helper::BoolToType<false>)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::create_from(" << getAddress(*this) << "; " << getAddress(mat) << ") [2&&]");
                d_data = StorageTraits::alloc_dontconstruct(implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                if (mat.get_coeff_alwayszero())
                {
                    size_type s = implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols();
                    for (size_type index = 0; index < s; ++index)
                    {
                        StorageTraits::construct(d_data[index]);
                        set_zero(d_data[index]);
                    }
                }
                else
                {
                    typename MatrixType::ConstEnumerator e = mat.enumerate();
                    for (size_type index = 0; e.has_current(); e.next(), ++index)
                        if (implementation::MatrixInfo<MatrixType>::coeffs_are_simple_expressions)
                            StorageTraits::move_construct(d_data[index], std::move(e.current()));
                        else
                        {
                            typename MatrixType::ConstEnumerator::GetCoeffSteps_Type coeff = e.get_current_steps();
                            StorageTraits::copy_construct(d_data[index], coeff.step1());
                            coeff.step2(d_data[index]);
                        }
                }
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
            }
#endif
            
        public:
            /**
               \brief Return type for the method `data() const`.
             */
            typedef typename StorageTraits::constpointer_type data_pointer_type;
            /**
               \brief Return type for the method `data(size_type) const`.
             */
            typedef typename StorageTraits::constref_type data_ref_type;
            
            /**
               \brief Creates a new matrix. The dimensions are set to zero if no compile-time dimensions
                      were specified. The elements are default-constructed.
             */
            base_matrix()
                : d_data(StorageTraits::alloc(implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols()))
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::base_matrix(" << getAddress(*this) << ")");
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
            }
            
            /**
               \brief Creates a copy of the matrix `mat`.
               
               \param mat The matrix to copy from.
             */
            base_matrix(const base_matrix<T, Rows, Cols, StorageTraits, MathObject> & mat)
                : implementation::row_count_storage<Rows>(mat), implementation::col_count_storage<Cols>(mat),
                  d_data(StorageTraits::clone(mat.d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols()))
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::base_matrix(" << getAddress(*this) << "; " << getAddress(mat) << ") [0]");
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
            }
            
            /**
               \brief Creates a copy of the matrix `mat`.
               
               \param mat The matrix to copy from.
             */
            template<int Rs, int Cs, bool MO>
            base_matrix(const base_matrix<T, Rs, Cs, StorageTraits, MO> & mat)
                : implementation::row_count_storage<Rows>(mat), implementation::col_count_storage<Cols>(mat),
                  d_data(StorageTraits::clone(mat.d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols()))
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::base_matrix(" << getAddress(*this) << "; " << getAddress(mat) << ") [1]");
                PLLL_INTERNAL_STATIC_CHECK((Rows < 0) || (Rs < 0) || (Rows == Rs), NotCompatible);
                PLLL_INTERNAL_STATIC_CHECK((Cols < 0) || (Cs < 0) || (Cols == Cs), NotCompatible);
                if (Rows >= 0) assert(static_cast<size_type>(Rows) == implementation::row_count_storage<Rows>::rows());
                if (Cols >= 0) assert(static_cast<size_type>(Cols) == implementation::col_count_storage<Cols>::cols());
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
            }
            
#if __cplusplus >= 201103L
            /**
               \brief Moves the content of the matrix `mat` to the current matrix.
               
               \param mat The matrix to move from.
             */
            template<bool MO>
            inline base_matrix(base_matrix<T, Rows, Cols, StorageTraits, MO> && mat) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : implementation::row_count_storage<Rows>(mat), implementation::col_count_storage<Cols>(mat), d_data(mat.d_data)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::base_matrix(" << getAddress(*this) << "; " << getAddress(mat) << ") [1&&a]");
                mat.set_rows(0);
                mat.set_cols(0);
                mat.d_data = NULL;
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
            }
            
            /**
               \brief Moves the content of the matrix `mat` to the current matrix.
               
               \param mat The matrix to move from.
             */
            template<int Rs, int Cs, bool MO>
            base_matrix(base_matrix<T, Rs, Cs, StorageTraits, MO> && mat) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : implementation::row_count_storage<Rows>(mat), implementation::col_count_storage<Cols>(mat), d_data(mat.d_data)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::base_matrix(" << getAddress(*this) << "; " << getAddress(mat) << ") [1&&]");
                PLLL_INTERNAL_STATIC_CHECK((Rows < 0) || (Rs < 0) || (Rows == Rs), NotCompatible);
                PLLL_INTERNAL_STATIC_CHECK((Cols < 0) || (Cs < 0) || (Cols == Cs), NotCompatible);
                if (Rows >= 0) assert(static_cast<size_type>(Rows) == implementation::row_count_storage<Rows>::rows());
                if (Cols >= 0) assert(static_cast<size_type>(Cols) == implementation::col_count_storage<Cols>::cols());
                mat.set_rows(0);
                mat.set_cols(0);
                mat.d_data = NULL;
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
            }
#endif
            
            /**
               \brief Creates a copy of the matrix `mat`.
               
               \param mat The matrix to copy from.
             */
            template<typename S, int Rs, int Cs, typename ST, bool MO>
            base_matrix(const base_matrix<S, Rs, Cs, ST, MO> & mat)
                : implementation::row_count_storage<Rows>(mat), implementation::col_count_storage<Cols>(mat),
                  d_data(StorageTraits::clone(mat.d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols()))
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::base_matrix(" << getAddress(*this) << "; " << getAddress(mat) << ") [2]");
                PLLL_INTERNAL_STATIC_CHECK((Rows < 0) || (Rs < 0) || (Rows == Rs), NotCompatible);
                PLLL_INTERNAL_STATIC_CHECK((Cols < 0) || (Cs < 0) || (Cols == Cs), NotCompatible);
                if (Rows >= 0) assert(static_cast<size_type>(Rows) == implementation::row_count_storage<Rows>::rows());
                if (Cols >= 0) assert(static_cast<size_type>(Cols) == implementation::col_count_storage<Cols>::cols());
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
            }
            
#if __cplusplus >= 201103L
            /**
               \brief Moves the content of the matrix `mat` to the current matrix.
               
               \param mat The matrix to move from.
             */
            template<typename S, int Rs, int Cs, typename ST, bool MO>
            base_matrix(base_matrix<S, Rs, Cs, ST, MO> && mat)
                : implementation::row_count_storage<Rows>(mat), implementation::col_count_storage<Cols>(mat),
                  d_data(StorageTraits::clone_move(mat.d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols()))
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::base_matrix(" << getAddress(*this) << "; " << getAddress(mat) << ") [2&&]");
                PLLL_INTERNAL_STATIC_CHECK((Rows < 0) || (Rs < 0) || (Rows == Rs), NotCompatible);
                PLLL_INTERNAL_STATIC_CHECK((Cols < 0) || (Cs < 0) || (Cols == Cs), NotCompatible);
                if (Rows >= 0) assert(static_cast<size_type>(Rows) == implementation::row_count_storage<Rows>::rows());
                if (Cols >= 0) assert(static_cast<size_type>(Cols) == implementation::col_count_storage<Cols>::cols());
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
            }
#endif
            
            /**
               \brief Creates a copy of the matrix `mat`.
               
               \param mat The matrix to copy from.
             */
            template<template<typename> class Op, typename Data>
            base_matrix(const implementation::expressions::expr<Op, Data> & mat)
                : implementation::row_count_storage<Rows>(mat.rows()), implementation::col_count_storage<Cols>(mat.cols()), d_data(NULL)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::base_matrix(" << getAddress(*this) << "; " << getAddress(mat) << ") [3]");
                PLLL_INTERNAL_STATIC_CHECK((Rows < 0) || (implementation::MatrixInfo<implementation::expressions::expr<Op, Data> >::rows < 0) ||
                                           (static_cast<size_type>(Rows) == static_cast<size_type>(implementation::MatrixInfo<implementation::expressions::expr<Op, Data> >::rows)), NotCompatible);
                PLLL_INTERNAL_STATIC_CHECK((Cols < 0) || (implementation::MatrixInfo<implementation::expressions::expr<Op, Data> >::cols < 0) ||
                                           (static_cast<size_type>(Cols) == static_cast<size_type>(implementation::MatrixInfo<implementation::expressions::expr<Op, Data> >::cols)), NotCompatible);
                if (Rows >= 0) assert(static_cast<size_type>(Rows) == implementation::row_count_storage<Rows>::rows());
                if (Cols >= 0) assert(static_cast<size_type>(Cols) == implementation::col_count_storage<Cols>::cols());
                create_from(mat, helper::BoolToType<implementation::expressions::expr<Op, Data>::has_direct_access>());
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
            }
            
#if __cplusplus >= 201103L
            /**
               \brief Moves the content of the matrix `mat` to the current matrix.
               
               \param mat The matrix to move from.
             */
            template<template<typename> class Op, typename Data>
            base_matrix(implementation::expressions::expr<Op, Data> && mat)
                : implementation::row_count_storage<Rows>(mat.rows()), implementation::col_count_storage<Cols>(mat.cols()), d_data(NULL)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::base_matrix(" << getAddress(*this) << "; " << getAddress(mat) << ") [3&&]");
                PLLL_INTERNAL_STATIC_CHECK((Rows < 0) || (implementation::MatrixInfo<implementation::expressions::expr<Op, Data> >::rows < 0) ||
                                           (static_cast<size_type>(Rows) == static_cast<size_type>(implementation::MatrixInfo<implementation::expressions::expr<Op, Data> >::rows)), NotCompatible);
                PLLL_INTERNAL_STATIC_CHECK((Cols < 0) || (implementation::MatrixInfo<implementation::expressions::expr<Op, Data> >::cols < 0) ||
                                           (static_cast<size_type>(Cols) == static_cast<size_type>(implementation::MatrixInfo<implementation::expressions::expr<Op, Data> >::cols)), NotCompatible);
                if (Rows >= 0) assert(static_cast<size_type>(Rows) == implementation::row_count_storage<Rows>::rows());
                if (Cols >= 0) assert(static_cast<size_type>(Cols) == implementation::col_count_storage<Cols>::cols());
                if (implementation::MatrixInfo<implementation::expressions::expr<Op, Data> >::can_move_from)
                    move_from(std::move(mat), helper::BoolToType<implementation::expressions::expr<Op, Data>::has_direct_access>());
                else
                    create_from(mat, helper::BoolToType<implementation::expressions::expr<Op, Data>::has_direct_access>());
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
            }
#endif
            
            /**
               \brief Creates a vector with `entries` entries. Its coefficients are
                      default-constructed.
               
               Note that the number of rows or columns must be forced to be 1 during compile time;
               otherwise, this constructor will not compile.
               
               \param entries The number of entries.
             */
            base_matrix(size_type entries)
                : implementation::row_count_storage<Rows>(entries), implementation::col_count_storage<Cols>(entries),
                  d_data(StorageTraits::alloc(implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols()))
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::base_matrix(" << getAddress(*this) << "; " << entries << ") [4]");
                PLLL_INTERNAL_STATIC_CHECK((Cols == 1) || (Rows == 1), ConstructorOnlyAvailableForVectors);
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
            }
            
            /**
               \brief Creates a vector with `entries` entries. Its coefficients are
                      default-constructed.
               
               Note that the number of rows or columns must be forced to be 1 during compile time;
               otherwise, this constructor will not compile.
               
               \param entries The number of entries.
             */
            base_matrix(int entries)
                : implementation::row_count_storage<Rows>(entries), implementation::col_count_storage<Cols>(entries),
                  d_data(StorageTraits::alloc(implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols()))
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::base_matrix(" << getAddress(*this) << "; " << entries << ") [4b]");
                PLLL_INTERNAL_STATIC_CHECK((Cols == 1) || (Rows == 1), ConstructorOnlyAvailableForVectors);
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
            }
            
            /**
               \brief Creates a vector with `entries` entries. Its coefficients are copy-constructed
                      from `i.ref()`.
               
               Note that the number of rows or columns must be forced to be 1 during compile time;
               otherwise, this constructor will not compile.
               
               \param entries The number of entries.
               \param i       The element from which to copy-construct the coefficients.
             */
            template<typename S, bool def>
            base_matrix(size_type entries, const implementation::Initialize_Impl<S, def> & i)
                : implementation::row_count_storage<Rows>(entries), implementation::col_count_storage<Cols>(entries),
                  d_data(StorageTraits::alloc(i.ref(), implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols()))
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::base_matrix(" << getAddress(*this) << "; " << entries << " " << getAddress(i.ref()) << ") [5]");
                PLLL_INTERNAL_STATIC_CHECK((Cols == 1) || (Rows == 1), ConstructorOnlyAvailableForVectors);
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
            }
            
            /**
               \brief Creates a matrix with `rows` times `cols` entries. Its coefficients are
                      default-constructed.
               
               \param rows The number of rows of the matrix.
               \param cols The number of columns of the matrix.
             */
            base_matrix(size_type rows, size_type cols)
                : implementation::row_count_storage<Rows>(rows), implementation::col_count_storage<Cols>(cols),
                  d_data(StorageTraits::alloc(implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols()))
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::base_matrix(" << getAddress(*this) << "; " << rows << " " << cols << ") [6]");
                if (Rows >= 0) assert(rows == static_cast<size_type>(Rows));
                if (Cols >= 0) assert(cols == static_cast<size_type>(Cols));
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
            }
            
            /**
               \brief Creates a matrix with `rows` times `cols` entries. Its coefficients are
                      copy-constructed from `i.ref()`.
               
               \param rows The number of rows of the matrix.
               \param cols The number of columns of the matrix.
               \param i       The element from which to copy-construct the coefficients.
             */
            template<typename S, bool def>
            base_matrix(size_type rows, size_type cols, const implementation::Initialize_Impl<S, def> & i)
                : implementation::row_count_storage<Rows>(rows), implementation::col_count_storage<Cols>(cols),
                  d_data(StorageTraits::alloc(i.ref(), implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols()))
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::base_matrix(" << getAddress(*this) << "; " << rows << " " << cols << " " << getAddress(i.ref()) << ") [7]");
                if (Rows >= 0) assert(rows == static_cast<size_type>(Rows));
                if (Cols >= 0) assert(cols == static_cast<size_type>(Cols));
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
            }

            /**
               \brief Releases the memory allocated for the matrix.
             */
            inline ~base_matrix() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::~base_matrix(" << getAddress(*this) << ")");
                StorageTraits::free(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
            }
            
            /**
               \brief Copies the contents of the matrix `m` into `*this`.
               
               \param m The matrix whose contents to copy into the current matrix.
             */
            base_matrix<T, Rows, Cols, StorageTraits, MathObject> & operator = (const base_matrix<T, Rows, Cols, StorageTraits, MathObject> & m)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::operator = (" << getAddress(*this) << "; " << getAddress(m) << ") [0]");
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                assign(*this, m);
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                return *this;
            }
            
            /**
               \brief Copies the contents of the matrix `m` into `*this`.
               
               \param m The matrix whose contents to copy into the current matrix.
             */
            template<typename T_, int Rows_, int Cols_, typename StorageTraits_, bool MathObject_>
            base_matrix<T, Rows, Cols, StorageTraits, MathObject> & operator = (const base_matrix<T_, Rows_, Cols_, StorageTraits_, MathObject_> & m)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::operator = (" << getAddress(*this) << "; " << getAddress(m) << ") [1]");
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                assign(*this, m);
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                return *this;
            }
            
            /**
               \brief Copies the contents of the matrix `m` into `*this`.
               
               \param m The matrix whose contents to copy into the current matrix.
             */
            template<template<typename> class Op, typename Data>
            base_matrix<T, Rows, Cols, StorageTraits, MathObject> & operator = (const implementation::expressions::expr<Op, Data> & m)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::operator = (" << getAddress(*this) << "; " << getAddress(m) << ") [2]");
                PLLL_INTERNAL_STATIC_CHECK((implementation::MatrixInfo<implementation::expressions::expr<Op, Data> >::is_matrix), RequiresMatrixType);
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                assign(*this, m);
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                return *this;
            }
            
#if __cplusplus >= 201103L
            /**
               \brief Moves the contents of the matrix `m` into `*this`.
               
               \param m The matrix whose contents to move into the current matrix.
             */
            base_matrix<T, Rows, Cols, StorageTraits, MathObject> & operator = (base_matrix<T, Rows, Cols, StorageTraits, MathObject> && m)
                PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::operator = (" << getAddress(*this) << "; " << getAddress(m) << ") [0&&]");
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                StorageTraits::free(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                implementation::row_count_storage<Rows>::set_rows(m.rows());
                implementation::col_count_storage<Cols>::set_cols(m.cols());
                m.implementation::row_count_storage<Rows>::set_rows(0);
                m.implementation::col_count_storage<Cols>::set_cols(0);
                d_data = m.d_data;
                m.d_data = NULL;
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                return *this;
            }
            
            /**
               \brief Moves the contents of the matrix `m` into `*this`.
               
               \param m The matrix whose contents to move into the current matrix.
             */
            template<typename T_, int Rows_, int Cols_, typename StorageTraits_, bool MathObject_>
            base_matrix<T, Rows, Cols, StorageTraits, MathObject> & operator = (base_matrix<T_, Rows_, Cols_, StorageTraits_, MathObject_> && m)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::operator = (" << getAddress(*this) << "; " << getAddress(m) << ") [1&&]");
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                assign(*this, std::move(m));
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                return *this;
            }
            
            /**
               \brief Moves the contents of the matrix `m` into `*this`.
               
               \param m The matrix whose contents to move into the current matrix.
             */
            template<template<typename> class Op, typename Data>
            base_matrix<T, Rows, Cols, StorageTraits, MathObject> & operator = (implementation::expressions::expr<Op, Data> && m)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::operator = (" << getAddress(*this) << "; " << getAddress(m) << ") [2&&]");
                PLLL_INTERNAL_STATIC_CHECK((implementation::MatrixInfo<implementation::expressions::expr<Op, Data> >::is_matrix), RequiresMatrixType);
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                assign(*this, std::move(m));
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                return *this;
            }
#endif
            
            /**
               \brief Retrieves the coefficient `(i, j)` and stores its value into `result`.
               
               \param result Where to store the value of coefficient `(i, j)`.
               \param i      The row of the coefficient.
               \param j      The column of the coefficient.
             */
            template<typename ResultType>
            inline void get_coeff(ResultType & result, size_type i, size_type j) const
                PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(result = d_data[i * implementation::col_count_storage<Cols>::cols() + j]))
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::get_coeff(" << getAddress(*this) << "; " << getAddress(result) << " " << i << " " << j << ") const");
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                assert(i < implementation::row_count_storage<Rows>::rows());
                assert(j < implementation::col_count_storage<Cols>::cols());
                result = d_data[i * implementation::col_count_storage<Cols>::cols() + j];
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
            }
            
            /**
               \brief Queries whether all coefficients are definitely zero.
               
               \return `true` if all coefficients will be zero, and `false` if we are not sure.
             */
            inline bool get_coeff_alwayszero() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::get_coeff_alwayszero(" << getAddress(*this) << ") const");
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                return false;
            }
            
            /**
               \brief Allows to retrieve a coefficient of the matrix in two steps.
             */
            class GetCoeffSteps_Type
            {
                friend class base_matrix<T, Rows, Cols, StorageTraits, MathObject>;
                
            protected:
                typename StorageTraits::constref_type d_value;
                
                inline GetCoeffSteps_Type(typename StorageTraits::constref_type value) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    : d_value(value)
                {
                }
                
            public:
                /** \brief The return type for `step1()`. */
                typedef typename StorageTraits::constref_type Step1_Type;
                
                /** \brief Retrieves the coefficient as a simple expression. The result should be
                           stored in a variable of the final coefficient type.
                    
                    \return A simple expression to begin the construction of the coefficient.
                */
                inline Step1_Type step1() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    return d_value;
                }
                
                /** \brief Modifies the result returned by `step1()` to contain the full
                           coefficient.
                    
                    \param result Where the result of `step1()` was stored.
                */
                template<typename ResultType>
                inline void step2(ResultType & result) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                }
            };
            
            /**
               \brief Returns a object of type `GetCoeffSteps_Type` which constructs the coefficient
                      at `(i, j)`.
               
               This should _only_ be called if `get_coeff_alwayszero()` returns `false`. Otherwise,
               the behavior is undefined.
               
               \param i The row of the coefficient.
               \param j The column of the coefficient.
               \return A `GetCoeffSteps_Type` object for coefficient `(i, j)`.
             */
            inline GetCoeffSteps_Type get_coeff_steps(size_type i, size_type j) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::get_coeff_steps(" << getAddress(*this) << "; " << i << " " << j << ") const");
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                assert(i < implementation::row_count_storage<Rows>::rows());
                assert(j < implementation::col_count_storage<Cols>::cols());
                return GetCoeffSteps_Type(d_data[i * implementation::col_count_storage<Cols>::cols() + j]);
            }
            
            /**
               \brief Retrieves and returns the coefficient `(i, j)`.
               
               \param i The row of the coefficient.
               \param j The column of the coefficient.
               \return The value of coefficient `(i, j)`.
             */
            inline typename StorageTraits::constref_type operator () (size_type i, size_type j) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::operator () (" << getAddress(*this) << "; " << i << " " << j << ") const");
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                assert(i < implementation::row_count_storage<Rows>::rows());
                assert(j < implementation::col_count_storage<Cols>::cols());
                return d_data[i * implementation::col_count_storage<Cols>::cols() + j];
            }
            
            /**
               \brief Retrieves and returns the coefficient `(i, j)`.
               
               \param i The row of the coefficient.
               \param j The column of the coefficient.
               \return The value of coefficient `(i, j)`.
             */
            inline typename StorageTraits::ref_type operator () (size_type i, size_type j) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::operator () (" << getAddress(*this) << "; " << i << " " << j << ")");
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                assert(i < implementation::row_count_storage<Rows>::rows());
                assert(j < implementation::col_count_storage<Cols>::cols());
                return d_data[i * implementation::col_count_storage<Cols>::cols() + j];
            }
            
            /**
               \brief Retrieves and returns the `i`-th entry of the vector.
               
               Note that the number of rows or columns must be forced to be 1 during compile time;
               otherwise, this method will not compile.
               
               \param i The entry of the vector to retrieve.
               \return The value of the `i`-th entry.
             */
            inline typename StorageTraits::constref_type operator [] (size_type i) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::operator [] (" << getAddress(*this) << "; " << i << ") const");
                PLLL_INTERNAL_STATIC_CHECK((Cols == 1) || (Rows == 1), OnlyDefinedForVectors);
                StorageTraits::check(d_data, (Cols == 1) ? implementation::row_count_storage<Rows>::rows() : implementation::col_count_storage<Cols>::cols());
                assert((Cols == 1) ? (i < implementation::row_count_storage<Rows>::rows()) : (i < implementation::col_count_storage<Cols>::cols()));
                return d_data[i];
            }
            
            /**
               \brief Retrieves and returns the `i`-th entry of the vector.
               
               Note that the number of rows or columns must be forced to be 1 during compile time;
               otherwise, this method will not compile.
               
               \param i The entry of the vector to retrieve.
               \return The value of the `i`-th entry.
             */
            inline typename StorageTraits::ref_type operator [] (size_type i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::operator [] (" << getAddress(*this) << "; " << i << ")");
                PLLL_INTERNAL_STATIC_CHECK((Cols == 1) || (Rows == 1), OnlyDefinedForVectors);
                StorageTraits::check(d_data, (Cols == 1) ? implementation::row_count_storage<Rows>::rows() : implementation::col_count_storage<Cols>::cols());
                assert((Cols == 1) ? (i < implementation::row_count_storage<Rows>::rows()) : (i < implementation::col_count_storage<Cols>::cols()));
                return d_data[i];
            }
            
            /**
               \brief Returns a pointer to the raw data of the matrix.
               
               \return A pointer to a block of memory containing `size()` elements of type `T`.
             */
            inline typename StorageTraits::constpointer_type data() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::data(" << getAddress(*this) << ") const");
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                return d_data;
            }
            
            /**
               \brief Returns a pointer to the raw data of the matrix.
               
               \return A pointer to a block of memory containing `size()` elements of type `T`.
             */
            inline typename StorageTraits::pointer_type data() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::data(" << getAddress(*this) << ")");
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                return d_data;
            }
            
            /**
               \brief Returns a (raw) reference to the `i`-th entry.
               
               \param i The index of the entry to retrieve.
               \return A (raw) refrence into a block of memory containing `size()` elements of type `T`.
             */
            inline typename StorageTraits::constref_type data(size_type i) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::data(" << getAddress(*this) << "; " << i << ") const");
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                return d_data[i];
            }
            
            /**
               \brief Returns a (raw) reference to the `i`-th entry.
               
               \param i The index of the entry to retrieve.
               \return A (raw) refrence into a block of memory containing `size()` elements of type `T`.
             */
            inline typename StorageTraits::ref_type data(size_type i) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::data(" << getAddress(*this) << "; " << i << ")");
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                return d_data[i];
            }
            
            /**
               \brief Returns the number of rows of the matrix.
               
               \return The number of rows.
             */
            inline size_type rows() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::rows(" << getAddress(*this) << ") const");
                return implementation::row_count_storage<Rows>::rows();
            }
            
            /**
               \brief Returns the number of columns of the matrix.
               
               \return The number of columns.
             */
            inline size_type cols() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::cols(" << getAddress(*this) << ") const");
                return implementation::col_count_storage<Cols>::cols();
            }
            
            /**
               \brief Returns the number of entries of the matrix. This equals `rows() * cols()`.
               
               \return The number of entries of the matrix.
             */
            inline size_type size() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::size(" << getAddress(*this) << ") const");
                return implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols();
            }

            /**
               \brief Returns a submatrix of `SRows` rows and `SCols` columns starting at `(r_ofs,
                      c_ofs)`.

               \param r_ofs The row to begin the submatrix in.
               \param c_ofs The column to begin the submatrix in.
               \return An object representing a submatrix of `*this` of a number of rows and columns
                       fixed at compile-time.
             */
            template<unsigned SRows, unsigned SCols>
            inline typename implementation::expressions::expr<implementation::expressions::sub<SRows, SCols>::template operation_generic,
                                                              implementation::expressions::MatrixWrapper<base_matrix<T, Rows, Cols, StorageTraits, MathObject> > >
            block(size_type r_ofs, size_type c_ofs)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::block<" << SRows << " " << SCols << ">(" << getAddress(*this) << "; " << r_ofs << " " << c_ofs << ")");
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                return typename implementation::expressions::expr<implementation::expressions::sub<SRows, SCols>::template operation_generic, implementation::expressions::MatrixWrapper<base_matrix<T, Rows, Cols, StorageTraits, MathObject> > >
                    (typename implementation::expressions::sub<SRows, SCols>::template operation_generic<implementation::expressions::MatrixWrapper<base_matrix<T, Rows, Cols, StorageTraits, MathObject> > >(SRows, SCols,
                                                                                                                                                                              r_ofs, c_ofs, *this),
                     implementation::expressions::make_matrix_wrapper(*this));
            }
            
            /**
               \brief Returns a submatrix of `SRows` rows and `SCols` columns starting at `(r_ofs,
                      c_ofs)`.

               The number of rows and columns must be specified at compile-time.

               \param r_ofs The row to begin the submatrix in.
               \param c_ofs The column to begin the submatrix in.
               \return An object representing a submatrix of `*this` of a number of rows and columns
                       fixed at compile-time.
             */
            template<unsigned SRows, unsigned SCols>
            inline typename implementation::expressions::expr<implementation::expressions::sub<SRows, SCols>::template operation_generic,
                                                              implementation::expressions::ConstMatrixWrapper<base_matrix<T, Rows, Cols, StorageTraits, MathObject> > >
            block(size_type r_ofs, size_type c_ofs) const
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::block<" << SRows << " " << SCols << ">(" << getAddress(*this) << "; " << r_ofs << " " << c_ofs << ") const");
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                return typename implementation::expressions::expr<implementation::expressions::sub<SRows, SCols>::template operation_generic, implementation::expressions::ConstMatrixWrapper<base_matrix<T, Rows, Cols, StorageTraits, MathObject> > >
                    (typename implementation::expressions::sub<SRows, SCols>::template operation_generic<implementation::expressions::ConstMatrixWrapper<base_matrix<T, Rows, Cols, StorageTraits, MathObject> > >(SRows, SCols,
                                                                                                                                                                                   r_ofs, c_ofs, *this),
                     implementation::expressions::make_matrix_wrapper(*this));
            }
            
            /**
               \brief Returns a submatrix of `rows` rows and `cols` columns starting at `(r_ofs,
                      c_ofs)`.

               The number of rows and columns must be specified at runtime.
               
               \param r_ofs The row to begin the submatrix in.
               \param c_ofs The column to begin the submatrix in.
               \param rows The number of rows of the submatrix.
               \param cols The number of columns of the submatrix.
               \return An object representing a submatrix of `*this`.
             */
            inline typename implementation::expressions::expr<implementation::expressions::sub<Flexible, Flexible>::template operation_generic,
                                                              implementation::expressions::MatrixWrapper<base_matrix<T, Rows, Cols, StorageTraits, MathObject> > >
            block(size_type r_ofs, size_type c_ofs, size_type rows, size_type cols)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::block(" << getAddress(*this) << "; " << r_ofs << " " << c_ofs << " " << rows << " " << cols << ")");
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                return typename implementation::expressions::expr<implementation::expressions::sub<Flexible, Flexible>::template operation_generic, implementation::expressions::MatrixWrapper<base_matrix<T, Rows, Cols, StorageTraits, MathObject> > >
                    (typename implementation::expressions::sub<Flexible, Flexible>::template operation_generic<implementation::expressions::MatrixWrapper<base_matrix<T, Rows, Cols, StorageTraits, MathObject> > >(rows, cols,
                                                                                                                                                                                    r_ofs, c_ofs, *this),
                     implementation::expressions::make_matrix_wrapper(*this));
            }
            
            /**
               \brief Returns a submatrix of `rows` rows and `cols` columns starting at `(r_ofs,
                      c_ofs)`.

               The number of rows and columns must be specified at runtime.
               
               \param r_ofs The row to begin the submatrix in.
               \param c_ofs The column to begin the submatrix in.
               \param rows The number of rows of the submatrix.
               \param cols The number of columns of the submatrix.
               \return An object representing a submatrix of `*this`.
             */
            inline typename implementation::expressions::expr<implementation::expressions::sub<Flexible, Flexible>::template operation_generic,
                                                              implementation::expressions::ConstMatrixWrapper<base_matrix<T, Rows, Cols, StorageTraits, MathObject> > >
            block(size_type r_ofs, size_type c_ofs, size_type rows, size_type cols) const
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::block(" << getAddress(*this) << "; " << r_ofs << " " << c_ofs << " " << rows << " " << cols << ") const");
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                return typename implementation::expressions::expr<implementation::expressions::sub<Flexible, Flexible>::template operation_generic, implementation::expressions::ConstMatrixWrapper<base_matrix<T, Rows, Cols, StorageTraits, MathObject> > >
                    (typename implementation::expressions::sub<Flexible, Flexible>::template operation_generic<implementation::expressions::ConstMatrixWrapper<base_matrix<T, Rows, Cols, StorageTraits, MathObject> > >(rows, cols,
                                                                                                                                                                                         r_ofs, c_ofs, *this),
                     implementation::expressions::make_matrix_wrapper(*this));
            }
            
            /**
               \brief Returns the `row_index`-th row of this matrix.
               
               \param row_index The index of the row.
               \return An object representing a row of `*this`.
             */
            inline typename implementation::expressions::expr<implementation::expressions::sub_1d<Cols>::template operation_row,
                                                              implementation::expressions::MatrixWrapper<base_matrix<T, Rows, Cols, StorageTraits, MathObject> > >
            row(size_type row_index)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::row(" << getAddress(*this) << "; " << row_index << ")");
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                return typename implementation::expressions::expr<implementation::expressions::sub_1d<Cols>::template operation_row, implementation::expressions::MatrixWrapper<base_matrix<T, Rows, Cols, StorageTraits, MathObject> > >
                    (typename implementation::expressions::sub_1d<Cols>::template operation_row<implementation::expressions::MatrixWrapper<base_matrix<T, Rows, Cols, StorageTraits, MathObject> > >(implementation::col_count_storage<Cols>::cols(),
                                                                                                                                                                     row_index, *this),
                     implementation::expressions::make_matrix_wrapper(*this));
            }
            
            /**
               \brief Returns the `row_index`-th row of this matrix.
               
               \param row_index The index of the row.
               \return An object representing a row of `*this`.
             */
            inline typename implementation::expressions::expr<implementation::expressions::sub_1d<Cols>::template operation_row,
                                                              implementation::expressions::ConstMatrixWrapper<base_matrix<T, Rows, Cols, StorageTraits, MathObject> > >
            row(size_type row_index) const
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::row(" << getAddress(*this) << "; " << row_index << ") const");
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                return typename implementation::expressions::expr<implementation::expressions::sub_1d<Cols>::template operation_row, implementation::expressions::ConstMatrixWrapper<base_matrix<T, Rows, Cols, StorageTraits, MathObject> > >
                    (typename implementation::expressions::sub_1d<Cols>::template operation_row<implementation::expressions::ConstMatrixWrapper<base_matrix<T, Rows, Cols, StorageTraits, MathObject> > >(implementation::col_count_storage<Cols>::cols(),
                                                                                                                                                                          row_index, *this),
                     implementation::expressions::make_matrix_wrapper(*this));
            }
            
            /**
               \brief Returns the `col_index`-th column of this matrix.
               
               \param col_index The index of the column.
               \return An object representing a column of `*this`.
             */
            inline typename implementation::expressions::expr<implementation::expressions::sub_1d<Rows>::template operation_col,
                                                              implementation::expressions::MatrixWrapper<base_matrix<T, Rows, Cols, StorageTraits, MathObject> > >
            col(size_type col_index)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::col(" << getAddress(*this) << "; " << col_index << ")");
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                return typename implementation::expressions::expr<implementation::expressions::sub_1d<Rows>::template operation_col, implementation::expressions::MatrixWrapper<base_matrix<T, Rows, Cols, StorageTraits, MathObject> > >
                    (typename implementation::expressions::sub_1d<Rows>::template operation_col<implementation::expressions::MatrixWrapper<base_matrix<T, Rows, Cols, StorageTraits, MathObject> > >(implementation::row_count_storage<Rows>::rows(),
                                                                                                                                                                     col_index, *this),
                     implementation::expressions::make_matrix_wrapper(*this));
            }
            
            /**
               \brief Returns the `col_index`-th column of this matrix.
               
               \param col_index The index of the column.
               \return An object representing a column of `*this`.
             */
            inline typename implementation::expressions::expr<implementation::expressions::sub_1d<Rows>::template operation_col,
                                                              implementation::expressions::ConstMatrixWrapper<base_matrix<T, Rows, Cols, StorageTraits, MathObject> > >
            col(size_type col_index) const
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::col(" << getAddress(*this) << "; " << col_index << ") const");
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                return typename implementation::expressions::expr<implementation::expressions::sub_1d<Rows>::template operation_col, implementation::expressions::ConstMatrixWrapper<base_matrix<T, Rows, Cols, StorageTraits, MathObject> > >
                    (typename implementation::expressions::sub_1d<Rows>::template operation_col<implementation::expressions::ConstMatrixWrapper<base_matrix<T, Rows, Cols, StorageTraits, MathObject> > >(implementation::row_count_storage<Rows>::rows(),
                                                                                                                                                                          col_index, *this),
                     implementation::expressions::make_matrix_wrapper(*this));
            }
            
            /**
               \brief Returns the transpose of this matrix.
               
               \return An object representing the transpose of `*this`.
             */
            inline implementation::expressions::expr<implementation::expressions::transpose,
                                                     implementation::expressions::MatrixWrapper<base_matrix<T, Rows, Cols, StorageTraits, MathObject> > >
            transpose() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::transpose(" << getAddress(*this) << ")");
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                return implementation::expressions::expr<implementation::expressions::transpose, implementation::expressions::MatrixWrapper<base_matrix<T, Rows, Cols, StorageTraits, MathObject> > >
                    (implementation::expressions::make_matrix_wrapper(*this));
            }
            
            /**
               \brief Returns the transpose of this matrix.
               
               \return An object representing the transpose of `*this`.
             */
            inline implementation::expressions::expr<implementation::expressions::transpose,
                                                     implementation::expressions::ConstMatrixWrapper<base_matrix<T, Rows, Cols, StorageTraits, MathObject> > >
            transpose() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::transpose(" << getAddress(*this) << ") const");
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                return implementation::expressions::expr<implementation::expressions::transpose, implementation::expressions::ConstMatrixWrapper<base_matrix<T, Rows, Cols, StorageTraits, MathObject> > >
                    (implementation::expressions::make_matrix_wrapper(*this));
            }
            
            /**
               \brief Resizes the vector to `new_size` entries. If `new_size` exceeds the current
                      size, new elements are default-constructed.
               
               Note that the number of rows or columns must be forced to be 1 during compile time;
               otherwise, this method will not compile.
               
               \param new_size The new number of entries of the vector.
             */
            void resize(size_type new_size)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::resize(" << getAddress(*this) << "; " << new_size << ") [1]");
                // Check
                PLLL_INTERNAL_STATIC_CHECK((Cols == 1) || (Rows == 1), OnlyDefinedForVectors);
                assert(((Cols == 1) && ((Rows < 0) || (new_size == static_cast<size_type>(Rows)))) ||
                       ((Rows == 1) && ((Cols < 0) || (new_size == static_cast<size_type>(Cols)))));
                if ((Cols == 1) ? (new_size == implementation::row_count_storage<Rows>::rows()) : (new_size == implementation::col_count_storage<Cols>::cols()))
                    return;
                // Allocate new space (don't construct yet!)
                StorageTraits::check(d_data, (Cols == 1) ? implementation::row_count_storage<Rows>::rows() : implementation::col_count_storage<Cols>::cols());
#if __cplusplus >= 201103L
                typename StorageTraits::pointer_type new_data = StorageTraits::clone_move(d_data, new_size,
                                                                                          ((Cols == 1) ? implementation::row_count_storage<Rows>::rows() : implementation::col_count_storage<Cols>::cols()) > new_size
                                                                                          ? new_size : ((Cols == 1) ? implementation::row_count_storage<Rows>::rows() : implementation::col_count_storage<Cols>::cols()));
#else
                typename StorageTraits::pointer_type new_data = StorageTraits::clone(d_data, new_size,
                                                                                     ((Cols == 1) ? implementation::row_count_storage<Rows>::rows() : implementation::col_count_storage<Cols>::cols()) > new_size
                                                                                     ? new_size : ((Cols == 1) ? implementation::row_count_storage<Rows>::rows() : implementation::col_count_storage<Cols>::cols()));
#endif
                // Copy over new data
                StorageTraits::free(d_data, (Cols == 1) ? implementation::row_count_storage<Rows>::rows() : implementation::col_count_storage<Cols>::cols());
                d_data = new_data;
                if (Cols == 1)
                    implementation::row_count_storage<Rows>::set_rows(new_size);
                else
                    implementation::col_count_storage<Cols>::set_cols(new_size);
                StorageTraits::check(d_data, (Cols == 1) ? implementation::row_count_storage<Rows>::rows() : implementation::col_count_storage<Cols>::cols());
            }
            
            /**
               \brief Resizes the vector to `new_size` entries. If `new_size` exceeds the current
                      size, new elements are copy-constructed from `init.ref()`.
               
               Note that the number of rows or columns must be forced to be 1 during compile time;
               otherwise, this method will not compile.
               
               \param new_size The new number of entries of the vector.
               \param init An initializer object carrying a reference to an instance of `S`.
             */
            template<typename S, bool def>
            void resize(size_type new_size, const implementation::Initialize_Impl<S, def> & init)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::resize(" << getAddress(*this) << "; " << new_size << " " << getAddress(init.ref()) << ") [1]");
                // Check
                PLLL_INTERNAL_STATIC_CHECK((Cols == 1) || (Rows == 1), OnlyDefinedForVectors);
                assert(((Cols == 1) && ((Rows < 0) || (new_size == static_cast<size_type>(Rows)))) ||
                       ((Rows == 1) && ((Cols < 0) || (new_size == static_cast<size_type>(Cols)))));
                if ((Cols == 1) ? (new_size == implementation::row_count_storage<Rows>::rows()) : (new_size == implementation::col_count_storage<Cols>::cols()))
                    return;
                // Allocate new space (don't construct yet!)
                StorageTraits::check(d_data, (Cols == 1) ? implementation::row_count_storage<Rows>::rows() : implementation::col_count_storage<Cols>::cols());
#if __cplusplus >= 201103L
                typename StorageTraits::pointer_type new_data = StorageTraits::clone_move(d_data, new_size,
                                                                                          ((Cols == 1) ? implementation::row_count_storage<Rows>::rows() : implementation::col_count_storage<Cols>::cols()) > new_size
                                                                                          ? new_size : ((Cols == 1) ? implementation::row_count_storage<Rows>::rows() : implementation::col_count_storage<Cols>::cols()), init.ref());
#else
                typename StorageTraits::pointer_type new_data = StorageTraits::clone(d_data, new_size,
                                                                                     ((Cols == 1) ? implementation::row_count_storage<Rows>::rows() : implementation::col_count_storage<Cols>::cols()) > new_size
                                                                                     ? new_size : ((Cols == 1) ? implementation::row_count_storage<Rows>::rows() : implementation::col_count_storage<Cols>::cols()), init.ref());
#endif
                // Copy over new data
                StorageTraits::free(d_data, (Cols == 1) ? implementation::row_count_storage<Rows>::rows() : implementation::col_count_storage<Cols>::cols());
                d_data = new_data;
                if (Cols == 1)
                    implementation::row_count_storage<Rows>::set_rows(new_size);
                else
                    implementation::col_count_storage<Cols>::set_cols(new_size);
                StorageTraits::check(d_data, (Cols == 1) ? implementation::row_count_storage<Rows>::rows() : implementation::col_count_storage<Cols>::cols());
            }
            
            /**
               \brief Resizes the matrix to `new_rows` times `new_cols` entries. If the new matrix
                      has entries which cannot be taken from the old matrix, they are
                      default-constructed.
               
               \param new_rows The new number of rows for the matrix.
               \param new_cols The new number of columns for the matrix.
             */
            void resize(size_type new_rows, size_type new_cols)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::resize(" << getAddress(*this) << "; " << new_rows << " " << new_cols << ") [2]");
                // Check
                if ((new_rows == implementation::row_count_storage<Rows>::rows()) && (new_cols == implementation::col_count_storage<Cols>::cols())) return;
                assert((Rows < 0) || (new_rows == static_cast<size_type>(Rows)));
                assert((Cols < 0) || (new_cols == static_cast<size_type>(Cols)));
                // Allocate new space (don't construct yet!)
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                typename StorageTraits::pointer_type new_data = StorageTraits::alloc_dontconstruct(new_rows * new_cols);
                // Copy existing and construct the rest
                size_type min_cols = std::min(new_cols, implementation::col_count_storage<Cols>::cols());
                size_type source_index = 0, dest_index = 0;
                for (size_type i = (Rows < 0 ? std::min(new_rows, implementation::row_count_storage<Rows>::rows()) : Rows); i > 0; --i)
                {
                    for (size_type j = (Cols < 0 ? min_cols : Cols); j > 0; --j, ++source_index, ++dest_index)
#if __cplusplus >= 201103L
                        StorageTraits::move_construct(new_data[dest_index], std::move(d_data[source_index]));
#else
                        StorageTraits::copy_construct(new_data[dest_index], d_data[source_index]);
#endif
                    if (implementation::col_count_storage<Cols>::cols() < new_cols)
                    {
                        for (size_type j = new_cols - implementation::col_count_storage<Cols>::cols(); j > 0; --j, ++dest_index)
                            StorageTraits::construct(new_data[dest_index]);
                    }
                    else
                        source_index += implementation::col_count_storage<Cols>::cols() - new_cols;
                }
                // Construct the rest
                if ((Rows < 0) && (implementation::row_count_storage<Rows>::rows() < new_rows))
                {
                    size_type s = (new_rows - implementation::row_count_storage<Rows>::rows()) * (Cols < 0 ? new_cols : Cols);
                    for (; s > 0; --s, ++dest_index)
                        StorageTraits::construct(new_data[dest_index]);
                }
                // Copy over new data
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                StorageTraits::free(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                d_data = new_data;
                implementation::row_count_storage<Rows>::set_rows(new_rows);
                implementation::col_count_storage<Cols>::set_cols(new_cols);
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
            }
            
            /**
               \brief Resizes the matrix to `new_rows` times `new_cols` entries. If the new matrix
                      has entries which cannot be taken from the old matrix, they are
                      copy-constructed from `init.ref()`.
               
               \param new_rows The new number of rows for the matrix.
               \param new_cols The new number of columns for the matrix.
               \param init An initializer object carrying a reference to an instance of `S`.
             */
            template<typename S, bool def>
            void resize(size_type new_rows, size_type new_cols, const implementation::Initialize_Impl<S, def> & init)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::resize(" << getAddress(*this) << "; " << new_rows << " " << new_cols << " " << getAddress(init.ref()) << ") [3]");
                // Check
                if ((new_rows == implementation::row_count_storage<Rows>::rows()) && (new_cols == implementation::col_count_storage<Cols>::cols())) return;
                assert((Rows < 0) || (new_rows == static_cast<size_type>(Rows)));
                assert((Cols < 0) || (new_cols == static_cast<size_type>(Cols)));
                // Allocate new space (don't construct yet!)
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                typename StorageTraits::pointer_type new_data = StorageTraits::alloc_dontconstruct(new_rows * new_cols);
                // Copy existing and construct the rest
                size_type min_cols = std::min(new_cols, implementation::col_count_storage<Cols>::cols());
                size_type source_index = 0, dest_index = 0;
                for (size_type i = (Rows < 0 ? std::min(new_rows, implementation::row_count_storage<Rows>::rows()) : Rows); i > 0; --i)
                {
                    for (size_type j = (Cols < 0 ? min_cols : Cols); j > 0; --j, ++source_index, ++dest_index)
#if __cplusplus >= 201103L
                        StorageTraits::move_construct(new_data[dest_index], std::move(d_data[source_index]));
#else
                        StorageTraits::copy_construct(new_data[dest_index], d_data[source_index]);
#endif
                    if (implementation::col_count_storage<Cols>::cols() < new_cols)
                    {
                        for (size_type j = new_cols - implementation::col_count_storage<Cols>::cols(); j > 0; --j, ++dest_index)
                            StorageTraits::copy_construct(new_data[dest_index], init.ref());
                    }
                    else
                        source_index += implementation::col_count_storage<Cols>::cols() - new_cols;
                }
                // Construct the rest
                if ((Rows < 0) && (implementation::row_count_storage<Rows>::rows() < new_rows))
                {
                    size_type s = (new_rows - implementation::row_count_storage<Rows>::rows()) * (Cols < 0 ? new_cols : Cols);
                    for (; s > 0; --s, ++dest_index)
                        StorageTraits::copy_construct(new_data[dest_index], init.ref());
                }
                // Copy over new data
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                StorageTraits::free(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                d_data = new_data;
                implementation::row_count_storage<Rows>::set_rows(new_rows);
                implementation::col_count_storage<Cols>::set_cols(new_cols);
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
            }

            /**
               \brief Resizes the matrix to fit the format of `mat` and assigns the contents of
                      `mat` to the matrix.
               
               Note that if `mat` is an expression involving `*this`, there will be no problems as
               the old elements will only be released after the new elements are created.
               
               \param mat The matrix to copy from.
             */
            template<typename MatrixType>
            void assign_resize(const MatrixType & mat)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::assign_resize(" << getAddress(*this) << "; " << getAddress(mat) << ")");
                // Check
                assert((Rows < 0) || (mat.rows() == static_cast<size_type>(Rows)));
                assert((Cols < 0) || (mat.cols() == static_cast<size_type>(Cols)));
                // Allocate new space (don't construct yet!)
                size_type new_rows = mat.rows(), new_cols = mat.cols();
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                typename StorageTraits::pointer_type new_data = StorageTraits::alloc_dontconstruct(new_rows * new_cols);
                // Construct entries from old matrix
                if (mat.get_coeff_alwayszero())
                {
                    size_type s = implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols();
                    for (size_type index = 0; index < s; ++index)
                    {
                        StorageTraits::construct(new_data[index]);
                        set_zero(new_data[index]);
                    }
                }
                else
                {
                    typename MatrixType::ConstEnumerator e = mat.enumerate();
                    for (size_type index = 0; e.has_current(); e.next(), ++index)
                        if (implementation::MatrixInfo<MatrixType>::coeffs_are_simple_expressions)
                            StorageTraits::copy_construct(new_data[index], e.current());
                        else
                        {
                            typename MatrixType::ConstEnumerator::GetCoeffSteps_Type coeff = e.get_current_steps();
                            StorageTraits::copy_construct(new_data[index], coeff.step1());
                            coeff.step2(new_data[index]);
                        }
                }
                // Free old matrix
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                StorageTraits::free(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                // Copy over new data
                d_data = new_data;
                implementation::row_count_storage<Rows>::set_rows(new_rows);
                implementation::col_count_storage<Cols>::set_cols(new_cols);
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
            }
            
#if __cplusplus >= 201103L
            template<bool MO>
            inline void move_resize(base_matrix<T, Rows, Cols, StorageTraits, MO> && mat) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                *this = std::move(mat);
            }
            
            /**
               \brief Resizes the matrix to fit the format of `mat` and moves the contents of
                      `mat` to the matrix.
               
               Note that if `mat` is an expression involving `*this`, there will be no problems as
               the old elements will only be released after the new elements are created.
               
               \param mat The matrix to move from.
             */
            template<typename MatrixType>
            void move_resize(MatrixType && mat)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::assign_resize(" << getAddress(*this) << "; " << getAddress(mat) << ")");
                PLLL_INTERNAL_STATIC_CHECK(implementation::MatrixInfo<MatrixType>::can_move_from, CannotMoveFromMatrix);
                // Check
                assert((Rows < 0) || (mat.rows() == static_cast<size_type>(Rows)));
                assert((Cols < 0) || (mat.cols() == static_cast<size_type>(Cols)));
                // Allocate new space (don't construct yet!)
                size_type new_rows = mat.rows(), new_cols = mat.cols();
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                typename StorageTraits::pointer_type new_data = StorageTraits::alloc_dontconstruct(new_rows * new_cols);
                // Construct entries from old matrix
                if (mat.get_coeff_alwayszero())
                {
                    size_type s = implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols();
                    for (size_type index = 0; index < s; ++index)
                    {
                        StorageTraits::construct(new_data[index]);
                        set_zero(new_data[index]);
                    }
                }
                else
                {
                    typename MatrixType::ConstEnumerator e = mat.enumerate();
                    for (size_type index = 0; e.has_current(); e.next(), ++index)
                        if (implementation::MatrixInfo<MatrixType>::coeffs_are_simple_expressions)
                            StorageTraits::move_construct(new_data[index], std::move(e.current()));
                        else
                        {
                            typename MatrixType::ConstEnumerator::GetCoeffSteps_Type coeff = e.get_current_steps();
                            StorageTraits::copy_construct(new_data[index], coeff.step1());
                            coeff.step2(new_data[index]);
                        }
                }
                // Free old matrix
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                StorageTraits::free(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                // Copy over new data
                d_data = new_data;
                implementation::row_count_storage<Rows>::set_rows(new_rows);
                implementation::col_count_storage<Cols>::set_cols(new_cols);
            }
#else
            /**
               \brief Resizes the matrix to fit the format of `mat` and moves the contents of
                      `mat` to the matrix.
               
               Note that if `mat` is an expression involving `*this`, there will be no problems as
               the old elements will only be released after the new elements are created.
               
               \param mat The matrix to move from.
             */
            template<typename MatrixType>
            void move_resize(const MatrixType & mat)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::assign_resize(" << getAddress(*this) << "; " << getAddress(mat) << ")");
                PLLL_INTERNAL_STATIC_CHECK(implementation::MatrixInfo<MatrixType>::can_move_from, CannotMoveFromMatrix);
                // Check
                assert((Rows < 0) || (mat.rows() == static_cast<size_type>(Rows)));
                assert((Cols < 0) || (mat.cols() == static_cast<size_type>(Cols)));
                // Allocate new space (don't construct yet!)
                size_type new_rows = mat.rows(), new_cols = mat.cols();
                typename StorageTraits::pointer_type new_data = StorageTraits::alloc_dontconstruct(new_rows * new_cols);
                // Construct entries from old matrix
                if (mat.get_coeff_alwayszero())
                {
                    size_type s = implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols();
                    for (size_type index = 0; index < s; ++index)
                    {
                        StorageTraits::construct(new_data[index]);
                        set_zero(new_data[index]);
                    }
                }
                else
                {
                    typename MatrixType::ConstEnumerator e = mat.enumerate();
                    for (size_type index = 0; e.has_current(); e.next(), ++index)
                        if (implementation::MatrixInfo<MatrixType>::coeffs_are_simple_expressions)
                            StorageTraits::copy_construct(new_data[index], e.current());
                        else
                        {
                            typename MatrixType::ConstEnumerator::GetCoeffSteps_Type coeff = e.get_current_steps();
                            StorageTraits::copy_construct(new_data[index], coeff.step1());
                            coeff.step2(new_data[index]);
                        }
                }
                // Free old matrix
                StorageTraits::free(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                // Copy over new data
                d_data = new_data;
                implementation::row_count_storage<Rows>::set_rows(new_rows);
                implementation::col_count_storage<Cols>::set_cols(new_cols);
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
            }
#endif
            
            /**
               \brief Swaps the contents of this matrix with the matrix `B`.
               
               \param B Another matrix to swap the contents with.
             */
            template<typename MatrixType>
            inline void swap(MatrixType & B)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::swap(" << getAddress(*this) << "; " << getAddress(B) << ")");
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
                linalg::swap(*this, B);
                StorageTraits::check(d_data, implementation::row_count_storage<Rows>::rows() * implementation::col_count_storage<Cols>::cols());
            }
            
            /**
               \brief Tests whether this matrix (or expression) involves `A`.
               
               This will effectively compare the addresses of `A` and `*this`.
               
               \param A The other matrix to compare this one to.
             */
            inline bool involves_this_matrix(const base_matrix<T, Rows, Cols, StorageTraits, MathObject> & A) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::involves_this_matrix(" << getAddress(*this) << "; " << getAddress(A) << ") [1]");
                return &A == this;
            }
            
            /**
               \brief Tests whether this matrix (or expression) involves `A`.
               
               This will always return `false`, since possible equality is already handled by the
               other overload.
               
               \param A The other matrix to compare this one to.
             */
            template<typename MatrixType>
            inline bool involves_this_matrix(const MatrixType & A) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::involves_this_matrix(" << getAddress(*this) << "; " << getAddress(A) << ") [2]");
                return false;
            }
            
            /**
               \brief For the underlying matrices `U` of the current expression, calls
                      `A.involves_this_matrix(U)`. Returns `true` if any of these returns `true`.
               
               Since this is already a matrix, will call `A.involves_this_matrix(*this)` directly.
               
               \param A A matrix to call `involves_this_matrix()` of.
             */
            template<typename MatrixType>
            inline bool test_involvement(const MatrixType & A) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_matrix::test_involvement(" << getAddress(*this) << "; " << getAddress(A) << ")");
                return A.involves_this_matrix(*this);
            }
            
            class ConstEnumerator;

            /**
               \brief An enumerator for matrix elements.
               
               It enumerates first through the rows and then through the columns. For example, a
              \f$2 \times 3\f$ matrix would be enumerated with indices `(0, 0)`, `(0, 1)`, `(0, 2)`,
              `(1, 0)`, `(1, 1)`, `(1, 2)`.
             */
            class Enumerator
            {
                friend class ConstEnumerator;
                
            private:
                T * d_entry;
                size_type d_left;
                
            public:
                /**
                   \brief The final type of the coefficients.
                 */
                typedef T & Type;
                /**
                   \brief The intermediate type of the coefficients, which is returned by `current()`.
                 */
                typedef base_matrix<T, Rows, Cols, StorageTraits, MathObject>::GetCoeffSteps_Type GetCoeffSteps_Type;
                
                inline Enumerator(T * entry, size_type left) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    : d_entry(entry), d_left(left)
                {
                }

                /**
                   \brief Returns the current matrix entry.
                   
                   \return A reference to the current entry.
                 */
                inline Type current() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    return *d_entry;
                }
                
                /**
                   \brief Creates a copy of the current matrix entry in `result`.
                   
                   \param result An object to store the current element in.
                 */
                template<typename Result>
                inline void get_current(Result & result) const PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(result = *d_entry))
                {
                    result = *d_entry;
                }
                
                /**
                   \brief Returns a `GetCoeffSteps_Type` object for the current entry.
                   
                   \return A `GetCoeffSteps_Type` object for the current element.
                 */
                inline GetCoeffSteps_Type get_current_steps() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    return GetCoeffSteps_Type(*d_entry);
                }
                
                /**
                   \brief Returns `true` if a current element exists and can be queried by the
                          `current()`, `get_current()` and `get_current_steps()` methods.
                   
                   \return Returns `true` if an element is available.
                 */
                inline bool has_current() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    return d_left > 0;
                }

                /**
                   \brief Goes to the next element.
                   
                   Must only be called if `has_current()` is true.
                 */
                inline void next() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    --d_left;
                    ++d_entry;
                }
            };
            
            /**
               \brief An enumerator for matrix elements.
               
               It enumerates first through the rows and then through the columns. For example, a
              \f$2 \times 3\f$ matrix would be enumerated with indices `(0, 0)`, `(0, 1)`, `(0, 2)`,
              `(1, 0)`, `(1, 1)`, `(1, 2)`.
             */
            class ConstEnumerator
            {
            private:
                const T * d_entry;
                size_type d_left;
                
            public:
                /**
                   \brief The final type of the coefficients.
                 */
                typedef const T & Type;
                /**
                   \brief The intermediate type of the coefficients, which is returned by `current()`.
                 */
                typedef base_matrix<T, Rows, Cols, StorageTraits, MathObject>::GetCoeffSteps_Type GetCoeffSteps_Type;
                
                inline ConstEnumerator(const Enumerator & e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    : d_entry(e.d_entry), d_left(e.d_left)
                {
                }
                
#if __cplusplus >= 201103L
                inline ConstEnumerator(Enumerator && e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    : d_entry(e.d_entry), d_left(e.d_left)
                {
                }
#endif
                
                inline ConstEnumerator(const T * entry, size_type left) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    : d_entry(entry), d_left(left)
                {
                }
                
                /**
                   \brief Returns the current matrix entry.
                   
                   \return A reference to the current entry.
                 */
                inline Type current() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    return *d_entry;
                }
                
                /**
                   \brief Creates a copy of the current matrix entry in `result`.
                   
                   \param result An object to store the current element in.
                 */
                template<typename Result>
                inline void get_current(Result & result) const PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(result = *d_entry))
                {
                    result = *d_entry;
                }
                
                /**
                   \brief Returns a `GetCoeffSteps_Type` object for the current entry.
                   
                   \return A `GetCoeffSteps_Type` object for the current element.
                 */
                inline GetCoeffSteps_Type get_current_steps() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    return GetCoeffSteps_Type(*d_entry);
                }
                
                /**
                   \brief Returns `true` if a current element exists and can be queried by the
                          `current()`, `get_current()` and `get_current_steps()` methods.
                   
                   \return Returns `true` if an element is available.
                 */
                inline bool has_current() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    return d_left > 0;
                }
                
                /**
                   \brief Goes to the next element.
                   
                   Must only be called if `has_current()` is true.
                 */
                inline void next() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    --d_left;
                    ++d_entry;
                }
            };
            
            typedef Enumerator DefaultEnumerator;
            
            typedef Enumerator RowEnumerator;
            typedef ConstEnumerator ConstRowEnumerator;
            typedef DefaultEnumerator DefaultRowEnumerator;
            
            class ConstColEnumerator;
            
            /**
               \brief An enumerator for a matrix' column.
               
               It sequentially enumerates the elements of a matrix' column
             */
            class ColEnumerator : private implementation::col_count_storage<Cols>
            {
                friend class ConstColEnumerator;
                
            private:
                size_type d_left;
                T * d_entry;
                
            public:
                /**
                   \brief The final type of the coefficients.
                 */
                typedef T & Type;
                /**
                   \brief The intermediate type of the coefficients, which is returned by `current()`.
                 */
                typedef base_matrix<T, Rows, Cols, StorageTraits, MathObject>::GetCoeffSteps_Type GetCoeffSteps_Type;
                
                inline ColEnumerator(T * entry, size_type left, size_type cols) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    : implementation::col_count_storage<Cols>(cols), d_left(left), d_entry(entry)
                {
                }
                
                /**
                   \brief Returns the current matrix entry.
                   
                   \return A reference to the current entry.
                 */
                inline Type current() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    return *d_entry;
                }
                
                /**
                   \brief Creates a copy of the current matrix entry in `result`.
                   
                   \param result An object to store the current element in.
                 */
                template<typename Result>
                inline void get_current(Result & result) const PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(result = *d_entry))
                {
                    result = *d_entry;
                }
                
                /**
                   \brief Returns a `GetCoeffSteps_Type` object for the current entry.
                   
                   \return A `GetCoeffSteps_Type` object for the current element.
                 */
                inline GetCoeffSteps_Type get_current_steps() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    return GetCoeffSteps_Type(*d_entry);
                }
                
                /**
                   \brief Returns `true` if a current element exists and can be queried by the
                          `current()`, `get_current()` and `get_current_steps()` methods.
                   
                   \return Returns `true` if an element is available.
                 */
                inline bool has_current() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    return d_left > 0;
                }
                
                /**
                   \brief Goes to the next element.
                   
                   Must only be called if `has_current()` is true.
                 */
                inline void next() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    --d_left;
                    d_entry += implementation::col_count_storage<Cols>::cols();
                }
            };
            
            /**
               \brief An enumerator for a matrix' column.
               
               It sequentially enumerates the elements of a matrix' column
             */
            class ConstColEnumerator : private implementation::col_count_storage<Cols>
            {
            private:
                size_type d_left;
                const T * d_entry;
                
            public:
                /**
                   \brief The final type of the coefficients.
                 */
                typedef const T & Type;
                /**
                   \brief The intermediate type of the coefficients, which is returned by `current()`.
                 */
                typedef base_matrix<T, Rows, Cols, StorageTraits, MathObject>::GetCoeffSteps_Type GetCoeffSteps_Type;
                
                inline ConstColEnumerator(const ColEnumerator & e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    : implementation::col_count_storage<Cols>(e), d_left(e.d_left), d_entry(e.d_entry)
                {
                }
                
#if __cplusplus >= 201103L
                inline ConstColEnumerator(ColEnumerator && e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    : implementation::col_count_storage<Cols>(e), d_left(e.d_left), d_entry(e.d_entry)
                {
                }
#endif
                
                inline ConstColEnumerator(const T * entry, size_type left, size_type cols) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    : implementation::col_count_storage<Cols>(cols), d_left(left), d_entry(entry)
                {
                }
                
                /**
                   \brief Returns the current matrix entry.
                   
                   \return A reference to the current entry.
                 */
                inline Type current() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    return *d_entry;
                }
                
                /**
                   \brief Creates a copy of the current matrix entry in `result`.
                   
                   \param result An object to store the current element in.
                 */
                template<typename Result>
                inline void get_current(Result & result) const PLLL_INTERNAL_NOTHROW_POSTFIX_CONDITIONAL(noexcept(result = *d_entry))
                {
                    result = *d_entry;
                }
                
                /**
                   \brief Returns a `GetCoeffSteps_Type` object for the current entry.
                   
                   \return A `GetCoeffSteps_Type` object for the current element.
                 */
                inline GetCoeffSteps_Type get_current_steps() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    return GetCoeffSteps_Type(*d_entry);
                }
                
                /**
                   \brief Returns `true` if a current element exists and can be queried by the
                          `current()`, `get_current()` and `get_current_steps()` methods.
                   
                   \return Returns `true` if an element is available.
                 */
                inline bool has_current() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    return d_left > 0;
                }
                
                /**
                   \brief Goes to the next element.
                   
                   Must only be called if `has_current()` is true.
                 */
                inline void next() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    --d_left;
                    d_entry += implementation::col_count_storage<Cols>::cols();
                }
            };
            
            typedef ColEnumerator DefaultColEnumerator;

            /**
               \brief Enumerates all entries of the current matrix.
               
               \return Returns a `Enumerator` object for the current matrix.
             */
            inline Enumerator enumerate() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return Enumerator(d_data, size());
            }
            
            /**
               \brief Enumerates all entries of the current matrix.
               
               \return Returns a `ConstEnumerator` object for the current matrix.
             */
            inline ConstEnumerator enumerate() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return ConstEnumerator(d_data, size());
            }
            
            /**
               \brief Enumerates all entries of the given row.
               
               \param row The row to enumerate.
               \return Returns a `RowEnumerator` object for the specified row.
             */
            inline RowEnumerator enumerate_row(size_type row) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                assert(row < implementation::row_count_storage<Rows>::rows());
                return RowEnumerator(d_data + row * cols(), cols());
            }
            
            /**
               \brief Enumerates all entries of the given row.
               
               \param row The row to enumerate.
               \return Returns a `ConstRowEnumerator` object for the specified row.
             */
            inline ConstRowEnumerator enumerate_row(size_type row) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                assert(row < implementation::row_count_storage<Rows>::rows());
                return ConstRowEnumerator(d_data + row * cols(), cols());
            }
            
            /**
               \brief Enumerates all entries of the given column.
               
               \param col The column to enumerate.
               \return Returns a `ColumnEnumerator` object for the specified column.
             */
            inline ColEnumerator enumerate_col(size_type col) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                assert(col < implementation::col_count_storage<Cols>::cols());
                return ColEnumerator(d_data + col, rows(), cols());
            }
            
            /**
               \brief Enumerates all entries of the given column.
               
               \param col The column to enumerate.
               \return Returns a `ConstColumnEnumerator` object for the specified column.
             */
            inline ConstColEnumerator enumerate_col(size_type col) const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                assert(col < implementation::col_count_storage<Cols>::cols());
                return ConstColEnumerator(d_data + col, rows(), cols());
            }
            
            class ConstRowsEnumerator;
            
            /**
               \brief Enumerates over all rows. For every row, a `RowsEnumerator::RowEnumerator`
                      object is returned which allows to enumerate over that row.
             */
            class RowsEnumerator : private implementation::col_count_storage<Cols>
            {
                friend class ConstRowsEnumerator;
                
            private:
                size_type d_left;
                T * d_entry;
                
            public:
                /**
                   \brief The row enumerator type.
                 */
                typedef RowEnumerator Type;
                
                inline RowsEnumerator(T * entry, size_type rows, size_type cols) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    : implementation::col_count_storage<Cols>(cols), d_left(rows), d_entry(entry)
                {
                }
                
                /**
                   \brief Returns a current row's enumerator.
                   
                   \return An enumerator for the current row.
                 */
                inline Type current() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    return Enumerator(d_entry, implementation::col_count_storage<Cols>::cols());
                }
                
                /**
                   \brief Returns `true` if a current row exists and can be queried by the
                          `current()` method.
                   
                   \return Returns `true` if a row is available.
                 */
                inline bool has_current() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    return d_left > 0;
                }
                
                /**
                   \brief Goes to the next element.
                   
                   Must only be called if `has_current()` is true.
                 */
                inline void next() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    --d_left;
                    d_entry += implementation::col_count_storage<Cols>::cols();
                }
            };
            
            /**
               \brief Enumerates over all rows. For every row, a
                      `ConstRowsEnumerator::RowEnumerator` object is returned which allows to
                      enumerate over that row.
             */
            class ConstRowsEnumerator : private implementation::col_count_storage<Cols>
            {
            private:
                size_type d_left;
                const T * d_entry;
                
            public:
                /**
                   \brief The row enumerator type.
                 */
                typedef ConstRowEnumerator Type;
                
                inline ConstRowsEnumerator(const RowsEnumerator & e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    : implementation::col_count_storage<Cols>(e), d_left(e.d_left), d_entry(e.d_entry)
                {
                }
                
#if __cplusplus >= 201103L
                inline ConstRowsEnumerator(RowsEnumerator && e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    : implementation::col_count_storage<Cols>(e), d_left(e.d_left), d_entry(e.d_entry)
                {
                }
#endif
                
                inline ConstRowsEnumerator(const T * entry, size_type rows, size_type cols) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    : implementation::col_count_storage<Cols>(cols), d_left(rows), d_entry(entry)
                {
                }
                
                /**
                   \brief Returns a current row's enumerator.
                   
                   \return An enumerator for the current row.
                 */
                inline Type current() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    return ConstEnumerator(d_entry, implementation::col_count_storage<Cols>::cols());
                }
                
                /**
                   \brief Returns `true` if a current row exists and can be queried by the
                          `current()` method.
                   
                   \return Returns `true` if a row is available.
                 */
                inline bool has_current() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    return d_left > 0;
                }
                
                /**
                   \brief Goes to the next element.
                   
                   Must only be called if `has_current()` is true.
                 */
                inline void next() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    --d_left;
                    d_entry += implementation::col_count_storage<Cols>::cols();
                }
            };
            
            class ConstColsEnumerator;
            
            /**
               \brief Enumerates over all columns. For every column, a
                     `ColsEnumerator::ColEnumerator` object is returned which allows to enumerate
                     over that column.
             */
            class ColsEnumerator : private implementation::row_count_storage<Rows>, implementation::col_count_storage<Cols>
            {
                friend class ConstColsEnumerator;
                
            private:
                size_type d_left;
                T * d_entry;
                
            public:
                /**
                   \brief The row enumerator type.
                 */
                typedef ColEnumerator Type;
                
                inline ColsEnumerator(T * entry, size_type rows, size_type cols) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    : implementation::row_count_storage<Rows>(rows), implementation::col_count_storage<Cols>(cols), d_left(cols), d_entry(entry)
                {
                }
                
                /**
                   \brief Returns a current column's enumerator.
                   
                   \return An enumerator for the current column.
                 */
                inline Type current() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    return ColEnumerator(d_entry, implementation::row_count_storage<Rows>::rows(), implementation::col_count_storage<Cols>::cols());
                }
                
                /**
                   \brief Returns `true` if a current column exists and can be queried by the
                          `current()` method.
                   
                   \return Returns `true` if a column is available.
                 */
                inline bool has_current() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    return d_left > 0;
                }
                
                /**
                   \brief Goes to the next element.
                   
                   Must only be called if `has_current()` is true.
                 */
                inline void next() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    --d_left;
                    ++d_entry;
                }
            };
            
            /**
               \brief Enumerates over all columns. For every column, a
                     `ConstColsEnumerator::ColEnumerator` object is returned which allows to
                     enumerate over that column.
             */
            class ConstColsEnumerator : private implementation::row_count_storage<Rows>, implementation::col_count_storage<Cols>
            {
            private:
                size_type d_left;
                const T * d_entry;
                
            public:
                /**
                   \brief The row enumerator type.
                 */
                typedef ConstColEnumerator Type;
                
                inline ConstColsEnumerator(const ColsEnumerator & e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    : implementation::row_count_storage<Rows>(e), d_left(e.d_left), d_entry(e.d_entry)
                {
                }
                
#if __cplusplus >= 201103L
                inline ConstColsEnumerator(ColsEnumerator && e) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    : implementation::row_count_storage<Rows>(e), d_left(e.d_left), d_entry(e.d_entry)
                {
                }
#endif
                
                inline ConstColsEnumerator(const T * entry, size_type rows, size_type cols) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                    : implementation::row_count_storage<Rows>(rows), implementation::col_count_storage<Cols>(cols), d_left(cols), d_entry(entry)
                {
                }
                
                /**
                   \brief Returns a current column's enumerator.
                   
                   \return An enumerator for the current column.
                 */
                inline Type current() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    return ConstColEnumerator(d_entry, implementation::row_count_storage<Rows>::rows(), implementation::col_count_storage<Cols>::cols());
                }
                
                /**
                   \brief Returns `true` if a current column exists and can be queried by the
                          `current()` method.
                   
                   \return Returns `true` if a column is available.
                 */
                inline bool has_current() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    return d_left > 0;
                }
                
                /**
                   \brief Goes to the next element.
                   
                   Must only be called if `has_current()` is true.
                 */
                inline void next() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                {
                    --d_left;
                    ++d_entry;
                }
            };
            
            typedef RowsEnumerator DefaultRowsEnumerator;
            typedef ColsEnumerator DefaultColsEnumerator;
            
            /**
               \brief Enumerates all rows. For each row, a row enumerator is given.
               
               \return Returns a `RowsEnumerator` object for this matrix.
             */
            inline RowsEnumerator enumerate_rows() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return RowsEnumerator(d_data, rows(), cols());
            }
            
            /**
               \brief Enumerates all rows. For each row, a row enumerator is given.
               
               \return Returns a `ConstRowsEnumerator` object for this matrix.
             */
            inline ConstRowsEnumerator enumerate_rows() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return ConstRowsEnumerator(d_data, rows(), cols());
            }
            
            /**
               \brief Enumerates all columns. For each column, a column enumerator is given.
               
               \return Returns a `ColsEnumerator` object for this matrix.
             */
            inline ColsEnumerator enumerate_cols() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return ColsEnumerator(d_data, rows(), cols());
            }
            
            /**
               \brief Enumerates all columns. For each column, a column enumerator is given.
               
               \return Returns a `ConstColsEnumerator` object for this matrix.
             */
            inline ConstColsEnumerator enumerate_cols() const PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                return ConstColsEnumerator(d_data, rows(), cols());
            }
        };
        
#if __cplusplus >= 201103L
        template<typename T, int Rows = Flexible, int Cols = Flexible, typename StorageTraits = storage_traits<T> >
        /**
           \brief Represents a math matrix with coefficients in `T`.
           
           This is in fact the same as `base_matrix<T, Rows, Cols, StorageTraits, true>`.
           
           \tparam T              The coefficient type.
           \tparam Rows           If set to `plll::linalg::Flexible`, allows to flexibly choose
                                  and change the number of rows. Otherwise, enforces fixed number
                                  of rows during compile time. Default value is
                                  `plll::linalg::Flexible`.
           \tparam Cols           If set to `plll::linalg::Flexible`, allows to flexibly choose
                                  and change the number of columns. Otherwise, enforces fixed number
                                  of columns during compile time. Default value is
                                  `plll::linalg::Flexible`.
           \tparam StorageTraits  Allows to determine how the matrix's contents are stored. The
                                  default value is `plll::linalg::storage_traits<T>`, which should
                                  usually not be changed.
           \tparam MathObject     Allows to determine whether the matrix is a math object (`true`)
                                  and thus operations can be performed on it, or whether it is a
                                  non-math object (`false`) which allows no operations except
                                  reading and storing coefficients and resizing.
         */
        using math_matrix = base_matrix<T, Rows, Cols, StorageTraits, true>;
        
        template<typename T, int Cols = Flexible, typename StorageTraits = storage_traits<T>, bool MathObject = false>
        /**
           \brief Represents a row vector with coefficients in `T`.
           
           This is in fact the same as `base_matrix<T, 1, Cols, StorageTraits, MathObject>`.
           
           \tparam T              The coefficient type.
           \tparam Cols           If set to `plll::linalg::Flexible`, allows to flexibly choose
                                  and change the number of columns. Otherwise, enforces fixed number
                                  of columns during compile time. Default value is
                                  `plll::linalg::Flexible`.
           \tparam StorageTraits  Allows to determine how the vector's contents are stored. The
                                  default value is `plll::linalg::storage_traits<T>`, which should
                                  usually not be changed.
           \tparam MathObject     Allows to determine whether the vector is a math object (`true`)
                                  and thus operations can be performed on it, or whether it is a
                                  non-math object (`false`) which allows no operations except
                                  reading and storing coefficients and resizing.
         */
        using base_rowvector = base_matrix<T, 1, Cols, StorageTraits, MathObject>;

        template<typename T, int Rows = Flexible, typename StorageTraits = storage_traits<T>, bool MathObject = false>
        /**
           \brief Represents a column vector with coefficients in `T`.
           
           This is in fact the same as `base_matrix<T, Rows, 1, StorageTraits, MathObject>`.
           
           \tparam T              The coefficient type.
           \tparam Rows           If set to `plll::linalg::Flexible`, allows to flexibly choose
                                  and change the number of rows. Otherwise, enforces fixed number
                                  of rows during compile time. Default value is
                                  `plll::linalg::Flexible`.
           \tparam StorageTraits  Allows to determine how the vector's contents are stored. The
                                  default value is `plll::linalg::storage_traits<T>`, which should
                                  usually not be changed.
           \tparam MathObject     Allows to determine whether the vector is a math object (`true`)
                                  and thus operations can be performed on it, or whether it is a
                                  non-math object (`false`) which allows no operations except
                                  reading and storing coefficients and resizing.
         */
        using base_colvector = base_matrix<T, Rows, 1, StorageTraits, MathObject>;
        
        template<typename T, int Cols = Flexible, typename StorageTraits = storage_traits<T> >
        /**
           \brief Represents a math row vector with coefficients in `T`.
           
           This is in fact the same as `base_matrix<T, 1, Cols, StorageTraits, true>`.
           
           \tparam T              The coefficient type.
           \tparam Cols           If set to `plll::linalg::Flexible`, allows to flexibly choose
                                  and change the number of columns. Otherwise, enforces fixed number
                                  of columns during compile time. Default value is
                                  `plll::linalg::Flexible`.
           \tparam StorageTraits  Allows to determine how the vector's contents are stored. The
                                  default value is `plll::linalg::storage_traits<T>`, which should
                                  usually not be changed.
         */
        using math_rowvector = base_rowvector<T, Cols, StorageTraits, true>;
        
        template<typename T, int Rows = Flexible, typename StorageTraits = storage_traits<T> >
        /**
           \brief Represents a math column vector with coefficients in `T`.
           
           This is in fact the same as `base_matrix<T, Rows, 1, StorageTraits, true>`.
           
           \tparam T              The coefficient type.
           \tparam Rows           If set to `plll::linalg::Flexible`, allows to flexibly choose
                                  and change the number of rows. Otherwise, enforces fixed number
                                  of rows during compile time. Default value is
                                  `plll::linalg::Flexible`.
           \tparam StorageTraits  Allows to determine how the vector's contents are stored. The
                                  default value is `plll::linalg::storage_traits<T>`, which should
                                  usually not be changed.
         */
        using math_colvector = base_colvector<T, Rows, StorageTraits, true>;
#else
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        //// math_matrix
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        
        template<typename T, int Rows = Flexible, int Cols = Flexible, typename StorageTraits = storage_traits<T> >
        /**
           \brief Represents a math matrix with coefficients in `T`.
           
           This is in fact the same as `base_matrix<T, Rows, Cols, StorageTraits, true>`.
           
           \tparam T              The coefficient type.
           \tparam Rows           If set to `plll::linalg::Flexible`, allows to flexibly choose
                                  and change the number of rows. Otherwise, enforces fixed number
                                  of rows during compile time. Default value is
                                  `plll::linalg::Flexible`.
           \tparam Cols           If set to `plll::linalg::Flexible`, allows to flexibly choose
                                  and change the number of columns. Otherwise, enforces fixed number
                                  of columns during compile time. Default value is
                                  `plll::linalg::Flexible`.
           \tparam StorageTraits  Allows to determine how the matrix's contents are stored. The
                                  default value is `plll::linalg::storage_traits<T>`, which should
                                  usually not be changed.
           \tparam MathObject     Allows to determine whether the matrix is a math object (`true`)
                                  and thus operations can be performed on it, or whether it is a
                                  non-math object (`false`) which allows no operations except
                                  reading and storing coefficients and resizing.
         */
        class math_matrix : public base_matrix<T, Rows, Cols, StorageTraits, true>
        {
        public:
            enum { has_direct_access = base_matrix<T, Rows, Cols, StorageTraits, true>::has_direct_access };

            /**
               \brief Creates a new matrix. The dimensions are set to zero if no compile-time dimensions
                      were specified. The elements are default-constructed.
             */
            math_matrix()
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("math_matrix::math_matrix(" << getAddress(*this) << ") [1]");
            }

            /**
               \brief Creates a copy of the matrix `mat`.
               
               \param mat The matrix to copy from.
             */
            template<bool MO>
            math_matrix(const base_matrix<T, Rows, Cols, StorageTraits, MO> & mat)
                : base_matrix<T, Rows, Cols, StorageTraits, true>(mat)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("math_matrix::math_matrix(" << getAddress(*this) << "; " << getAddress(mat) << ") [2]");
            }
            
            /**
               \brief Creates a copy of the matrix `mat`.
               
               \param mat The matrix to copy from.
             */
            math_matrix(const math_matrix<T, Rows, Cols, StorageTraits> & mat)
                : base_matrix<T, Rows, Cols, StorageTraits, true>(mat)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("math_matrix::math_matrix(" << getAddress(*this) << "; " << getAddress(mat) << ") [3]");
            }
            
            /**
               \brief Creates a copy of the matrix `mat`.
               
               \param mat The matrix to copy from.
             */
            math_matrix(const base_matrix<T, Rows, Cols, StorageTraits, false> & mat)
                : base_matrix<T, Rows, Cols, StorageTraits, true>(mat)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("math_matrix::math_matrix(" << getAddress(*this) << "; " << getAddress(mat) << ") [3.b]");
            }
            
            /**
               \brief Creates a copy of the matrix `mat`.
               
               \param mat The matrix to copy from.
             */
            template<typename S, int Rs, int Cs, typename ST, bool MO>
            math_matrix(const base_matrix<S, Rs, Cs, ST, MO> & mat)
                : base_matrix<T, Rows, Cols, StorageTraits, true>(mat)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("math_matrix::math_matrix(" << getAddress(*this) << "; " << getAddress(mat) << ") [4]");
            }
            
            /**
               \brief Creates a copy of the matrix `mat`.
               
               \param mat The matrix to copy from.
             */
            template<typename MatrixType>
            math_matrix(const MatrixType & mat)
                : base_matrix<T, Rows, Cols, StorageTraits, true>(mat)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("math_matrix::math_matrix(" << getAddress(*this) << "; " << getAddress(mat) << ") [5]");
            }
            
            /**
               \brief Creates a vector with `entries` entries. Its coefficients are
                      default-constructed.
               
               Note that the number of rows or columns must be forced to be 1 during compile time;
               otherwise, this constructor will not compile.
               
               \param entries The number of entries.
             */
            math_matrix(size_type entries)
                : base_matrix<T, Rows, Cols, StorageTraits, true>(entries)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("math_matrix::math_matrix(" << getAddress(*this) << "; " << entries << ") [6]");
            }
            
            /**
               \brief Creates a vector with `entries` entries. Its coefficients are
                      default-constructed.
               
               Note that the number of rows or columns must be forced to be 1 during compile time;
               otherwise, this constructor will not compile.
               
               \param entries The number of entries.
             */
            math_matrix(int entries)
                : base_matrix<T, Rows, Cols, StorageTraits, true>(entries)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("math_matrix::math_matrix(" << getAddress(*this) << "; " << entries << ") [6b]");
            }
            
            /**
               \brief Creates a vector with `entries` entries. Its coefficients are copy-constructed
                      from `i.ref()`.
               
               Note that the number of rows or columns must be forced to be 1 during compile time;
               otherwise, this constructor will not compile.
               
               \param entries The number of entries.
               \param i       The element from which to copy-construct the coefficients.
             */
            template<typename S, bool def>
            math_matrix(size_type entries, const implementation::Initialize_Impl<S, def> & i)
                : base_matrix<T, Rows, Cols, StorageTraits, true>(entries, i)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("math_matrix::math_matrix(" << getAddress(*this) << "; " << entries << " " << getAddress(i.ref()) << ") [7]");
            }
            
            /**
               \brief Creates a matrix with `rows` times `cols` entries. Its coefficients are
                      default-constructed.
               
               \param rows The number of rows of the matrix.
               \param cols The number of columns of the matrix.
             */
            math_matrix(size_type rows, size_type cols)
                : base_matrix<T, Rows, Cols, StorageTraits, true>(rows, cols)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("math_matrix::math_matrix(" << getAddress(*this) << "; " << rows << " " << cols << ") [8]");
            }
            
            /**
               \brief Creates a matrix with `rows` times `cols` entries. Its coefficients are
                      copy-constructed from `i.ref()`.
               
               \param rows The number of rows of the matrix.
               \param cols The number of columns of the matrix.
               \param i    The element from which to copy-construct the coefficients.
             */
            template<typename S, bool def>
            math_matrix(size_type rows, size_type cols, const implementation::Initialize_Impl<S, def> & i)
                : base_matrix<T, Rows, Cols, StorageTraits, true>(rows, cols, i)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("math_matrix::math_matrix(" << getAddress(*this) << "; " << rows << " " << cols << " " << getAddress(i.ref()) << ") [9]");
            }
            
            /**
               \brief Releases the memory allocated for the matrix.
             */
            ~math_matrix() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("math_matrix::~math_matrix(" << getAddress(*this) << ")");
            }
            
            /**
               \brief Copies the contents of the matrix `m` into `*this`.
               
               \param m The matrix whose contents to copy into the current matrix.
             */
            template<typename MatrixType>
            math_matrix<T, Rows, Cols, StorageTraits> & operator = (const MatrixType & m)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("math_matrix::operator = (" << getAddress(*this) << "; " << getAddress(m) << ")");
                base_matrix<T, Rows, Cols, StorageTraits, true>::operator = (m);
                return *this;
            }
        };
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        //// base_rowvector
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        
        /**
           \brief Represents a row vector with coefficients in `T`.
           
           This is in fact the same as `base_matrix<T, 1, Cols, StorageTraits, MathObject>`.
           
           \tparam T              The coefficient type.
           \tparam Cols           If set to `plll::linalg::Flexible`, allows to flexibly choose
                                  and change the number of columns. Otherwise, enforces fixed number
                                  of columns during compile time. Default value is
                                  `plll::linalg::Flexible`.
           \tparam StorageTraits  Allows to determine how the vector's contents are stored. The
                                  default value is `plll::linalg::storage_traits<T>`, which should
                                  usually not be changed.
           \tparam MathObject     Allows to determine whether the vector is a math object (`true`)
                                  and thus operations can be performed on it, or whether it is a
                                  non-math object (`false`) which allows no operations except
                                  reading and storing coefficients and resizing.
         */
        template<typename T, int Cols = Flexible, typename StorageTraits = storage_traits<T>, bool MathObject = false>
        class base_rowvector : public base_matrix<T, 1, Cols, StorageTraits, MathObject>
        {
        public:
            enum { has_direct_access = base_matrix<T, 1, Cols, StorageTraits, MathObject>::has_direct_access };
            
            /**
               \brief Creates a new row vector. The dimensions are set to zero if no compile-time
                      dimensions were specified. The elements are default-constructed.
             */
            base_rowvector()
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_rowvector::base_rowvector(" << getAddress(*this) << ") [1]");
            }

            /**
               \brief Creates a copy of the row vector `mat`.
               
               \param mat The row vector to copy from.
             */
            template<bool MO>
            base_rowvector(const base_matrix<T, 1, Cols, StorageTraits, MO> & mat)
                : base_matrix<T, 1, Cols, StorageTraits, MathObject>(mat)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_rowvector::base_rowvector(" << getAddress(*this) << "; " << getAddress(mat) << ") [2]");
            }
            
            /**
               \brief Creates a copy of the row vector `mat`.
               
               \param mat The row vector to copy from.
             */
            base_rowvector(const base_rowvector<T, Cols, StorageTraits> & mat)
                : base_matrix<T, 1, Cols, StorageTraits, MathObject>(mat)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_rowvector::base_rowvector(" << getAddress(*this) << "; " << getAddress(mat) << ") [3]");
            }
            
            /**
               \brief Creates a copy of the row vector `mat`.
               
               \param mat The row vector to copy from.
             */
            template<typename S, int Rs, int Cs, typename ST, bool MO>
            base_rowvector(const base_matrix<S, Rs, Cs, ST, MO> & mat)
                : base_matrix<T, 1, Cols, StorageTraits, MathObject>(mat)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_rowvector::base_rowvector(" << getAddress(*this) << "; " << getAddress(mat) << ") [4]");
            }
            
            /**
               \brief Creates a copy of the row vector `mat`.
               
               \param mat The row vector to copy from.
             */
            template<typename MatrixType>
            base_rowvector(const MatrixType & mat)
                : base_matrix<T, 1, Cols, StorageTraits, MathObject>(mat)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_rowvector::base_rowvector(" << getAddress(*this) << "; " << getAddress(mat) << ") [5]");
            }
            
            /**
               \brief Creates a row vector with `entries` entries. Its coefficients are
                      default-constructed.
               
               \param entries The number of entries.
             */
            base_rowvector(size_type entries)
                : base_matrix<T, 1, Cols, StorageTraits, MathObject>(entries)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_rowvector::base_rowvector(" << getAddress(*this) << "; " << entries << ") [6]");
            }
            
            /**
               \brief Creates a row vector with `entries` entries. Its coefficients are
                      default-constructed.
               
               \param entries The number of entries.
             */
            base_rowvector(int entries)
                : base_matrix<T, 1, Cols, StorageTraits, MathObject>(entries)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_rowvector::base_rowvector(" << getAddress(*this) << "; " << entries << ") [6b]");
            }
            
            /**
               \brief Creates a row vector with `entries` entries. Its coefficients are
                      copy-constructed from `i.ref()`.
               
               \param entries The number of entries.
               \param i       The element from which to copy-construct the coefficients.
             */
            template<typename S, bool def>
            base_rowvector(size_type entries, const implementation::Initialize_Impl<S, def> & i)
                : base_matrix<T, 1, Cols, StorageTraits, MathObject>(entries, i)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_rowvector::base_rowvector(" << getAddress(*this) << "; " << entries << " " << getAddress(i.ref()) << ") [7]");
            }
            
            /**
               \brief Creates a row vector with `cols` entries. Its coefficients are
                      default-constructed.
               
               \param rows Must always be 1.
               \param cols The number of columns of the row vector.
             */
            base_rowvector(size_type rows, size_type cols)
                : base_matrix<T, 1, Cols, StorageTraits, MathObject>(rows, cols)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_rowvector::base_rowvector(" << getAddress(*this) << "; " << rows << " " << cols << ") [8]");
            }
            
            /**
               \brief Creates a row vector with `cols` entries. Its coefficients are
                      copy-constructed from `i.ref()`.
               
               \param rows Must always be 1.
               \param cols The number of columns of the row vector.
               \param i    The element from which to copy-construct the coefficients.
             */
            template<typename S, bool def>
            base_rowvector(size_type rows, size_type cols, const implementation::Initialize_Impl<S, def> & i)
                : base_matrix<T, 1, Cols, StorageTraits, MathObject>(rows, cols, i)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_rowvector::base_rowvector(" << getAddress(*this) << "; " << rows << " " << cols << " " << getAddress(i.ref()) << ") [9]");
            }
            
            /**
               \brief Releases the memory allocated for the row vector.
             */
            ~base_rowvector() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_rowvector::~base_rowvector(" << getAddress(*this) << ")");
            }
            
            /**
               \brief Copies the contents of the row vector `m` into `*this`.
               
               \param m The row vector whose contents to copy into the current row vector.
             */
            template<typename MatrixType>
            base_rowvector<T, Cols, StorageTraits, MathObject> & operator = (const MatrixType & m)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_rowvector::operator = (" << getAddress(*this) << "; " << getAddress(m) << ")");
                base_matrix<T, 1, Cols, StorageTraits, MathObject>::operator = (m);
                return *this;
            }
        };
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        //// base_colvector
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        
        /**
           \brief Represents a column vector with coefficients in `T`.
           
           This is in fact the same as `base_matrix<T, Rows, 1, StorageTraits, MathObject>`.
           
           \tparam T              The coefficient type.
           \tparam Rows           If set to `plll::linalg::Flexible`, allows to flexibly choose
                                  and change the number of rows. Otherwise, enforces fixed number
                                  of rows during compile time. Default value is
                                  `plll::linalg::Flexible`.
           \tparam StorageTraits  Allows to determine how the vector's contents are stored. The
                                  default value is `plll::linalg::storage_traits<T>`, which should
                                  usually not be changed.
           \tparam MathObject     Allows to determine whether the vector is a math object (`true`)
                                  and thus operations can be performed on it, or whether it is a
                                  non-math object (`false`) which allows no operations except
                                  reading and storing coefficients and resizing.
         */
        template<typename T, int Rows = Flexible, typename StorageTraits = storage_traits<T>, bool MathObject = false>
        class base_colvector : public base_matrix<T, Rows, 1, StorageTraits, MathObject>
        {
        public:
            enum { has_direct_access = base_matrix<T, Rows, 1, StorageTraits, MathObject>::has_direct_access };
            
            /**
               \brief Creates a new column vector. The dimensions are set to zero if no compile-time
                      dimensions were specified. The elements are default-constructed.
             */
            base_colvector()
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_colvector::base_colvector(" << getAddress(*this) << ") [1]");
            }

            /**
               \brief Creates a copy of the column vector `mat`.
               
               \param mat The column vector to copy from.
             */
            template<bool MO>
            base_colvector(const base_matrix<T, Rows, 1, StorageTraits, MO> & mat)
                : base_matrix<T, Rows, 1, StorageTraits, MathObject>(mat)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_colvector::base_colvector(" << getAddress(*this) << "; " << getAddress(mat) << ") [2]");
            }
            
            /**
               \brief Creates a copy of the column vector `mat`.
               
               \param mat The column vector to copy from.
             */
            base_colvector(const base_colvector<T, Rows, StorageTraits> & mat)
                : base_matrix<T, Rows, 1, StorageTraits, MathObject>(mat)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_colvector::base_colvector(" << getAddress(*this) << "; " << getAddress(mat) << ") [3]");
            }
            
            /**
               \brief Creates a copy of the column vector `mat`.
               
               \param mat The column vector to copy from.
             */
            template<typename S, int Rs, int Cs, typename ST, bool MO>
            base_colvector(const base_matrix<S, Rs, Cs, ST, MO> & mat)
                : base_matrix<T, Rows, 1, StorageTraits, MathObject>(mat)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_colvector::base_colvector(" << getAddress(*this) << "; " << getAddress(mat) << ") [4]");
            }
            
            /**
               \brief Creates a copy of the column vector `mat`.
               
               \param mat The column vector to copy from.
             */
            template<typename MatrixType>
            base_colvector(const MatrixType & mat)
                : base_matrix<T, Rows, 1, StorageTraits, MathObject>(mat)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_colvector::base_colvector(" << getAddress(*this) << "; " << getAddress(mat) << ") [5]");
            }
            
            /**
               \brief Creates a column vector with `entries` entries. Its coefficients are
                      default-constructed.
               
               \param entries The number of entries.
             */
            base_colvector(size_type entries)
                : base_matrix<T, Rows, 1, StorageTraits, MathObject>(entries)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_colvector::base_colvector(" << getAddress(*this) << "; " << entries << ") [6]");
            }
            
            /**
               \brief Creates a column vector with `entries` entries. Its coefficients are
                      default-constructed.
               
               \param entries The number of entries.
             */
            base_colvector(int entries)
                : base_matrix<T, Rows, 1, StorageTraits, MathObject>(entries)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_colvector::base_colvector(" << getAddress(*this) << "; " << entries << ") [6b]");
            }
            
            /**
               \brief Creates a column vector with `entries` entries. Its coefficients are
                      copy-constructed from `i.ref()`.
               
               \param entries The number of entries.
               \param i       The element from which to copy-construct the coefficients.
             */
            template<typename S, bool def>
            base_colvector(size_type entries, const implementation::Initialize_Impl<S, def> & i)
                : base_matrix<T, Rows, 1, StorageTraits, MathObject>(entries, i)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_colvector::base_colvector(" << getAddress(*this) << "; " << entries << " " << getAddress(i.ref()) << ") [7]");
            }
            
            /**
               \brief Creates a column vector with `cols` entries. Its coefficients are
                      default-constructed.
               
               \param rows Must always be 1.
               \param cols The number of columns of the column vector.
             */
            base_colvector(size_type rows, size_type cols)
                : base_matrix<T, Rows, 1, StorageTraits, MathObject>(rows, cols)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_colvector::base_colvector(" << getAddress(*this) << "; " << rows << " " << cols << ") [8]");
            }
            
            /**
               \brief Creates a column vector with `cols` entries. Its coefficients are
                      copy-constructed from `i.ref()`.
               
               \param rows Must always be 1.
               \param cols The number of columns of the column vector.
               \param i    The element from which to copy-construct the coefficients.
             */
            template<typename S, bool def>
            base_colvector(size_type rows, size_type cols, const implementation::Initialize_Impl<S, def> & i)
                : base_matrix<T, Rows, 1, StorageTraits, MathObject>(rows, cols, i)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_colvector::base_colvector(" << getAddress(*this) << "; " << rows << " " << cols << " " << getAddress(i.ref()) << ") [9]");
            }
            
            /**
               \brief Releases the memory allocated for the column vector.
             */
            ~base_colvector() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_colvector::~base_colvector(" << getAddress(*this) << ")");
            }
            
            /**
               \brief Copies the contents of the column vector `m` into `*this`.
               
               \param m The column vector whose contents to copy into the current column vector.
             */
            template<typename MatrixType>
            base_colvector<T, Rows, StorageTraits, MathObject> & operator = (const MatrixType & m)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("base_colvector::operator = (" << getAddress(*this) << "; " << getAddress(m) << ")");
                base_matrix<T, Rows, 1, StorageTraits, MathObject>::operator = (m);
                return *this;
            }
        };
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        //// math_rowvector
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        
        /**
           \brief Represents a math row vector with coefficients in `T`.
           
           This is in fact the same as `base_matrix<T, 1, Cols, StorageTraits, true>`.
           
           \tparam T              The coefficient type.
           \tparam Cols           If set to `plll::linalg::Flexible`, allows to flexibly choose
                                  and change the number of columns. Otherwise, enforces fixed number
                                  of columns during compile time. Default value is
                                  `plll::linalg::Flexible`.
           \tparam StorageTraits  Allows to determine how the vector's contents are stored. The
                                  default value is `plll::linalg::storage_traits<T>`, which should
                                  usually not be changed.
         */
        template<typename T, int Cols = Flexible, typename StorageTraits = storage_traits<T> >
        class math_rowvector : public base_rowvector<T, Cols, StorageTraits, true>
        {
        public:
            enum { has_direct_access = base_rowvector<T, Cols, StorageTraits, true>::has_direct_access };

            /**
               \brief Creates a new row vector. The dimensions are set to zero if no compile-time
                      dimensions were specified. The elements are default-constructed.
             */
            math_rowvector()
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("math_rowvector::math_rowvector(" << getAddress(*this) << ") [1]");
            }

            /**
               \brief Creates a copy of the row vector `mat`.
               
               \param mat The row vector to copy from.
             */
            template<bool MO>
            math_rowvector(const base_matrix<T, 1, Cols, StorageTraits, MO> & mat)
                : base_rowvector<T, Cols, StorageTraits, true>(mat)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("math_rowvector::math_rowvector(" << getAddress(*this) << "; " << getAddress(mat) << ") [2]");
            }
            
            /**
               \brief Creates a copy of the row vector `mat`.
               
               \param mat The row vector to copy from.
             */
            math_rowvector(const math_rowvector<T, Cols, StorageTraits> & mat)
                : base_rowvector<T, Cols, StorageTraits, true>(mat)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("math_rowvector::math_rowvector(" << getAddress(*this) << "; " << getAddress(mat) << ") [3]");
            }
            
            /**
               \brief Creates a copy of the row vector `mat`.
               
               \param mat The row vector to copy from.
             */
            template<typename S, int Rs, int Cs, typename ST, bool MO>
            math_rowvector(const base_matrix<S, Rs, Cs, ST, MO> & mat)
                : base_rowvector<T, Cols, StorageTraits, true>(mat)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("math_rowvector::math_rowvector(" << getAddress(*this) << "; " << getAddress(mat) << ") [4]");
            }
            
            /**
               \brief Creates a copy of the row vector `mat`.
               
               \param mat The row vector to copy from.
             */
            template<typename MatrixType>
            math_rowvector(const MatrixType & mat)
                : base_rowvector<T, Cols, StorageTraits, true>(mat)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("math_rowvector::math_rowvector(" << getAddress(*this) << "; " << getAddress(mat) << ") [5]");
            }
            
            /**
               \brief Creates a row vector with `entries` entries. Its coefficients are
                      default-constructed.
               
               \param entries The number of entries.
             */
            math_rowvector(size_type entries)
                : base_rowvector<T, Cols, StorageTraits, true>(entries)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("math_rowvector::math_rowvector(" << getAddress(*this) << "; " << entries << ") [6]");
            }
            
            /**
               \brief Creates a row vector with `entries` entries. Its coefficients are
                      default-constructed.
               
               \param entries The number of entries.
             */
            math_rowvector(int entries)
                : base_rowvector<T, Cols, StorageTraits, true>(entries)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("math_rowvector::math_rowvector(" << getAddress(*this) << "; " << entries << ") [6b]");
            }
            
            /**
               \brief Creates a row vector with `entries` entries. Its coefficients are
                      copy-constructed from `i.ref()`.
               
               \param entries The number of entries.
               \param i       The element from which to copy-construct the coefficients.
             */
            template<typename S, bool def>
            math_rowvector(size_type entries, const implementation::Initialize_Impl<S, def> & i)
                : base_rowvector<T, Cols, StorageTraits, true>(entries, i)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("math_rowvector::math_rowvector(" << getAddress(*this) << "; " << entries << " " << getAddress(i.ref()) << ") [7]");
            }
            
            /**
               \brief Creates a row vector with `cols` entries. Its coefficients are
                      default-constructed.
               
               \param rows Must always be 1.
               \param cols The number of columns of the row vector.
             */
            math_rowvector(size_type rows, size_type cols)
                : base_rowvector<T, Cols, StorageTraits, true>(rows, cols)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("math_rowvector::math_rowvector(" << getAddress(*this) << "; " << rows << " " << cols << ") [8]");
            }
            
            /**
               \brief Creates a row vector with `cols` entries. Its coefficients are
                      copy-constructed from `i.ref()`.
               
               \param rows Must always be 1.
               \param cols The number of columns of the row vector.
               \param i    The element from which to copy-construct the coefficients.
             */
            template<typename S, bool def>
            math_rowvector(size_type rows, size_type cols, const implementation::Initialize_Impl<S, def> & i)
                : base_rowvector<T, Cols, StorageTraits, true>(rows, cols, i)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("math_rowvector::math_rowvector(" << getAddress(*this) << "; " << rows << " " << cols << " " << getAddress(i.ref()) << ") [9]");
            }
            
            /**
               \brief Releases the memory allocated for the row vector.
             */
            ~math_rowvector() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("math_rowvector::~math_rowvector(" << getAddress(*this) << ")");
            }
            
            /**
               \brief Copies the contents of the row vector `m` into `*this`.
               
               \param m The row vector whose contents to copy into the current row vector.
             */
            template<typename MatrixType>
            math_rowvector<T, Cols, StorageTraits> & operator = (const MatrixType & m)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("math_rowvector::operator = (" << getAddress(*this) << "; " << getAddress(m) << ")");
                base_rowvector<T, Cols, StorageTraits, true>::operator = (m);
                return *this;
            }
        };
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        //// math_colvector
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        
        /**
           \brief Represents a math column vector with coefficients in `T`.
           
           This is in fact the same as `base_matrix<T, Rows, 1, StorageTraits, true>`.
           
           \tparam T              The coefficient type.
           \tparam Rows           If set to `plll::linalg::Flexible`, allows to flexibly choose
                                  and change the number of rows. Otherwise, enforces fixed number
                                  of rows during compile time. Default value is
                                  `plll::linalg::Flexible`.
           \tparam StorageTraits  Allows to determine how the vector's contents are stored. The
                                  default value is `plll::linalg::storage_traits<T>`, which should
                                  usually not be changed.
         */
        template<typename T, int Rows = Flexible, typename StorageTraits = storage_traits<T> >
        class math_colvector : public base_colvector<T, Rows, StorageTraits, true>
        {
        public:
            enum { has_direct_access = base_colvector<T, Rows, StorageTraits, true>::has_direct_access };

            /**
               \brief Creates a new column vector. The dimensions are set to zero if no compile-time
                      dimensions were specified. The elements are default-constructed.
             */
            math_colvector()
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("math_colvector::math_colvector(" << getAddress(*this) << ") [1]");
            }

            /**
               \brief Creates a copy of the column vector `mat`.
               
               \param mat The column vector to copy from.
             */
            template<bool MO>
            math_colvector(const base_matrix<T, Rows, 1, StorageTraits, MO> & mat)
                : base_colvector<T, Rows, StorageTraits, true>(mat)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("math_colvector::math_colvector(" << getAddress(*this) << "; " << getAddress(mat) << ") [2]");
            }
            
            /**
               \brief Creates a copy of the column vector `mat`.
               
               \param mat The column vector to copy from.
             */
            math_colvector(const math_colvector<T, Rows, StorageTraits> & mat)
                : base_colvector<T, Rows, StorageTraits, true>(mat)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("math_colvector::math_colvector(" << getAddress(*this) << "; " << getAddress(mat) << ") [3]");
            }
            
            /**
               \brief Creates a copy of the column vector `mat`.
               
               \param mat The column vector to copy from.
             */
            template<typename S, int Rs, int Cs, typename ST, bool MO>
            math_colvector(const base_matrix<S, Rs, Cs, ST, MO> & mat)
                : base_colvector<T, Rows, StorageTraits, true>(mat)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("math_colvector::math_colvector(" << getAddress(*this) << "; " << getAddress(mat) << ") [4]");
            }
            
            /**
               \brief Creates a copy of the column vector `mat`.
               
               \param mat The column vector to copy from.
             */
            template<typename MatrixType>
            math_colvector(const MatrixType & mat)
                : base_colvector<T, Rows, StorageTraits, true>(mat)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("math_colvector::math_colvector(" << getAddress(*this) << "; " << getAddress(mat) << ") [5]");
            }
            
            /**
               \brief Creates a column vector with `entries` entries. Its coefficients are
                      default-constructed.
               
               \param entries The number of entries.
             */
            math_colvector(size_type entries)
                : base_colvector<T, Rows, StorageTraits, true>(entries)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("math_colvector::math_colvector(" << getAddress(*this) << "; " << entries << ") [6]");
            }
            
            /**
               \brief Creates a column vector with `entries` entries. Its coefficients are
                      default-constructed.
               
               \param entries The number of entries.
             */
            math_colvector(int entries)
                : base_colvector<T, Rows, StorageTraits, true>(entries)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("math_colvector::math_colvector(" << getAddress(*this) << "; " << entries << ") [6b]");
            }
            
            /**
               \brief Creates a column vector with `entries` entries. Its coefficients are
                      copy-constructed from `i.ref()`.
               
               \param entries The number of entries.
               \param i       The element from which to copy-construct the coefficients.
             */
            template<typename S, bool def>
            math_colvector(size_type entries, const implementation::Initialize_Impl<S, def> & i)
                : base_colvector<T, Rows, StorageTraits, true>(entries, i)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("math_colvector::math_colvector(" << getAddress(*this) << "; " << entries << " " << getAddress(i.ref()) << ") [7]");
            }
            
            /**
               \brief Creates a column vector with `cols` entries. Its coefficients are
                      default-constructed.
               
               \param rows Must always be 1.
               \param cols The number of columns of the column vector.
             */
            math_colvector(size_type rows, size_type cols)
                : base_colvector<T, Rows, StorageTraits, true>(rows, cols)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("math_colvector::math_colvector(" << getAddress(*this) << "; " << rows << " " << cols << ") [8]");
            }
            
            /**
               \brief Creates a column vector with `cols` entries. Its coefficients are
                      copy-constructed from `i.ref()`.
               
               \param rows Must always be 1.
               \param cols The number of columns of the column vector.
               \param i    The element from which to copy-construct the coefficients.
             */
            template<typename S, bool def>
            math_colvector(size_type rows, size_type cols, const implementation::Initialize_Impl<S, def> & i)
                : base_colvector<T, Rows, StorageTraits, true>(rows, cols, i)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("math_colvector::math_colvector(" << getAddress(*this) << "; " << rows << " " << cols << " " << getAddress(i.ref()) << ") [9]");
            }
            
            /**
               \brief Releases the memory allocated for the column vector.
             */
            ~math_colvector() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("math_colvector::~math_colvector(" << getAddress(*this) << ")");
            }
            
            /**
               \brief Copies the contents of the column vector `m` into `*this`.
               
               \param m The column vector whose contents to copy into the current column vector.
             */
            template<typename MatrixType>
            math_colvector<T, Rows, StorageTraits> & operator = (const MatrixType & m)
            {
                PLLL_DEBUG_OUTPUT_MESSAGE("math_colvector::operator = (" << getAddress(*this) << "; " << getAddress(m) << ")");
                base_colvector<T, Rows, StorageTraits, true>::operator = (m);
                return *this;
            }
        };
#endif
        
        /**@{
           \name Stream input/output
        */
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        //// output
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        
        namespace implementation
        {
            template<typename MatrixType, int Rows, int Cols, bool forceGeneric>
            class MatrixPrinter
            {
            public:
                static void print(std::ostream & s, const MatrixType & mat)
                {
                    PLLL_DEBUG_OUTPUT_MESSAGE("MatrixPrinter::operator () (" << getAddress(s) << " " << getAddress(mat) << ") [0]");
                    s << "mat<" << mat.rows() << "," << mat.cols() << ">[";
                    bool rfirst = true;
                    for (typename MatrixType::ConstRowsEnumerator re = mat.enumerate_rows(); re.has_current(); re.next())
                    {
                        if (rfirst)
                            rfirst = false;
                        else
                            s << ", ";
                        s << "[";
                        bool first = true;
                        for (typename MatrixType::ConstRowsEnumerator::Type e = re.current(); e.has_current(); e.next())
                        {
                            if (first)
                                first = false;
                            else
                                s << ", ";
                            s << e.current();
                        }
                        s << "]";
                    }
                    s << "]";
                }
            };
            
            template<typename MatrixType, int Rows>
            class MatrixPrinter<MatrixType, Rows, 1, false>
            {
            public:
                static void print(std::ostream & s, const MatrixType & mat)
                {
                    PLLL_DEBUG_OUTPUT_MESSAGE("MatrixPrinter::operator () (" << getAddress(s) << " " << getAddress(mat) << ") [1]");
                    s << "vec<" << mat.rows() << ">[";
                    bool first = true;
                    for (typename MatrixType::ConstEnumerator e = mat.enumerate(); e.has_current(); e.next())
                    {
                        if (first)
                            first = false;
                        else
                            s << ", ";
                        s << e.current();
                    }
                    s << "]";
                }
            };
        }

        /**
           \brief Prints the given matrix `mat` to the output stream `s`.
           
           \tparam forceGeneric If `true`, always does the output in the `mat<r,c>[...]` format. If
                                `false`, uses the `vec<r>[...]` format for column vectors.
           \param s   An output stream.
           \param mat The matrix to print.
           
           \sa \ref matrixvector_io
          */
        template<typename MatrixType, bool forceGeneric>
        inline void print_matrix(std::ostream & s, const MatrixType & mat)
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("print_matrix(" << getAddress(s) << " " << getAddress(mat) << ")");
            PLLL_INTERNAL_STATIC_CHECK(implementation::MatrixInfo<MatrixType>::is_matrix, RequiresMatrix);
            implementation::MatrixPrinter<MatrixType, implementation::MatrixInfo<MatrixType>::rows, implementation::MatrixInfo<MatrixType>::cols, forceGeneric>::print(s, mat);
        }
        
        /**
           \brief Prints the given matrix `mat` to the output stream `s`.
           
           \param s            An output stream.
           \param mat          The matrix to print.
           \param forceGeneric If `true`, always does the output in the `mat<r,c>[...]` format. If
                               `false`, uses the `vec<r>[...]` format for column vectors. The
                               default value is `false`.
           
           \sa \ref matrixvector_io
          */
        template<typename MatrixType>
        inline void print_matrix(std::ostream & s, const MatrixType & mat, bool forceGeneric = false)
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("print_matrix(" << getAddress(s) << " " << getAddress(mat) << " " << forceGeneric << ")");
            if (forceGeneric)
                print_matrix<MatrixType, true>(s, mat);
            else
                print_matrix<MatrixType, false>(s, mat);
        }
        
        /**
           \brief Prints the given matrix `mat` to the output stream `s`.
           
           \param s            An output stream.
           \param mat          The matrix to print.
           
           \sa \ref matrixvector_io
          */
        template<typename T, int Rows, int Cols, typename ST, bool MO>
        inline std::ostream & operator << (std::ostream & s, const base_matrix<T, Rows, Cols, ST, MO> & mat)
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator << (" << getAddress(s) << " " << getAddress(mat) << ") [1]");
            print_matrix(s, mat, false);
            return s;
        }
        
        /**
           \brief Prints the given matrix `mat` to the output stream `s`.
           
           \param s            An output stream.
           \param mat          The matrix to print.
           
           \sa \ref matrixvector_io
          */
        template<template<typename DataType> class Operator, typename Data>
        inline std::ostream & operator << (std::ostream & s, const implementation::expressions::expr<Operator, Data> & mat)
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator << (" << getAddress(s) << " " << getAddress(mat) << ") [4]");
            print_matrix(s, mat, false);
            return s;
        }
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        //// input
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        
        namespace implementation
        {
            template<typename MatrixType, int Rows, int Cols>
                class MatrixScanner
            {
            public:
                template<typename S>
                static inline bool scan(std::istream & s, MatrixType & mat, const S & DefaultObject)
                {
                    PLLL_DEBUG_OUTPUT_MESSAGE("MatrixScanner::operator () (" << getAddress(s) << " " << getAddress(mat) << ") [0]");
                    if (!s) return false;
                    // Skip whitespace
                    char c = s.get();
                    while (std::isspace(c) && s)
                        c = s.get();
                    if (!s) return false;
                    // Parse matrix/vector
                    bool is_vector = false;
                    int mrows = -1, mcols = -1;
                    if (c == 'm')
                    {
                        // Scan "mat<"
                        c = s.get();
                        if (c != 'a')
                            return false;
                        c = s.get();
                        if (c != 't')
                            return false;
                        c = s.get();
                        if (c != '<')
                            return false;
                        c = s.get();
                        while (std::isspace(c) && s)
                            c = s.get();
                        if (!isdigit(c))
                            return false;
                        // Read number of rows
                        s.unget();
                        s >> mrows;
                        if (!s || (mrows < 0))
                            return false;
                        if ((Rows >= 0) && (mrows != Rows))
                            return false;
                        if (!implementation::MatrixInfo<MatrixType>::can_resize_rows && (mat.rows() != static_cast<size_type>(mrows)))
                            return false;
                        // Skip "," with surrounding whitespaces
                        c = s.get();
                        while (std::isspace(c) && s)
                            c = s.get();
                        if (c != ',')
                            return false;
                        c = s.get();
                        while (std::isspace(c) && s)
                            c = s.get();
                        // Read number of columns
                        s.unget();
                        s >> mcols;
                        if (!s || (mcols < 0))
                            return false;
                        if ((Cols >= 0) && (mcols != Cols))
                            return false;
                        if (!implementation::MatrixInfo<MatrixType>::can_resize_cols && (mat.cols() != static_cast<size_type>(mcols)))
                            return false;
                        // Read ">"
                        c = s.get();
                        while (std::isspace(c) && s)
                            c = s.get();
                        if (c != '>')
                            return false;
                        c = s.get();
                    }
                    else if (c == 'v')
                    {
                        is_vector = true;
                        // Is the given matrix object possibly a vector?
                        mcols = 1;
                        if ((Cols >= 0) && (mcols != Cols))
                            return false;
                        if (!implementation::MatrixInfo<MatrixType>::can_resize_cols && (mat.cols() != static_cast<size_type>(mcols)))
                            return false;
                        // Scan "vec<"
                        c = s.get();
                        if (c != 'e')
                            return false;
                        c = s.get();
                        if (c != 'c')
                            return false;
                        c = s.get();
                        if (c != '<')
                            return false;
                        c = s.get();
                        while (std::isspace(c) && s)
                            c = s.get();
                        // Scan number of elements
                        s.unget();
                        s >> mrows;
                        if (!s || (mrows < 0))
                            return false;
                        if ((Rows >= 0) && (mrows != Rows))
                            return false;
                        if (!implementation::MatrixInfo<MatrixType>::can_resize_rows && (mat.rows() != static_cast<size_type>(mrows)))
                            return false;
                        // Read ">"
                        c = s.get();
                        while (std::isspace(c) && s)
                            c = s.get();
                        if (c != '>')
                            return false;
                        c = s.get();
                    }
                    // We now scan (size-agnosticly) the matrix from the given input. We later compare
                    // the findings with the information extracted above.
                    size_type rows = 0, cols = 0;
                    // Now proceed with '['
                    while (std::isspace(c) && s)
                        c = s.get();
                    if ((c != '[') || !s)
                        return false;
                    c = s.get(); // remove "["
                    std::list<typename implementation::MatrixInfo<MatrixType>::Type> vals;
                    // Scan matrix elements
                    while (s)
                    {
                        // Find something non-trivial (or stop when finding ']')
                        while (std::isspace(c) && s)
                            c = s.get();
                        if (!s)
                            return false;
                        if (c == ']')
                            // last line?
                            break;
                        // Skip "," if avaiable
                        if (c == ',')
                            c = s.get();
                        while (std::isspace(c) && s)
                            c = s.get();
                        // Remove '[' (if not a vector!)
                        if (!is_vector)
                        {
                            if (c != '[')
                                return false;
                            c = s.get();
                        }
                        // remove numbers
                        unsigned i = 0;
                        ++rows;
                        while (s)
                        {
                            // Skip whitespace
                            while (std::isspace(c) && s)
                                c = s.get();
                            if (!s)
                                return false;
                            // Are we done? (if not a vector!)
                            if (!is_vector)
                                if ((c == ']') && (i == cols))
                                {
                                    c = s.get();
                                    break;
                                }
                            // Read value
                            s.unget();
                            vals.push_back(DefaultObject);
                            s >> vals.back();
                            if (!s)
                                return false;
                            c = s.get();
                            // Increase position/size
                            ++i;
                            if (rows == 1)
                                ++cols;
                            // Continue with next entry (if not a vector!)
                            if (is_vector)
                                break;
                            else
                            {
                                // Skip whitespace
                                while (std::isspace(c) && s)
                                    c = s.get();
                                // Scan "," if available
                                if (c == ',')
                                    c = s.get();
                                // Skip more whitespace
                                while (std::isspace(c) && s)
                                    c = s.get();
                                // Error?
                                if (!s)
                                    return false;
                            }
                        }
                        // Error?
                        if (!s)
                            return false;
                    }
                    // Error?
                    if (!s)
                        return false;
                    // Compare whether scanned data is consistent with known data
                    if (mrows >= 0)
                        if ((static_cast<size_type>(mrows) != rows) || (static_cast<size_type>(mcols) != cols))
                            return false;
                    if (!implementation::MatrixInfo<MatrixType>::can_resize_rows && (mat.rows() != rows))
                        return false;
                    if (!implementation::MatrixInfo<MatrixType>::can_resize_cols && (mat.cols() != cols))
                        return false;
                    // Assign matrix
                    typename std::list<typename implementation::MatrixInfo<MatrixType>::Type>::iterator it = vals.begin();
                    if ((mat.rows() != rows) || (mat.cols() != cols))
                        mat.resize(rows, cols, Initialize(DefaultObject));
                    for (typename MatrixType::Enumerator e = mat.enumerate(); e.has_current(); e.next(), ++it)
#if __cplusplus >= 201103L
                        e.current() = std::move(*it);
#else
                    e.current() = *it;
#endif
                    return true;
                }
            
                static inline bool scan(std::istream & s, MatrixType & mat)
                {
                    typename implementation::MatrixInfo<MatrixType>::Type def;
                    return scan(s, mat, def);
                }
            };
        }
        
        /**
           \brief Reads a matrix from the input stream `s` and stores it in `mat`.
           
           \param s   An input stream.
           \param mat Where to store the scanned matrix in.
           \return    Returns `true` if the matrix could be successfully read,
                      or `false if an error occured.
           
           \sa \ref matrixvector_io
          */
        template<typename MatrixType>
        inline bool scan_matrix(std::istream & s, MatrixType & mat)
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("scan_matrix(" << getAddress(s) << " " << getAddress(mat) << ")");
            PLLL_INTERNAL_STATIC_CHECK(implementation::MatrixInfo<MatrixType>::is_matrix, RequiresMatrix);
            PLLL_INTERNAL_STATIC_CHECK(!implementation::MatrixInfo<MatrixType>::is_const, RequiresNonconstMatrix);
            return implementation::MatrixScanner<MatrixType, implementation::MatrixInfo<MatrixType>::rows, implementation::MatrixInfo<MatrixType>::cols>::scan(s, mat);
        }
        
        /**
           \brief Reads a matrix from the input stream `s` and stores it in `mat`.
           
           \param s   An input stream.
           \param mat Where to store the scanned matrix in.
           \return    Returns `true` if the matrix could be successfully read,
                      or `false if an error occured.
           
           \sa \ref matrixvector_io
          */
        template<typename MatrixType>
        inline bool scan_matrix(std::istream & s, const MatrixType & mat)
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("scan_matrix(" << getAddress(s) << " " << getAddress(mat) << ")");
            PLLL_INTERNAL_STATIC_CHECK(implementation::MatrixInfo<const MatrixType>::is_matrix, RequiresMatrix);
            PLLL_INTERNAL_STATIC_CHECK(!implementation::MatrixInfo<const MatrixType>::is_const, RequiresNonconstMatrix);
            return implementation::MatrixScanner<const MatrixType, implementation::MatrixInfo<const MatrixType>::rows, implementation::MatrixInfo<const MatrixType>::cols>::scan(s, mat);
        }
        
        /**
           \brief Reads a matrix from the input stream `s` and stores it in `mat`.
           
           \param s              An input stream.
           \param mat            Where to store the scanned matrix in.
           \param default_object The default object to use to initialize new elements. Could be a
                                 context.
           \return               Returns `true` if the matrix could be successfully read, or
                                 `false` if an error occured.
           
           \sa \ref matrixvector_io
          */
        template<typename MatrixType, typename T>
        inline bool scan_matrix(std::istream & s, MatrixType & mat, const T & default_object)
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("scan_matrix(" << getAddress(s) << " " << getAddress(mat) << " " << getAddress(default_object) << ")");
            PLLL_INTERNAL_STATIC_CHECK(!implementation::MatrixInfo<MatrixType>::is_const, RequiresNonconstMatrix);
            PLLL_INTERNAL_STATIC_CHECK(implementation::MatrixInfo<MatrixType>::is_matrix, RequiresMatrix);
            return implementation::MatrixScanner<MatrixType, implementation::MatrixInfo<MatrixType>::rows, implementation::MatrixInfo<MatrixType>::cols>::scan(s, mat, default_object);
        }
        
        /**
           \brief Reads a matrix from the input stream `s` and stores it in `mat`.
           
           \param s              An input stream.
           \param mat            Where to store the scanned matrix in.
           \param default_object The default object to use to initialize new elements. Could be a
                                 context.
           \return               Returns `true` if the matrix could be successfully read, or
                                 `false` if an error occured.
           
           \sa \ref matrixvector_io
          */
        template<typename MatrixType, typename T>
        inline bool scan_matrix(std::istream & s, const MatrixType & mat, const T & default_object)
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("scan_matrix(" << getAddress(s) << " " << getAddress(mat) << " " << getAddress(default_object) << ")");
            PLLL_INTERNAL_STATIC_CHECK(!implementation::MatrixInfo<const MatrixType>::is_const, RequiresNonconstMatrix);
            PLLL_INTERNAL_STATIC_CHECK(implementation::MatrixInfo<const MatrixType>::is_matrix, RequiresMatrix);
            return implementation::MatrixScanner<const MatrixType, implementation::MatrixInfo<const MatrixType>::rows, implementation::MatrixInfo<const MatrixType>::cols>::scan(s, mat, default_object);
        }
        
        /**
           \brief Reads a matrix from the input stream `s` and stores it in `mat`.
           
           In case of an error, the `std::ios_base::bad_bit` flag of `s` is set.
           
           \param s   An input stream.
           \param mat Where to store the scanned matrix in.
           
           \sa \ref matrixvector_io
          */
        template<typename T, int Rows, int Cols, typename ST, bool MO>
        inline std::istream & operator >> (std::istream & s, base_matrix<T, Rows, Cols, ST, MO> & mat)
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator >> (" << getAddress(s) << " " << getAddress(mat) << ") [1]");
            if (!scan_matrix(s, mat))
                s.setstate(std::ios_base::failbit);
            return s;
        }
        
        /**
           \brief Reads a matrix from the input stream `s` and stores it in `mat`.
           
           In case of an error, the `std::ios_base::bad_bit` flag of `s` is set.
           
           \param s   An input stream.
           \param mat Where to store the scanned matrix in.
           
           \sa \ref matrixvector_io
          */
        template<template<typename DataType> class Op, typename Data>
        inline std::istream & operator >> (std::istream & s, const implementation::expressions::expr<Op, Data> & mat)
        {
            PLLL_DEBUG_OUTPUT_MESSAGE("operator >> (" << getAddress(s) << " " << getAddress(mat) << ") [2]");
            PLLL_INTERNAL_STATIC_CHECK((implementation::MatrixInfo<const implementation::expressions::expr<Op, Data> >::can_assign_to), RequiresAssignableMatrix);
            if (!scan_matrix(s, mat))
                s.setstate(std::ios_base::failbit);
            return s;
        }
        ///@}
    }
}

#include "matrix-ops.hpp"
#include "matrix-ops2.hpp"

#endif

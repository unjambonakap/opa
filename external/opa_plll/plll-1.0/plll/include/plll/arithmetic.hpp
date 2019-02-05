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

#ifndef PLLL_INCLUDE_GUARD__ARITHMETIC_HPP
#define PLLL_INCLUDE_GUARD__ARITHMETIC_HPP

#include <plll/config.hpp>
#include <plll/helper.hpp>
#include <sstream>
#include <iomanip>
#include <cassert>

/**
   \file
   \brief Main arithmetic header.
   
   This is the main arithmetic header for the `plll` lattice reduction library. It provides the
   basic definitions and functions for arbitrary precision integer arithmetic and arbitrary
   precision floating point arithmetic. Moreover, it provides a framework for converting between
   types, for determining the result types of arithmetic operations, and for handling of "global"
   invariants such as default precision for floating point values.
 */

namespace plll
{
    /**
       \brief Contains the arithmetic backend of `plll`.
       
       The `plll::arithmetic` namespace contains all arithmetic code of `plll`, such as contexts,
       expression templates, conversion facilities and traits.
    */
    namespace arithmetic
    {
        class Integer;
        class IntegerContext;
        class Real;
        class RealContext;
        
        namespace op
        {
            struct addition { };
            struct subtraction { };
            struct multiplication { };
            struct division { };
            struct modulo { };
        }
        
        namespace implementation
        {
            template<typename A, typename B, typename Op>
            /**
               \brief Provides information on the result of a binary arithmetic operation.
               
               This template is used to implement `plll::arithmetic::binary_operation<A, B, Op>`. It
               should not be used directly, and should be defined only for the types
               themselves. `plll::arithmetic::binary_operation<A, B, Op>` takes care to strip
               decorations from the types `A` and `B` (using the
               `plll::helper::remove_decorations<>` template) before using this template.
             */
            struct binary_operation_impl;
        }
        
        template<typename A, typename B, typename Op>
        /**
           \brief Provides information on the result of a binary arithmetic operation.
           
           See \ref arithctxts_resvals for documentation.
           
           \tparam A The first operand.
           \tparam B The second operand.
           \tparam Op The operator. Should be one of `op::addition`, `op::subtraction`,
                      `op::multiplication`, `op::division` and `op::modulo`.
         */
        struct binary_operation
        {
            enum
            {
                supported = implementation::binary_operation_impl<typename helper::remove_decorations<A>::Result,
                                                                  typename helper::remove_decorations<B>::Result, Op>::supported,
                /**< A boolean telling whether the current operation is supported, i.e. whether the
                     other fields and properties are defined. This should always be `true` except if
                     the instantiation is not defined. */
                intermediate_expression = implementation::binary_operation_impl<typename helper::remove_decorations<A>::Result,
                                                                                typename helper::remove_decorations<B>::Result, Op>::intermediate_expression
                /**< A boolean telling whether the result of the expression is an intermediate
                     expression, for example in case of expression templates, or already a final
                     value. In case this is `true`, `ResultType` and `IntermediateType` are
                     different types, and in case this is `false`, they are the same type. */
            };
            typedef typename implementation::binary_operation_impl<typename helper::remove_decorations<A>::Result,
                                                                   typename helper::remove_decorations<B>::Result, Op>::ResultType ResultType;
            /**< The final result type of the operation `op` involving two operands of type `A` and
                 `B`. */
            typedef typename implementation::binary_operation_impl<typename helper::remove_decorations<A>::Result,
                                                                   typename helper::remove_decorations<B>::Result, Op>::IntermediateType IntermediateType;
            /**< The intermediate result of the operation `op` involving the two operands of type
                 `A` and `B`. This will be the direct result type of `A op B`, and usually is a type
                 which contains almost no informations and which can be copied efficiently (as
                 opposed to `ResultType`, which can be expensive to copy). */
        };
        
        namespace op
        {
            struct negation { };
        }
        
        namespace implementation
        {
            template<typename A, typename Op>
            /**
               \brief Provides information on the result of an unary arithmetic operation.
               
               This template is used to implement `plll::arithmetic::unary_operation<A, Op>`. It
               should not be used directly, and should be defined only for the type
               itself. `plll::arithmetic::unary_operation<A, Op>` takes care to strip decorations
               from the type `A` (using the `plll::helper::remove_decorations<>` template) before
               using this template.
             */
            struct unary_operation_impl;
        }
        
        template<typename A, typename Op>
        /**
           \brief Provides information on the result of an unary arithmetic operation.
           
           See \ref arithctxts_resvals for documentation.
           
           \tparam A The operand.
           \tparam Op The operator. Should be `op::negation`.
         */
        struct unary_operation
        {
            enum
            {
                supported = implementation::unary_operation_impl<typename helper::remove_decorations<A>::Result, Op>::supported,
                /**< A boolean telling whether the current operation is supported, i.e. whether the
                     other fields and properties are defined. This should always be `true` except if
                     the instantiation is not defined. */
                intermediate_expression = implementation::unary_operation_impl<typename helper::remove_decorations<A>::Result, Op>::intermediate_expression
                /**< A boolean telling whether the result of the expression is an intermediate
                     expression, for example in case of expression templates, or already a final
                     value. In case this is `true`, `ResultType` and `IntermediateType` are
                     different types, and in case this is `false`, they are the same type. */
            };
            typedef typename implementation::unary_operation_impl<typename helper::remove_decorations<A>::Result, Op>::ResultType ResultType;
            /**< The final result type of the operation `op` involving one operand of type `A`. */
            typedef typename implementation::unary_operation_impl<typename helper::remove_decorations<A>::Result, Op>::IntermediateType IntermediateType;
            /**< The intermediate result of the operation `op` involving one operand of type
                 `A`. This will be the direct result type of `op A`, and usually is a type which
                 contains almost no informations and which can be copied efficiently (as opposed to
                 `ResultType`, which can be expensive to copy). */
        };
        
        class StringContext;
        
        namespace traits
        {
            template<class Type>
            /** \brief Provides information on arithmetic (and string) types.
                
                See \ref arithctxts_traits for documentation of the `type_traits<>` template.
                
                \tparam Type The arithmetic or string type we want to query information about.
             */
            struct type_traits;
            
            template<>
            struct type_traits<Integer>
            {
                enum
                {
                    is_number = true,
                    is_cputype = false,
                    is_realtype = false,
                    is_inttype = true,
                    is_string = false,
                    is_cpp_string = false,
                    is_c_string = false,
                    is_exact = true,
                    has_uniform_rng = true,
                    has_context = true,
                    is_native = false,
                    is_modulo = false,
                    has_infinity = false,
                    is_variable_precision = false,
                    has_squareroot = false,
                    has_full_power = false,
                    has_special_fns = false,
                    has_huge_exponent = false,
                    has_constants = false,
                    has_trigonometric = false
                };
                
                typedef IntegerContext Context;
                typedef Integer PromoteType;
                typedef const Integer & ConstReferenceType;
            };
            
            template<>
            struct type_traits<Real>
            {
                enum
                {
                    is_number = true,
                    is_cputype = false,
                    is_realtype = true,
                    is_inttype = false,
                    is_string = false,
                    is_cpp_string = false,
                    is_c_string = false,
                    is_exact = false,
                    has_uniform_rng = true,
                    has_context = true,
                    is_native = false,
                    is_modulo = false,
                    is_variable_precision = true,
                    has_squareroot = true,
                    has_full_power = true,
                    has_special_fns = true,
                    has_huge_exponent = true,
                    has_constants = true,
                    has_trigonometric = true
                };
                
                typedef IntegerContext Context;
                typedef Real PromoteType;
                typedef const Real & ConstReferenceType;
            };
            
            template<>
            struct type_traits<signed int>
            {
                enum
                {
                    is_number = true,
                    is_cputype = true,
                    is_realtype = false,
                    is_inttype = true,
                    is_string = false,
                    is_cpp_string = false,
                    is_c_string = false,
                    is_exact = true,
                    has_uniform_rng = false,
                    has_context = false,
                    is_native = true,
                    is_modulo = true,
                    has_infinity = false,
                    is_variable_precision = false,
                    has_squareroot = false,
                    has_full_power = false,
                    has_special_fns = false,
                    has_huge_exponent = false,
                    has_constants = false,
                    has_trigonometric = false
                };
                
                typedef signed long PromoteType;
                typedef signed int ConstReferenceType;
            };
            
            template<>
            struct type_traits<unsigned int>
            {
                enum
                {
                    is_number = true,
                    is_cputype = true,
                    is_realtype = false,
                    is_inttype = true,
                    is_string = false,
                    is_cpp_string = false,
                    is_c_string = false,
                    is_exact = true,
                    has_uniform_rng = false,
                    has_context = false,
                    is_native = true,
                    is_modulo = true,
                    has_infinity = false,
                    is_variable_precision = false,
                    has_squareroot = false,
                    has_full_power = false,
                    has_special_fns = false,
                    has_huge_exponent = false,
                    has_constants = false,
                    has_trigonometric = false
                };
                
                typedef unsigned long PromoteType;
                typedef unsigned int ConstReferenceType;
            };
            
            template<>
            struct type_traits<signed long>
            {
                enum
                {
                    is_number = true,
                    is_cputype = true,
                    is_realtype = false,
                    is_inttype = true,
                    is_string = false,
                    is_cpp_string = false,
                    is_c_string = false,
                    is_exact = true,
                    has_uniform_rng = false,
                    has_context = false,
                    is_native = true,
                    is_modulo = true,
                    has_infinity = false,
                    is_variable_precision = false,
                    has_squareroot = false,
                    has_full_power = false,
                    has_special_fns = false,
                    has_huge_exponent = false,
                    has_constants = false,
                    has_trigonometric = false
                };
                
                typedef signed long PromoteType;
                typedef signed long ConstReferenceType;
            };
            
            template<>
            struct type_traits<unsigned long>
            {
                enum
                {
                    is_number = true,
                    is_cputype = true,
                    is_realtype = false,
                    is_inttype = true,
                    is_string = false,
                    is_cpp_string = false,
                    is_c_string = false,
                    is_exact = true,
                    has_uniform_rng = false,
                    has_context = false,
                    is_native = true,
                    is_modulo = true,
                    has_infinity = false,
                    is_variable_precision = false,
                    has_squareroot = false,
                    has_full_power = false,
                    has_special_fns = false,
                    has_huge_exponent = false,
                    has_constants = false,
                    has_trigonometric = false
                };
                
                typedef unsigned long PromoteType;
                typedef unsigned long ConstReferenceType;
            };
            
            template<>
            struct type_traits<long long>
            {
                enum
                {
                    is_number = true,
                    is_cputype = true,
                    is_realtype = false,
                    is_inttype = true,
                    is_string = false,
                    is_cpp_string = false,
                    is_c_string = false,
                    is_exact = true,
                    has_uniform_rng = false,
                    has_context = false,
                    is_native = true,
                    is_modulo = true,
                    has_infinity = false,
                    is_variable_precision = false,
                    has_squareroot = false,
                    has_full_power = false,
                    has_special_fns = false,
                    has_huge_exponent = false,
                    has_constants = false,
                    has_trigonometric = false
                };
                
                typedef long long PromoteType;
                typedef long long ConstReferenceType;
            };
            
            template<>
            struct type_traits<float>
            {
                enum
                {
                    is_number = true,
                    is_cputype = true,
                    is_realtype = true,
                    is_inttype = false,
                    is_string = false,
                    is_cpp_string = false,
                    is_c_string = false,
                    is_exact = false,
                    has_uniform_rng = false,
                    has_context = false,
                    is_native = true,
                    is_modulo = false,
                    is_variable_precision = false,
                    has_squareroot = true,
                    has_full_power = true,
                    has_special_fns = false,
                    has_huge_exponent = false,
                    has_constants = false,
                    has_trigonometric = true
                };
                
                typedef double PromoteType;
                typedef float ConstReferenceType;
            };
            
            template<>
            struct type_traits<double>
            {
                enum
                {
                    is_number = true,
                    is_cputype = true,
                    is_realtype = true,
                    is_inttype = false,
                    is_string = false,
                    is_cpp_string = false,
                    is_c_string = false,
                    is_exact = false,
                    has_uniform_rng = false,
                    has_context = false,
                    is_native = true,
                    is_modulo = false,
                    is_variable_precision = false,
                    has_squareroot = true,
                    has_full_power = true,
                    has_special_fns = false,
                    has_huge_exponent = false,
                    has_constants = false,
                    has_trigonometric = true
                };
                
                typedef double PromoteType;
                typedef double ConstReferenceType;
            };
            
            template<>
            struct type_traits<long double>
            {
                enum
                {
                    is_number = true,
                    is_cputype = true,
                    is_realtype = true,
                    is_inttype = false,
                    is_string = false,
                    is_cpp_string = false,
                    is_c_string = false,
                    is_exact = false,
                    has_uniform_rng = false,
                    has_context = false,
                    is_native = true,
                    is_modulo = false,
                    is_variable_precision = false,
                    has_squareroot = true,
                    has_full_power = true,
                    has_special_fns = false,
                    has_huge_exponent = false,
                    has_constants = false,
                    has_trigonometric = true
                };
                
                typedef long double PromoteType;
                typedef long double ConstReferenceType;
            };
            
            template<>
            struct type_traits<std::string>
            // TODO: what about unicode strings (wchar_t in C++98, and new character types in
            // C++11)? ??? !!! ...
            {
                typedef StringContext Context;
                
                enum
                {
                    is_number = false,
                    is_realtype = false,
                    is_inttype = false,
                    is_string = true,
                    is_cpp_string = true,
                    is_c_string = false,
                    
                    has_context = true
                };
            };
            
            template<>
            struct type_traits<const char *>
            {
                enum
                {
                    is_number = false,
                    is_realtype = false,
                    is_inttype = false,
                    is_string = true,
                    is_cpp_string = false,
                    is_c_string = true
                };
            };
            
            template<>
            struct type_traits<char *>
            {
                enum
                {
                    is_number = false,
                    is_realtype = false,
                    is_inttype = false,
                    is_string = true,
                    is_cpp_string = false,
                    is_c_string = true
                };
            };
            
            template<unsigned n>
            struct type_traits<const char[n]>
            {
                enum
                {
                    is_number = false,
                    is_realtype = false,
                    is_inttype = false,
                    is_string = true,
                    is_cpp_string = false,
                    is_c_string = true
                };
            };
            
            template<unsigned n>
            struct type_traits<char[n]>
            {
                enum
                {
                    is_number = false,
                    is_realtype = false,
                    is_inttype = false,
                    is_string = true,
                    is_cpp_string = false,
                    is_c_string = true
                };
            };
        }
        
        using namespace traits;

        namespace implementation
        {
            #define PLLL_INTERNAL_DEFINE_OPERATION(A, B, R) template<typename Op> struct binary_operation_impl<A, B, Op> \
                { enum { supported = true, intermediate_expression = false }; typedef R ResultType; typedef R IntermediateType; }
            // signed int
            PLLL_INTERNAL_DEFINE_OPERATION(signed int, signed int,    signed int   );
            PLLL_INTERNAL_DEFINE_OPERATION(signed int, unsigned int,  unsigned int );
            PLLL_INTERNAL_DEFINE_OPERATION(signed int, signed long,   signed long  );
            PLLL_INTERNAL_DEFINE_OPERATION(signed int, unsigned long, unsigned long);
            PLLL_INTERNAL_DEFINE_OPERATION(signed int, long long,     long long    );
            PLLL_INTERNAL_DEFINE_OPERATION(signed int, float,         float        );
            PLLL_INTERNAL_DEFINE_OPERATION(signed int, double,        double       );
            PLLL_INTERNAL_DEFINE_OPERATION(signed int, long double,   long double  );
            // unsigned int
            PLLL_INTERNAL_DEFINE_OPERATION(unsigned int, signed int,    unsigned int );
            PLLL_INTERNAL_DEFINE_OPERATION(unsigned int, unsigned int,  unsigned int );
            PLLL_INTERNAL_DEFINE_OPERATION(unsigned int, signed long,   unsigned long);
            PLLL_INTERNAL_DEFINE_OPERATION(unsigned int, unsigned long, unsigned long);
            PLLL_INTERNAL_DEFINE_OPERATION(unsigned int, long long,     long long    );
            PLLL_INTERNAL_DEFINE_OPERATION(unsigned int, float,         float        );
            PLLL_INTERNAL_DEFINE_OPERATION(unsigned int, double,        double       );
            PLLL_INTERNAL_DEFINE_OPERATION(unsigned int, long double,   long double  );
            // signed long
            PLLL_INTERNAL_DEFINE_OPERATION(signed long, signed int,    signed long  );
            PLLL_INTERNAL_DEFINE_OPERATION(signed long, unsigned int,  unsigned long);
            PLLL_INTERNAL_DEFINE_OPERATION(signed long, signed long,   signed long  );
            PLLL_INTERNAL_DEFINE_OPERATION(signed long, unsigned long, unsigned long);
            PLLL_INTERNAL_DEFINE_OPERATION(signed long, long long,     long long    );
            PLLL_INTERNAL_DEFINE_OPERATION(signed long, float,         float        );
            PLLL_INTERNAL_DEFINE_OPERATION(signed long, double,        double       );
            PLLL_INTERNAL_DEFINE_OPERATION(signed long, long double,   long double  );
            // unsigned long
            PLLL_INTERNAL_DEFINE_OPERATION(unsigned long, signed int,    unsigned long);
            PLLL_INTERNAL_DEFINE_OPERATION(unsigned long, unsigned int,  unsigned long);
            PLLL_INTERNAL_DEFINE_OPERATION(unsigned long, signed long,   unsigned long);
            PLLL_INTERNAL_DEFINE_OPERATION(unsigned long, unsigned long, unsigned long);
            PLLL_INTERNAL_DEFINE_OPERATION(unsigned long, long long,     long long    );
            PLLL_INTERNAL_DEFINE_OPERATION(unsigned long, float,         float        );
            PLLL_INTERNAL_DEFINE_OPERATION(unsigned long, double,        double       );
            PLLL_INTERNAL_DEFINE_OPERATION(unsigned long, long double,   long double  );
            // long long
            PLLL_INTERNAL_DEFINE_OPERATION(long long, signed int,    long long  );
            PLLL_INTERNAL_DEFINE_OPERATION(long long, unsigned int,  long long  );
            PLLL_INTERNAL_DEFINE_OPERATION(long long, signed long,   long long  );
            PLLL_INTERNAL_DEFINE_OPERATION(long long, unsigned long, long long  );
            PLLL_INTERNAL_DEFINE_OPERATION(long long, long long,     long long  );
            PLLL_INTERNAL_DEFINE_OPERATION(long long, float,         float      );
            PLLL_INTERNAL_DEFINE_OPERATION(long long, double,        double     );
            PLLL_INTERNAL_DEFINE_OPERATION(long long, long double,   long double);
            // float
            PLLL_INTERNAL_DEFINE_OPERATION(float, signed int,      float      );
            PLLL_INTERNAL_DEFINE_OPERATION(float, unsigned int,    float      );
            PLLL_INTERNAL_DEFINE_OPERATION(float, signed long,     float      );
            PLLL_INTERNAL_DEFINE_OPERATION(float, unsigned long,   float      );
            PLLL_INTERNAL_DEFINE_OPERATION(float, long long,       float      );
            PLLL_INTERNAL_DEFINE_OPERATION(float, float,           float      );
            PLLL_INTERNAL_DEFINE_OPERATION(float, double,          double     );
            PLLL_INTERNAL_DEFINE_OPERATION(float, long double,     long double);
            // double
            PLLL_INTERNAL_DEFINE_OPERATION(double, signed int,     double     );
            PLLL_INTERNAL_DEFINE_OPERATION(double, unsigned int,   double     );
            PLLL_INTERNAL_DEFINE_OPERATION(double, signed long,    double     );
            PLLL_INTERNAL_DEFINE_OPERATION(double, unsigned long,  double     );
            PLLL_INTERNAL_DEFINE_OPERATION(double, long long,      double     );
            PLLL_INTERNAL_DEFINE_OPERATION(double, float,          double     );
            PLLL_INTERNAL_DEFINE_OPERATION(double, double,         double     );
            PLLL_INTERNAL_DEFINE_OPERATION(double, long double,    long double);
            // long double
            PLLL_INTERNAL_DEFINE_OPERATION(long double, signed int,    long double);
            PLLL_INTERNAL_DEFINE_OPERATION(long double, unsigned int,  long double);
            PLLL_INTERNAL_DEFINE_OPERATION(long double, signed long,   long double);
            PLLL_INTERNAL_DEFINE_OPERATION(long double, unsigned long, long double);
            PLLL_INTERNAL_DEFINE_OPERATION(long double, long long,     long double);
            PLLL_INTERNAL_DEFINE_OPERATION(long double, float,         long double);
            PLLL_INTERNAL_DEFINE_OPERATION(long double, double,        long double);
            PLLL_INTERNAL_DEFINE_OPERATION(long double, long double,   long double);
            #undef PLLL_INTERNAL_DEFINE_OPERATION
            
            #define PLLL_INTERNAL_DEFINE_OPERATION(A, R) template<> struct unary_operation_impl<A, op::negation> \
                { enum { supported = true, intermediate_expression = false }; typedef R ResultType; typedef R IntermediateType; }
            PLLL_INTERNAL_DEFINE_OPERATION(signed int,    signed int   );
            PLLL_INTERNAL_DEFINE_OPERATION(unsigned int,  unsigned int );
            PLLL_INTERNAL_DEFINE_OPERATION(signed long,   signed long  );
            PLLL_INTERNAL_DEFINE_OPERATION(unsigned long, unsigned long);
            PLLL_INTERNAL_DEFINE_OPERATION(long long,     long long    );
            PLLL_INTERNAL_DEFINE_OPERATION(float,         float        );
            PLLL_INTERNAL_DEFINE_OPERATION(double,        double       );
            PLLL_INTERNAL_DEFINE_OPERATION(long double,   long double  );
            #undef PLLL_INTERNAL_DEFINE_OPERATION
            
            template<class SourceType, class DestContext>
            /**
               \brief Provides conversion implementation from type `SourceType` to
                      `DestContext::Type`.
               
               Every instantiation should provide the following member types:
               
               - `RetVal`: return value for the `convert()` function.
               - `RetVal_Frac`: return value for the `convert_frac()` function. Should be `void` if
                 `DestContext::is_real` is `false`.
               - `RetVal_Floor`: return value for the `convert_floor()` function. Should be `void`
                 if `DestContext::is_int` is `false` or `type_traits<SourceType>::is_real` is
                 `false`
               - `RetVal_Round`: return value for the `convert_round()` function. Should be `void`
                 if `DestContext::is_int` is `false` or `type_traits<SourceType>::is_real` is
                 `false`
               - `RetVal_Round2`: return value for the `convert_round()` function which also
                 retrieves rounding information. Should be `void` if `DestContext::is_int` is
                 `false` or `type_traits<SourceType>::is_real` is `false`
               - `RetVal_Ceil`: return value for the `convert_ceil()` function. Should be `void` if
                 `DestContext::is_int` is `false` or `type_traits<SourceType>::is_real` is `false`
               
               The following static functions should always be provided:
               
               - `void convert(typename Context::Type &, const typename Context::Type &, const Context &)`
               - `RetVal convert(const typename Context::Type &, const Context &)`
               
               The following static functions should be provided if `DestContext::is_real` is true:
               
               - `void convert_frac(typename Context::Type &, const typename Context::Type &, const
                 typename Context::Type &, const Context &)`
               - `RetVal_Frac convert_frac(const typename Context::Type &, const typename
                 Context::Type &, const Context &)`
               
               The following static functions should be provided if both `DestContext::is_int` and
               `type_traits<SourceType>::is_real` are `true`:
               
               - `void floor(typename Context::Type &, const typename Context::Type &, const Context
                 &)`
               - `RetVal_Floor floor(const typename Context::Type &, const Context &)`
               - `void round(typename Context::Type &, const typename Context::Type &, const Context
                  &)`
               - `RetVal_Round round(const typename Context::Type &, const Context &)`
               - `void round(typename Context::Type &, const typename Context::Type &, bool &, const
                  Context &)`
               - `RetVal_Round2 round(const typename Context::Type &, bool &, const Context &)`
               - `void ceil(typename Context::Type &, const typename Context::Type &, const Context
                  &)`
               - `RetVal_Ceil ceil(const typename Context::Type &, const Context &)`
             */
            class conversion_impl;
        }
        
        /**@{
           \name Conversion functions.
        */
        
        template<class SourceType, class DestContext>
        /**
           \brief Converts the value in `src` to the type described by `context` and stores the
                  result in `dst`.
           
           \tparam SourceType  The source type.
           \tparam DestContext The type of the arithmetic context for the destination.
           \param dst     The destination variable.
           \param src     The source variable. Can also be a C or C++ string.
           \param context The arithmetic context for the destination.
           
           \sa \ref arithctxts_convs_w_context
           \sa \ref arithctxts_convs_string
         */
        void convert(typename DestContext::Type & dst, const SourceType & src, const DestContext & context)
        {
            implementation::conversion_impl<SourceType, DestContext>::convert(dst, src, context);
        }
        
        template<class SourceType, class DestContext>
        /**
           \brief Converts the value in `src` to the type described by `context` and returns the
                  result.
           
           \tparam SourceType  The source type.
           \tparam DestContext The type of the arithmetic context for the destination.
           \return        The result of the conversion.
           \param src     The source variable. Can also be a C or C++ string.
           \param context The arithmetic context for the destination.
           
           \sa \ref arithctxts_convs_w_context
           \sa \ref arithctxts_convs_string
         */
        typename implementation::conversion_impl<SourceType, DestContext>::RetVal convert(const SourceType & src, const DestContext & context)
        {
            return implementation::conversion_impl<SourceType, DestContext>::convert(src, context);
        }
        
        template<class SourceType, class DestContext>
        /**
           \brief Converts the values in `src1` and `scr2` to the type described by `context` and
                  stores their quotient in `dst`.
           
           This method is only allowed in case `context` a real context.
           
           \tparam SourceType  The source type.
           \tparam DestContext The type of the arithmetic context for the destination.
           \param dst     The destination variable.
           \param src1    The source variable for the numerator.
           \param src2    The source variable for the denominator.
           \param context The arithmetic context for the destination.
           
           \sa \ref arithctxts_convs_w_context
         */
        void convert_fraction(typename DestContext::Type & dst, const SourceType & src1, const SourceType & src2,
                              const DestContext & context)
        {
            PLLL_INTERNAL_STATIC_CHECK(DestContext::is_realtype, NeedsRealTypeDestination);
            implementation::conversion_impl<SourceType, DestContext>::convert_frac(dst, src1, src2, context);
        }
        
        template<class SourceType, class DestContext>
        /**
           \brief Converts the values in `src1` and `scr2` to the type described by `context` and
                  returns their quotient.
           
           This method is only allowed in case `context` a real context.
           
           \tparam SourceType  The source type.
           \tparam DestContext The type of the arithmetic context for the destination.
           \return        The result of the conversion.
           \param src1    The source variable for the numerator.
           \param src2    The source variable for the denominator.
           \param context The arithmetic context for the destination.
           
           \sa \ref arithctxts_convs_w_context
         */
        typename implementation::conversion_impl<SourceType, DestContext>::RetVal_Frac convert_fraction(const SourceType & src1, const SourceType & src2,
                                                                                                        const DestContext & context)
        {
            PLLL_INTERNAL_STATIC_CHECK(DestContext::is_realtype, NeedsRealTypeDestination);
            return implementation::conversion_impl<SourceType, DestContext>::convert_frac(src1, src2, context);
        }
        
        template<class SourceType, class DestContext>
        /**
           \brief Converts the floor of the value in `src` to the type described by `context` and
                  stores the result in `dst`.
           
           This method is only allowed in case `src` is a real type and `dst` an integer type.
           
           \tparam SourceType  The source type.
           \tparam DestContext The type of the arithmetic context for the destination.
           \param dst     The destination variable.
           \param src     The source variable.
           \param context The arithmetic context for the destination.
           
           \sa \ref arithctxts_convs_w_context
         */
        void convert_floor(typename DestContext::Type & dst, const SourceType & src, const DestContext & context)
        {
            PLLL_INTERNAL_STATIC_CHECK(traits::type_traits<SourceType>::is_realtype, NeedsRealTypeSource);
            PLLL_INTERNAL_STATIC_CHECK(DestContext::is_inttype, NeedsIntegerTypeDestination);
            implementation::conversion_impl<SourceType, DestContext>::floor(dst, src, context);
        }
        
        template<class SourceType, class DestContext>
        /**
           \brief Converts the floor of the value in `src` to the type described by `context` and
                  returns the result.
           
           This method is only allowed in case `src` is a real type and `context` an integer context.
           
           \tparam SourceType  The source type.
           \tparam DestContext The type of the arithmetic context for the destination.
           \return        The result of the conversion.
           \param src     The source variable.
           \param context The arithmetic context for the destination.
           
           \sa \ref arithctxts_convs_w_context
         */
        typename implementation::conversion_impl<SourceType, DestContext>::RetVal_Floor convert_floor(const SourceType & src, const DestContext & context)
        {
            PLLL_INTERNAL_STATIC_CHECK(traits::type_traits<SourceType>::is_realtype, NeedsRealTypeSource);
            PLLL_INTERNAL_STATIC_CHECK(DestContext::is_inttype, NeedsIntegerTypeDestination);
            return implementation::conversion_impl<SourceType, DestContext>::floor(src, context);
        }
        
        template<class SourceType, class DestContext>
        /**
           \brief Converts the rounded value in `src` to the type described by `context` and stores
                  the result in `dst`.
           
           This method is only allowed in case `src` is a real type and `dst` an integer type.
           
           \tparam SourceType  The source type.
           \tparam DestContext The type of the arithmetic context for the destination.
           \param dst     The destination variable.
           \param src     The source variable.
           \param context The arithmetic context for the destination.
           
           \sa \ref arithctxts_convs_w_context
         */
        void convert_round(typename DestContext::Type & dst, const SourceType & src, const DestContext & context)
        {
            PLLL_INTERNAL_STATIC_CHECK(traits::type_traits<SourceType>::is_realtype, NeedsRealTypeSource);
            PLLL_INTERNAL_STATIC_CHECK(DestContext::is_inttype, NeedsIntegerTypeDestination);
            implementation::conversion_impl<SourceType, DestContext>::round(dst, src, context);
        }
        
        template<class SourceType, class DestContext>
        /**
           \brief Converts the rounded value in `src` to the type described by `context` and returns
                  the result.
           
           This method is only allowed in case `src` is a real type and `context` an integer context.
           
           \tparam SourceType  The source type.
           \tparam DestContext The type of the arithmetic context for the destination.
           \return        The result of the conversion.
           \param src     The source variable.
           \param context The arithmetic context for the destination.
           
           \sa \ref arithctxts_convs_w_context
         */
        typename implementation::conversion_impl<SourceType, DestContext>::RetVal_Round convert_round(const SourceType & src, const DestContext & context)
        {
            PLLL_INTERNAL_STATIC_CHECK(traits::type_traits<SourceType>::is_realtype, NeedsRealTypeSource);
            PLLL_INTERNAL_STATIC_CHECK(DestContext::is_inttype, NeedsIntegerTypeDestination);
            return implementation::conversion_impl<SourceType, DestContext>::round(src, context);
        }
        
        template<class SourceType, class DestContext>
        /**
           \brief Converts the rounded value in `src` to the type described by `context` and stores
                  the result in `dst`. Stores in `up` whether the value was rounded up or down.
           
           This method is only allowed in case `src` is a real type and `dst` an integer type.
           
           \tparam SourceType  The source type.
           \tparam DestContext The type of the arithmetic context for the destination.
           \param dst     The destination variable.
           \param src     The source variable.
           \param up      Will be set to `true` if `src` was rounded up, or `false` if it was
                          rounded down.
           \param context The arithmetic context for the destination.
           
           \sa \ref arithctxts_convs_w_context
         */
        void convert_round(typename DestContext::Type & dst, const SourceType & src, bool & up, const DestContext & context)
        {
            PLLL_INTERNAL_STATIC_CHECK(traits::type_traits<SourceType>::is_realtype, NeedsRealTypeSource);
            PLLL_INTERNAL_STATIC_CHECK(DestContext::is_inttype, NeedsIntegerTypeDestination);
            implementation::conversion_impl<SourceType, DestContext>::round(dst, src, up, context);
        }
        
        template<class SourceType, class DestContext>
        /**
           \brief Converts the rounded value in `src` to the type described by `context` and returns
                  the result. Stores in `up` whether the value was rounded up or down.
           
           This method is only allowed in case `src` is a real type and `context` an integer context.
           
           \tparam SourceType  The source type.
           \tparam DestContext The type of the arithmetic context for the destination.
           \return        The result of the conversion.
           \param src     The source variable.
           \param up      Will be set to `true` if `src` was rounded up, or `false` if it was
                          rounded down.
           \param context The arithmetic context for the destination.
           
           \sa \ref arithctxts_convs_w_context
         */
        typename implementation::conversion_impl<SourceType, DestContext>::RetVal_Round2 convert_round(const SourceType & src, bool & up, const DestContext & context)
        {
            PLLL_INTERNAL_STATIC_CHECK(traits::type_traits<SourceType>::is_realtype, NeedsRealTypeSource);
            PLLL_INTERNAL_STATIC_CHECK(DestContext::is_inttype, NeedsIntegerTypeDestination);
            return implementation::conversion_impl<SourceType, DestContext>::round(src, up, context);
        }
        
        template<class SourceType, class DestContext>
        /**
           \brief Converts the ceil of the value in `src` to the type described by `context` and
                  stores the result in `dst`.
           
           This method is only allowed in case `src` is a real type and `dst` an integer type.
           
           \tparam SourceType  The source type.
           \tparam DestContext The type of the arithmetic context for the destination.
           \param dst     The destination variable.
           \param src     The source variable.
           \param context The arithmetic context for the destination.
           
           \sa \ref arithctxts_convs_w_context
         */
        void convert_ceil(typename DestContext::Type & dst, const SourceType & src, const DestContext & context)
        {
            PLLL_INTERNAL_STATIC_CHECK(traits::type_traits<SourceType>::is_realtype, NeedsRealTypeSource);
            PLLL_INTERNAL_STATIC_CHECK(DestContext::is_inttype, NeedsIntegerTypeDestination);
            implementation::conversion_impl<SourceType, DestContext>::ceil(dst, src, context);
        }
        
        template<class SourceType, class DestContext>
        /**
           \brief Converts the ceil of the value in `src` to the type described by `context` and
                  returns the result in `dst`.
           
           This method is only allowed in case `src` is a real type and `context` an integer context.
           
           \tparam SourceType  The source type.
           \tparam DestContext The type of the arithmetic context for the destination.
           \return        The result of the conversion.
           \param src     The source variable.
           \param context The arithmetic context for the destination.
           
           \sa \ref arithctxts_convs_w_context
         */
        typename implementation::conversion_impl<SourceType, DestContext>::RetVal_Ceil convert_ceil(const SourceType & src, const DestContext & context)
        {
            PLLL_INTERNAL_STATIC_CHECK(traits::type_traits<SourceType>::is_realtype, NeedsRealTypeSource);
            PLLL_INTERNAL_STATIC_CHECK(DestContext::is_inttype, NeedsIntegerTypeDestination);
            return implementation::conversion_impl<SourceType, DestContext>::ceil(src, context);
        }
        
        namespace implementation
        {
            template<class SourceType>
            /**
               \brief Implementation backend for conversions from context types to native types.
               
               Provides functions to convert the type `SourceType` to all native types. This
               template is used by `plll::arithmetic::implementation::nativeconversion_impl2<>` to
               do the conversions.
             */
            class nativeconversion_impl;
            
            template<class Native, class SourceType, bool NativeInt, bool SourceInt>
            /**
               \brief Implementation backend for native conversion functions.
               
               Provides functions to convert the type `SourceType` to the given native type
               `Native`. The booleans `NativeInt` and `SourceInt` indicate whether the native
               respectively source type are integer types.
               
               It uses `nativeconversion_impl<SourceType>`, which provides functions for all
               supported native types.
               
               Note that for `Native == Integer`, it supports conversions of any type to arbitrary
               precision integers (using
               `plll::arithmetic::implementation::nativeconversion_impl2<Integer, Source, true,
               SourceInt>`). Moreover, for `Native == std::string`, it supports type to string
               conversions (using
               `plll::arithmetic::implementation::nativeconversion_impl2<std::string, SourceType,
               true, SourceInt>`).
             */
            class nativeconversion_impl2;
        }

        template<class Native, class SourceType>
        typename implementation::nativeconversion_impl2<Native, SourceType, traits::type_traits<Native>::is_inttype || !traits::type_traits<Native>::is_number,
                                                        traits::type_traits<SourceType>::is_inttype || !traits::type_traits<SourceType>::is_number>::RetVal
        /**
           \brief Converts the value in `src` to the native type `Native` and returns the result.
           
           \tparam Native     The native type to convert to. Can also be `std::string`.
           \tparam SourceType The source type.
           \param src The source variable.
           \return The converted result of type `Native`.
           
           \sa \ref arithctxts_convs_native
           \sa \ref arithctxts_convs_string
         */
        convert(const SourceType & src)
        {
            return implementation::nativeconversion_impl2<Native, SourceType, traits::type_traits<Native>::is_inttype || !traits::type_traits<Native>::is_number,
                                                          traits::type_traits<SourceType>::is_inttype || !traits::type_traits<SourceType>::is_number>::convert(src);
        }
        
        template<class Native, class SourceType>
        typename implementation::nativeconversion_impl2<Native, SourceType, traits::type_traits<Native>::is_inttype || !traits::type_traits<Native>::is_number,
                                                        traits::type_traits<SourceType>::is_inttype || !traits::type_traits<SourceType>::is_number>::RetVal_Floor
        /**
           \brief Converts the floor of the value in `src` to the native type `Native` and returns
                  the result.
           
           \tparam Native     The native type to convert to.
           \tparam SourceType The source type.
           \param src The source variable.
           \return The converted result.
           
           \sa \ref arithctxts_convs_native
         */
        convert_floor(const SourceType & src)
        {
            PLLL_INTERNAL_STATIC_CHECK(traits::type_traits<SourceType>::is_realtype, NeedsRealTypeSource);
            PLLL_INTERNAL_STATIC_CHECK(traits::type_traits<Native>::is_inttype, NeedsIntegerTypeDestination);
            return implementation::nativeconversion_impl2<Native, SourceType, traits::type_traits<Native>::is_inttype || !traits::type_traits<Native>::is_number,
                                                          traits::type_traits<SourceType>::is_inttype || !traits::type_traits<SourceType>::is_number>::convert_floor(src);
        }
        
        template<class Native, class SourceType>
        /**
           \brief Converts the rounded value of `src` to the native type `Native` and returns the
                  result.
           
           \tparam Native     The native type to convert to.
           \tparam SourceType The source type.
           \param src The source variable.
           \return The converted result.
           
           \sa \ref arithctxts_convs_native
         */
        typename implementation::nativeconversion_impl2<Native, SourceType, traits::type_traits<Native>::is_inttype || !traits::type_traits<Native>::is_number,
                                                        traits::type_traits<SourceType>::is_inttype || !traits::type_traits<SourceType>::is_number>::RetVal_Round
        convert_round(const SourceType & src)
        {
            PLLL_INTERNAL_STATIC_CHECK(traits::type_traits<SourceType>::is_realtype, NeedsRealTypeSource);
            PLLL_INTERNAL_STATIC_CHECK(traits::type_traits<Native>::is_inttype, NeedsIntegerTypeDestination);
            return implementation::nativeconversion_impl2<Native, SourceType, traits::type_traits<Native>::is_inttype || !traits::type_traits<Native>::is_number,
                                                          traits::type_traits<SourceType>::is_inttype || !traits::type_traits<SourceType>::is_number>::convert_round(src);
        }
        
        template<class Native, class SourceType>
        /**
           \brief Converts the rounded value of `src` to the native type `Native` and returns the
                  result. Returns whether the value was rounded up or not in `up`.
           
           \tparam Native     The native type to convert to.
           \tparam SourceType The source type.
           \param src The source variable.
           \param up  Will be set to `true` if `src` was rounded up, or `false` if it was
                      rounded down.
           \return The converted result.
           
           \sa \ref arithctxts_convs_native
         */
        typename implementation::nativeconversion_impl2<Native, SourceType, traits::type_traits<Native>::is_inttype || !traits::type_traits<Native>::is_number,
                                                        traits::type_traits<SourceType>::is_inttype || !traits::type_traits<SourceType>::is_number>::RetVal_Round2
        convert_round(const SourceType & src, bool & up)
        {
            PLLL_INTERNAL_STATIC_CHECK(traits::type_traits<SourceType>::is_realtype, NeedsRealTypeSource);
            PLLL_INTERNAL_STATIC_CHECK(traits::type_traits<Native>::is_inttype, NeedsIntegerTypeDestination);
            return implementation::nativeconversion_impl2<Native, SourceType, traits::type_traits<Native>::is_inttype || !traits::type_traits<Native>::is_number,
                                                          traits::type_traits<SourceType>::is_inttype || !traits::type_traits<SourceType>::is_number>::convert_round(src, up);
        }
        
        template<class Native, class SourceType>
        typename implementation::nativeconversion_impl2<Native, SourceType, traits::type_traits<Native>::is_inttype || !traits::type_traits<Native>::is_number,
                                                        traits::type_traits<SourceType>::is_inttype || !traits::type_traits<SourceType>::is_number>::RetVal_Ceil
        /**
           \brief Converts the ceil of the value in `src` to the native type `Native` and returns
                  the result.
           
           \tparam Native     The native type to convert to.
           \tparam SourceType The source type.
           \param src The source variable.
           \return The converted result.
           
           \sa \ref arithctxts_convs_native
         */
        convert_ceil(const SourceType & src)
        {
            PLLL_INTERNAL_STATIC_CHECK(traits::type_traits<SourceType>::is_realtype, NeedsRealTypeSource);
            PLLL_INTERNAL_STATIC_CHECK(traits::type_traits<Native>::is_inttype, NeedsIntegerTypeDestination);
            return implementation::nativeconversion_impl2<Native, SourceType, traits::type_traits<Native>::is_inttype || !traits::type_traits<Native>::is_number,
                                                          traits::type_traits<SourceType>::is_inttype || !traits::type_traits<SourceType>::is_number>::convert_ceil(src);
        }
        
        ///@}
        
        class StringContext
        /**
           \brief An "arithmetic" context for conversions to `std::string`s.
           
           Note that most integer types only support bases 8, 10 and 16, and most floating point
           types only base 10. The exception is `plll::arithmetic::Integer`, which should support at
           least bases 2, 3, 4, ..., up to 36. The digits used should be '0' up to '9' followed by
           'a' up to 'z' (lower case!).
          */
        {
        private:
            unsigned d_base;
            
        public:
            /**
               \brief The type.
             */
            typedef std::string Type;
            
            /**
               \brief The properties of strings.
             */
            enum { is_cputype = false, is_realtype = false, is_inttype = false };
            
            /**
               \brief Creates a string context for base 10 (decimal).
             */
            inline StringContext()
                : d_base(10)
            {
            }
            
            /**
               \brief Creates a string context for an arbitrary base.
               
               \param base The base. Must be at least 2.
             */
            inline StringContext(unsigned base)
                : d_base(base)
            {
                assert(base >= 2);
            }
            
            /**
               \brief Queries the base of this context.
               
               \return The base.
             */
            inline unsigned base() const
            {
                return d_base;
            }
        };
        
        class HexStringContext : public StringContext
        /**
           \brief An "arithmetic" context for hexadecimal (base 16) conversions to `std::string`s.
           
           Note that floating point types (`type_traits<Type>::is_real` is `true`) might not
           support htis.
          */
        {
        public:
            inline HexStringContext()
                : StringContext(16)
            {
            }
        };
        
        class OctalStringContext : public StringContext
        /**
           \brief An "arithmetic" context for octal (base 8) conversions to `std::string`s.
           
           Note that floating point types (`type_traits<Type>::is_real` is `true`) might not
           support htis.
          */
        {
        public:
            inline OctalStringContext()
                : StringContext(8)
            {
            }
        };
        
        namespace implementation
        {
            template<class SourceType, bool SourceInt>
            /**
               \brief Provides an implementation of the
                      `plll::arithmetic::implementation::nativeconversion_impl2<>` template for
                      conversions to `signed int`.
             */
            class nativeconversion_impl2<signed int, SourceType, true, SourceInt>
            {
            public:
                typedef signed int RetVal;
                typedef typename helper::SelectFirstType<SourceInt, void, signed int>::result RetVal_Floor;
                typedef typename helper::SelectFirstType<SourceInt, void, signed int>::result RetVal_Round;
                typedef typename helper::SelectFirstType<SourceInt, void, signed int>::result RetVal_Round2;
                typedef typename helper::SelectFirstType<SourceInt, void, signed int>::result RetVal_Ceil;
                
                static RetVal convert(const SourceType & src)
                {
                    return nativeconversion_impl<SourceType>::toInt(src);
                }
            
                static RetVal_Floor convert_floor(const SourceType & src)
                {
                    return nativeconversion_impl<SourceType>::toInt_Floor(src);
                }
            
                static RetVal_Round convert_round(const SourceType & src)
                {
                    return nativeconversion_impl<SourceType>::toInt_Round(src);
                }
            
                static RetVal_Round2 convert_round(const SourceType & src, bool & up)
                {
                    return nativeconversion_impl<SourceType>::toInt_Round(src, up);
                }
            
                static RetVal_Ceil convert_ceil(const SourceType & src)
                {
                    return nativeconversion_impl<SourceType>::toInt_Ceil(src);
                }
            };
        
            template<class SourceType, bool SourceInt>
            /**
               \brief Provides an implementation of the
                      `plll::arithmetic::implementation::nativeconversion_impl2<>` template for
                      conversions to `unsigned int`.
               
               Note that `plll::arithmetic::implementation::nativeconversion_impl<SourceType>` only
               provides the `_Floor`, `_Round` and `_Ceil` functions in case `SourceInt` is `true`.
             */
            class nativeconversion_impl2<unsigned int, SourceType, true, SourceInt>
            {
            public:
                typedef unsigned int RetVal;
                typedef typename helper::SelectFirstType<SourceInt, void, unsigned int>::result RetVal_Floor;
                typedef typename helper::SelectFirstType<SourceInt, void, unsigned int>::result RetVal_Round;
                typedef typename helper::SelectFirstType<SourceInt, void, unsigned int>::result RetVal_Round2;
                typedef typename helper::SelectFirstType<SourceInt, void, unsigned int>::result RetVal_Ceil;
            
                static RetVal convert(const SourceType & src)
                {
                    return nativeconversion_impl<SourceType>::toUInt(src);
                }
            
                static RetVal_Floor convert_floor(const SourceType & src)
                {
                    return nativeconversion_impl<SourceType>::toUInt_Floor(src);
                }
            
                static RetVal_Round convert_round(const SourceType & src)
                {
                    return nativeconversion_impl<SourceType>::toUInt_Round(src);
                }
            
                static RetVal_Round2 convert_round(const SourceType & src, bool & up)
                {
                    return nativeconversion_impl<SourceType>::toUInt_Round(src, up);
                }
            
                static RetVal_Ceil convert_ceil(const SourceType & src)
                {
                    return nativeconversion_impl<SourceType>::toUInt_Ceil(src);
                }
            };
        
            template<class SourceType, bool SourceInt>
            /**
               \brief Provides an implementation of the
                      `plll::arithmetic::implementation::nativeconversion_impl2<>` template for
                      conversions to `signed long`.
               
               Note that `plll::arithmetic::implementation::nativeconversion_impl<SourceType>` only
               provides the `_Floor`, `_Round` and `_Ceil` functions in case `SourceInt` is `true`.
             */
            class nativeconversion_impl2<signed long, SourceType, true, SourceInt>
            {
            public:
                typedef signed long RetVal;
                typedef typename helper::SelectFirstType<SourceInt, void, signed long>::result RetVal_Floor;
                typedef typename helper::SelectFirstType<SourceInt, void, signed long>::result RetVal_Round;
                typedef typename helper::SelectFirstType<SourceInt, void, signed long>::result RetVal_Round2;
                typedef typename helper::SelectFirstType<SourceInt, void, signed long>::result RetVal_Ceil;
                
                static RetVal convert(const SourceType & src)
                {
                    return nativeconversion_impl<SourceType>::toLong(src);
                }
            
                static RetVal_Floor convert_floor(const SourceType & src)
                {
                    return nativeconversion_impl<SourceType>::toLong_Floor(src);
                }
            
                static RetVal_Round convert_round(const SourceType & src)
                {
                    return nativeconversion_impl<SourceType>::toLong_Round(src);
                }
            
                static RetVal_Round2 convert_round(const SourceType & src, bool & up)
                {
                    return nativeconversion_impl<SourceType>::toLong_Round(src, up);
                }
            
                static RetVal_Ceil convert_ceil(const SourceType & src)
                {
                    return nativeconversion_impl<SourceType>::toLong_Ceil(src);
                }
            };
        
            template<class SourceType, bool SourceInt>
            /**
               \brief Provides an implementation of the
                      `plll::arithmetic::implementation::nativeconversion_impl2<>` template for
                      conversions to `unsigned long`.
               
               Note that `plll::arithmetic::implementation::nativeconversion_impl<SourceType>` only
               provides the `_Floor`, `_Round` and `_Ceil` functions in case `SourceInt` is `true`.
             */
            class nativeconversion_impl2<unsigned long, SourceType, true, SourceInt>
            {
            public:
                typedef unsigned long RetVal;
                typedef typename helper::SelectFirstType<SourceInt, void, unsigned long>::result RetVal_Floor;
                typedef typename helper::SelectFirstType<SourceInt, void, unsigned long>::result RetVal_Round;
                typedef typename helper::SelectFirstType<SourceInt, void, unsigned long>::result RetVal_Round2;
                typedef typename helper::SelectFirstType<SourceInt, void, unsigned long>::result RetVal_Ceil;
            
                static unsigned long convert(const SourceType & src)
                {
                    return nativeconversion_impl<SourceType>::toULong(src);
                }
            
                static RetVal_Floor convert_floor(const SourceType & src)
                {
                    return nativeconversion_impl<SourceType>::toULong_Floor(src);
                }
            
                static RetVal_Round convert_round(const SourceType & src)
                {
                    return nativeconversion_impl<SourceType>::toULong_Round(src);
                }
            
                static RetVal_Round2 convert_round(const SourceType & src, bool & up)
                {
                    return nativeconversion_impl<SourceType>::toULong_Round(src, up);
                }
            
                static RetVal_Ceil convert_ceil(const SourceType & src)
                {
                    return nativeconversion_impl<SourceType>::toULong_Ceil(src);
                }
            };
            
            template<class SourceType, bool SourceInt>
            /**
               \brief Provides an implementation of the
                      `plll::arithmetic::implementation::nativeconversion_impl2<>` template for
                      conversions to `long long`.
               
               Note that `plll::arithmetic::implementation::nativeconversion_impl<SourceType>` only
               provides the `_Floor`, `_Round` and `_Ceil` functions in case `SourceInt` is `true`.
             */
            class nativeconversion_impl2<long long, SourceType, true, SourceInt>
            {
            public:
                typedef long long RetVal;
                typedef typename helper::SelectFirstType<SourceInt, void, long long>::result RetVal_Floor;
                typedef typename helper::SelectFirstType<SourceInt, void, long long>::result RetVal_Round;
                typedef typename helper::SelectFirstType<SourceInt, void, long long>::result RetVal_Round2;
                typedef typename helper::SelectFirstType<SourceInt, void, long long>::result RetVal_Ceil;
                
                static RetVal convert(const SourceType & src)
                {
                    return nativeconversion_impl<SourceType>::toLongLong(src);
                }
                
                static RetVal_Floor convert_floor(const SourceType & src)
                {
                    return nativeconversion_impl<SourceType>::toLongLong_Floor(src);
                }
                
                static RetVal_Round convert_round(const SourceType & src)
                {
                    return nativeconversion_impl<SourceType>::toLongLong_Round(src);
                }
                
                static RetVal_Round2 convert_round(const SourceType & src, bool & up)
                {
                    return nativeconversion_impl<SourceType>::toLongLong_Round(src, up);
                }
                
                static RetVal_Ceil convert_ceil(const SourceType & src)
                {
                    return nativeconversion_impl<SourceType>::toLongLong_Ceil(src);
                }
            };
            
            template<class SourceType, bool SourceInt>
            /**
               \brief Provides an implementation of the
                      `plll::arithmetic::implementation::nativeconversion_impl2<>` template for
                      conversions to `float`.
               
               Note that `plll::arithmetic::implementation::nativeconversion_impl<SourceType>` only
               provides the `_Floor`, `_Round` and `_Ceil` functions in case `SourceInt` is `true`.
             */
            class nativeconversion_impl2<float, SourceType, false, SourceInt>
            {
            public:
                typedef float RetVal;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static RetVal convert(const SourceType & src)
                {
                    return nativeconversion_impl<SourceType>::toFloat(src);
                }
            };
            
            template<class SourceType, bool SourceInt>
            /**
               \brief Provides an implementation of the
                      `plll::arithmetic::implementation::nativeconversion_impl2<>` template for
                      conversions to `double`.
             */
            class nativeconversion_impl2<double, SourceType, false, SourceInt>
            {
            public:
                typedef double RetVal;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static RetVal convert(const SourceType & src)
                {
                    return nativeconversion_impl<SourceType>::toDouble(src);
                }
            };
            
            template<class SourceType, bool SourceInt>
            /**
               \brief Provides an implementation of the
                      `plll::arithmetic::implementation::nativeconversion_impl2<>` template for
                      conversions to `long double`.
             */
            class nativeconversion_impl2<long double, SourceType, false, SourceInt>
            {
            public:
                typedef long double RetVal;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static RetVal convert(const SourceType & src)
                {
                    return nativeconversion_impl<SourceType>::toLongDouble(src);
                }
            };
            
            template<class Context> class conversion_impl<typename Context::Type, Context>
            /**
               \brief Provides an identity conversion from a type to itself. Note that the context
                      is ignored.
               
               Note that types where the context matters are supposed to specialize this!
             */
            {
            public:
                typedef const typename Context::Type & RetVal;
                typedef typename helper::SelectFirstType<Context::is_realtype, typename Context::Type, void>::result RetVal_Frac;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static void convert(typename Context::Type & d, const typename Context::Type & v, const Context & c)
                {
                    d = v;
                }
                
                static RetVal convert(const typename Context::Type & v, const Context & c)
                {
                    return v;
                }
                
                static void convert_frac(typename Context::Type & d, const typename Context::Type & v1, const typename Context::Type & v2, const Context & c)
                {
                    d = v1 / v2;
                }
                
                static RetVal_Frac convert_frac(const typename Context::Type & v1, const typename Context::Type & v2, const Context & c)
                {
                    return v1 / v2;
                }
            };
            
            template<class Context>
            /**
               \brief Provides facilities to convert strings to types.
               
               It should define two member types `RetVal1` and `RetVal2`, and four static member
               functions:
               
               - `bool convert(Context::Type &, const std::string &, const Context &)`: converts the
                 given `std::string` string to the type using the given context's information. The
                 result is stored in the first argument. The return value is `true` if the
                 conversion succeeded, and `false` if it failed.
               
               - `bool convert(Context::Type &, const char *, const Context &)`: converts the given
                 C-style string to the type using the given context's information. The result is
                 stored in the first argument. The return value is `true` if the conversion
                 succeeded, and `false` if it failed.
               
               - `RetVal1 convert(const std::string &, const Context &)`: converts the given
                 `std::string` string to the type using the given context's information. The result
                 is returned. The behavior is undefined in case the string cannot be converted.
               
               - `RetVal2 convert(const char *, const Context &)`: converts the given C-style string
                 to the type using the given context's information. The result is returned. The
                 behavior is undefined in case the string cannot be converted.
             */
            class from_string_conversion;
            
            template<class Type>
            /**
               \brief Provides facilities to convert types to `std::string`s.
               
               It should define a member type `RetVal`, and four static member functions:
               
               - `static RetVal convert(const Type &)`: converts the given value to a string
                 representation. For numerical types, decimal representation should be used.
               
               - `static RetVal convert(const Type &, unsigned basis)`: converts the given value to
                 a string representation. For numerical types, the given basis should be used.
               
               - `static void convert(std::string &, const Type &)`: converts the given value to a
                 string representation and stores the result the given string. For numerical types,
                 decimal representation should be used.
               
               - `static void convert(std::string &, const Type &, unsigned basis)`: converts the
                 given value to a string representation and stores the result the given string. For
                 numerical types, the given basis should be used.
             */
            class to_string_conversion;
            
            template<class Context>
            /**
               \brief Default implementation for type to string conversion. Uses
                      `std::istringstream`.
             */
            class from_string_conversion
            {
            public:
                typedef typename Context::Type RetVal1;
                typedef typename Context::Type RetVal2;
            
                static bool convert(typename Context::Type & res, const std::string & s, const Context &)
                {
                    std::istringstream ss(s);
                    ss >> res;
                    return ss;
                }
            
                static bool convert(typename Context::Type & res, const char * s, const Context & c)
                {
                    std::istringstream ss(s);
                    ss >> res;
                    return ss;
                }
            
                static RetVal1 convert(const std::string & s, const Context & c)
                {
                    RetVal1 res(c);
                    convert(res, s, c);
                    return res;
                }
            
                static RetVal2 convert(const char * s, const Context & c)
                {
                    RetVal2 res(c);
                    convert(res, s, c);
                    return res;
                }
            };
            
            template<class Type>
            class to_string_conversion
            /**
               \brief Default implementation for string to type conversion. Uses
                      `std::ostringstream`. Only supports bases 8, 10 and 16.
             */
            {
            public:
                typedef std::string RetVal;
                
                static RetVal convert(const Type & t)
                {
                    std::ostringstream ss;
                    ss << t;
                    return ss.str();
                }
                
                static RetVal convert(const Type & t, unsigned base)
                {
                    std::ostringstream ss;
                    ss << std::setbase(base) << t;
                    return ss.str();
                }
                
                static void convert(std::string & str, const Type & t)
                {
                    std::ostringstream ss;
                    ss << t;
                    str = ss.str();
                }
                
                static void convert(std::string & str, const Type & t, unsigned base)
                {
                    std::ostringstream ss;
                    ss << std::setbase(base) << t;
                    str = ss.str();
                }
            };
            
            template<class Type>
            class conversion_impl<Type, StringContext>
            /**
               \brief Provides type to string conversions using the usual type conversion facility.
               
               Internally uses `plll::arithmetic::implementation::to_string_conversion<Type>`.
             */
            {
            public:
                typedef typename to_string_conversion<Type>::RetVal RetVal;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static inline void convert(std::string & str, const Type & v, const StringContext & c)
                {
                    to_string_conversion<Type>::convert(str, v, c.base());
                }
                
                static inline RetVal convert(const Type & v, const StringContext & c)
                {
                    return to_string_conversion<Type>::convert(v, c.base());
                }
            };
            
            template<class Type>
            /**
               \brief Provides type to string conversions using the usual type conversion facility
                      (for base 16).
               
               Internally uses `plll::arithmetic::implementation::to_string_conversion<Type>`.
             */
            class conversion_impl<Type, HexStringContext>
            {
            public:
                typedef typename to_string_conversion<Type>::RetVal RetVal;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static inline void convert(std::string & str, const Type & v, const HexStringContext & c)
                {
                    to_string_conversion<Type>::convert(str, v, c.base());
                }
                
                static inline RetVal convert(const Type & v, const HexStringContext & c)
                {
                    return to_string_conversion<Type>::convert(v, c.base());
                }
            };
            
            template<class Type>
            /**
               \brief Provides type to string conversions using the usual type conversion facility
                      (for base 8).
               
               Internally uses `plll::arithmetic::implementation::to_string_conversion<Type>`.
             */
            class conversion_impl<Type, OctalStringContext>
            {
            public:
                typedef typename to_string_conversion<Type>::RetVal RetVal;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static inline void convert(std::string & str, const Type & v, const OctalStringContext & c)
                {
                    to_string_conversion<Type>::convert(str, v, c.base());
                }
                
                static inline RetVal convert(const Type & v, const OctalStringContext & c)
                {
                    return to_string_conversion<Type>::convert(v, c.base());
                }
            };
            
            template<class Context>
            /**
               \brief Provides `std::string` to type conversions using the usual type conversion
                      facility.
               
               Internally uses `plll::arithmetic::implementation::from_string_conversion<Context>`.
             */
            class conversion_impl<std::string, Context>
            {
            public:
                typedef typename from_string_conversion<Context>::RetVal1 RetVal;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
            
                static void convert(typename Context::Type & d, const std::string & s, const Context & c)
                {
                    from_string_conversion<Context>::convert(d, s, c);
                }
            
                static RetVal convert(const std::string & s, const Context & c)
                {
                    return from_string_conversion<Context>::convert(s, c);
                }
            };
            
            template<class Context>
            /**
               \brief Provides C-style string to type conversions using the usual type conversion
                      facility.
               
               Internally uses `plll::arithmetic::implementation::from_string_conversion<Context>`.
             */
            class conversion_impl<const char *, Context>
            {
            public:
                typedef typename from_string_conversion<Context>::RetVal2 RetVal;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
            
                static inline void convert(typename Context::Type & d, const char * s, const Context & c)
                {
                    from_string_conversion<Context>::convert(d, s, c);
                }
            
                static inline RetVal convert(const char * s, const Context & c)
                {
                    return from_string_conversion<Context>::convert(s, c);
                }
            };
            
            template<class Context>
            /**
               \brief Provides C-style string to type conversions using the usual type conversion
                      facility.
               
               Internally uses `plll::arithmetic::implementation::from_string_conversion<Context>`.
             */
            class conversion_impl<char *, Context>
            {
            public:
                typedef typename from_string_conversion<Context>::RetVal2 RetVal;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
            
                static inline void convert(typename Context::Type & d, char * s, const Context & c)
                {
                    from_string_conversion<Context>::convert(d, static_cast<const char *>(s), c);
                }
            
                static inline RetVal convert(char * s, const Context & c)
                {
                    return from_string_conversion<Context>::convert(static_cast<const char *>(s), c);
                }
            };
            
            template<unsigned n, class Context>
            /**
               \brief Provides C-style string to type conversions using the usual type conversion
                      facility.
               
               Internally uses `plll::arithmetic::implementation::from_string_conversion<Context>`.
             */
            class conversion_impl<const char[n], Context>
            {
            public:
                typedef typename from_string_conversion<Context>::RetVal2 RetVal;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static inline void convert(typename Context::Type & d, const char * s, const Context & c)
                {
                    from_string_conversion<Context>::convert(d, s, c);
                }
                
                static inline RetVal convert(const char * s, const Context & c)
                {
                    return from_string_conversion<Context>::convert(s, c);
                }
            };
            
            template<unsigned n, class Context>
            /**
               \brief Provides C-style string to type conversions using the usual type conversion
                      facility.
               
               Internally uses `plll::arithmetic::implementation::from_string_conversion<Context>`.
             */
            class conversion_impl<char[n], Context>
            {
            public:
                typedef typename from_string_conversion<Context>::RetVal2 RetVal;
                typedef void RetVal_Floor;
                typedef void RetVal_Ceil;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                
                static inline void convert(typename Context::Type & d, const char * s, const Context & c)
                {
                    from_string_conversion<Context>::convert(d, static_cast<const char *>(s), c);
                }
                
                static inline RetVal convert(const char * s, const Context & c)
                {
                    return from_string_conversion<Context>::convert(static_cast<const char *>(s), c);
                }
            };
            
            template<class SourceType, bool SourceInt>
            /**
               \brief Provides type to `std::string` conversion using the native type conversion
                      facility.
               
               Internally uses `plll::arithmetic::implementation::to_string_conversion<Context>`.
             */
            class nativeconversion_impl2<std::string, SourceType, true, SourceInt>
            {
            public:
                typedef typename to_string_conversion<SourceType>::RetVal RetVal;
                typedef void RetVal_Floor;
                typedef void RetVal_Round;
                typedef void RetVal_Round2;
                typedef void RetVal_Ceil;
                
                static inline RetVal convert(const SourceType & v)
                {
                    return to_string_conversion<SourceType>::convert(v);
                }
            };
            
            template<> class nativeconversion_impl<std::string>
            /**
               \brief Provides `std::string` to native type conversion using the native type
                      conversion facility.
               
               Internally uses `std::istringstream`.
             */
            {
            public:
                static int toInt(const std::string & v)
                {
                    std::istringstream s(v);
                    int res;
                    s >> res;
                    return res;
                }
                
                static unsigned int toUInt(const std::string & v)
                {
                    std::istringstream s(v);
                    unsigned int res;
                    s >> res;
                    return res;
                }
                
                static long toLong(const std::string & v)
                {
                    std::istringstream s(v);
                    long res;
                    s >> res;
                    return res;
                }
                
                static long toULong(const std::string & v)
                {
                    std::istringstream s(v);
                    unsigned long res;
                    s >> res;
                    return res;
                }
                
                static long long toLongLong(const std::string & v)
                {
                    std::istringstream s(v);
                    long long res;
                    s >> res;
                    return res;
                }
                
                static float toFloat(const std::string & v)
                {
                    std::istringstream s(v);
                    float res;
                    s >> res;
                    return res;
                }
                
                static double toDouble(const std::string & v)
                {
                    std::istringstream s(v);
                    double res;
                    s >> res;
                    return res;
                }
                
                static long double toLongDouble(const std::string & v)
                {
                    std::istringstream s(v);
                    long double res;
                    s >> res;
                    return res;
                }
            };
            
            template<> class nativeconversion_impl<const char *>
            /**
               \brief Provides C-style string to native type conversion using the native type
                      conversion facility.
               
               Internally uses `std::istringstream`.
             */
            {
            public:
                static int toInt(const char * v)
                {
                    std::istringstream s(v);
                    int res;
                    s >> res;
                    return res;
                }
                
                static unsigned int toUInt(const char * v)
                {
                    std::istringstream s(v);
                    unsigned int res;
                    s >> res;
                    return res;
                }
                
                static long toLong(const char * v)
                {
                    std::istringstream s(v);
                    long res;
                    s >> res;
                    return res;
                }
                
                static unsigned long toULong(const char * v)
                {
                    std::istringstream s(v);
                    unsigned long res;
                    s >> res;
                    return res;
                }
                
                static long long toLongLong(const char * v)
                {
                    std::istringstream s(v);
                    long long res;
                    s >> res;
                    return res;
                }
                
                static float toFloat(const char * v)
                {
                    std::istringstream s(v);
                    float res;
                    s >> res;
                    return res;
                }
                
                static double toDouble(const char * v)
                {
                    std::istringstream s(v);
                    double res;
                    s >> res;
                    return res;
                }
                
                static long double toLongDouble(const char * v)
                {
                    std::istringstream s(v);
                    long double res;
                    s >> res;
                    return res;
                }
            };
            
            template<> class nativeconversion_impl<char *>
            /**
               \brief Provides C-style string to native type conversion using the native type
                      conversion facility.
               
               Internally uses `std::istringstream`.
             */
            {
            public:
                static inline int toInt(char * v)
                {
                    return nativeconversion_impl<const char *>::toInt(v);
                }
                
                static inline unsigned int toUInt(char * v)
                {
                    return nativeconversion_impl<const char *>::toInt(v);
                }
                
                static inline long toLong(char * v)
                {
                    return nativeconversion_impl<const char *>::toLong(v);
                }
                
                static inline unsigned long toULong(char * v)
                {
                    return nativeconversion_impl<const char *>::toLong(v);
                }
                
                static inline long long toLongLong(char * v)
                {
                    return nativeconversion_impl<const char *>::toLongLong(v);
                }
                
                static inline float toFloat(char * v)
                {
                    return nativeconversion_impl<const char *>::toFloat(v);
                }
                
                static inline double toDouble(char * v)
                {
                    return nativeconversion_impl<const char *>::toDouble(v);
                }
                
                static inline long double toLongDouble(char * v)
                {
                    return nativeconversion_impl<const char *>::toLongDouble(v);
                }
            };
            
            template<unsigned n> class nativeconversion_impl<char [n]>
            /**
               \brief Provides C-style string to native type conversion using the native type
                      conversion facility.
               
               Internally uses `std::istringstream`.
             */
            {
            public:
                static inline int toInt(const char * v)
                {
                    return nativeconversion_impl<const char *>::toInt(v);
                }
                
                static inline unsigned int toUInt(const char * v)
                {
                    return nativeconversion_impl<const char *>::toInt(v);
                }
                
                static inline long toLong(const char * v)
                {
                    return nativeconversion_impl<const char *>::toLong(v);
                }
                
                static inline unsigned long toULong(const char * v)
                {
                    return nativeconversion_impl<const char *>::toLong(v);
                }
                
                static inline long long toLongLong(const char * v)
                {
                    return nativeconversion_impl<const char *>::toLongLong(v);
                }
                
                static inline float toFloat(const char * v)
                {
                    return nativeconversion_impl<const char *>::toFloat(v);
                }
                
                static inline double toDouble(const char * v)
                {
                    return nativeconversion_impl<const char *>::toDouble(v);
                }
                
                static inline long double toLongDouble(const char * v)
                {
                    return nativeconversion_impl<const char *>::toLongDouble(v);
                }
            };
            
            template<unsigned n> class nativeconversion_impl<const char [n]>
            /**
               \brief Provides C-style string to native type conversion using the native type
                      conversion facility.
               
               Internally uses `std::istringstream`.
             */
            {
            public:
                static inline int toInt(const char * v)
                {
                    return nativeconversion_impl<const char *>::toInt(v);
                }
                
                static inline unsigned int toUInt(const char * v)
                {
                    return nativeconversion_impl<const char *>::toInt(v);
                }
                
                static inline long toLong(const char * v)
                {
                    return nativeconversion_impl<const char *>::toLong(v);
                }
                
                static inline unsigned long toULong(const char * v)
                {
                    return nativeconversion_impl<const char *>::toLong(v);
                }
                
                static inline long long toLongLong(const char * v)
                {
                    return nativeconversion_impl<const char *>::toLongLong(v);
                }
                
                static inline float toFloat(const char * v)
                {
                    return nativeconversion_impl<const char *>::toFloat(v);
                }
                
                static inline double toDouble(const char * v)
                {
                    return nativeconversion_impl<const char *>::toDouble(v);
                }
                
                static inline long double toLongDouble(const char * v)
                {
                    return nativeconversion_impl<const char *>::toLongDouble(v);
                }
            };
        }
    }
}

#include <plll/arithmetic-expressions.hpp>

namespace plll
{
    namespace arithmetic
    {
        namespace implementation
        {
            extern IntegerContext g_intcontext;
            
            template<class Source, bool SourceInt>
            class nativeconversion_impl2<Integer, Source, true, SourceInt>
            /**
               \brief This class allows a shortcut: converting types to Integer objects by writing
                      convert<Integer>(source).
            */
            {
            public:
                typedef typename conversion_impl<Source, IntegerContext>::RetVal RetVal;
                typedef typename helper::SelectFirstType<SourceInt, void, typename conversion_impl<Source, IntegerContext>::RetVal_Floor>::result RetVal_Floor;
                typedef typename helper::SelectFirstType<SourceInt, void, typename conversion_impl<Source, IntegerContext>::RetVal_Round>::result RetVal_Round;
                typedef typename helper::SelectFirstType<SourceInt, void, typename conversion_impl<Source, IntegerContext>::RetVal_Round2>::result RetVal_Round2;
                typedef typename helper::SelectFirstType<SourceInt, void, typename conversion_impl<Source, IntegerContext>::RetVal_Ceil>::result RetVal_Ceil;
                
                static RetVal convert(const Source & v)
                {
                    return conversion_impl<Source, IntegerContext>::convert(v, g_intcontext);
                }
                
                static RetVal_Floor convert_floor(const Source & v)
                {
                    return conversion_impl<Source, IntegerContext>::floor(v, g_intcontext);
                }
                
                static RetVal_Round convert_round(const Source & v)
                {
                    return conversion_impl<Source, IntegerContext>::round(v, g_intcontext);
                }
                
                static RetVal_Round2 convert_round(const Source & v, bool & up)
                {
                    return conversion_impl<Source, IntegerContext>::round(v, up, g_intcontext);
                }
                
                static RetVal_Ceil convert_ceil(const Source & v)
                {
                    return conversion_impl<Source, IntegerContext>::ceil(v, g_intcontext);
                }
            };
        }
    }
}

#include <plll/arithmetic-gmp.hpp>

#endif

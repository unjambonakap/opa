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

#ifndef PLLL_INCLUDE_GUARD__HELPER_HPP
#define PLLL_INCLUDE_GUARD__HELPER_HPP

#include <plll/config.hpp>

/**
   \file
   \brief Helper templates.
   
   This header provides helper templates, such as static assertions (`PLLL_INTERNAL_STATIC_CHECK`
   macro), type selection (`plll::SelectFirstType<>`), int-to-type and bool-to-type facilities
   (`IntToType<>` and `BoolToType<>`) as well as a decoration remover.
*/
namespace plll
{
    /**
       \brief Contains some helper templates and classes.
       
       The `plll::helper` namespace contains helper templates used throughout the `plll`
       library. Examples are a static assertion macro (`PLLL_INTERNAL_STATIC_CHECK`), a type
       selector (`plll::helper::SelectFirstType<>`), int/bool to type converters
       (`plll::helper::IntToType<>` and `plll::helper::BoolToType<>`) and a decoration remover
       (`plll::helper::remove_decorations<>`).
       
       Another important member is the `plll::helper::ArgumentParser` class which allows to parse
       command line arguments.
    */
    namespace helper
    {
        namespace implementation
        {
            template<bool no_error> struct CompileTimeError;
            template<> struct CompileTimeError<true> {};
        }
        
        /**
           Allows to do compile-time assertions. By writing
           
               PLLL_INTERNAL_STATIC_CHECK(condition, IdentifierWhichIsAMessage)
               
           the compiler will evaluate the resulting expression to nothing if `condition` evaluates
           to `true`, and produce a compilation error that it cannot instantiate
           `ERROR_IdentifierWhichIsAMessage` (or something like that) in case `condition` evaluates
           to `false`.
           
           \param condition                  An expression which evaluates to an integer or `bool`
                                             at compile time.
           \param IdentifierWhichIsAMessage  A C++ identifier which will be prefixed with `ERROR_`
                                             and hopefully included in the compiler's error message.
           
           The code was essentially taken from the <a
           href="http://loki-lib.cvs.sourceforge.net/loki-lib/loki/include/loki/static_check.h?view=markup">Loki
           library</a> by Andrei Alexandrescu. Note that in C++11, <a
           href="http://en.cppreference.com/w/cpp/language/static_assert">`static_assert()<>`</a>
           provides the same functionality.
        */
#if __cplusplus >= 201103L
        #define PLLL_INTERNAL_STATIC_CHECK(condition, IdentifierWhichIsAMessage) \
            static_assert(condition, #IdentifierWhichIsAMessage);
#else
        #define PLLL_INTERNAL_STATIC_CHECK(condition, IdentifierWhichIsAMessage) \
            { plll::helper::implementation::CompileTimeError<((condition) != 0)> ERROR_##IdentifierWhichIsAMessage; (void)ERROR_##IdentifierWhichIsAMessage; }
#endif
        
        namespace implementation
        {
            template<bool cond, typename A, typename B>
            class select_first_type_impl
            {
            public:
                typedef A result;
            };
            
            template<typename A, typename B>
            class select_first_type_impl<false, A, B>
            {
            public:
                typedef B result;
            };
        }
        
        template<bool cond, typename A, typename B>
        /**
           \brief A type selector template.
           
           Defines `SelectFirstType<cond, A, B>::result` as `A` in case `cond` is `true`, and `B` in
           case `cond` is false. This is used in template metaprogramming to realize certain compile
           time `if` clauses.
           
           \tparam cond An expression which evaluates to a `bool` at compile time.
           \tparam A    A type name; `result` will have this type in case `cond` evalutes to `true`.
           \tparam B    A type name; `result` will have this type in case `cond` evalutes to `false`.
           
           This was taken from "Modern C++ Design" by Andrei Alexandrescu. Note that in C++11, <a
           href="http://en.cppreference.com/w/cpp/types/conditional">`std::conditional<>`</a>
           provides the same functionality.
        */
        class SelectFirstType
        {
        public:
            /** \brief The result type. Either `A` or `B`, depending on the value of `cond`. */
            typedef typename implementation::select_first_type_impl<cond, A, B>::result result;
        };
        
        namespace implementation
        {
            template<typename A, typename B>
            struct type_equality
            {
                enum { result = false };
            };
            
            template<typename A>
            struct type_equality<A, A>
            {
                enum { result = true };
            };
        }
        
        template<typename A, typename B>
        /**
           \brief A type comparison template.
           
           Defines `CompareTypes<A, B>::equal` as `true` in case `A` and `B` are the same type, and
           `false` if they are different types.
           
           \tparam A    A type name.
           \tparam B    A type name.
        */
        class CompareTypes
        {
        public:
            /** \brief The result. Either `true` if types `A` and `B` are the same, or `false` otherwise. */
            enum { equal = implementation::type_equality<A, B>::result, not_equal = !implementation::type_equality<A, B>::result };
        };
        
        namespace implementation
        {
            template<bool cond, typename T>
            class enable_if_impl
            {
            public:
            };
            
            template<typename T>
            class enable_if_impl<true, T>
            {
            public:
            };
        }
        
        template<bool cond, typename Type = void>
        /**
           \brief Allows to remove functions from overload resolution.
           
           This can be used as <a
           href="http://en.cppreference.com/w/cpp/types/enable_if.html">`std::enable_if<>`</a> in
           C++11.
         */
        class EnableIf
        {
        public:
            typedef typename implementation::enable_if_impl<cond, Type>::result result;
        };
        
        template<int val>
        /**
           \brief Conversion from compile-time known ints to types.
           
           Provides a unique type for every `int` by writing `IntToType<specific_integer>`. This is
           used in template metaprogramming to realize compile time selections.
           
           \tparam val The "value" of this type.
           
           This was taken from "Modern C++ Design" by Andrei Alexandrescu. Note that in C++11, <a
           href="http://en.cppreference.com/w/cpp/types/integral_constant">`std::integral_constant<>`</a>
           provides the same functionality.
        */
        struct IntToType
        {
            enum { value = val }; /**< \brief The type's value. */
        };
        
        template<bool val>
        /**
           \brief Conversion from compile-time known bools to types.
           
           Provides a unique type for `true` and `false` writing `IntToType<true>` respectively
           `IntToType<false>`. This is used in template metaprogramming to realize compile time `if`
           clauses.
           
           \tparam val The "value" of this type.
           
           This was taken from "Modern C++ Design" by Andrei Alexandrescu. Note that in C++11, <a
           href="http://en.cppreference.com/w/cpp/types/integral_constant">`std::integral_constant<>`</a>
           provides the same functionality.
        */
        struct BoolToType
        {
            enum { value = val }; /**< \brief The type's value. */
        };
        
        namespace implementation
        {
            template<typename A>
            struct remove_decorations_impl
            {
                typedef A result; /**< \brief The stripped type. */
            };
            
            template<typename A>
            struct remove_decorations_impl<volatile A>
            {
                typedef typename remove_decorations_impl<A>::result result;
            };
            
            template<typename A>
            struct remove_decorations_impl<const A>
            {
                typedef typename remove_decorations_impl<A>::result result;
            };
            
            template<typename A>
            struct remove_decorations_impl<A &>
            {
                typedef typename remove_decorations_impl<A>::result result;
            };
            
#if __cplusplus >= 201103L
            template<typename A>
            struct remove_decorations_impl<A &&>
            {
                typedef typename remove_decorations_impl<A>::result result;
            };
#endif
        }
        
        template<typename A>
        /**
           \brief Strips decorations `const`, `&`, `volatile` and `&&` from types.
           
           Strips the decorations `const`, `&`, `volatile` and `&&` (in case of C++11) from
           types. More precisely, the four types
           
           - `plll::remove_decorations<T>::Result`
           - `plll::remove_decorations<T &>::Result`
           - `plll::remove_decorations<const T>::Result`
           - `plll::remove_decorations<const T &>::Result`
           - `plll::remove_decorations<volatile T>::Result`
           - `plll::remove_decorations<volatile T &>::Result`
           - `plll::remove_decorations<volatile const T>::Result`
           - `plll::remove_decorations<volatile const T &>::Result`
           
           will be all equal to T in case T is a class, native type or pointer. (In C++11, one could
           add eight more examples with `&&`.)
           
           \tparam A The type to strip decorations from.
           
           Note that in C++11, <a
           href="http://en.cppreference.com/w/cpp/types/remove_reference">`std::remove_reference<>`</a>
           and <a href="http://en.cppreference.com/w/cpp/types/remove_cv">`std::remove_cv<>`</a>
           combined would provide the same functionality.
        */
        struct remove_decorations
        {
            /** \brief The stripped type. */
            typedef typename implementation::remove_decorations_impl<A>::result Result;
        };
        
        /**
           \brief This is a pseudo-template which should only be used in expressions which are never
                  evaluated, such as `sizeof()` and `noexcept()` (in C++11).
           
           There exists no implementation (and it would be hard to actually create one) without
           wasting ressources. Note that this function is guaranteed to not throw an exception by
           using `throw()` in C++98/03 and `noexcept` in C++11 and later.
         */
        template<typename T>
        T & make_type_lvalue() PLLL_INTERNAL_NOTHROW_POSTFIX_ENFORCE;
    }
}

#endif

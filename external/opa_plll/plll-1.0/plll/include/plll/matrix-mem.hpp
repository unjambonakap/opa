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

#ifndef PLLL_INCLUDE_GUARD__MATRIX_MEM_HPP
#define PLLL_INCLUDE_GUARD__MATRIX_MEM_HPP

#include <cassert>
#include <cstddef>
#include <cstdlib>
//#include <iostream>

#if __cplusplus >= 201103L
    #include <utility>
#endif

/**
   \file
   \brief Memory management for matrix and vector operations.
   
   This header provides memory management traits (`plll::linalg::storage_traits<>` template) for the
   matrix and vector templates in `matrix.hpp`.
   
   Note that by modifying this header, one can switch between different memory managers for the
   default `plll::linalg::storage_traits<>` traits. One can choose between

   - a debug variant, which pads blocks of memory and periodically checks them if out-of-bound write
     accesses happened or if something was freed twice;

   - a failsafe variant, which uses no tricks and classical new and delete for memory
     management. Will call superfluous constructors and use operator= instead of copy-constructing;

   - the "release" mode variant, which is the default and should always be used except when
     debugging.
*/
namespace plll
{
    namespace linalg
    {
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Memory management
        
        //#define __VECMAT_MM_DEBUG
        //#define __VECMAT_MM_FAILSAFE
        #define __VECMAT_MM_RELEASE
        
        #ifdef __VECMAT_MM_DEBUG
        
        #define __VECMAT_MM_DEBUG__BEGIN_OK  0xAFFE1234
        #define __VECMAT_MM_DEBUG__END_OK    0x1234AFFE
        #define __VECMAT_MM_DEBUG__BEGIN_REL 0xF8E8F8E8
        #define __VECMAT_MM_DEBUG__END_REL   0x1E1F1E1F
        #define __VECMAT_MM_DEBUG__DEAD      0xDEADBEEF
          
        template<class T>
        class storage_traits
        {
        private:
            enum { DEBUG_PAD_BEGIN = 4,
                   DEBUG_PAD_END = 4 };
      
            static T * create(unsigned int n, bool dontconstruct)
            {
                unsigned SIZE = ((n * sizeof(T) + 3) & (unsigned)(~3));
                assert(SIZE >= n * sizeof(T));
                char * realptr = reinterpret_cast<char *>(malloc(SIZE + (DEBUG_PAD_BEGIN + DEBUG_PAD_END) * 4));
                char * realptrend = realptr + SIZE + DEBUG_PAD_BEGIN * 4;
                for (int i = 0; i < DEBUG_PAD_BEGIN; ++i)
                    ((unsigned int*)realptr)[i] = __VECMAT_MM_DEBUG__BEGIN_OK;
                for (int i = 0; i < DEBUG_PAD_END; ++i)
                    ((unsigned int*)realptrend)[i] = __VECMAT_MM_DEBUG__END_OK;
                T * ptr = reinterpret_cast<T*>(realptr + DEBUG_PAD_BEGIN * 4);
//          std::cout << "ALC: " << (void*)realptr << " " << (void*)ptr << " " << sizeof(T) << " x " << n << " " << (void*)realptrend << "\n";
                if (!dontconstruct)
                    for (unsigned i = 0; i < n; ++i)
                        ::new(static_cast<void *>(ptr + i)) T();
                for (int i = 0; i < DEBUG_PAD_BEGIN; ++i)
                    assert(((unsigned int*)realptr)[i] == __VECMAT_MM_DEBUG__BEGIN_OK);
                for (int i = 0; i < DEBUG_PAD_END; ++i)
                    assert(((unsigned int*)realptrend)[i] == __VECMAT_MM_DEBUG__END_OK);
                return ptr;
            }
      
            static void release(T * ptr, unsigned int n)
            {
                if (ptr == 0)
                {
                    assert(n == 0);
                    return;
                }
                unsigned SIZE = ((n * sizeof(T) + 3) & (unsigned)(~3));
                char * realptr = reinterpret_cast<char *>(ptr) - DEBUG_PAD_BEGIN * 4;
                char * realptrend = realptr + SIZE + DEBUG_PAD_BEGIN * 4;
//          std::cout << "REL: " << (void*)realptr << " " << (void*)ptr << " " << sizeof(T) << " x " << n << " " << (void*)realptrend << "\n";
                for (int i = 0; i < DEBUG_PAD_BEGIN; ++i)
                    assert(((unsigned int*)realptr)[i] == __VECMAT_MM_DEBUG__BEGIN_OK);
                for (int i = 0; i < DEBUG_PAD_END; ++i)
                    assert(((unsigned int*)realptrend)[i] == __VECMAT_MM_DEBUG__END_OK);
                // Deconstruct
                for (unsigned i = 0; i < n; ++i)
                    (ptr + i)->~T();
                // Fill with __VECMAT_MM_DEBUG__DEAD
                for (unsigned i = 0; i < n * sizeof(T) / 4; ++i)
                    (reinterpret_cast<unsigned int*>(ptr))[i] = __VECMAT_MM_DEBUG__DEAD;
                // Check padding and stop
                for (int i = 0; i < DEBUG_PAD_BEGIN; ++i)
                {
                    assert((reinterpret_cast<unsigned int*>(realptr))[i] == __VECMAT_MM_DEBUG__BEGIN_OK);
                    (reinterpret_cast<unsigned int*>(realptr))[i] = __VECMAT_MM_DEBUG__BEGIN_REL;
                }
                for (int i = 0; i < DEBUG_PAD_END; ++i)
                {
                    assert((reinterpret_cast<unsigned int*>(realptrend))[i] == __VECMAT_MM_DEBUG__END_OK);
                    (reinterpret_cast<unsigned int*>(realptrend))[i] = __VECMAT_MM_DEBUG__END_REL;
                }
                ::free(realptr);
            }
      
        public:
            static void check(T * ptr, unsigned int n)
            {
                if (ptr == 0)
                {
                    assert(n == 0);
                    return;
                }
                unsigned SIZE = ((n * sizeof(T) + 3) & (unsigned)(~3));
                char * realptr = reinterpret_cast<char *>(ptr) - DEBUG_PAD_BEGIN * 4;
                char * realptrend = realptr + SIZE + DEBUG_PAD_BEGIN * 4;
//          std::cout << "CHK: " << (void*)realptr << " " << (void*)ptr << " " << sizeof(T) << "x" << n << "=" << SIZE << " " << (void*)realptrend << "\n";
                for (int i = 0; i < DEBUG_PAD_BEGIN; ++i)
                    assert((reinterpret_cast<unsigned int*>(realptr))[i] == __VECMAT_MM_DEBUG__BEGIN_OK);
                for (int i = 0; i < DEBUG_PAD_END; ++i)
                    assert((reinterpret_cast<unsigned int*>(realptrend))[i] == __VECMAT_MM_DEBUG__END_OK);
            }
      
            typedef T * pointer_type;
            typedef const T * constpointer_type;
            typedef T & ref_type;
            typedef const T & constref_type;
      
            static void construct(ref_type ref)
            // Constructs an object via default constructor.
            {
                ::new(static_cast<void *>(&ref)) T();
            }

            template<class S>
            static void copy_construct(ref_type ref, const S & obj)
            // Constructs an object via copy constructor.
            {
                ::new(static_cast<void *>(&ref)) T(obj);
            }
      
#if __cplusplus >= 201103L
            template<class S>
            static void move_construct(ref_type ref, S && obj)
            // Constructs an object via move constructor.
            {
                ::new(static_cast<void *>(&ref)) T(std::move(obj));
            }
#endif
            
            static pointer_type alloc(unsigned int n)
            // Creates an array of n elements of type T. The entries are created using the default
            // constructor.
            {
                pointer_type newptr = create(n, false);
//          std::cout << "[#:" << newptr << "] Creating array of " << n << " elements of type " << typeid(T).name() << "\n";
                return newptr;
            }
      
            static pointer_type alloc_dontconstruct(unsigned int n)
            // Creates an array of n elements of type T. The entries are not constructed.
            {
                pointer_type newptr = create(n, true);
//          std::cout << "[#:" << newptr << "] Creating array of " << n << " elements of type " << typeid(T).name() << "\n";
                return newptr;
            }
      
            template<class S>
            static pointer_type alloc(const S & obj, unsigned int n)
            // Creates an array of n elements of type T, all being copies of obj (of type S).
            {
                pointer_type newptr = create(n, true);
//          std::cout << "[#:" << newptr << "] Creating array of " << n << " elements of type " << typeid(T).name() << ", constructed from object [" << &obj << "].\n";
                for (unsigned i = 0; i < n; ++i)
                    copy_construct(newptr[i], obj);
                return newptr;
            }
    
            template<class S>
            static pointer_type clone(S * ptr, unsigned int n)
            // Creates an array of n elements of type T. The entries are copy-constructed from the given
            // array ptr (of type S).
            {
                pointer_type newptr = create(n, true);
//          std::cout << "[#:" << newptr << "] Cloning array of " << n << " elements of type " << typeid(T).name() << " [" << ptr << "].\n";
                for (unsigned i = 0; i < n; ++i)
                    copy_construct(newptr[i], ptr[i]);
                return newptr;
            }
      
            template<class S>
            static pointer_type clone(S * ptr, unsigned int n, unsigned int ncopy)
            // Creates an array of n elements of type T. The first ncopy elements are copy-constructed from
            // the given array ptr, the last ones are constructed with the default constructor.
            {
                pointer_type newptr = create(n, true);
//          std::cout << "[#:" << newptr << "] Cloning array of " << n << " elements of type " << typeid(T).name() << " (only the " << ncopy << " first) [" << ptr << "].\n";
                for (unsigned i = 0; i < ncopy; ++i)
                    copy_construct(newptr[i], ptr[i]);
                for (unsigned i = ncopy; i < n; ++i)
                    construct(newptr[i]);
                return newptr;
            }
      
            template<class S, class SS>
            static pointer_type clone(S * ptr, unsigned int n, unsigned int ncopy, const SS & copyobj)
            // Creates an array of n elements of type T. The first ncopy elements are copy-constructed from
            // the given array ptr, the last ones are constructed with the default constructor.
            {
                pointer_type newptr = create(n, true);
//          std::cout << "[#:" << newptr << "] Cloning array of " << n << " elements of type " << typeid(T).name() << " (only the " << ncopy << " first) [" << ptr << "].\n";
                for (unsigned i = 0; i < ncopy; ++i)
                    copy_construct(newptr[i], ptr[i]);
                for (unsigned i = ncopy; i < n; ++i)
                    construct(newptr[i]);
                return newptr;
            }
      
#if __cplusplus >= 201103L
            template<class S>
            static pointer_type clone_move(S * ptr, unsigned int n)
            // Creates an array of n elements of type T. The entries are copy-constructed from the given
            // array ptr (of type S).
            {
                pointer_type newptr = create(n, true);
//          std::cout << "[#:" << newptr << "] Cloning array of " << n << " elements of type " << typeid(T).name() << " [" << ptr << "].\n";
                for (unsigned i = 0; i < n; ++i)
                    copy_construct(newptr[i], std::move(ptr[i]));
                return newptr;
            }
      
            template<class S>
            static pointer_type clone_move(S * ptr, unsigned int n, unsigned int ncopy)
            // Creates an array of n elements of type T. The first ncopy elements are copy-constructed from
            // the given array ptr, the last ones are constructed with the default constructor.
            {
                pointer_type newptr = create(n, true);
//          std::cout << "[#:" << newptr << "] Cloning array of " << n << " elements of type " << typeid(T).name() << " (only the " << ncopy << " first) [" << ptr << "].\n";
                for (unsigned i = 0; i < ncopy; ++i)
                    copy_construct(newptr[i], std::move(ptr[i]));
                for (unsigned i = ncopy; i < n; ++i)
                    construct(newptr[i]);
                return newptr;
            }
      
            template<class S, class SS>
            static pointer_type clone_move(S * ptr, unsigned int n, unsigned int ncopy, const & SS copyobj)
            // Creates an array of n elements of type T. The first ncopy elements are copy-constructed from
            // the given array ptr, the last ones are constructed with the default constructor.
            {
                pointer_type newptr = create(n, true);
//          std::cout << "[#:" << newptr << "] Cloning array of " << n << " elements of type " << typeid(T).name() << " (only the " << ncopy << " first) [" << ptr << "].\n";
                for (unsigned i = 0; i < ncopy; ++i)
                    copy_construct(newptr[i], std::move(ptr[i]));
                for (unsigned i = ncopy; i < n; ++i)
                    construct(newptr[i]);
                return newptr;
            }
#endif
            
            static void free(pointer_type ptr, unsigned int n)
            // Destroys an  array of n elements of type T.
            {
//          std::cout << "[#:" << ptr << "] Releasing array of " << n << " elements of type " << typeid(T).name() << ".\n";
                release(ptr, n);
            }
        };
        
        #endif
        
        #ifdef __VECMAT_MM_FAILSAFE
        
        template<class T>
        class storage_traits
        {
        public:
            typedef T * pointer_type;
            typedef T & ref_type;
            typedef const T * constpointer_type;
            typedef const T & constref_type;
      
            static void construct(ref_type ref)
            // Constructs an object via default constructor.
            {
                ::new(static_cast<void *>(&ref)) T();
            }

            template<class S>
            static void copy_construct(ref_type ref, const S & obj)
            // Constructs an object via copy constructor.
            {
                ::new(static_cast<void *>(&ref)) T(obj);
            }
      
#if __cplusplus >= 201103L
            template<class S>
            static void move_construct(ref_type ref, S && obj)
            // Constructs an object via move constructor.
            {
                ::new(static_cast<void *>(&ref)) T(std::move(obj));
            }
#endif
            
            static pointer_type alloc(unsigned int n)
            // Creates an array of n elements of type T. The entries are created using the default
            // constructor.
            {
                return new T[n];
            }
      
            static pointer_type alloc_dontconstruct(unsigned int n)
            // Creates an array of n elements of type T. The entries are not constructed.
            {
                return reinterpret_cast<pointer_type>(::operator new(sizeof(T) * n));
            }
      
            template<class S>
            static pointer_type alloc(const S & obj, unsigned int n)
            // Creates an array of n elements of type T, all being copies of obj (of type S).
            {
                pointer_type newptr = new T[n];
                for (unsigned i = 0; i < n; ++i)
                    newptr[i] = obj;
                return newptr;
            }
      
            template<class S>
            static pointer_type clone(S * ptr, unsigned int n)
            // Creates an array of n elements of type T. The entries are copy-constructed from the given
            // array ptr (of type S).
            {
                pointer_type newptr = new T[n];
                for (unsigned i = 0; i < n; ++i)
                    newptr[i] = ptr[i];
                return newptr;
            }
      
            template<class S>
            static pointer_type clone(S * ptr, unsigned int n, unsigned int ncopy)
            // Creates an array of n elements of type T. The first ncopy elements are copy-constructed from
            // the given array ptr, the last ones are constructed with the default constructor.
            {
                pointer_type newptr = new T[n];
                for (unsigned i = 0; i < ncopy; ++i)
                    newptr[i] = ptr[i];
                return newptr;
            }
      
            template<class S, class SS>
            static pointer_type clone(S * ptr, unsigned int n, unsigned int ncopy, const SS & copyobj)
            // Creates an array of n elements of type T. The first ncopy elements are copy-constructed from
            // the given array ptr, the last ones are constructed with the default constructor.
            {
                pointer_type newptr = new T[n];
                for (unsigned i = 0; i < ncopy; ++i)
                    newptr[i] = ptr[i];
                for (unsigned i = ncopy; i < n; ++i)
                    newptr[i] = copyobj;
                return newptr;
            }
            
#if __cplusplus >= 201103L
            template<class S>
            static pointer_type clone_move(S * ptr, unsigned int n)
            // Creates an array of n elements of type T. The entries are copy-constructed from the given
            // array ptr (of type S).
            {
                pointer_type newptr = new T[n];
                for (unsigned i = 0; i < n; ++i)
                    newptr[i] = std::move(ptr[i]);
                return newptr;
            }
      
            template<class S>
            static pointer_type clone_move(S * ptr, unsigned int n, unsigned int ncopy)
            // Creates an array of n elements of type T. The first ncopy elements are copy-constructed from
            // the given array ptr, the last ones are constructed with the default constructor.
            {
                pointer_type newptr = new T[n];
                for (unsigned i = 0; i < ncopy; ++i)
                    newptr[i] = std::move(ptr[i]);
                return newptr;
            }
      
            template<class S, class SS>
            static pointer_type clone_move(S * ptr, unsigned int n, unsigned int ncopy, const SS & copyobj)
            // Creates an array of n elements of type T. The first ncopy elements are copy-constructed from
            // the given array ptr, the last ones are constructed with the default constructor.
            {
                pointer_type newptr = new T[n];
                for (unsigned i = 0; i < ncopy; ++i)
                    newptr[i] = std::move(ptr[i]);
                for (unsigned i = ncopy; i < n; ++i)
                    newptr[i] = copyobj;
                return newptr;
            }
#endif
      
            static void free(pointer_type ptr, unsigned int n)
            // Destroys an  array of n elements of type T.
            {
                delete[] ptr;
            }
      
            inline static void check(pointer_type ptr, unsigned int n)
            {
            }
        };
        
        #endif
        
        #ifdef __VECMAT_MM_RELEASE
        
        template<class T>
        /**
           \brief Provides means to allocate, reallocate and release (linear) arrays of objects of
                  type T.
         */
        class storage_traits
        {
        public:
            typedef T * pointer_type;
            typedef T & ref_type;
            typedef const T * constpointer_type;
            typedef const T & constref_type;
      
            static void construct(ref_type ref)
            /** \brief Constructs an object via default constructor. */
            {
                ::new(static_cast<void *>(&ref)) T();
            }

            template<class S>
            static void copy_construct(ref_type ref, const S & obj)
            /** \brief Constructs an object via copy constructor. */
            {
                ::new(static_cast<void *>(&ref)) T(obj);
            }
      
#if __cplusplus >= 201103L
            template<class S>
            static void move_construct(ref_type ref, S && obj)
            /** \brief Constructs an object via move constructor. Only available for C++11. */
            {
                ::new(static_cast<void *>(&ref)) T(std::move(obj));
            }
#endif
            
            static pointer_type alloc(unsigned int n)
            /** \brief Creates an array of `n` elements of type `T`. The entries are created using
                       the default constructor. */
            {
                if (n == 0)
                    return NULL;
                pointer_type newptr = reinterpret_cast<pointer_type>(::operator new(sizeof(T) * n));
                for (unsigned i = 0; i < n; ++i)
                    ::new(static_cast<void *>(newptr + i)) T();
                return newptr;
            }
      
            static pointer_type alloc_dontconstruct(unsigned int n)
            /** \brief Creates an array of `n` elements of type `T`. The entries are not
                       constructed. */
            {
                if (n == 0)
                    return NULL;
                pointer_type newptr = reinterpret_cast<pointer_type>(::operator new(sizeof(T) * n));
                return newptr;
            }
      
            template<class S>
            static pointer_type alloc(const S & obj, unsigned int n)
            /** \brief Creates an array of `n` elements of type `T`, all being copies of `obj` (of
                       type `S`). */
            {
                if (n == 0)
                    return NULL;
                pointer_type newptr = reinterpret_cast<pointer_type>(::operator new(sizeof(T) * n));
                for (unsigned i = 0; i < n; ++i)
                    ::new(static_cast<void *>(newptr + i)) T(obj);
                return newptr;
            }
      
            template<class S>
            static pointer_type clone(S * ptr, unsigned int n)
            /** \brief Creates an array of `n` elements of type `T`. The entries are
                       copy-constructed from the given array `ptr` (with elements of type `S`). */
            {
                if (n == 0)
                    return NULL;
                pointer_type newptr = reinterpret_cast<pointer_type>(::operator new(sizeof(T) * n));
                for (unsigned i = 0; i < n; ++i)
                    ::new(static_cast<void *>(newptr + i)) T(ptr[i]);
                return newptr;
            }
      
            template<class S>
            static pointer_type clone(S * ptr, unsigned int n, unsigned int ncopy)
            /** \brief Creates an array of `n` elements of type `T`. The first `ncopy` elements are
                       copy-constructed from the given array `ptr`, the last ones are constructed
                       with the default constructor. */
            {
                if (n == 0)
                    return NULL;
                assert(ncopy <= n);
                pointer_type newptr = reinterpret_cast<pointer_type>(::operator new(sizeof(T) * n));
                for (unsigned i = 0; i < ncopy; ++i)
                    ::new(static_cast<void *>(newptr + i)) T(ptr[i]);
                for (unsigned i = ncopy; i < n; ++i)
                    ::new(static_cast<void *>(newptr + i)) T();
                return newptr;
            }
      
            template<class S, class SS>
            static pointer_type clone(S * ptr, unsigned int n, unsigned int ncopy, const SS & copyobj)
            /** \brief Creates an array of `n` elements of type `T`. The first `ncopy` elements are
                       copy-constructed from the given array `ptr`, the last ones are constructed
                       with the default constructor. */
            {
                if (n == 0)
                    return NULL;
                assert(ncopy <= n);
                pointer_type newptr = reinterpret_cast<pointer_type>(::operator new(sizeof(T) * n));
                for (unsigned i = 0; i < ncopy; ++i)
                    ::new(static_cast<void *>(newptr + i)) T(ptr[i]);
                for (unsigned i = ncopy; i < n; ++i)
                    ::new(static_cast<void *>(newptr + i)) T(copyobj);
                return newptr;
            }
      
#if __cplusplus >= 201103L
            template<class S>
            static pointer_type clone_move(S * ptr, unsigned int n)
            /** \brief Creates an array of `n` elements of type `T`. The entries are
                       move-constructed from the given array `ptr` (of type `S`). */
            {
                if (n == 0)
                    return NULL;
                pointer_type newptr = reinterpret_cast<pointer_type>(::operator new(sizeof(T) * n));
                for (unsigned i = 0; i < n; ++i)
                    ::new(static_cast<void *>(newptr + i)) T(std::move(ptr[i]));
                return newptr;
            }
       
            template<class S>
            static pointer_type clone_move(S * ptr, unsigned int n, unsigned int ncopy)
            /** \brief Creates an array of `n` elements of type `T`. The first `ncopy` elements are
                       move-constructed from the given array `ptr`, the last ones are constructed
                       with the default constructor. */
            {
                if (n == 0)
                    return NULL;
                assert(ncopy <= n);
                pointer_type newptr = reinterpret_cast<pointer_type>(::operator new(sizeof(T) * n));
                for (unsigned i = 0; i < ncopy; ++i)
                    ::new(static_cast<void *>(newptr + i)) T(std::move(ptr[i]));
                for (unsigned i = ncopy; i < n; ++i)
                    ::new(static_cast<void *>(newptr + i)) T();
                return newptr;
            }
      
            template<class S, class SS>
            static pointer_type clone_move(S * ptr, unsigned int n, unsigned int ncopy, const SS & copyobj)
            /** \brief Creates an array of `n` elements of type `T`. The first `ncopy` elements are
                       move-constructed from the given array `ptr`, the last ones are constructed
                       with the default constructor. */
            {
                if (n == 0)
                    return NULL;
                assert(ncopy <= n);
                pointer_type newptr = reinterpret_cast<pointer_type>(::operator new(sizeof(T) * n));
                for (unsigned i = 0; i < ncopy; ++i)
                    ::new(static_cast<void *>(newptr + i)) T(std::move(ptr[i]));
                for (unsigned i = ncopy; i < n; ++i)
                    ::new(static_cast<void *>(newptr + i)) T(copyobj);
                return newptr;
            }
#endif
            
            static void free(pointer_type ptr, unsigned int n)
            /** \brief Destroys an array of `n` elements of type `T`. */
            {
                if (n == 0)
                    return;
                for (unsigned i = 0; i < n; ++i)
                    (ptr + i)->~T();
                ::operator delete(ptr);
            }
      
            inline static void check(pointer_type ptr, unsigned int n)
            /** \brief Verifies the integrity of an array of `n` elements of type `T`. Does nothing
                       for the release memory manager. */
            {
            }
        };
        
        #endif
    }
}

#endif

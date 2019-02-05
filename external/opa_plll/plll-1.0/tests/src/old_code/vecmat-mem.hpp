#ifndef PLLL_INCLUDE_GUARD__VECMAT_MEM_HPP
#define PLLL_INCLUDE_GUARD__VECMAT_MEM_HPP

#include <cassert>
#include <cstddef>
//#include <iostream>

namespace plll
{
    namespace linalgold
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
                char * realptr = (char *)malloc(SIZE + (DEBUG_PAD_BEGIN + DEBUG_PAD_END) * 4);
                char * realptrend = realptr + SIZE + DEBUG_PAD_BEGIN * 4;
                for (int i = 0; i < DEBUG_PAD_BEGIN; ++i)
                    ((unsigned int*)realptr)[i] = __VECMAT_MM_DEBUG__BEGIN_OK;
                for (int i = 0; i < DEBUG_PAD_END; ++i)
                    ((unsigned int*)realptrend)[i] = __VECMAT_MM_DEBUG__END_OK;
                T * ptr = (T*)(realptr + DEBUG_PAD_BEGIN * 4);
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
                char * realptr = (char *)ptr - DEBUG_PAD_BEGIN * 4;
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
                    ((unsigned int*)ptr)[i] = __VECMAT_MM_DEBUG__DEAD;
                // Check padding and stop
                for (int i = 0; i < DEBUG_PAD_BEGIN; ++i)
                {
                    assert(((unsigned int*)realptr)[i] == __VECMAT_MM_DEBUG__BEGIN_OK);
                    ((unsigned int*)realptr)[i] = __VECMAT_MM_DEBUG__BEGIN_REL;
                }
                for (int i = 0; i < DEBUG_PAD_END; ++i)
                {
                    assert(((unsigned int*)realptrend)[i] == __VECMAT_MM_DEBUG__END_OK);
                    ((unsigned int*)realptrend)[i] = __VECMAT_MM_DEBUG__END_REL;
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
                char * realptr = (char *)ptr - DEBUG_PAD_BEGIN * 4;
                char * realptrend = realptr + SIZE + DEBUG_PAD_BEGIN * 4;
//          std::cout << "CHK: " << (void*)realptr << " " << (void*)ptr << " " << sizeof(T) << "x" << n << "=" << SIZE << " " << (void*)realptrend << "\n";
                for (int i = 0; i < DEBUG_PAD_BEGIN; ++i)
                    assert(((unsigned int*)realptr)[i] == __VECMAT_MM_DEBUG__BEGIN_OK);
                for (int i = 0; i < DEBUG_PAD_END; ++i)
                    assert(((unsigned int*)realptrend)[i] == __VECMAT_MM_DEBUG__END_OK);
            }
      
            typedef T * pointer_type;
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
                pointer_type newptr = create(n);
//          std::cout << "[#:" << newptr << "] Creating array of " << n << " elements of type " << typeid(T).name() << ", constructed from object [" << &obj << "].\n";
                for (unsigned i = 0; i < n; ++i)
                    newptr[i] = obj;
                return newptr;
            }
    
            template<class S>
            static pointer_type clone(S * ptr, unsigned int n)
            // Creates an array of n elements of type T. The entries are copy-constructed from the given
            // array ptr (of type S).
            {
                pointer_type newptr = create(n);
//          std::cout << "[#:" << newptr << "] Cloning array of " << n << " elements of type " << typeid(T).name() << " [" << ptr << "].\n";
                for (unsigned i = 0; i < n; ++i)
                    newptr[i] = ptr[i];
                return newptr;
            }
      
            template<class S>
            static pointer_type clone(S * ptr, unsigned int n, unsigned int ncopy)
            // Creates an array of n elements of type T. The first ncopy elements are copy-constructed from
            // the given array ptr, the last ones are constructed with the default constructor.
            {
                pointer_type newptr = create(n);
//          std::cout << "[#:" << newptr << "] Cloning array of " << n << " elements of type " << typeid(T).name() << " (only the " << ncopy << " first) [" << ptr << "].\n";
                for (unsigned i = 0; i < ncopy; ++i)
                    newptr[i] = ptr[i];
                return newptr;
            }
      
            template<class S, class SS>
            static pointer_type clone(S * ptr, unsigned int n, unsigned int ncopy, const & SS copyobj)
            // Creates an array of n elements of type T. The first ncopy elements are copy-constructed from
            // the given array ptr, the last ones are constructed with the default constructor.
            {
                pointer_type newptr = create(n);
//          std::cout << "[#:" << newptr << "] Cloning array of " << n << " elements of type " << typeid(T).name() << " (only the " << ncopy << " first) [" << ptr << "].\n";
                for (unsigned i = 0; i < ncopy; ++i)
                    newptr[i] = ptr[i];
                for (unsigned i = ncopy; i < n; ++i)
                    newptr[i] = copyobj;
                return newptr;
            }
      
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
        class storage_traits
        {
        public:
            typedef T * pointer_type;
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
      
            static pointer_type alloc(unsigned int n)
            // Creates an array of n elements of type T. The entries are created using the default
            // constructor.
            {
                if (n == 0)
                    return NULL;
                pointer_type newptr = reinterpret_cast<pointer_type>(::operator new(sizeof(T) * n));
                for (unsigned i = 0; i < n; ++i)
                    ::new(static_cast<void *>(newptr + i)) T();
                return newptr;
            }
      
            static pointer_type alloc_dontconstruct(unsigned int n)
            // Creates an array of n elements of type T. The entries are not constructed.
            {
                if (n == 0)
                    return NULL;
                pointer_type newptr = reinterpret_cast<pointer_type>(::operator new(sizeof(T) * n));
                return newptr;
            }
      
            template<class S>
            static pointer_type alloc(const S & obj, unsigned int n)
            // Creates an array of n elements of type T, all being copies of obj (of type S).
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
            // Creates an array of n elements of type T. The entries are copy-constructed from the given
            // array ptr (of type S).
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
            // Creates an array of n elements of type T. The first ncopy elements are copy-constructed from
            // the given array ptr, the last ones are constructed with the default constructor.
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
            // Creates an array of n elements of type T. The first ncopy elements are copy-constructed from
            // the given array ptr, the last ones are constructed with the default constructor.
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
      
            static void free(pointer_type ptr, unsigned int n)
            // Destroys an  array of n elements of type T.
            {
                if (n == 0)
                    return;
                for (unsigned i = 0; i < n; ++i)
                    (ptr + i)->~T();
                ::operator delete(ptr);
            }
      
            inline static void check(pointer_type ptr, unsigned int n)
            {
            }
        };
        
        #endif
    }
}

#endif

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

#ifndef PLLL_INCLUDE_GUARD__MYALLOC_HPP
#define PLLL_INCLUDE_GUARD__MYALLOC_HPP

#include <memory>
#include <limits>

namespace plll
{
    void initTLAlloc(); // should be called once per thread

    void * TLalloc(std::size_t);
    void * TLrealloc(void *, std::size_t);
    void * TLrealloc(void *, std::size_t, std::size_t);
    void TLfree(void *);
    void TLfree(void *, std::size_t);

    void * TLSalloc(std::size_t);
    void * TLSrealloc(void *, std::size_t);
    void * TLSrealloc(void *, std::size_t, std::size_t);
    void TLSfree(void *);
    void TLSfree(void *, std::size_t);

    template<class T>
    class TLAlloc
    /*
      A thread-local allocator. If something has been allocated in one thread, it has to be released in
      that same thread. Moreover, when the thread ends, the memory is freed automatically.
  
      The allocator also has a reallocate method, which works as realloc():
      pointer reallocate(pointer, size_type);
    */
    {
    public:
        typedef T value_type;
        typedef value_type * pointer;
        typedef const value_type * const_pointer;
        typedef value_type & reference;
        typedef const value_type & const_reference;
        typedef std::size_t size_type;
        typedef std::ptrdiff_t difference_type;

        template<class U>
        struct rebind
        {
            typedef TLAlloc<U> other;
        };
    
    public: 
        inline explicit TLAlloc()
        {
        }
    
        inline ~TLAlloc()
        {
        }
    
        inline explicit TLAlloc(const TLAlloc &)
        {
        }
    
        template<class U>
        inline explicit TLAlloc(const TLAlloc<U> &)
        {
        }
    
        static inline pointer address(reference r)
        {
            return &r;
        }
    
        static inline const_pointer address(const_reference r)
        {
            return &r;
        }
    
        static inline pointer reallocate(pointer ptr, size_type cnt)
        {
            return reinterpret_cast<pointer>(TLrealloc(ptr, cnt * sizeof(T))); 
        }
    
        static inline pointer allocate(size_type cnt, typename std::allocator<void>::const_pointer = 0)
        {
            return reinterpret_cast<pointer>(TLalloc(cnt * sizeof(T))); 
        }
    
        static inline void deallocate(pointer p, size_type cnt)
        {
            TLfree(p, cnt);
        }
    
        static inline size_type max_size()
        { 
            return std::numeric_limits<size_type>::max() / sizeof(T);
        }
    
        static inline void construct(pointer p, const T & t)
        {
            new(p) T(t);
        }
    
        static inline void destroy(pointer p)
        {
            p->~T();
        }
    
        inline bool operator == (const TLAlloc &) const
        {
            return true;
        }
    
        inline bool operator != (const TLAlloc &) const
        {
            return false;
        }
    };

    template<class T>
    class TLSAlloc
    /*
      A "safe" thread-local allocator. Memory can be released in other threads as well, and when the
      thread is terminated the memory is kept.
  
      The allocator also has a reallocate method, which works as realloc():
      pointer reallocate(pointer, size_type);
    */
    {
    public:
        typedef T value_type;
        typedef value_type * pointer;
        typedef const value_type * const_pointer;
        typedef value_type & reference;
        typedef const value_type & const_reference;
        typedef std::size_t size_type;
        typedef std::ptrdiff_t difference_type;

        template<class U>
        struct rebind
        {
            typedef TLSAlloc<U> other;
        };
    
    public: 
        inline explicit TLSAlloc()
        {
        }
    
        inline ~TLSAlloc()
        {
        }
    
        inline explicit TLSAlloc(const TLSAlloc &)
        {
        }
    
        template<class U>
        inline explicit TLSAlloc(const TLSAlloc<U> &)
        {
        }
    
        static inline pointer address(reference r)
        {
            return &r;
        }
    
        static inline const_pointer address(const_reference r)
        {
            return &r;
        }
    
        static inline pointer reallocate(pointer ptr, size_type cnt)
        {
            return reinterpret_cast<pointer>(TLSrealloc(ptr, cnt * sizeof(T))); 
        }
    
        static inline pointer allocate(size_type cnt, typename std::allocator<void>::const_pointer = 0)
        {
            return reinterpret_cast<pointer>(TLSalloc(cnt * sizeof(T))); 
        }
    
        static inline void deallocate(pointer p, size_type cnt)
        {
            TLSfree(p, cnt);
        }
    
        static inline size_type max_size()
        { 
            return std::numeric_limits<size_type>::max() / sizeof(T);
        }
    
        static inline void construct(pointer p, const T & t)
        {
            new(p) T(t);
        }
    
        static inline void destroy(pointer p)
        {
            p->~T();
        }
    
        inline bool operator == (const TLSAlloc &) const
        {
            return true;
        }
    
        inline bool operator != (const TLSAlloc &) const
        {
            return false;
        }
    };
}

#endif

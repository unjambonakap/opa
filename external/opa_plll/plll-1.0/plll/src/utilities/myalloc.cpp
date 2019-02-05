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

#include <plll/config.hpp>
#include "myalloc.hpp"
#include <boost/thread.hpp>
#include <cstring>

#include "dlalloc.h"

namespace plll
{
    namespace
    {
        class NoCheck
        {
        public:
            static unsigned extra_before()
            {
                return 0;
            }
        
            static unsigned extra_after()
            {
                return 0;
            }
        
            static void alloc(void *, std::size_t)
            {
            }
        
            static void realloc_before(void *, std::size_t)
            {
            }
        
            static void realloc_after(void *, std::size_t)
            {
            }
        
            static void free(void *)
            {
            }
        };
    
        class SimpleCheck
        {
        private:
            enum { c_ok = 0xAFFEAFFE, c_bad = 0xDEADBEEF };
        
        public:
            static unsigned extra_before()
            {
                return sizeof(unsigned);
            }
        
            static unsigned extra_after()
            {
                return 0;
            }
        
            static void alloc(void * m, std::size_t)
            {
                *((unsigned*)m) = c_ok;
            }
        
            static void realloc_before(void * m, std::size_t)
            {
                assert(*(unsigned*)m == c_ok);
            }
        
            static void realloc_after(void * m, std::size_t)
            {
                assert(*(unsigned*)m == c_ok);
            }
        
            static void free(void * m)
            {
                assert(*(unsigned*)m == c_ok);
                *(unsigned*)m = c_bad;
            }
        };
    
        class ExtendedCheck
        {
        private:
            enum { c_ok = 0xAFFEAFFE, c_bad = 0xDEADBEEF, c_before = 16, c_after = 16 };
        
            static void check(unsigned * mm)
            {
                assert(mm[0] == c_ok);
                assert(mm[1] != c_bad);
                for (unsigned i = 0; i < c_before; ++i)
                    assert(mm[2 + i] == c_ok);
#ifndef NDEBUG
                unsigned * mmm = (unsigned*)((char*)mm + ((mm[1] + sizeof(unsigned) - 1) & ~(sizeof(unsigned) - 1)) + extra_before());
                for (unsigned i = 0; i < c_after; ++i)
                    assert(mmm[i] == c_ok);
#endif
            }
        
        public:
            static unsigned extra_before()
            {
                return (2 + c_before) * sizeof(unsigned);
            }
        
            static unsigned extra_after()
            {
                return (1 + c_after) * sizeof(unsigned) - 1;
            }
        
            static void alloc(void * m, std::size_t s)
            {
                unsigned * mm = (unsigned*)m;
                mm[0] = c_ok;
                mm[1] = s;
                for (unsigned i = 0; i < c_before; ++i)
                    mm[2 + i] = c_ok;
                assert(mm[1] == s);
                unsigned * mmm = (unsigned*)((char*)mm + ((mm[1] + sizeof(unsigned) - 1) & ~(sizeof(unsigned) - 1)) + extra_before());
                for (unsigned i = 0; i < c_after; ++i)
                    mmm[i] = c_ok;
            }
        
            static void realloc_before(void * m, std::size_t s)
            {
                unsigned * mm = (unsigned*)m;
                check(mm);
            }
        
            static void realloc_after(void * m, std::size_t s)
            {
                unsigned * mm = (unsigned*)m;
                mm[1] = s;
                unsigned * mmm = (unsigned*)((char*)mm + ((mm[1] + sizeof(unsigned) - 1) & ~(sizeof(unsigned) - 1)) + extra_before());
                for (unsigned i = 0; i < c_after; ++i)
                    mmm[i] = c_ok;
            }
        
            static void free(void * m)
            {
                unsigned * mm = (unsigned*)m;
                check(mm);
                unsigned * mmm = (unsigned*)((char*)mm + ((mm[1] + sizeof(unsigned) - 1) & ~(sizeof(unsigned) - 1)) + extra_before());
                for (unsigned i = 0; i < c_before + 2; ++i)
                    mm[i] = c_bad;
                for (unsigned i = 0; i < c_after; ++i)
                    mmm[i] = c_bad;
            }
        };
    
        class NoTS
        {
        public:
            class Lock
            {
            public:
                inline Lock(NoTS &)
                {
                }
            };
        };
    
        class TSMutex
        {
        private:
            boost::mutex d_mutex;
        
        public:
            class Lock
            {
            private:
                boost::mutex & d_mutex_ref;
            
            public:
                inline Lock(TSMutex & mutex)
                    : d_mutex_ref(mutex.d_mutex)
                {
                    d_mutex_ref.lock();
                }
            
                inline ~Lock()
                {
                    d_mutex_ref.unlock();
                }
            };
        };
    
        template<class ThreadSafe, class Checker>
        class RealAllocatorImpl
        {
        private:
            ThreadSafe d_ts;
            Checker d_check;
            mspace d_mspace;
        
        public:
            inline RealAllocatorImpl()
            {
                d_mspace = create_mspace(0, 0);
            }

            inline ~RealAllocatorImpl()
            {
                destroy_mspace(d_mspace);
            }
        
            inline void * alloc(std::size_t s)
            {
                typename ThreadSafe::Lock lock(d_ts);
                void * m = mspace_malloc(d_mspace, s + d_check.extra_before() + d_check.extra_after());
//            void * m = ::malloc(s + d_check.extra_before() + d_check.extra_after());
                d_check.alloc(m, s);
                return (char*)m + d_check.extra_before();
            }
        
            inline void * realloc(void * p, std::size_t s)
            {
                typename ThreadSafe::Lock lock(d_ts);
                void * m = (char*)p - d_check.extra_before();
                d_check.realloc_before(m, s);
                m = mspace_realloc(d_mspace, m, s + d_check.extra_before() + d_check.extra_after());
//            m = ::realloc(m, s + d_check.extra_before() + d_check.extra_after());
                d_check.realloc_after(m, s);
                return (char*)m + d_check.extra_before();
            }
        
            inline void free(void * p)
            {
                typename ThreadSafe::Lock lock(d_ts);
                void * m = (char*)p - d_check.extra_before();
                d_check.free(m);
                return mspace_free(d_mspace, m);
//            return ::free(m);
            }
        
            inline std::size_t size(void * p)
            // Gives an upper bound for the space allocated at p. This amount of bytes can actually be
            // used.
            {
                typename ThreadSafe::Lock lock(d_ts);
                return mspace_usable_size((char*)p - d_check.extra_before()) - d_check.extra_before() - d_check.extra_after();
            }
        };
    
        // The template RealAllocator<ThreadSafe, storeAsPointer, freeAllocator> essentially stores a
        // RealAllocatorImpl<ThreadSafe> object. It either does so by pointer or not (storeAsPointer),
        // and it allows to not free the allocator when using a pointer (freeAllocator). The combination
        // storeAsPointer = false and freeAllocator = false is not allowed.
        template<class ThreadSafe, class Checker, bool storeAsPointer, bool freeAllocator>
        class RealAllocator;
    
        template<class ThreadSafe, class Checker, bool freeAllocator>
        class RealAllocator<ThreadSafe, Checker, true, freeAllocator>
        {
        private:
            RealAllocatorImpl<ThreadSafe, Checker> * d_alloc;
        
        public:
            inline RealAllocator()
            {
                d_alloc = new RealAllocatorImpl<ThreadSafe, Checker>();
            }

            inline ~RealAllocator()
            {
                if (freeAllocator)
                    delete d_alloc;
            }
        
            inline void * alloc(std::size_t s)
            {
                return d_alloc->alloc(s);
            }
        
            inline void * realloc(void * p, std::size_t s)
            {
                return d_alloc->realloc(p, s);
            }
        
            inline void free(void * p)
            {
                d_alloc->free(p);
            }
        
            inline std::size_t size(void * p)
            // Gives an upper bound for the space allocated at p. This amount of bytes can actually be
            // used.
            {
                return d_alloc->size(p);
            }
        
            inline RealAllocatorImpl<ThreadSafe, Checker> * getImplAddress()
            {
                return d_alloc;
            }
        };
    
        template<class ThreadSafe, class Checker>
        class RealAllocator<ThreadSafe, Checker, false, true>
        {
        private:
            RealAllocatorImpl<ThreadSafe, Checker> d_alloc;
        
        public:
            inline void * alloc(std::size_t s)
            {
                return d_alloc.alloc(s);
            }
        
            inline void * realloc(void * p, std::size_t s)
            {
                return d_alloc.realloc(p, s);
            }
        
            inline void free(void * p)
            {
                d_alloc.free(p);
            }
        
            inline std::size_t size(void * p)
            // Gives an upper bound for the space allocated at p. This amount of bytes can actually be
            // used.
            {
                return d_alloc.size(p);
            }
        
            inline RealAllocatorImpl<ThreadSafe, Checker> * getImplAddress()
            {
                return &d_alloc;
            }
        };
    
        template<class ThreadSafe, class Checker, bool storeDeallocAddress, bool freeAllocator>
        class Allocator;
    
        template<class ThreadSafe, class Checker, bool freeAllocator>
        class Allocator<ThreadSafe, Checker, true, freeAllocator>
        {
        public:
            typedef RealAllocatorImpl<ThreadSafe, Checker> RealAllocImpl;
            typedef RealAllocator<ThreadSafe, Checker, true, freeAllocator> RealAlloc;
        
        private:
            RealAlloc d_alloc;
        
        public:
            inline void * alloc(std::size_t s)
            {
                char * p = static_cast<char*>(d_alloc.alloc(sizeof(RealAllocImpl*) + s));
                reinterpret_cast<RealAllocImpl**>(p)[0] = d_alloc.getImplAddress();
                return p + sizeof(RealAllocImpl*);
            }
        
            inline void * realloc(void * p, std::size_t s)
            {
                if (p != NULL)
                {
                    char * rp = static_cast<char*>(p) - sizeof(RealAllocImpl*);
                    RealAllocImpl * ra = reinterpret_cast<RealAllocImpl**>(rp)[0];
                    const bool reallocChangesOwnership = false; // don't change ownership, i.e. use
                    // second branch of the following `if'
                    // statement
                    if (reallocChangesOwnership && (ra != d_alloc.getImplAddress()))
                    {
                        char * np = static_cast<char*>(alloc(s));
                        std::memcpy(np, p, std::min(s, ra->size(p)));
                        free(p);
                        return np;
                    }
                    else
                    {
                        rp = static_cast<char*>(ra->realloc(rp, sizeof(RealAllocImpl*) + s));
                        return rp + sizeof(RealAllocImpl*);
                    }
                }
                else
                    return alloc(s);
            }
        
            inline void free(void * p)
            {
                if (p != NULL)
                {
                    char * rp = static_cast<char*>(p) - sizeof(RealAllocImpl*);
                    reinterpret_cast<RealAllocImpl**>(rp)[0]->free(rp);
                }
            }
        };
    
        template<class ThreadSafe, class Checker, bool freeAllocator>
        class Allocator<ThreadSafe, Checker, false, freeAllocator>
        {
        public:
            typedef RealAllocator<ThreadSafe, Checker, !freeAllocator, freeAllocator> RealAlloc;
        
        private:
            RealAlloc d_alloc;
        
        public:
            inline void * alloc(std::size_t s)
            {
                return d_alloc.alloc(s);
            }
        
            inline void * realloc(void * p, std::size_t s)
            {
                return d_alloc.realloc(p, s);
            }
        
            inline void free(void * p)
            {
                d_alloc.free(p);
            }
        };

//    typedef Allocator<NoTS, ExtendedCheck, false, true> TLAllocData_T;
//    typedef Allocator<TSMutex, ExtendedCheck, true, false> TLSAllocData_T;
    
//    typedef Allocator<NoTS, SimpleCheck, false, true> TLAllocData_T;
//    typedef Allocator<TSMutex, SimpleCheck, true, false> TLSAllocData_T;
    
        typedef Allocator<NoTS, NoCheck, false, true> TLAllocData_T;
        typedef Allocator<TSMutex, NoCheck, true, false> TLSAllocData_T;
    
        boost::thread_specific_ptr<TLAllocData_T> TLAllocData;
        boost::thread_specific_ptr<TLSAllocData_T> TLSAllocData;
    }

// Initialization

    void initTLAlloc()
    {
        if (TLAllocData.get() == NULL)
            TLAllocData.reset(new TLAllocData_T());
        if (TLSAllocData.get() == NULL)
            TLSAllocData.reset(new TLSAllocData_T());
    }

// Thread-Local Allocator (unsafe)

    void * TLalloc(std::size_t s)
    {
        return TLAllocData->alloc(s);
    }

    void * TLrealloc(void * p, std::size_t s)
    {
        return TLAllocData->realloc(p, s);
    }

    void * TLrealloc(void * p, std::size_t /*olds*/, std::size_t news)
    {
        return TLAllocData->realloc(p, news);
    }

    void TLfree(void * p)
    {
        TLAllocData->free(p);
    }

    void TLfree(void * p, std::size_t /*s*/)
    {
        TLAllocData->free(p);
    }

// Safe Thread Local Allocator

    void * TLSalloc(std::size_t s)
    {
        return TLSAllocData->alloc(s);
    }

    void * TLSrealloc(void * p, std::size_t s)
    {
        return TLSAllocData->realloc(p, s);
    }

    void * TLSrealloc(void * p, std::size_t /*olds*/, std::size_t news)
    {
        return TLSAllocData->realloc(p, news);
    }

    void TLSfree(void * p)
    {
        TLSAllocData->free(p);
    }

    void TLSfree(void * p, std::size_t /*s*/)
    {
        TLSAllocData->free(p);
    }
}

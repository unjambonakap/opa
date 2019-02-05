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

#ifndef PLLL_INCLUDE_GUARD__TIMER_HPP
#define PLLL_INCLUDE_GUARD__TIMER_HPP

#include <cassert>

//////////////////////////////////////////////////////////////////////
/////////////////////////////////////// SOLARIS High Resolution Timer
#if defined(PLLL_INTERNAL_HR_SOLARIS)

#include <sys/time.h>
#include <sys/times.h>
#include <stdint.h>
#include <unistd.h>

namespace plll
{
    class Timer
    {
    public:
        typedef hrtime_t CountType;
    
    private:
        mutable CountType d_counter;
        mutable bool d_running;
    
    public:
        static const CountType d_one_timeunit;
        static const char * d_timeunit_name;
        
        Timer()
            : d_counter(0),
              d_running(false)
        {
        }
        
        void start() const
        {
            assert(!d_running);
            d_counter -= gethrtime();
            d_running = true;
        }
        
        void stop() const
        {
            assert(d_running);
            d_counter += gethrtime();
            d_running = false;
        }
        
        CountType elapsed() const
        {
            return d_running ? d_counter + gethrtime() : d_counter;
        }
    };

    class CPUTimer
    {
    public:
        typedef uint64_t CountType;
        
    private:
        mutable tms d_tms;
        mutable CountType d_counter;
        mutable bool d_running;
        
    public:
        static const CountType d_one_timeunit;
        static const char * d_timeunit_name;
        
        CPUTimer()
            : d_counter(0), d_running(false)
        {
        }
        
        void start() const
        {
            assert(!d_running);
            times(&d_tms);
            d_counter -= d_tms.tms_utime;
            d_running = true;
        }
        
        void stop() const
        {
            assert(d_running);
            times(&d_tms);
            d_counter += d_tms.tms_utime;
            d_running = false;
        }
        
        CountType elapsed() const
        {
            if (d_running)
            {
                times(&d_tms);
                return d_counter + d_tms.tms_utime;
            }
            else
                return d_counter;
        }
    };
}

//////////////////////////////////////////////////////////////////////
/////////////////////////////////////// POSIX High Resolution Timer
#elif defined(PLLL_INTERNAL_HR_POSIX)

#include <time.h>

namespace plll
{
    class Timer
    {
    public:
        typedef unsigned long long CountType;
    
    private:
        mutable CountType d_counter;
        mutable bool d_running;
    
        static inline unsigned long long timestamp()
        {
            timespec t;
            clock_gettime(CLOCK_MONOTONIC, &t);
            return ((unsigned long long)t.tv_sec) * 1000000000ul + (unsigned long long)t.tv_nsec;
        }
        
    public:
        static const CountType d_one_timeunit;
        static const char * d_timeunit_name;
        
        Timer()
            : d_counter(0), d_running(false)
        {
        }
        
        void start() const
        {
            assert(!d_running);
            d_counter -= timestamp();
            d_running = true;
        }
        
        void stop() const
        {
            assert(d_running);
            d_counter += timestamp();
            d_running = false;
        }
        
        CountType elapsed() const
        {
            return d_running ? d_counter + timestamp() : d_counter;
        }
    };
    
    class CPUTimer
    {
    public:
        typedef unsigned long long CountType;
        
    private:
        mutable CountType d_counter;
        mutable bool d_running;
        
        static inline unsigned long long timestamp()
        {
            timespec t;
            clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t);
            return ((unsigned long long)t.tv_sec) * 1000000000ul + (unsigned long long)t.tv_nsec;
        }
        
    public:
        static const CountType d_one_timeunit;
        static const char * d_timeunit_name;
        
        CPUTimer()
            : d_counter(0), d_running(false)
        {
        }
        
        void start() const
        {
            assert(!d_running);
            d_counter -= timestamp();
            d_running = true;
        }
        
        void stop() const
        {
            assert(d_running);
            d_counter += timestamp();
            d_running = false;
        }
        
        CountType elapsed() const
        {
            return d_running ? d_counter + timestamp() : d_counter;
        }
    };
}

//////////////////////////////////////////////////////////////////////
/////////////////////////////////////// Fallback: RDTSC readout
#else

#include <stdint.h>
#include <sys/times.h>
#include <stdint.h>
#include <unistd.h>

namespace plll
{
    class Timer
    {
    public:
        typedef uint64_t CountType;
    
    private:
        mutable CountType d_counter;
        mutable bool d_running;
    
        static inline uint64_t rdtsc()
        {
            uint32_t low, high;
            __asm volatile ("rdtsc" : "=a" (low), "=d" (high));
            return ((uint64_t)high << 32) | low;
        }
    
    public:
        static const CountType d_one_timeunit;
        static const char * d_timeunit_name;
    
        Timer()
            : d_counter(0), d_running(false)
        {
        }
        
        void start() const
        {
            assert(!d_running);
            d_counter -= rdtsc();
            d_running = true;
        }
        
        void stop() const
        {
            assert(d_running);
            d_counter += rdtsc();
            d_running = false;
        }
        
        CountType elapsed() const
        {
            return d_running ? d_counter + rdtsc() : d_counter;
        }
    };
    
    class CPUTimer
    {
    public:
        typedef uint64_t CountType;
        
    private:
        mutable tms d_tms;
        mutable CountType d_counter;
        mutable bool d_running;
        
    public:
        static const CountType d_one_timeunit;
        static const char * d_timeunit_name;
        
        CPUTimer()
            : d_counter(0), d_running(false)
        {
        }
        
        void start() const
        {
            assert(!d_running);
            times(&d_tms);
            d_counter -= d_tms.tms_utime;
            d_running = true;
        }
        
        void stop() const
        {
            assert(d_running);
            times(&d_tms);
            d_counter += d_tms.tms_utime;
            d_running = false;
        }
        
        CountType elapsed() const
        {
            if (d_running)
            {
                times(&d_tms);
                return d_counter + d_tms.tms_utime;
            }
            else
                return d_counter;
        }
    };
}

#endif

#endif

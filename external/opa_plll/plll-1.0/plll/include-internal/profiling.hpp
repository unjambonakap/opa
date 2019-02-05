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

#ifndef PLLL_INCLUDE_GUARD__PROFILING_HPP
#define PLLL_INCLUDE_GUARD__PROFILING_HPP

#include "timer.hpp"
#include <iostream>
#include <iomanip>
#include <climits>
#include <unistd.h>
#include <sys/times.h>

namespace plll
{
    class CheckTime
    {
    private:
        tms d_tms;
        clock_t d_starttime;
        clock_t d_endtime;
        
    public:
        CheckTime()
        {
            d_starttime = times(&d_tms);
        }
    
        ~CheckTime()
        {
            d_endtime = times(&d_tms);
            output();
        }
    
        void output() const
        {
          #ifndef PLLL_INTERNAL_HR_SOLARIS
            clock_t t = sysconf(_SC_CLK_TCK);
          #else
            clock_t t = CLK_TCK;
          #endif
            std::cout << "Real:   " << std::setiosflags(std::ios::fixed) << (double)(d_endtime - d_starttime) / (double)t << " seconds\n";
            std::cout << "User:   " << (double)d_tms.tms_utime / (double)t << " seconds\n";
            std::cout << "System: " << (double)d_tms.tms_stime / (double)t << " seconds\n";
        }
    };
    
    class TimedDataCollector
    {
    private:
        Timer d_timer;
        int d_no_calls;
        const char * d_name;
        
    public:
        const Timer & getTimer() const
        {
            return d_timer;
        }
        
        int getCalls() const
        {
            return d_no_calls;
        }
        
        const char * getName() const
        {
            return d_name;
        }
        
        TimedDataCollector(const char * name = NULL)
            : d_no_calls(0), d_name(name)
        {
        }
        
        void setName(const char * name)
        {
            d_name = name;
        }
        
        void output() const
        {
            if (d_name)
                std::cout << d_name << ": ";
            if (d_no_calls)
                std::cout << "total " << (double)d_timer.elapsed() / (double)Timer::d_one_timeunit
                          << " " << Timer::d_timeunit_name << "s during " << d_no_calls << " calls, i.e. ~"
                          << ((double)d_timer.elapsed() / (double)d_no_calls / (double)Timer::d_one_timeunit)
                          << " " << Timer::d_timeunit_name << "s per call.\n";
            else
                std::cout << "never called.\n";
        }
        
        ~TimedDataCollector()
        {
            output();
        }
        
        void begin()
        {
            ++d_no_calls;
            d_timer.start();
        }
        
        void end()
        {
            d_timer.stop();
        }
        
        void restart()
        // Should only be called between begin() and end()
        {
            ++d_no_calls;
        }
    };
    
    class TFDataCollector
    {
    private:
        long d_cur_total;
        long d_cur_true;
        long d_sum_total;
        long d_sum_true;
        int d_no_calls;
        double d_sum_fract;
        const char * d_name;
        
    public:
        const char * getName() const
        {
            return d_name;
        }
        
        TFDataCollector(const char * name = NULL)
            : d_cur_total(0), d_cur_true(0), d_sum_total(0), d_sum_true(0), d_no_calls(0), d_sum_fract(0.0), d_name(name)
        {
        }
        
        void output()
        {
            if (d_name)
                std::cout << d_name << ": ";
            if (d_no_calls)
                std::cout << "Total " << d_sum_true << " out of " << d_sum_total << " hits ("
                          << ((double)d_sum_true / (double)d_sum_total * 100.0) << "%); we had "
                          << d_no_calls << " instances with an average rate of "
                          << (d_sum_fract * 100.0 / (double)d_no_calls) << "%.\n";
            else
                std::cout << "never called.\n";
        }
        
        ~TFDataCollector()
        {
            output();
        }
        
        void begin()
        {
            ++d_no_calls;
            d_cur_total = d_cur_true = 0;
        }
        
        void hit(bool state)
        {
            ++d_cur_total;
            if (state)
                ++d_cur_true;
        }
        
        void end()
        {
            d_sum_total += d_cur_total;
            d_sum_true += d_cur_true;
            d_sum_fract += (double)d_cur_true / (double)d_cur_total;
        }
        
        void restart()
        // Should only be called between begin() and end()
        {
            end();
            begin();
        }
    };
    
    template<class T>
    class Profiler
    {
    private:
        T & d_dc;
        
    public:
        Profiler(T & dc)
            : d_dc(dc)
        {
            d_dc.begin();
        }
        
        ~Profiler()
        {
            d_dc.end();
        }
        
        void restart()
        {
            d_dc.restart();
        }
    };
}

#endif

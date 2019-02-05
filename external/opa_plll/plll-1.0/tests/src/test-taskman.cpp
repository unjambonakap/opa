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

#include "profiling.hpp"
#include "taskmanager.hpp"
#include <plll/arithmetic.hpp>

using namespace plll;

class StupidJob : public TaskManager::Job
{
private:
    unsigned d_id;
    unsigned d_tmp;
    
public:
    StupidJob(unsigned id, arithmetic::RandomNumberGenerator & rng)
        : d_id(id), d_tmp(rng.random(100000000))
    {
    }
    
    StupidJob(unsigned id, unsigned tmp)
        : d_id(id), d_tmp(tmp)
    {
    }
    
    virtual ~StupidJob() { }
    
    virtual void run(unsigned thread_no, TaskManager * tm)
    // The parameter thread_no ranges from 0 to number of threads - 1. The parameter tm points to
    // the task manager.
    {
        unsigned tmp = d_tmp;
        for (unsigned i = 0; i < d_tmp; ++i)
        {
            if (tmp & 1)
                tmp += i;
            else
                tmp -= i;
            tmp %= d_tmp;
        }
        std::ostringstream s;
        s << "#" << d_id << ": running at thread #" << thread_no << " with result " << tmp << " and initialization parameter " << d_tmp << "\n";
//        std::cout << s.str();
        
        if (tmp > 10)
        {
            tm->enqueueJob(new StupidJob(d_id, tmp));
            tm->enqueueJob(new StupidJob(d_id * 10, sqrt(tmp)));
        }
    }
};

int main()
{
    arithmetic::initArithmeticThreadAllocators();
    arithmetic::RandomNumberGenerator rng;
    TaskManager tm(4);
    for (unsigned i = 0; i < 10; ++i)
        tm.enqueueJob(new StupidJob(i, rng));
    tm.waitForDone();
}

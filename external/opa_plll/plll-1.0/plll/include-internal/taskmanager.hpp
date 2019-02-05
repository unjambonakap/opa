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

#ifndef PLLL_INCLUDE_GUARD__TASKMANAGER_HPP
#define PLLL_INCLUDE_GUARD__TASKMANAGER_HPP

#include <deque>
#include <boost/thread.hpp>

namespace plll
{
    class TaskManager
    {
    public:
        class Job
        {
        public:
            virtual ~Job() { }
            
            virtual void run(unsigned thread_no, TaskManager * tm) = 0;
            // The parameter thread_no ranges from 0 to number of threads - 1. The parameter tm points to
            // the task manager.
        };
        
    private:
        unsigned d_nothreads;
        
        boost::condition_variable_any d_jobs_cond;
        mutable boost::shared_mutex d_jobs_mutex;
        std::list<Job *> d_jobs;
        unsigned d_running;
        
        boost::condition_variable d_waitfordone_cond;
        boost::mutex d_waitfordone_mutex;
        
        boost::condition_variable d_waitfordone_single_cond;
        bool d_done_single;
        boost::mutex d_waitfordone_single_mutex;
        
        class Thread;
        
        Job * getJob(Job * prev, Thread * thread);
        
        friend class Thread;
        
        class Thread
        {
            friend class TaskManager;
            
        private:
            unsigned d_id;
            boost::thread * d_thread;
            
            static void runThread(TaskManager * tm, Thread * thread);
            
        public:
            bool d_stop;
            
            Thread(unsigned id, TaskManager * tm);
            ~Thread();
        };
        
        std::deque<Thread *> d_threads;
        
        class RAIIAdding
        {
        private:
            TaskManager & d_tm;
            
        public:
            RAIIAdding(TaskManager & tm);
            ~RAIIAdding();
        };
        
    public:
        TaskManager(unsigned nothreads);
        ~TaskManager();
        
        static unsigned getNumberOfCores()
        {
            unsigned r = boost::thread::hardware_concurrency();
            if (r == 0) // always return something positive!
                r = 1;
            return r;
        }
        
        unsigned getNumberOfThreads() const
        {
            return d_nothreads;
        }
        
        // Add jobs
        void enqueueJob(Job * j)
        {
            RAIIAdding a(*this);
            d_jobs.push_back(j);
        }
        
        void enqueueJobs(Job ** j, unsigned n)
        {
            if (n == 0)
                return;
            RAIIAdding a(*this);
            for (unsigned i = 0; i < n; ++i)
                d_jobs.push_back(j[i]);
        }
        
        template<class It>
        void enqueueJobs(It begin, It end)
        {
            if (begin == end)
                return;
            RAIIAdding a(*this);
            d_jobs.insert(d_jobs.end(), begin, end);
        }
        
        void clearAllJobs()
        {
            RAIIAdding a(*this);
            for (std::list<Job *>::iterator i = d_jobs.begin(); i != d_jobs.end(); ++i)
            {
                delete *i;
                *i = NULL;
            }
            d_jobs.clear();
        }
        
        // Check if some tasks are waiting for jobs
        bool isPartiallyIdling() const;
        
        // How many jobs are currently running?
        unsigned howManyAreRunning() const;
        unsigned howManyJobsAreEnqueued() const;
        bool areJobsEnqueued() const;
        std::pair<unsigned, unsigned> getJobsInfo() const; // first will contain how many jobs are running, second how many jobs are in the queue
        
        // Check if done
        bool isDone() const;
        
        // Wait until everything is done
        void waitForDone();
        
        // Wait until at least one thread is done
        void waitForDone_Single();
    };
}

#endif

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
#include "taskmanager.hpp"
#include "myalloc.hpp"
#include <plll/arithmetic.hpp>
#include <iostream>

namespace plll
{
    void TaskManager::Thread::runThread(TaskManager * tm, Thread * thread)
    {
        initTLAlloc();
        TaskManager::Job * job = NULL;
        while (true)
        {
            job = tm->getJob(job, thread);
            if (job == NULL)
                break;
            job->run(thread->d_id, tm);
        }
    }
    
    TaskManager::Thread::Thread(unsigned id, TaskManager * tm)
        : d_id(id), d_thread(new boost::thread(runThread, tm, this)), d_stop(false)
    {
    }
    
    TaskManager::Thread::~Thread()
    {
        d_thread->join();
    }
    
    TaskManager::Job * TaskManager::getJob(TaskManager::Job * previous, TaskManager::Thread * thread)
    {
        bool dec = false;
        if (previous)
        {
            delete previous;
            dec = true;
        }
        boost::unique_lock<boost::shared_mutex> lock(d_jobs_mutex);
        // Wait until something interesting happens
        while (d_jobs.empty() && !thread->d_stop)
        {
            if (dec)
            {
                --d_running;
                {
                    boost::unique_lock<boost::mutex> lock(d_waitfordone_single_mutex);
                    d_done_single = true;
                }
                d_waitfordone_single_cond.notify_all();
                dec = false;
            }
            if (d_running == 0)
            {
                // Inform that we are done
                d_waitfordone_cond.notify_all();
            }
            d_jobs_cond.wait(lock);
        }
        
        if (dec)
        {
            {
                boost::unique_lock<boost::mutex> lock(d_waitfordone_single_mutex);
                d_done_single = true;
            }
            d_waitfordone_single_cond.notify_all();
        }
        
        // Are we done?
        if (thread->d_stop)
        {
            if (dec)
                --d_running;
            return NULL;
        }
        
        // Get job
        Job * res = d_jobs.front();
        d_jobs.pop_front();
        if (!dec)
            ++d_running;
        return res;
    }
    
    TaskManager::RAIIAdding::RAIIAdding(TaskManager & tm)
        : d_tm(tm)
    {
        d_tm.d_jobs_mutex.lock();
    }
    
    TaskManager::RAIIAdding::~RAIIAdding()
    {
        unsigned c = d_tm.d_jobs.size();
        d_tm.d_jobs_mutex.unlock();
        if (c == 1)
            d_tm.d_jobs_cond.notify_one();
        else
            d_tm.d_jobs_cond.notify_all();
    }
    
    TaskManager::TaskManager(unsigned nothreads)
        : d_nothreads(nothreads), d_running(0), d_done_single(false)
    {
        if (d_nothreads == 0)
            d_nothreads = getNumberOfCores();
        
        std::cout << "(starting " << d_nothreads << " threads)\n";
        try
        {
            for (unsigned i = 0; i < d_nothreads; ++i)
                d_threads.push_back(new Thread(i, this));
        }
        catch(...)
        {
            try
            {
                {
                    boost::lock_guard<boost::shared_mutex> guard(d_jobs_mutex);
                    for (std::deque<Thread *>::iterator i = d_threads.begin(); i != d_threads.end(); ++i)
                        (*i)->d_stop = true;
                }
                d_jobs_cond.notify_all();
                for (std::deque<Thread *>::iterator i = d_threads.begin(); i != d_threads.end(); ++i)
                {
                    delete *i;
                    *i = NULL;
                }
            }
            catch (...)
            {
            }
            throw;
        }
    }
    
    TaskManager::~TaskManager()
    {
        std::cout << "(stopping threads)\n";
        try
        {
            {
                boost::lock_guard<boost::shared_mutex> guard(d_jobs_mutex);
                for (std::deque<Thread *>::iterator i = d_threads.begin(); i != d_threads.end(); ++i)
                    (*i)->d_stop = true;
            }
            d_jobs_cond.notify_all();
        }
        catch (...)
        {
        }
        for (std::deque<Thread *>::iterator i = d_threads.begin(); i != d_threads.end(); ++i)
        {
            try
            {
                delete *i;
            }
            catch (...)
            {
            }
            *i = NULL;
        }
        d_threads.clear();
        std::cout << "(all threads have stopped)\n";
    }
    
    bool TaskManager::isPartiallyIdling() const
    // Check if some tasks are waiting for jobs
    {
        return d_running < d_nothreads;
    }
    
    unsigned TaskManager::howManyAreRunning() const
    {
        return d_running;
    }
    
    unsigned TaskManager::howManyJobsAreEnqueued() const
    {
        boost::shared_lock<boost::shared_mutex> lock(d_jobs_mutex);
        return d_jobs.size();
    }
    
    bool TaskManager::areJobsEnqueued() const
    {
        boost::shared_lock<boost::shared_mutex> lock(d_jobs_mutex);
        return !d_jobs.empty();
    }
    
    std::pair<unsigned, unsigned> TaskManager::getJobsInfo() const
    // first will contain how many jobs are running, second how many jobs are in the queue
    {
        boost::shared_lock<boost::shared_mutex> lock(d_jobs_mutex);
        return std::make_pair(d_running, d_jobs.size());
    }
    
    bool TaskManager::isDone() const
    // Check if done
    {
        boost::shared_lock<boost::shared_mutex> lock(d_jobs_mutex);
        return d_jobs.empty() && (d_running == 0);
    }
    
    void TaskManager::waitForDone()
    {
        boost::unique_lock<boost::mutex> lock(d_waitfordone_mutex);
        while ((d_running > 0) || !d_jobs.empty())
            d_waitfordone_cond.wait(lock);
    }
    
    void TaskManager::waitForDone_Single()
    {
        boost::unique_lock<boost::mutex> lock(d_waitfordone_single_mutex);
        while (!d_done_single || (d_running > 0))
            d_waitfordone_single_cond.wait(lock);
        d_done_single = false;
    }
}


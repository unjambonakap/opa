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

#ifndef PLLL_INCLUDE_GUARD__ENUMIMPL_PARALLELSERIAL_CPP
#define PLLL_INCLUDE_GUARD__ENUMIMPL_PARALLELSERIAL_CPP

#include "enumimpl-parallelserial-enum.cpp"
#include "taskmanager.hpp"

namespace plll
{
    template<class RealTypeContext, class IntTypeContext>
    class ParallelEnumerator
    {
    public:
        typedef boost::function<void(const linalg::math_matrix<typename IntTypeContext::Integer> & basis,
                                     int p,
                                     const linalg::math_rowvector<typename IntTypeContext::Integer> & vec)> CallbackFunction;
        
    private:
        TaskManager d_tm;
        unsigned d_dim; // d_dim is the real dimension
        unsigned d_begin, d_end;
        Updater<RealTypeContext, IntTypeContext> d_update;
        Lattice<RealTypeContext, IntTypeContext> * d_lattice;
        
        class JobManager;
        
        class TMJob : public TaskManager::Job
        {
        private:
            ParallelEnumerator & d_pe;
            typename JobManager::JobDataT * d_job;
            
        public:
            TMJob(ParallelEnumerator & pe, typename JobManager::JobDataT * job)
                : d_pe(pe), d_job(job)
            {
            }
            
            virtual ~TMJob()
            {
            }
            
            virtual void run(unsigned thread_no, TaskManager * tm);
        };
        
        class JobManager
        {
        public:
            typedef ThreadData<RealTypeContext, IntTypeContext, JobManager> ThreadDataT;
            typedef JobData<RealTypeContext, IntTypeContext, JobManager> JobDataT;
            
        private:
            std::list<JobDataT> d_jobs;
            ParallelEnumerator & d_pe;
            
            boost::mutex d_jobs_mutex;
            std::list<JobDataT*> d_jobs_to_run;
            
            class JobsMutexLock
            {
            private:
                JobManager & d_jm;
                
            public:
                inline JobsMutexLock(JobManager & jm)
                    : d_jm(jm)
                {
                    if (d_jm.d_parallel_mode) d_jm.d_jobs_mutex.lock();
                }
                
                inline ~JobsMutexLock()
                {
                    if (d_jm.d_parallel_mode) d_jm.d_jobs_mutex.unlock();
                }
            };
            
            linalg::math_rowvector<typename IntTypeContext::Integer> d_result; // is zero vector (of no length) until a solution was found
            volatile bool d_has_bound;
            typename RealTypeContext::Real d_bound;
            volatile bool d_bound_updated;
            CallbackFunction d_ecf;
            
            bool d_parallel_mode;
            
            long double computeHeuristic(const typename JobManager::JobDataT & job, const typename RealTypeContext::Real & bound)
            // Computes the Gauss volume heuristic for the given job
            {
                if (job.dimension() == 0)
                    return 0.0;
                unsigned stage = d_pe.d_dim - job.dimension();
                long double v = std::pow(arithmetic::convert<long double>(bound - job.ell()), 0.5 * stage);
                v *= d_pe.d_gaussianfactors.HPF(stage);
                long double d = arithmetic::convert<double>(d_pe.d_lattice->getNormSq(d_pe.d_begin));
                for (unsigned i = 1; i <= stage; ++i)
                    d *= arithmetic::convert<double>(d_pe.d_lattice->getNormSq(d_pe.d_begin + i));
                return v * std::pow(d, -0.5l / (long double)stage);
            }
            
            void updateBound(const linalg::math_rowvector<typename IntTypeContext::Integer> & result, const typename RealTypeContext::Real bound)
            {
                bool did_update_here = false;
                bool hasSolutionYet = !d_has_bound;
                if (d_pe.d_update(*d_pe.d_lattice, d_pe.d_begin, d_pe.d_end, d_result, d_bound, result, bound, hasSolutionYet))
                {
                    d_bound_updated = true;
                    did_update_here = true;
                    // Call EnumCallbackFunction!
                    if (!d_ecf.empty())
                    {
                        std::pair<linalg::math_matrix<typename IntTypeContext::Integer> *, unsigned> l = d_pe.d_lattice->getMatrixIndex(d_pe.d_begin);
                        if (l.first)
                            d_ecf(*l.first, l.second, result);
                    }
                }
                d_has_bound = !hasSolutionYet;
                if (did_update_here)
                    d_pe.updateBound(d_bound);
            }
            
            void addResult(JobDataT * data)
            {
                // Mark as done
                data->markDone();
                // Check out result
                if (data->d_result.size())
                    updateBound(data->d_result, data->d_bound);
                // All jobs done?
                if (!d_pe.d_tm.areJobsEnqueued())
                    d_pe.split();
            }
            
            JobDataT * getJobImpl()
            {
                JobDataT * res = NULL;
                while (res == NULL)
                {
                    if (d_jobs_to_run.empty())
                        break;
                    res = d_jobs_to_run.front();
                    d_jobs_to_run.pop_front();
                    if (d_has_bound)
                        if (!res->update(d_bound))
                            res = NULL; // job not needed anymore
                }
                return res;
            }
            
        public:
            PruningHelper<RealTypeContext, IntTypeContext> d_pruninghelper;
            
            JobManager(ParallelEnumerator & pe)
                : d_pe(pe), d_has_bound(false), d_bound_updated(false), d_ecf(NULL), d_parallel_mode(false), d_pruninghelper()
            {
            }
            
            void setCallback(CallbackFunction ecf)
            {
                d_ecf = ecf;
            }
            
            void restart(const typename RealTypeContext::Real & bound)
            {
                d_jobs.clear();
                d_jobs_to_run.clear();
                d_has_bound = false;
                d_bound = bound;
                d_result.resize(0);
                d_bound_updated = false;
            }
            
            class HSorter
            {
            public:
                bool operator() (const std::pair<JobDataT*, long double> & p1, const std::pair<JobDataT*, long double> & p2) const
                {
                    // > seems to minimize total CPU time, but to "maximize" total wallclock time
                    // < seems to "maximize" total CPU time, but to minimize total wallclock time
                    return p1.second > p2.second;
                }
            };
            
            void optimize()
            {
                JobsMutexLock jml(*this);
                
                // Store queue into std::vector including the heuristic's values
                std::vector<std::pair<JobDataT*, long double> > tmp;
                tmp.reserve(d_jobs_to_run.size());
                for (typename std::list<JobDataT*>::iterator i = d_jobs_to_run.begin(); i != d_jobs_to_run.end(); ++i)
                    tmp.push_back(std::make_pair(*i, computeHeuristic(**i, d_bound)));
                d_jobs_to_run.clear();
                
                // Sort
                HSorter sorter;
                std::sort(tmp.begin(), tmp.end(), sorter);
                
                // Add jobs back to queue
                for (unsigned i = 0; i < tmp.size(); ++i)
                    d_jobs_to_run.push_back(tmp[i].first);
            }
            
            // For scheduler
            
            void addJob(unsigned dim,
                        const linalg::math_rowvector<typename IntTypeContext::Integer> & x_end,
                        const typename RealTypeContext::Real & ell,
                        const typename RealTypeContext::Real & bound,
                        bool iszero)
            {
                JobsMutexLock jml(*this);
                d_jobs.push_back(JobDataT(dim, x_end, ell, d_has_bound ? (bound < d_bound ? bound : d_bound) : bound, iszero));
                if (d_parallel_mode)
                    d_pe.d_tm.enqueueJob(new TMJob(d_pe, &d_jobs.back()));
                else
                    d_jobs_to_run.push_back(&d_jobs.back());
            }
            
            template<class It>
            void addJobs(It begin, It end)
            {
                if (begin == end)
                    return;
                JobsMutexLock jml(*this);
                while (begin != end)
                {
                    bool skip = false;
                    if (d_has_bound)
                        if (!begin->update(d_bound))
                            skip = true;
                    if (!skip)
                    {
                        d_jobs.push_back(*begin);
                        if (d_parallel_mode)
                            d_pe.d_tm.enqueueJob(new TMJob(d_pe, &d_jobs.back()));
                        else
                            d_jobs_to_run.push_back(&d_jobs.back());
                    }
                    ++begin;
                }
            }
            
            void launchJobsInParallel()
            {
                boost::unique_lock<boost::mutex> lock(d_jobs_mutex);
                d_parallel_mode = true;
                for (typename std::list<JobDataT*>::iterator i = d_jobs_to_run.begin(); i != d_jobs_to_run.end(); ++i)
                    d_pe.d_tm.enqueueJob(new TMJob(d_pe, *i));
                d_jobs_to_run.clear();
            }
            
            void stopLaunchingJobsInParallel()
            {
                boost::unique_lock<boost::mutex> lock(d_jobs_mutex);
                d_parallel_mode = false;
            }
            
            void tryToUpdateBound(const linalg::math_rowvector<typename IntTypeContext::Integer> & result, const typename RealTypeContext::Real bound)
            {
                JobsMutexLock jml(*this);
                updateBound(result, bound);
            }
            
            bool isBoundUpdated() const
            {
                return d_bound_updated;
            }
            
            typename RealTypeContext::Real getUpdatedBound()
            {
                JobsMutexLock jml(*this);
                typename RealTypeContext::Real bound(d_bound);
                d_bound_updated = false;
                return bound;
            }
            
            const linalg::math_rowvector<typename IntTypeContext::Integer> & getResult() const
            {
                return d_result;
            }
            
            const typename RealTypeContext::Real & getBound() const
            {
                return d_bound;
            }
            
            // For JobRunner:
            
            bool getJob(JobDataT * job)
            // Returns false if no job is left
            {
                if (d_has_bound)
                {
                    boost::unique_lock<boost::mutex> lock(d_jobs_mutex);
                    if (!job->update(d_bound))
                        return false;
                }
                return true;
            }
            
            JobDataT * getJob()
            // Returns NULL if no job is left
            {
                if (d_parallel_mode)
                {
                    boost::unique_lock<boost::mutex> lock(d_jobs_mutex);
                    return getJobImpl();
                }
                else
                    return getJobImpl();
            }
            
            void sendResult(JobDataT * data)
            // Sends back result
            {
                if (d_parallel_mode)
                {
                    boost::unique_lock<boost::mutex> lock(d_jobs_mutex);
                    addResult(data);
                }
                else
                    addResult(data);
            }
            
            // Stats
            
            unsigned getNoToDo()
            {
                boost::unique_lock<boost::mutex> lock(d_jobs_mutex);
                return d_jobs_to_run.size();
            }
            
            void printStats()
            {
                std::pair<unsigned, unsigned> j = d_pe.d_tm.getJobsInfo();
                unsigned jsize = d_jobs.size();
                arithmetic::RealContext rc;
                arithmetic::Real bound = arithmetic::convert(d_bound, rc);
                d_pe.d_verbose(LatticeReduction::VL_Information) << jsize << " jobs, " << j.second << " to do, " << j.first << " running right now (" << d_pe.d_end - d_pe.d_begin + 1 << "-dimensional enumerate, indices " << d_pe.d_begin << "--" << d_pe.d_end << "); bound is " << bound;
            }
            
            // Accessors for JobDataT
            
            static const linalg::math_rowvector<typename IntTypeContext::Integer> & jobGetResult(const JobDataT & td)
            {
                return td.d_result;
            }
            
            static const typename RealTypeContext::Real & jobGetBound(const JobDataT & td)
            {
                return td.d_bound;
            }
        };
        
        JobManager d_jm;
        
        std::vector<typename JobManager::ThreadDataT *> d_thread_data;
        
        void setup(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned begin, unsigned end)
        {
            for (typename std::vector<typename JobManager::ThreadDataT *>::iterator i = d_thread_data.begin(); i != d_thread_data.end(); ++i)
                (*i)->flagToSetup(lattice, begin, end);
        }
        
        void split()
        {
            for (typename std::vector<typename JobManager::ThreadDataT *>::iterator i = d_thread_data.begin(); i != d_thread_data.end(); ++i)
                (*i)->flagToSplit();
        }
        
        void stop()
        {
            d_tm.clearAllJobs();
            for (typename std::vector<typename JobManager::ThreadDataT *>::iterator i = d_thread_data.begin(); i != d_thread_data.end(); ++i)
                (*i)->flagToStop();
        }
        
        void updateBound(const typename RealTypeContext::Real & bound)
        {
            for (typename std::vector<typename JobManager::ThreadDataT *>::iterator i = d_thread_data.begin(); i != d_thread_data.end(); ++i)
                (*i)->setUpdatedBound(bound);
        }
        
        volatile bool d_status_active, d_status_quit;
        boost::thread d_status_thread;
        
        static void statusThread(ParallelEnumerator * pe)
        {
            initTLAlloc(); // This might not seem necessary, but as soon as a temporary Real and/or
            // Integer object is created while outputting something, it is necessary.
//        boost::posix_time::milliseconds sleepTime(500);
            boost::posix_time::milliseconds sleepTime(5000);
            while (!pe->d_status_quit)
            {
                if (pe->d_status_active)
                    pe->d_jm.printStats();
                boost::this_thread::sleep(sleepTime);
            }
        }
        
        typename JobManager::ThreadDataT d_main_threaddata;
        StandardEnumSettings::PruningMethod d_pruning;
        double d_pruning_parameter; // needed for PM_Step, PM_Piecewise and PM_SchnorrHoerner2
        
        Verbose & d_verbose;
        GaussianFactorComputer & d_gaussianfactors;
        
        enum StopType { ST_Normal, ST_StopEnum, ST_StopReduction, ST_BadAlloc, ST_Unknown };
        
        StopType d_threads_stop;
        
    public:
        ParallelEnumerator(Verbose & v, GaussianFactorComputer & gf, unsigned enumdimension, unsigned nothreads)
            : d_tm(nothreads), d_dim(0), d_begin(0), d_end(0), d_update(), d_lattice(NULL),
              d_jm(*this), d_status_active(false), d_status_quit(false), d_status_thread(statusThread, this),
              d_main_threaddata(d_jm, d_update),
              d_pruning(StandardEnumSettings::PM_None), d_pruning_parameter(0),
              d_verbose(v), d_gaussianfactors(gf), d_threads_stop(ST_Normal)
        {
            d_thread_data.reserve(d_tm.getNumberOfThreads());
            for (unsigned i = 0; i < d_tm.getNumberOfThreads(); ++i)
                d_thread_data.push_back(new typename JobManager::ThreadDataT(d_jm, d_update));
        }
        
        ~ParallelEnumerator()
        {
            d_status_quit = true;
            d_status_thread.interrupt();
            d_status_thread.join();
            for (unsigned i = 0; i < d_thread_data.size(); ++i)
                delete d_thread_data[i];
        }
        
        bool enumerate(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned begin, unsigned end,
                       linalg::math_rowvector<typename IntTypeContext::Integer> & result,
                       typename RealTypeContext::Real & bound, const StandardEnumSettings & s)
        // Finds a shortest vector in the lattice generated by the orthogonal projections of the vectors
        // A.row(begin) to A.row(end) into the orthogonal complement of the vectors A.row(0) to
        // A.row(begin-1). Uses the Kannan-Schnorr-Euchner enumeration method.
        {
            d_pruning = s.d_pruning;
            d_pruning_parameter = s.d_pruning_parameter;
            
            // Initialize updater
            d_lattice = &lattice;
            d_update.initialize(d_verbose, lattice, begin, end, bound);
            d_jm.d_pruninghelper.setup(begin, end - begin + 1, lattice, d_pruning, d_pruning_parameter);
            
            d_dim = end - begin + 1;
            d_begin = begin;
            d_end = end;
            
            bool verbose = false;
            
            // Restart job manager
            d_jm.restart(bound);
            d_threads_stop = ST_Normal;
            
            // Insert beginning stub
            linalg::math_rowvector<typename IntTypeContext::Integer> x;
            typename RealTypeContext::Real ell(lattice.rc());
            setZero(ell);
            d_jm.addJob(0, x, ell, bound, true);
            
            // Refine queue
            d_main_threaddata.setup(lattice, d_begin, d_end);
            typename JobManager::JobDataT * job = d_jm.getJob();
            unsigned last_loD = d_dim;
            while (job)
            {
                // Run the job
                unsigned loD = job->leftoverDims(d_main_threaddata);
                if (loD < last_loD)
                {
                    // Sort
                    d_jm.optimize();
                    // Stats
                    if (verbose)
                        d_jm.printStats();
                    last_loD = loD;
                }
                if (loD < 15)
                    job->run(d_main_threaddata);
                else
                    job->refine(d_main_threaddata, 1);
                // Return job
                d_jm.sendResult(job);
                // Is there enough?
                if (d_jm.getNoToDo() >= d_tm.getNumberOfThreads() * 75)
                    break;
                // Get next (if available)
                job = d_jm.getJob();
            }
            
            if (d_jm.getNoToDo() > 0) // only if there's something to do!
            {
                verbose = true;
                if (verbose)
                {
                    d_verbose(LatticeReduction::VL_Information) << "(after refining)";
                    d_jm.printStats();
                }
                if (verbose)
                    d_verbose(LatticeReduction::VL_Information) << "(optimizing queue)";
                d_jm.optimize();
                setup(lattice, d_begin, d_end);
                if (verbose)
                    d_verbose(LatticeReduction::VL_Information) << "(continuing threads)";
                d_jm.launchJobsInParallel();
                d_status_active = true;
                d_tm.waitForDone();
                d_jm.stopLaunchingJobsInParallel();
                d_status_active = false;
                if (verbose)
                    d_verbose(LatticeReduction::VL_Information) << "(all threads done)";
                if (d_threads_stop != ST_Normal)
                    switch (d_threads_stop)
                    {
                    case ST_StopEnum: break;
                    case ST_StopReduction: throw LatticeReduction::stop_reduction();
                    case ST_BadAlloc: throw std::bad_alloc();
                    default:
                    case ST_Unknown: throw std::runtime_error("Unknown exception caught in enumeration thread");
                    }
            }
            
            // Fetch solutions
            if (d_jm.getResult().size() > 0)
            {
                for (unsigned i = 0; i < result.size(); ++i)
                    result[i] = d_jm.getResult()[i];
                bound = d_jm.getBound();
                if (verbose)
                    d_verbose(LatticeReduction::VL_Information) << "[done]";
                return true;
            }
            else
            {
                if (verbose)
                {
                    d_verbose(LatticeReduction::VL_Information) << "(no solution found!)";
                    d_verbose(LatticeReduction::VL_Information) << "[done]";
                }
                return false;
            }
        }
        
        void setCallback(CallbackFunction cf)
        {
            d_jm.setCallback(cf);
        }
    };
    
    template<class RealTypeContext, class IntTypeContext>
    class SerialEnumerator
    {
    public:
        typedef boost::function<void(const linalg::math_matrix<typename IntTypeContext::Integer> & basis,
                                     int p,
                                     const linalg::math_rowvector<typename IntTypeContext::Integer> & vec)> CallbackFunction;
        
    private:
        Lattice<RealTypeContext, IntTypeContext> * d_lattice;
        Updater<RealTypeContext, IntTypeContext> d_update;
        unsigned d_begin;
        
        class JobManager;
        
        typedef ThreadData<RealTypeContext, IntTypeContext, JobManager> ThreadDataT;
        typedef JobData<RealTypeContext, IntTypeContext, JobManager> JobDataT;
        
        class JobManager
        {
        private:
            SerialEnumerator & d_se;
            CallbackFunction d_ecf;
            
        public:
            PruningHelper<RealTypeContext, IntTypeContext> d_pruninghelper;
            
            JobManager(SerialEnumerator & se)
                : d_se(se), d_ecf(NULL), d_pruninghelper()
            {
            }
            
            static const linalg::math_rowvector<typename IntTypeContext::Integer> & jobGetResult(const JobDataT & td)
            {
                return td.d_result;
            }
            
            static const typename RealTypeContext::Real & jobGetBound(const JobDataT & td)
            {
                return td.d_bound;
            }
            
            template<class It>
            static inline void addJobs(It, It)
            // Dummy function which just ignores the given data
            {
            }
            
            void tryToUpdateBound(const linalg::math_rowvector<typename IntTypeContext::Integer> & result, const typename RealTypeContext::Real bound)
            {
                if (!d_ecf.empty())
                {
                    std::pair<linalg::math_matrix<typename IntTypeContext::Integer> *, unsigned> l = d_se.d_lattice->getMatrixIndex(d_se.d_begin);
                    if (l.first)
                        d_ecf(*l.first, l.second, result);
                }
            }
            
            inline bool hasCallback() const
            {
                return !d_ecf.empty();
            }
            
            inline void setCallback(CallbackFunction ecf)
            {
                d_ecf = ecf;
            }
        };
        
        JobManager d_jm;
        ThreadDataT d_threaddata;
        StandardEnumSettings::PruningMethod d_pruning;
        double d_pruning_parameter; // needed for PM_Step, PM_Piecewise and PM_SchnorrHoerner2
        
        Verbose & d_verbose;
        GaussianFactorComputer & d_gaussianfactors;
        
    public:
        SerialEnumerator(Verbose & v, GaussianFactorComputer & gf, unsigned enumdimension)
            : d_lattice(NULL), d_update(), d_begin(0), d_jm(*this), d_threaddata(d_jm, d_update),
              d_pruning(StandardEnumSettings::PM_None), d_pruning_parameter(0),
              d_verbose(v), d_gaussianfactors(gf)
        {
        }
        
        bool enumerate(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned begin, unsigned end,
                       linalg::math_rowvector<typename IntTypeContext::Integer> & result,
                       typename RealTypeContext::Real & bound, const StandardEnumSettings & s)
        // Finds a shortest vector in the lattice generated by the orthogonal projections of the vectors
        // A.row(begin) to A.row(end) into the orthogonal complement of the vectors A.row(0) to
        // A.row(begin-1). Uses the Kannan-Schnorr-Euchner enumeration method.
        {
            d_pruning = s.d_pruning;
            d_pruning_parameter = s.d_pruning_parameter;
            d_begin = begin;
            
            // Initialize updater
            d_lattice = &lattice;
            d_update.initialize(d_verbose, lattice, begin, end, bound);
            d_jm.d_pruninghelper.setup(begin, end - begin + 1, lattice, d_pruning, d_pruning_parameter);
            
            // Beginning stub
            linalg::math_rowvector<typename IntTypeContext::Integer> x;
            typename RealTypeContext::Real ell(lattice.rc());
            setZero(ell);
            JobDataT job(0, x, ell, bound, true);
            // Process
            d_threaddata.setup(lattice, begin, end);
            if (d_jm.hasCallback())
            {
                try
                {
                    // Use two-stage enum
                    job.run2(d_threaddata);
                }
                catch (LatticeReduction::stop_enumeration &)
                {
                    d_verbose(LatticeReduction::VL_Chatter) << "Stopping enumeration";
                }
            }
            else
                // Do just one stage which covers everything
                job.run(d_threaddata);
            
            // Fetch solution
            if (d_jm.jobGetResult(job).size() > 0)
            {
                for (unsigned i = 0; i < d_jm.jobGetResult(job).size(); ++i)
                    result[i] = d_jm.jobGetResult(job)[i];
                bound = d_jm.jobGetBound(job);
                return true;
            }
            else
                return false;
        }
        
        void setCallback(CallbackFunction cf)
        {
            d_jm.setCallback(cf);
        }
    };
    
    template<class RealTypeContext, class IntTypeContext>
    class ParallelSerialEnumerator
    {
    private:
        bool d_parallel;
        STD_AUTO_PTR<ParallelEnumerator<RealTypeContext, IntTypeContext> > d_enum_parallel;
        STD_AUTO_PTR<SerialEnumerator<RealTypeContext, IntTypeContext> > d_enum_serial;
        
        // Disable copying and copy-constructing
        ParallelSerialEnumerator(const ParallelSerialEnumerator &);
        ParallelSerialEnumerator & operator = (const ParallelSerialEnumerator &);
        
        typedef boost::function<void(const linalg::math_matrix<typename IntTypeContext::Integer> & basis,
                                     int p,
                                     const linalg::math_rowvector<typename IntTypeContext::Integer> & vec)> CallbackFunction;
        
    public:
        ParallelSerialEnumerator(Verbose & v, GaussianFactorComputer & gf, unsigned enumdimension, unsigned max_threads, LatticeReduction::Statistics & stats)
            : d_enum_parallel(), d_enum_serial()
        {
            unsigned nothreads = boost::thread::hardware_concurrency();
            if (max_threads)
            {
                if (nothreads > max_threads)
                    nothreads = max_threads;
            }
            if (nothreads == 1)
            {
                d_parallel = false;
                d_enum_serial.reset(new SerialEnumerator<RealTypeContext, IntTypeContext>(v, gf, enumdimension));
            }
            else
            {
                d_parallel = true;
                d_enum_parallel.reset(new ParallelEnumerator<RealTypeContext, IntTypeContext>(v, gf, enumdimension, nothreads));
            }
        }
        
        bool enumerate(unsigned begin, unsigned end, linalg::math_rowvector<typename IntTypeContext::Integer> & result,
                       typename RealTypeContext::Real & bound, const StandardEnumSettings & s)
        // Finds a shortest vector in the lattice generated by the orthogonal projections of the vectors
        // A.row(begin) to A.row(end) into the orthogonal complement of the vectors A.row(0) to
        // A.row(begin-1). Uses the Kannan-Schnorr-Euchner enumeration method.
        {
            return d_parallel ? d_enum_parallel->enumerate(begin, end, result, bound, s)
                : d_enum_serial->enumerate(begin, end, result, bound, s);
        }
        
        void setCallback(CallbackFunction cf)
        {
            if (d_parallel)
                d_enum_parallel->setCallback(cf);
            else
                d_enum_serial->setCallback(cf);
        }
    };
    
    template<class RealTypeContext, class IntTypeContext>
    void ParallelEnumerator<RealTypeContext, IntTypeContext>::TMJob::run(unsigned thread_no, TaskManager * tm)
    {
        try
        {
            if (d_pe.d_jm.getJob(d_job))
            {
                // Run the job
                d_job->run2(*d_pe.d_thread_data[thread_no]);
                // Return
                d_pe.d_jm.sendResult(d_job);
            }
        }
        catch (LatticeReduction::stop_enumeration &)
        {
            d_pe.stop();
            d_pe.d_threads_stop = ST_StopEnum;
            d_pe.d_verbose(LatticeReduction::VL_Chatter) << "Stopping enumeration";
        }
        catch (LatticeReduction::stop_reduction &)
        {
            d_pe.stop();
            d_pe.d_threads_stop = ST_StopReduction;
            d_pe.d_verbose(LatticeReduction::VL_Chatter) << "Stopping enumeration (stop reduction)";
        }
        catch (std::bad_alloc &)
        {
            d_pe.stop();
            d_pe.d_threads_stop = ST_BadAlloc;
            d_pe.d_verbose(LatticeReduction::VL_Chatter) << "Stopping enumeration (bad alloc)";
        }
        catch(std::exception & e)
        {
            d_pe.stop();
            d_pe.d_threads_stop = ST_Unknown;
            d_pe.d_verbose(LatticeReduction::VL_Chatter) << "Stopping enumeration (unknown exception: " << e.what() << ")";
        }
        catch(...)
        {
            d_pe.stop();
            d_pe.d_threads_stop = ST_Unknown;
            d_pe.d_verbose(LatticeReduction::VL_Chatter) << "Stopping enumeration (unknown exception)";
        }
    }
}

#endif

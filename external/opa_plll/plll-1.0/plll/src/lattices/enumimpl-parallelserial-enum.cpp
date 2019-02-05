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

#ifndef PLLL_INCLUDE_GUARD__ENUMIMPL_PARALLELSERIAL_ENUM_CPP
#define PLLL_INCLUDE_GUARD__ENUMIMPL_PARALLELSERIAL_ENUM_CPP

#include <boost/thread.hpp>

namespace plll
{
    template<class RealTypeContext, class IntTypeContext, class JobManager>
    /*
      The JobManager class must have a method
          void addJob(unsigned dim,
                      const linalg::math_rowvector<typename IntTypeContext::Integer> & x_end,
                      const typename RealTypeContext::Real & ell,
                      const typename RealTypeContext::Real & bound,
                      bool iszero);
    */
    class JobData;
    
    template<class RealTypeContext, class IntTypeContext, class JobManager>
    class ThreadData
    // Stores per thread data, like certain arrays
    {
        friend class JobData<RealTypeContext, IntTypeContext, JobManager>;
        
    private:
        Lattice<RealTypeContext, IntTypeContext> * d_lattice;
        JobManager & d_jm;
        Updater<RealTypeContext, IntTypeContext> & d_update;
        unsigned d_begin;
        unsigned d_end;
        unsigned d_dim;
        
        volatile bool d_need_to_setup;
        volatile unsigned d_need_to_setup_begin;
        volatile unsigned d_need_to_setup_end;
        volatile Lattice<RealTypeContext, IntTypeContext> * d_need_to_setup_lattice;
        
        linalg::math_rowvector<typename IntTypeContext::Integer> d_result;
        linalg::math_rowvector<typename IntTypeContext::Integer> d_x;
        linalg::math_rowvector<typename RealTypeContext::Real> d_x_real;
        linalg::base_rowvector<long> d_delta;
        linalg::base_rowvector<int> d_delta2;
        linalg::math_rowvector<int> d_r; /* Used for improvement sketched in Appendix B of Gama,
        Nguyen, Regev: "Lattice Enumeration Using Extreme Pruning", EUROCRYPT 2010. Stores
        information on how much of the cache in d_sigma can be used. */
        linalg::math_matrix<typename RealTypeContext::Real> d_sigma; /* used for improvement
        sketched in Appendix B of Gama, Nguyen, Regev: "Lattice Enumeration Using Extreme Pruning",
        EUROCRYPT 2010. The i-th row caches intermediate values used during the computation of
        d_c[i]. */
        linalg::math_rowvector<typename RealTypeContext::Real> d_c;
        linalg::math_rowvector<typename RealTypeContext::Real> d_ell;
        linalg::math_rowvector<typename RealTypeContext::Real> d_bounds; /* For pruning. This is
        updated by JobData, not by ThreadData, and any relation to d_bound is managed by JobData. */
        std::vector<bool> d_iszero;
        
        std::list<JobData<RealTypeContext, IntTypeContext, JobManager> > d_job_list;
        
        boost::mutex d_bound_update_mutex;
        volatile bool d_bound_updated;
        typename RealTypeContext::Real d_bound;
        
        volatile bool d_split;
        volatile bool d_stop;
        
    public:
        ThreadData(JobManager & jm, Updater<RealTypeContext, IntTypeContext> & update)
            : d_lattice(NULL), d_jm(jm), d_update(update), d_begin(0), d_end(0), d_dim(0),
              d_need_to_setup(false), d_need_to_setup_begin(0), d_need_to_setup_end(0),
              d_need_to_setup_lattice(NULL), d_bound_updated(false), d_split(false),
              d_stop(false)
        {
        }
        
        void setup(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned begin, unsigned end)
        // Called at beginning of specific enumerate
        {
            d_need_to_setup = false;
            d_lattice = &lattice;
            d_begin = begin;
            d_end = end;
            d_dim = d_end - d_begin + 1;
            // The updater functions rely on the size of d_result and d_x, hence always resize them!
            d_result.resize(d_dim);
            d_x.resize(d_dim);
            // Only resize the others when necessary
            if (d_x_real.size() < d_dim)
            {
                d_x_real.resize(d_dim, linalg::Initialize(lattice.rc()));
                d_delta.resize(d_dim);
                d_delta2.resize(d_dim);
                d_r.resize(d_dim + 1);
                d_sigma.resize(d_dim, d_dim, linalg::Initialize(lattice.rc()));
                d_c.resize(d_dim, linalg::Initialize(lattice.rc()));
                d_ell.resize(d_dim + 1, linalg::Initialize(lattice.rc()));
                d_bounds.resize(d_dim, linalg::Initialize(lattice.rc()));
                d_iszero.resize(d_dim + 1);
            }
            for (unsigned i = 0; i < d_dim; ++i)
            {
                setZero(d_x_real[i]);
                setZero(d_c[i]);
                setZero(d_ell[i]);
                setZero(d_bounds[i]);
//            for (unsigned j = 0; j < d_dim; ++j)
//                setZero(d_sigma(i, j));
            }
            setZero(d_ell[d_dim]);
            d_bound_updated = false;
            setZero(d_bound);
            d_split = false;
            d_stop = false;
            updateContext();
        }
        
        void updateContext()
        {
            RealTypeContext & rc = d_lattice->rc();
            for (unsigned i = 0; i < d_sigma.rows(); ++i)
                for (unsigned j = 0; j < d_sigma.cols(); ++j)
                    d_sigma(i, j).setContext(rc);
            for (unsigned i = 0; i < d_x_real.size(); ++i)
                d_x_real[i].setContext(rc);
            for (unsigned i = 0; i < d_c.size(); ++i)
                d_c[i].setContext(rc);
            for (unsigned i = 0; i < d_ell.size(); ++i)
                d_ell[i].setContext(rc);
            for (unsigned i = 0; i < d_bounds.size(); ++i)
                d_bounds[i].setContext(rc);
            d_bound.setContext(rc);
        }
        
        inline void flagToSetup(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned begin, unsigned end)
        {
            d_need_to_setup_lattice = &lattice;
            d_need_to_setup_begin = begin;
            d_need_to_setup_end = end;
            d_need_to_setup = true;
        }
        
        inline void doSetup()
        {
            if (d_need_to_setup)
            {
                setup((Lattice<RealTypeContext, IntTypeContext>&)*d_need_to_setup_lattice, d_need_to_setup_begin, d_need_to_setup_end);
                d_need_to_setup = false;
            }
        }
        
        inline void setUpdatedBound(const typename RealTypeContext::Real & bound)
        {
            boost::lock_guard<boost::mutex> guard(d_bound_update_mutex);
            d_bound_updated = true;
            d_bound = bound;
        }
        
        inline bool hasUpdatedBound() const
        {
            return d_bound_updated;
        }
        
        inline bool updateBound(typename RealTypeContext::Real & bound)
        {
            bool updated = false;
            if (d_bound_updated)
            {
                boost::lock_guard<boost::mutex> guard(d_bound_update_mutex);
                if (d_bound_updated)
                {
                    d_bound_updated = false;
                    bound = d_bound;
                    updated = true;
                }
            }
            return updated;
        }
        
        inline void flagToStop()
        {
            d_stop = true;
        }
        
        inline bool shouldIStop()
        {
            return d_stop;
        }
        
        inline void flagToSplit()
        {
            d_split = true;
        }
        
        inline bool shouldISplit()
        {
            if (d_split)
            {
                d_split = false;
                return true;
            }
            else
                return false;
        }
    };
    
    template<class RealTypeContext, class IntTypeContext, class JobManager>
    class JobData
    // The job itself. Contains info on job parameters, local parameters while running, as well as contains result
    {
    private:
        linalg::math_rowvector<typename IntTypeContext::Integer> d_x_end; // td.d_dim - d_x_end.size() is the dimension which is
        // still to be enumerated
        typename RealTypeContext::Real d_ell;
        bool d_iszero;
        
    public:
        linalg::math_rowvector<typename IntTypeContext::Integer> d_result; // is zero vector (of no length) until a solution was found
        typename RealTypeContext::Real d_bound;
        
    private:
        void propagate_new_bound(ThreadData<RealTypeContext, IntTypeContext, JobManager> & td)
        {
            td.d_jm.d_pruninghelper.computeBounds(td.d_bounds, td.d_dim, d_bound);
        }
        
        inline void run_init(ThreadData<RealTypeContext, IntTypeContext, JobManager> & td, unsigned realdim, unsigned enumbegin, unsigned enumdim)
        {
            // Fill in x, delta, delta2, ell
            for (unsigned i = 0; i < d_x_end.size(); ++i)
            {
                td.d_x[enumdim + i] = d_x_end[i];
                arithmetic::convert(td.d_x_real[enumdim + i], d_x_end[i], td.d_lattice->rc());
            }
            td.d_ell[enumdim] = d_ell;
            td.d_iszero[enumdim] = d_iszero;
            
            // Prepare enumeration
            for (unsigned i = 0; i < enumdim; ++i)
            {
                td.d_r[i] = realdim - 1;
                setZero(td.d_sigma(i, realdim - 1));
            }
            typename RealTypeContext::Real tmp(td.d_lattice->rc());
            
            for (unsigned j = td.d_r[enumdim - 1]; j >= enumdim; --j) // this is the summation order considered in Pujol-Stehle
            {
                tmp = td.d_x_real[j] * td.d_lattice->getCoeff(enumbegin + j, enumbegin + enumdim - 1);
                td.d_sigma(enumdim - 1, j - 1) = td.d_sigma(enumdim - 1, j) + tmp;
            }
            td.d_c[enumdim - 1] = -td.d_sigma(enumdim - 1, enumdim - 1);
            bool ru;
            arithmetic::convert_round(td.d_x[enumdim - 1], td.d_c[enumdim - 1], ru, td.d_lattice->ic());
            arithmetic::convert(td.d_x_real[enumdim - 1], td.d_x[enumdim - 1], td.d_lattice->rc());
            td.d_iszero[enumdim - 1] = td.d_iszero[enumdim] && isZero(td.d_x[enumdim - 1]);
            td.d_delta[enumdim - 1] = 0;
            td.d_delta2[enumdim - 1] = ru ? 1 : -1;
            propagate_new_bound(td);
        }
        
        class CallUpdater
        {
        public:
            inline bool operator() (ThreadData<RealTypeContext, IntTypeContext, JobManager> & td, JobData & jd,
                                    unsigned realdim, unsigned enumbegin, unsigned enumdim, unsigned stage, bool & noSolutionYet,
                                    typename RealTypeContext::Real & tmp) const
            {
                if (td.d_update(*td.d_lattice, td.d_begin, td.d_end, jd.d_result, jd.d_bound, td.d_x, td.d_ell[0], noSolutionYet))
                    jd.propagate_new_bound(td);
                return true;
            }
        };
        
        class CallUpdaterNotify
        {
        private:
            mutable bool d_was_updated;
            
        public:
            inline CallUpdaterNotify()
                : d_was_updated(false)
            {
            }
            
            inline bool wasUpdated() const
            {
                return d_was_updated;
            }
            
            inline bool operator() (ThreadData<RealTypeContext, IntTypeContext, JobManager> & td, JobData & jd,
                                    unsigned realdim, unsigned enumbegin, unsigned enumdim, unsigned stage, bool & noSolutionYet,
                                    typename RealTypeContext::Real & tmp) const
            {
                if (td.d_update(*td.d_lattice, td.d_begin, td.d_end, jd.d_result, jd.d_bound, td.d_x, td.d_ell[0], noSolutionYet))
                {
                    jd.propagate_new_bound(td);
                    d_was_updated = true;
                }
                return true;
            }
        };
        
        class CallRefiner
        {
        public:
            inline bool operator() (ThreadData<RealTypeContext, IntTypeContext, JobManager> & td, JobData & jd,
                                    unsigned realdim, unsigned enumbegin, unsigned enumdim, unsigned stage, bool & noSolutionYet,
                                    typename RealTypeContext::Real & tmp) const
            {
                td.d_jm.addJob(realdim, td.d_x, td.d_ell[0], jd.d_bound, td.d_iszero[0]);
                // update bound (improves pruning)
                if (td.d_jm.isBoundUpdated())
                {
                    jd.d_bound = td.d_jm.getUpdatedBound();
                    jd.propagate_new_bound(td);
                }
                return true;
            }
        };
        
        class CallRunner;
        friend class CallRunner;
        
        class CallRunner
        {
        private:
            static void enqueue(ThreadData<RealTypeContext, IntTypeContext, JobManager> & td, JobData & jd,
                                unsigned realdim, unsigned enumbegin, unsigned enumdim, unsigned stage,
                                typename RealTypeContext::Real & tmp)
            {
                while (true)
                {
                    tmp = td.d_x_real[stage] - td.d_c[stage];
                    square(td.d_ell[stage], tmp);
                    td.d_ell[stage] *= td.d_lattice->getNormSq(enumbegin + stage);
                    td.d_ell[stage] += td.d_ell[stage + 1];
                    if (td.d_ell[stage] <= td.d_bounds[stage])
                        td.d_job_list.push_back(JobData(td.d_dim - stage, td.d_x, td.d_ell[stage], jd.d_bound, td.d_iszero[stage], stage));
                    else
                    {
                        if (++stage >= enumdim)
                            break;
                    }
                    // Modify d_x[stage]
                    if (td.d_iszero[stage + 1])
                    {
                        ++td.d_x[stage];
                    }
                    else
                    {
                        td.d_delta2[stage] = -td.d_delta2[stage];
                        td.d_delta[stage] = -td.d_delta[stage] + td.d_delta2[stage];
                        td.d_x[stage] += arithmetic::convert(td.d_delta[stage], td.d_lattice->ic());
                    }
                    arithmetic::convert(td.d_x_real[stage], td.d_x[stage], td.d_lattice->rc());
                    td.d_iszero[stage] = false;
                }
                
                td.d_jm.addJobs(td.d_job_list.begin(), td.d_job_list.end());
                td.d_job_list.clear();
            }
            
        public:
            inline bool operator() (ThreadData<RealTypeContext, IntTypeContext, JobManager> & td, JobData & jd,
                                    unsigned realdim, unsigned enumbegin, unsigned enumdim, unsigned stage, bool & noSolutionYet,
                                    typename RealTypeContext::Real & tmp) const
            {
                if (td.shouldIStop())
                    return false;
                if (td.shouldISplit())
                {
                    enqueue(td, jd, realdim, enumbegin, enumdim, stage, tmp);
                    return false;
                }
                if (td.updateBound(jd.d_bound)) // update bound
                    jd.propagate_new_bound(td);
                jd.go_stage_down<false>(td, realdim, enumbegin, enumdim, stage - 1, 0, tmp);
                CallUpdaterNotify caller;
                jd.do_the_run<true, true>(td, realdim, enumbegin, stage, stage - 1, noSolutionYet, caller, 0);
                if (caller.wasUpdated())
                {
                    // Inform other threads about update
                    td.d_jm.tryToUpdateBound(jd.d_result, jd.d_bound);
                }
                return true;
            }
        };
        
        template<bool avoid_bottom>
        static inline void go_stage_down(ThreadData<RealTypeContext, IntTypeContext, JobManager> & td,
                                  unsigned realdim, unsigned enumbegin, unsigned enumdim, unsigned stage, unsigned zerostage,
                                  typename RealTypeContext::Real & tmp)
        {
            if (stage > 0)
            {
                if (td.d_r[stage - 1] < td.d_r[stage])
                    td.d_r[stage - 1] = td.d_r[stage];
            }
            for (unsigned j = td.d_r[stage]; j > stage; --j) // this is the summation order considered in Pujol-Stehle
            {
                tmp = td.d_x_real[j] * td.d_lattice->getCoeff(enumbegin + j, enumbegin + stage);
                td.d_sigma(stage, j - 1) = td.d_sigma(stage, j) + tmp;
            }
            td.d_c[stage] = -td.d_sigma(stage, stage);
            bool ru;
            arithmetic::convert_round(td.d_x[stage], td.d_c[stage], ru, td.d_lattice->ic());
            td.d_iszero[stage] = td.d_iszero[stage + 1] && isZero(td.d_x[stage]);
            if (avoid_bottom && (stage == zerostage) && td.d_iszero[stage])
                // fuck it. this should _not_ be zero! (in case avoidbottom is set!)
                ++td.d_x[stage];
            arithmetic::convert(td.d_x_real[stage], td.d_x[stage], td.d_lattice->rc());
            td.d_delta[stage] = 0;
            td.d_delta2[stage] = ru ? 1 : -1;
        }
        
        template<bool avoid_bottom, bool inc_stage_after_call, class Caller>
        inline bool do_the_run(ThreadData<RealTypeContext, IntTypeContext, JobManager> & td,
                               unsigned realdim, unsigned enumbegin, unsigned enumdim, unsigned stage,
                               bool & noSolutionYet, const Caller & caller, unsigned zerostage = 0)
        {
            typename RealTypeContext::Real tmp(td.d_lattice->rc());
            
            while (true)
            {
                tmp = td.d_x_real[stage] - td.d_c[stage];
                square(td.d_ell[stage], tmp);
                td.d_ell[stage] *= td.d_lattice->getNormSq(enumbegin + stage);
                td.d_ell[stage] += td.d_ell[stage + 1];
                if (td.d_ell[stage] <= td.d_bounds[stage])
                {
                    if (stage == zerostage)
                    {
                        if (!caller(td, *this, realdim, enumbegin, enumdim, stage, noSolutionYet, tmp))
                            return false;
                        if (inc_stage_after_call)
                        {
                            if (++stage >= enumdim)
                                break;
                            // Update r for the lower stage
                            td.d_r[stage - 1] = stage;
                        }
                    }
                    else
                    {
                        --stage;
                        go_stage_down<avoid_bottom>(td, realdim, enumbegin, enumdim, stage, zerostage, tmp);
                        continue;
                    }
                }
                else
                {
                    if (++stage >= enumdim)
                        break;
                    // Update r for the lower stage
                    td.d_r[stage - 1] = stage;
                }
                // Modify d_x[stage]
                if (td.d_iszero[stage + 1])
                {
                    ++td.d_x[stage];
                }
                else
                {
                    td.d_delta2[stage] = -td.d_delta2[stage];
                    td.d_delta[stage] = -td.d_delta[stage] + td.d_delta2[stage];
                    td.d_x[stage] += arithmetic::convert(td.d_delta[stage], td.d_lattice->ic());
                }
                arithmetic::convert(td.d_x_real[stage], td.d_x[stage], td.d_lattice->rc());
                td.d_iszero[stage] = false;
            }
            return true;
        }
        
        bool d_done;
        
    public:
        JobData(unsigned dim, const linalg::math_rowvector<typename IntTypeContext::Integer> & x_end,
                const typename RealTypeContext::Real & ell,
                const typename RealTypeContext::Real & bound,
                bool iszero, unsigned ofs = 0)
            : d_x_end(dim), d_ell(ell), d_iszero(iszero), d_bound(bound), d_done(false)
        {
            for (unsigned i = 0; i < dim; ++i)
                d_x_end[i] = x_end[ofs + i];
        }
        
        void markDone()
        {
            d_done = true;
        }
        
        bool isDone() const
        {
            return d_done;
        }
        
        const typename RealTypeContext::Real & ell() const
        {
            return d_ell;
        }
        
        unsigned dimension() const
        {
            return d_x_end.size();
        }
        
        bool update(typename RealTypeContext::Real & bound)
        // Returns false if job is superfluous
        {
            d_bound = bound;
            return d_ell <= d_bound;
        }
        
        unsigned leftoverDims(const ThreadData<RealTypeContext, IntTypeContext, JobManager> & td) const
        {
            return td.d_dim - d_x_end.size();
        }
        
        void print()
        {
            for (unsigned i = 0; i < d_x_end.size(); ++i)
                std::cout << " " << d_x_end[i];
            std::cout << " -- " << d_ell << " -- " << d_bound << " -- " << (d_iszero ? "true" : "false") << "\n";
        }
        
        void run(ThreadData<RealTypeContext, IntTypeContext, JobManager> & td)
        {
            // Prepare enumeration
            td.doSetup();
            unsigned enumdim = td.d_dim - d_x_end.size(); // td.d_enumdim;
            bool noSolutionYet = true;
            run_init(td, td.d_dim, td.d_begin, enumdim);
            // Do enumeration
            CallUpdater caller;
            do_the_run<true, true>(td, td.d_dim, td.d_begin, enumdim, enumdim - 1, noSolutionYet, caller, 0);
        }
        
        void run2(ThreadData<RealTypeContext, IntTypeContext, JobManager> & td)
        {
            enum { THRESHOLD = 40 };
            
            // Prepare enumeration
            td.doSetup();
            unsigned enumdim = td.d_dim - d_x_end.size(); // td.d_enumdim;
            bool noSolutionYet = true;
            run_init(td, td.d_dim, td.d_begin, enumdim);
            if (enumdim > THRESHOLD)
            {
                // Do enumeration in two (main) stages
                CallRunner caller;
                do_the_run<false, false>(td, td.d_dim, td.d_begin, enumdim, enumdim - 1, noSolutionYet, caller, THRESHOLD);
            }
            else
            {
                // Do enumeration
                CallUpdater caller;
                do_the_run<true, true>(td, td.d_dim, td.d_begin, enumdim, enumdim - 1, noSolutionYet, caller, 0);
            }
        }
        
        void refine(ThreadData<RealTypeContext, IntTypeContext, JobManager> & td, unsigned dim)
        // Computes the next dim dimensions. We must have dim < leftoverDims(td).
        {
            // Prepare refinement
            td.doSetup();
            assert(dim < leftoverDims(td));
            unsigned realdim = dim + d_x_end.size();
            bool noSolutionYet; // not needed
            unsigned enumbegin = td.d_begin + td.d_dim - realdim;
            run_init(td, realdim, enumbegin, dim);
            
            // Do refinement
            CallRefiner caller;
            do_the_run<false, false>(td, realdim, enumbegin, dim, dim - 1, noSolutionYet, caller, 0);
        }
    };
}

#endif

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

#ifndef PLLL_INCLUDE_GUARD__LLL2_CALLBACK_CPP
#define PLLL_INCLUDE_GUARD__LLL2_CALLBACK_CPP

namespace plll
{
    template<class RealTypeContext, class IntTypeContext>
    class Callback
    {
    public:
        typedef boost::function<void(const linalg::math_matrix<typename IntTypeContext::Integer> &)> CallbackFunction;
        typedef boost::function<void(const linalg::math_matrix<typename IntTypeContext::Integer> &, unsigned,
                                     const typename IntTypeContext::Integer &)> MinCallbackFunction;
        
    private:
        CallbackFunction d_cf;
        unsigned d_begin;
        bool d_begin_note;
        unsigned long long d_prev_timestamp;
        unsigned long long d_interval;
        
        MinCallbackFunction d_mcf;
        TransformDataChangeLogger<IntTypeContext> d_T;
        
        MaxBitsCallbackFunction d_mbcf;
        unsigned d_last_max_bits;
        TransformDataMaxBitsCounter<IntTypeContext> d_Tb;
        
        bool d_min_init;
        typename IntTypeContext::Integer d_min;
        
        static inline unsigned long long timestamp()
        {
            timespec t;
            clock_gettime(CLOCK_MONOTONIC, &t);
            return ((unsigned long long)t.tv_sec) * 1000000000ul + (unsigned long long)t.tv_nsec;
        }
        
    public:
        Callback(unsigned dimension)
            : d_cf(NULL), d_begin(0), d_begin_note(false), d_prev_timestamp(0), d_interval(1),
              d_mcf(NULL), d_T(dimension), d_last_max_bits(0), d_min_init(false)
        {
        }
        
        void setCallback(CallbackFunction cf, double interval = 60*5, unsigned begin = 0, bool begin_note = false)
        {
            assert(interval >= 0);
            d_cf = cf;
            d_begin = begin;
            d_begin_note = begin_note;
            d_interval = interval * 1000000000ul;
            d_prev_timestamp = timestamp();
        }
        
        void setMinCallback(Lattice<RealTypeContext, IntTypeContext> & lattice, MinCallbackFunction mcf)
        {
            d_mcf = mcf;
            d_min_init = true;
            lattice.addNotifier(&d_T);
        }
        
        void setMaxBitsCallback(Lattice<RealTypeContext, IntTypeContext> & lattice, MaxBitsCallbackFunction mbcf)
        {
            d_mbcf = mbcf;
            d_Tb.addTo(lattice);
        }
        
        inline void operator() (Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned k)
        {
            // Take care of min callback
            if (!d_mcf.empty())
            {
                if (lattice.canGetIntegerVector())
                {
                    if (d_T.wasSomethingChanged())
                    {
                        bool changed = false;
                        unsigned idx;
                        typename IntTypeContext::Integer t;
                        // Initialization needed?
                        if (d_min_init)
                        {
                            lattice.getNormSqUP(d_min, 0);
                            idx = 0;
                            changed = true;
                            d_min_init = false;
                        }
                        // Seek for new minimum
                        for (unsigned i = d_T.getMinChanged(); i <= d_T.getMaxChanged(); ++i)
                            if (d_T.wasChanged(i))
                            {
                                lattice.getNormSqUP(t, i);
                                if ((t < d_min) && !isZero(t))
                                {
                                    d_min = t;
                                    idx = i;
                                    changed = true;
                                }
                            }
                        // Reset changes in logger
                        d_T.resetChangeLogger();
                    
                        // New minimum found?
                        if (changed)
                        {
                            std::pair<linalg::math_matrix<typename IntTypeContext::Integer> *, unsigned> ll = lattice.getMatrixIndex(idx);
                            if (ll.first)
                                d_mcf(*ll.first, ll.second, d_min);
                        }
                    }
                }
            }
            
            // Take care of usual callback
            if (!d_cf.empty())
            {
                bool do_now;
                if (d_interval > 0)
                {
                    // Measure time since last callback
                    unsigned long long ts = timestamp();
                    unsigned long long elapsed = ts - d_prev_timestamp;
                    // Check whether we want to call the callback function
                    do_now = ((k <= d_begin + 1) && d_begin_note);
                    if (elapsed >= d_interval)
                    {
                        do_now = true;
                        d_prev_timestamp = ts;
                    }
                }
                else
                    do_now = true;
                // Execute callback if wanted
                if (do_now)
                {
                    linalg::math_matrix<typename IntTypeContext::Integer> * ll = lattice.getMatrix();
                    if (ll)
                        d_cf(*ll);
                }
            }
            
            // Take care of max bits callback
            if (!d_mbcf.empty())
            {
                if (d_Tb.changed())
                {
                    if (d_last_max_bits != d_Tb.maxBits())
                    {
                        d_last_max_bits = d_Tb.maxBits();
                        d_mbcf(lattice.dimension(), d_Tb.maxBits());
                    }
                    // Reset changes in logger
                    d_Tb.resetChange();
                }
            }
        }
    };
}

#endif

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
#include "lll2-internal.hpp"
#include "transform.cpp"

#include <sstream>
#include <boost/bind.hpp>

#include <plll/arithmetic.hpp>
#include <plll/rational.hpp>
#if !defined(PLLL_CONFIG_NO_ARITHMETIC_LONGDOUBLE) || !defined(PLLL_CONFIG_NO_ARITHMETIC_DOUBLE)
  #include "nfp-wrapper.hpp"
#endif
#if !defined(PLLL_CONFIG_NO_ARITHMETIC_QUADDOUBLE) || !defined(PLLL_CONFIG_NO_ARITHMETIC_DOUBLEDOUBLE)
  #include "ddqd-wrapper.hpp"
#endif

namespace plll
{
    class LRIISelector : public LRIInterface
    {
    private:
        class Interface
        {
        private:
            long d_max_bits;
            boost::function<LRIInterface*(linalg::math_matrix<arithmetic::Integer> &)> d_interface_generator;
            mutable const LatticeReduction::GramSchmidtInformer * d_gsi;
            
        public:
            Interface(long max_bits, boost::function<LRIInterface*(linalg::math_matrix<arithmetic::Integer> &)> interfacegen)
                : d_max_bits(max_bits), d_interface_generator(interfacegen), d_gsi(NULL)
            {
            }
            
            ~Interface()
            {
            }
            
            LRIInterface * generateInterface(linalg::math_matrix<arithmetic::Integer> & matrix) const
            {
                return d_interface_generator(matrix);
            }
            
            long maxSupportedBits() const
            {
                return d_max_bits;
            }
            
            bool isFine(unsigned bits, unsigned dimension) const
            {
                long bits_needed = 2 * bits + arithmetic::approxLog2(dimension) + 4;
                return bits_needed <= d_max_bits;
            }
        };
        
        std::list<Interface> d_interfaces;
        mutable std::list<Interface>::const_iterator d_interface_gen, d_better_interface_gen;
        mutable STD_AUTO_PTR<LRIInterface> d_interface;
        mutable bool d_continue, d_reselect;
        mutable MaxBitsCallbackFunction d_maxbitscf;
        
        inline bool isBetter(unsigned bits, unsigned dimension) const
        {
            return (d_better_interface_gen == d_interfaces.end()) ? false : d_better_interface_gen->isFine(bits, dimension);
        }        
        
        void maxBitsCallbackFunction(unsigned dimension, unsigned maxbits)
        {
            if (!d_maxbitscf.empty())
                d_maxbitscf(dimension, maxbits);
            if (d_verbose_function && Verbose::yieldsOutput(d_verbose_outputlevel, LatticeReduction::VL_Chatter))
            {
                std::ostringstream s;
                s << "...base has " << maxbits << " bits...";
                d_verbose_function(LatticeReduction::VL_Chatter, s.str());
            }
            if (!d_interface_gen->isFine(maxbits, dimension) || isBetter(maxbits, dimension))
            {
                d_reselect = true;
                throw change_interface_exception();
            }
        }
        
        bool isContinuing() const
        {
            if (d_reselect)
                selectInterface(true);
            bool v = d_continue;
            d_continue = false;
            return v;
        }
        
        unsigned maxBits() const
        {
            unsigned maxbits = 0;
            for (unsigned i = 0; i < d_lattice.rows(); ++i)
                for (unsigned j = 0; j < d_lattice.cols(); ++j)
                {
                    long b = arithmetic::approxLog2(d_lattice(i, j));
                    if (b > maxbits)
                        maxbits = b;
                }
            return maxbits;
        }
        
        void selectInterface(bool withContinue = false) const
        {
            if (d_interface.get())
            {
                if (d_verbose_function && Verbose::yieldsOutput(d_verbose_outputlevel, LatticeReduction::VL_Information))
                    d_verbose_function(LatticeReduction::VL_Information, "Applying interface change");
                d_interface.reset();
            }
            
            d_continue = withContinue;
            d_reselect = false;
            unsigned bits = maxBits();
            d_interface_gen = d_interfaces.end();
            for (std::list<Interface>::const_iterator i = d_interfaces.begin(); i != d_interfaces.end(); ++i)
            {
                if (i->isFine(bits, d_lattice.rows()))
                    d_interface_gen = i;
            }
            assert(d_interface_gen != d_interfaces.end());
            
            if (d_verbose_function && Verbose::yieldsOutput(d_verbose_outputlevel, LatticeReduction::VL_Information))
            {
                std::ostringstream s;
                s << "Selecting interface which supports up to " << d_interface_gen->maxSupportedBits() << " bits";
                d_verbose_function(LatticeReduction::VL_Information, s.str());
            }
            
            // Generate interface
            d_interface.reset(d_interface_gen->generateInterface(d_lattice));
            if (d_ensure_min_prec)
                d_interface->ensureMinimumPrecision(d_min_prec);
            if (d_set_verbose)
                d_interface->setupVerbose(d_verbose_outputlevel, d_verbose_function);
            
            // Select next interface
            d_better_interface_gen = d_interface_gen;
            ++d_better_interface_gen;
        }
        
        const LRIInterface & getInterface() const
        {
            return *d_interface;
        }
        
        LRIInterface & getInterface()
        {
            return *d_interface;
        }        
        
        void unselectInterface() const
        {
            d_interface.reset();
        }
        
        class GSI : public LatticeReduction::GramSchmidtInformer
        {
        private:
            LRIISelector & d_selector;
            
        public:
            GSI(LRIISelector & selector)
                : d_selector(selector)
            {
            }
            
            virtual ~GSI()
            {
            }
            
            virtual double getGSCoefficientD(unsigned i, unsigned j) const
            {
                d_selector.selectInterface();
                double r = d_selector.getInterface().getInformer()->getGSCoefficientD(i, j);
                d_selector.unselectInterface();
                return r;
            }
            
            virtual long double getGSCoefficientLD(unsigned i, unsigned j) const
            {
                d_selector.selectInterface();
                long double r = d_selector.getInterface().getInformer()->getGSCoefficientLD(i, j);
                d_selector.unselectInterface();
                return r;
            }
            
            virtual arithmetic::Real getGSCoefficientR(unsigned i, unsigned j, const arithmetic::RealContext & rc) const
            {
                d_selector.selectInterface();
                arithmetic::Real r = d_selector.getInterface().getInformer()->getGSCoefficientR(i, j, rc);
                d_selector.unselectInterface();
                return r;
            }
            
            virtual double getGSSqNormD(unsigned i) const
            {
                d_selector.selectInterface();
                double r = d_selector.getInterface().getInformer()->getGSSqNormD(i);
                d_selector.unselectInterface();
                return r;
            }
            
            virtual long double getGSSqNormLD(unsigned i) const
            {
                d_selector.selectInterface();
                long double r = d_selector.getInterface().getInformer()->getGSSqNormLD(i);
                d_selector.unselectInterface();
                return r;
            }
            
            virtual arithmetic::Real getGSSqNormR(unsigned i, const arithmetic::RealContext & rc) const
            {
                d_selector.selectInterface();
                arithmetic::Real r = d_selector.getInterface().getInformer()->getGSSqNormR(i, rc);
                d_selector.unselectInterface();
                return r;
            }
            
            virtual double computeProjectionLengthD(unsigned k, unsigned b, const linalg::math_rowvector<arithmetic::Integer> & vec) const
            {
                d_selector.selectInterface();
                double r = d_selector.getInterface().getInformer()->computeProjectionLengthD(k, b, vec);
                d_selector.unselectInterface();
                return r;
            }
            
            virtual long double computeProjectionLengthLD(unsigned k, unsigned b, const linalg::math_rowvector<arithmetic::Integer> & vec) const
            {
                d_selector.selectInterface();
                long double r = d_selector.getInterface().getInformer()->computeProjectionLengthLD(k, b, vec);
                d_selector.unselectInterface();
                return r;
            }
            
            virtual arithmetic::Real computeProjectionLengthR(unsigned k, unsigned b, const linalg::math_rowvector<arithmetic::Integer> & vec, const arithmetic::RealContext & rc) const
            {
                d_selector.selectInterface();
                arithmetic::Real r = d_selector.getInterface().getInformer()->computeProjectionLengthR(k, b, vec, rc);
                d_selector.unselectInterface();
                return r;
            }
        };
        
        friend class GSI;
        
        GSI d_gsi;
        
        bool d_ensure_min_prec;
        unsigned long d_min_prec;
        
        bool d_set_verbose;
        LatticeReduction::VerboseOutputLevel d_verbose_outputlevel;
        LatticeReduction::VerboseFunction d_verbose_function;
        
    public:
        LRIISelector(linalg::math_matrix<arithmetic::Integer> & lattice, LatticeReduction::VerboseOutputLevel vol, LatticeReduction::VerboseFunction vf)
            : LRIInterface(lattice), d_interface(), d_continue(false), d_reselect(false), d_gsi(*this),
              d_ensure_min_prec(false), d_min_prec(0), d_set_verbose(false), d_verbose_outputlevel(vol), d_verbose_function(vf)
        {
        }
        
        virtual ~LRIISelector()
        {
        }
        
        void addInterface(long max_bits, boost::function<LRIInterface*(linalg::math_matrix<arithmetic::Integer> &)> interfacegen)
        {
            d_interfaces.push_back(Interface(max_bits, interfacegen));
        }
        
        virtual void setupVerbose(LatticeReduction::VerboseOutputLevel vol, LatticeReduction::VerboseFunction vf)
        {
            d_set_verbose = true;
            d_verbose_outputlevel = vol;
            d_verbose_function = vf;
        }
        
        virtual const LatticeReduction::GramSchmidtInformer * getInformer() const
        {
            return &d_gsi;
        }
        
        virtual void ensureMinimumPrecision(unsigned long prec)
        {
            d_ensure_min_prec = prec;
            d_min_prec = prec;
        }
        
        virtual void forceGSRebuild(bool b)
        {
        }
        
        virtual double getGSCoefficientD(unsigned i, unsigned j) const
        {
            selectInterface();
            double r = getInterface().getGSCoefficientD(i, j);
            unselectInterface();
            return r;
        }
        
        virtual long double getGSCoefficientLD(unsigned i, unsigned j) const
        {
            selectInterface();
            long double r = getInterface().getGSCoefficientLD(i, j);
            unselectInterface();
            return r;
        }
        
        virtual arithmetic::Real getGSCoefficientR(unsigned i, unsigned j, const arithmetic::RealContext & rc) const
        {
            selectInterface();
            arithmetic::Real r = getInterface().getGSCoefficientR(i, j, rc);
            unselectInterface();
            return r;
        }
        
        virtual double getGSSqNormD(unsigned i) const
        {
            selectInterface();
            double r = getInterface().getGSSqNormD(i);
            unselectInterface();
            return r;
        }
        
        virtual long double getGSSqNormLD(unsigned i) const
        {
            selectInterface();
            long double r = getInterface().getGSSqNormLD(i);
            unselectInterface();
            return r;
        }
        
        virtual arithmetic::Real getGSSqNormR(unsigned i, const arithmetic::RealContext & rc) const
        {
            selectInterface();
            arithmetic::Real r = getInterface().getGSSqNormR(i, rc);
            unselectInterface();
            return r;
        }
        
        virtual void modFlip(Transform & transform, unsigned i)
        {
            selectInterface();
            getInterface().modFlip(transform, i);
            unselectInterface();
        }
        
        virtual void modSwap(Transform & transform, unsigned i, unsigned j)
        {
            selectInterface();
            getInterface().modSwap(transform, i, j);
            unselectInterface();
        }
        
        virtual void modAdd(Transform & transform, unsigned i, unsigned j, const arithmetic::Integer & m)
        {
            selectInterface();
            getInterface().modAdd(transform, i, j, m);
            unselectInterface();
        }
        
        virtual void sortProjected(unsigned & begin, unsigned & end, Transform & transform)
        {
            selectInterface();
            getInterface().sortProjected(begin, end, transform);
            unselectInterface();
        }
        
        virtual void sizereduction(unsigned & begin, unsigned & end, Transform & transform)
        {
            selectInterface();
            getInterface().sizereduction(begin, end, transform);
            unselectInterface();
        }
        
        virtual void lll(unsigned & begin, unsigned & end, Transform & transform, double alpha, LatticeReduction::LLLMode mode,
                         LatticeReduction::CallbackFunction cf, LatticeReduction::CallbackFunction_LI cf2, double cf_int,
                         LatticeReduction::MinCallbackFunction mcf, LatticeReduction::MinCallbackFunction_LI mcf2,
                         MaxBitsCallbackFunction mbcf, bool anneal, LatticeReduction::AnnealCallbackFunction acf,
                         LatticeReduction::LLL_AnnealFunction af, LatticeReduction::DIMethod di, LatticeReduction::DIMode di_mode,
                         LatticeReduction::DIChoice di_choice, unsigned di_bs)
        {
            selectInterface(true);
            d_maxbitscf = mbcf;
            while (isContinuing())
                getInterface().lll(begin, end, transform, alpha, mode, cf, cf2, cf_int, mcf, mcf2,
                                   boost::bind(&LRIISelector::maxBitsCallbackFunction, this, _1, _2),
                                   anneal, acf, af, di, di_mode, di_choice, di_bs);
            unselectInterface();
        }
        
        virtual void bkz(unsigned & begin, unsigned & end, Transform & transform, double alpha, unsigned blocksize, LatticeReduction::BKZMode mode,
                         LatticeReduction::CallbackFunction cf, LatticeReduction::CallbackFunction_LI cf2, double cf_int,
                         LatticeReduction::MinCallbackFunction mcf, LatticeReduction::MinCallbackFunction_LI mcf2, MaxBitsCallbackFunction mbcf,
                         LatticeReduction::EnumCallbackFunction ecf, LatticeReduction::EnumCallbackFunction_LI ecf2, bool anneal,
                         LatticeReduction::AnnealCallbackFunction acf, LatticeReduction::BKZ_AnnealFunction af,
                         LatticeReduction::DIMethod di, LatticeReduction::DIMode di_mode, LatticeReduction::DIChoice di_choice, unsigned di_bs)
        {
            selectInterface(true);
            d_maxbitscf = mbcf;
            while (isContinuing())
                getInterface().bkz(begin, end, transform, alpha, blocksize, mode, cf, cf2, cf_int, mcf, mcf2,
                                   boost::bind(&LRIISelector::maxBitsCallbackFunction, this, _1, _2),
                                   ecf, ecf2, anneal, acf, af, di, di_mode, di_choice, di_bs);
            unselectInterface();
        }
        
        virtual void hkz(unsigned & begin, unsigned & end, Transform & transform, bool dual,
                         LatticeReduction::CallbackFunction cf, LatticeReduction::CallbackFunction_LI cf2, double cf_int,
                         LatticeReduction::MinCallbackFunction mcf, LatticeReduction::MinCallbackFunction_LI mcf2,
                         MaxBitsCallbackFunction mbcf,
                         LatticeReduction::EnumCallbackFunction ecf, LatticeReduction::EnumCallbackFunction_LI ecf2)
        {
            selectInterface(true);
            d_maxbitscf = mbcf;
            while (isContinuing())
                getInterface().hkz(begin, end, transform, dual, cf, cf2, cf_int, mcf, mcf2,
                                   boost::bind(&LRIISelector::maxBitsCallbackFunction, this, _1, _2),
                                   ecf, ecf2);
            unselectInterface();
        }
        
        virtual void svp(unsigned & begin, unsigned & end, Transform & transform, bool make_basis, bool extreme, bool dual,
                         LatticeReduction::CallbackFunction cf, LatticeReduction::CallbackFunction_LI cf2, double cf_int,
                         LatticeReduction::MinCallbackFunction mcf, LatticeReduction::MinCallbackFunction_LI mcf2,
                         MaxBitsCallbackFunction mbcf,
                         LatticeReduction::EnumCallbackFunction ecf, LatticeReduction::EnumCallbackFunction_LI ecf2)
        {
            selectInterface(true);
            d_maxbitscf = mbcf;
            while (isContinuing())
                getInterface().svp(begin, end, transform, make_basis, extreme, dual, cf, cf2, cf_int, mcf, mcf2,
                                   boost::bind(&LRIISelector::maxBitsCallbackFunction, this, _1, _2),
                                   ecf, ecf2);
            unselectInterface();
        }
        
        virtual bool isSizeReduced(unsigned begin, unsigned end) const
        {
            selectInterface();
            bool ret = getInterface().isSizeReduced(begin, end);
            unselectInterface();
            return ret;
        }
        
        virtual bool isLLLBasis(unsigned begin, unsigned end, double alpha, LatticeReduction::LLLMode mode,
                                LatticeReduction::DIMethod di, LatticeReduction::DIChoice di_choice, unsigned di_bs) const
        {
            selectInterface();
            bool ret = getInterface().isLLLBasis(begin, end, alpha, mode, di, di_choice, di_bs);
            unselectInterface();
            return ret;
        }
        
        virtual bool isBKZBasis(unsigned begin, unsigned end, double alpha, unsigned blocksize, LatticeReduction::BKZMode mode,
                                LatticeReduction::DIMethod di, LatticeReduction::DIChoice di_choice, unsigned di_bs) const
        {
            selectInterface();
            bool ret = getInterface().isBKZBasis(begin, end, alpha, blocksize, mode, di, di_choice, di_bs);
            unselectInterface();
            return ret;
        }
        
        virtual bool isHKZBasis(unsigned begin, unsigned end, bool dual) const
        {
            selectInterface();
            bool ret = getInterface().isHKZBasis(begin, end, dual);
            unselectInterface();
            return ret;
        }
        
        virtual bool isSVPBasis(unsigned begin, unsigned end, bool dual) const
        {
            selectInterface();
            bool ret = getInterface().isSVPBasis(begin, end, dual);
            unselectInterface();
            return ret;
        }
    };
    
    template<class RealTypeContext, class IntTypeContext>
    LRIInterface * CreateLRIInterfaceWithContexts(LatticeReduction::VerboseOutputLevel vol, LatticeReduction::VerboseFunction vf,
                                                  linalg::math_matrix<arithmetic::Integer> & lattice, LatticeReduction::GramSchmidt gs, bool gsr,
                                                  LatticeReduction::SVPMode svp, unsigned max_cores,
                                                  const RealTypeContext & rc, const IntTypeContext & ic);
    
    template<class RealTypeContext>
    LRIInterface * CreateLRIInterface(LatticeReduction::VerboseOutputLevel vol, LatticeReduction::VerboseFunction vf,
                                      linalg::math_matrix<arithmetic::Integer> & lattice, LatticeReduction::GramSchmidt gs, bool gsr,
                                      LatticeReduction::SVPMode svp, unsigned max_cores,
                                      LatticeReduction::Integers ints, const RealTypeContext & rc)
    {
        switch (ints)
        {
#if !defined(PLLL_CONFIG_NO_ARITHMETIC_BIGINT) && !defined(PLLL_CONFIG_NO_ARITHMETIC_LONGINT)
        case LatticeReduction::I_Auto:
        {
            LRIISelector * sel = new LRIISelector(lattice, vol, vf);
            // Add in decreasing order of supported precisions
            sel->addInterface(std::numeric_limits<long>::max(),
                              boost::bind(CreateLRIInterfaceWithContexts<RealTypeContext, arithmetic::IntegerContext>,
                                          vol, vf, _1, gs, gsr, svp, max_cores, rc, arithmetic::IntegerContext()));
            sel->addInterface(std::numeric_limits<long int>::digits - 1,
                              boost::bind(CreateLRIInterfaceWithContexts<RealTypeContext, arithmetic::NIntContext<long int> >,
                                          vol, vf, _1, gs, gsr, svp, max_cores, rc, arithmetic::NIntContext<long int>()));
            return sel;
        }
#endif
#ifndef PLLL_CONFIG_NO_ARITHMETIC_BIGINT
        case LatticeReduction::I_ArbitraryPrecision:
            return CreateLRIInterfaceWithContexts(vol, vf, lattice, gs, gsr, svp, max_cores, rc, arithmetic::IntegerContext());
#endif
#ifndef PLLL_CONFIG_NO_ARITHMETIC_LONGINT
        case LatticeReduction::I_LongInt:
            return CreateLRIInterfaceWithContexts(vol, vf, lattice, gs, gsr, svp, max_cores, rc, arithmetic::NIntContext<long int>());
#endif
        }
        assert(!"Integer arithmetic not supported!");
        return NULL;
    }
    
    LRIInterface * CreateLRIInterface(LatticeReduction::VerboseOutputLevel vol, LatticeReduction::VerboseFunction vf,
                                      linalg::math_matrix<arithmetic::Integer> & lattice,
                                      LatticeReduction::Arithmetic arith, LatticeReduction::Integers ints,
                                      LatticeReduction::GramSchmidt gs, bool gsr,
                                      LatticeReduction::SVPMode svp, unsigned max_cores)
    {
        switch(arith)
        {
        default:
#ifndef PLLL_CONFIG_NO_ARITHMETIC_LONGDOUBLE
        case LatticeReduction::A_LongDouble:   return CreateLRIInterface(vol, vf, lattice, gs, gsr, svp, max_cores, ints, arithmetic::NFPContext<long double>());
#endif
#ifndef PLLL_CONFIG_NO_ARITHMETIC_REAL
        case LatticeReduction::A_Real:         return CreateLRIInterface(vol, vf, lattice, gs, gsr, svp, max_cores, ints, arithmetic::RealContext());
#endif
#ifndef PLLL_CONFIG_NO_ARITHMETIC_RATIONAL
        case LatticeReduction::A_Rational:     return CreateLRIInterface(vol, vf, lattice, gs, gsr, svp, max_cores, ints, arithmetic::RationalContext());
#endif
#ifndef PLLL_CONFIG_NO_ARITHMETIC_DOUBLE
        case LatticeReduction::A_Double:       return CreateLRIInterface(vol, vf, lattice, gs, gsr, svp, max_cores, ints, arithmetic::NFPContext<double>());
#endif
#ifndef PLLL_CONFIG_NO_ARITHMETIC_DOUBLEDOUBLE
        case LatticeReduction::A_DoubleDouble: return CreateLRIInterface(vol, vf, lattice, gs, gsr, svp, max_cores, ints, arithmetic::DDQDContext<dd_real>());
#endif
#ifndef PLLL_CONFIG_NO_ARITHMETIC_QUADDOUBLE
        case LatticeReduction::A_QuadDouble:   return CreateLRIInterface(vol, vf, lattice, gs, gsr, svp, max_cores, ints, arithmetic::DDQDContext<qd_real>());
#endif
        }
        assert(!"Real arithmetic not supported!");
        return NULL;
    }
    
    // Explicit instantiation of templates
    
#ifndef PLLL_CONFIG_NO_ARITHMETIC_BIGINT
    template void setUnit<arithmetic::IntegerContext>(linalg::math_matrix<arithmetic::Integer> &);
#endif
    
#ifndef PLLL_CONFIG_NO_ARITHMETIC_LONGINT
    template void setUnit<arithmetic::NIntContext<long> >(linalg::math_matrix<arithmetic::NInt<long> > &);
#endif
}

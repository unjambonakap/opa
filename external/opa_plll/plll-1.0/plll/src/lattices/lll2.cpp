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
#include "lattice.cpp"
#include "buffer.hpp"
#include <istream>
#include <list>
#include <cmath>

#include <limits>
#ifndef PLLL_CONFIG_NO_ARITHMETIC_DOUBLEDOUBLE
  #include <qd/dd_real.h>
#endif
#ifndef PLLL_CONFIG_NO_ARITHMETIC_QUADDOUBLE
  #include <qd/qd_real.h>
#endif

#include "matrixconversion.hpp"

namespace plll
{
    void CallbackAdaptor(const linalg::math_matrix<arithmetic::Integer> & lattice, const LatticeReduction::CallbackFunction_LI & cf)
    {
        arithmetic::NIntContext<long int> ic;
        MatrixConversionConst<arithmetic::NIntContext<long int> > m(ic, lattice);
        cf(m.matrix());
    }
    
    void CallbackAdaptor_LI(const linalg::math_matrix<arithmetic::NInt<long int> > & lattice, const LatticeReduction::CallbackFunction & cf)
    {
        arithmetic::NIntContext<long int> ic;
        MatrixConversionConst2<arithmetic::NIntContext<long int> > m(ic, lattice);
        cf(m.matrix());
    }
    
    void MinCallbackAdaptor(const linalg::math_matrix<arithmetic::Integer> & lattice, unsigned idx,
                            const arithmetic::Integer & len, const LatticeReduction::MinCallbackFunction_LI & cf)
    {
        arithmetic::NIntContext<long int> ic;
        MatrixConversionConst<arithmetic::NIntContext<long int> > m(ic, lattice);
        cf(m.matrix(), idx, arithmetic::convert(len, ic));
    }
    
    void MinCallbackAdaptor_LI(const linalg::math_matrix<arithmetic::NInt<long int> > & lattice, unsigned idx,
                               const arithmetic::NInt<long int> & len, const LatticeReduction::MinCallbackFunction & cf)
    {
        arithmetic::NIntContext<long int> ic;
        MatrixConversionConst2<arithmetic::NIntContext<long int> > m(ic, lattice);
        cf(m.matrix(), idx, arithmetic::convert<arithmetic::Integer>(len));
    }
    
    void EnumCallbackAdaptor(const linalg::math_matrix<arithmetic::Integer> & lattice, int idx,
                             const linalg::math_rowvector<arithmetic::Integer> & vec, const LatticeReduction::EnumCallbackFunction_LI & cf)
    {
        arithmetic::NIntContext<long int> ic;
        MatrixConversionConst<arithmetic::NIntContext<long int> > m(ic, lattice);
        VectorConversionConst<arithmetic::NIntContext<long int> > v(ic, vec);
        cf(m.matrix(), idx, v.vector());
    }
    
    void EnumCallbackAdaptor_LI(const linalg::math_matrix<arithmetic::NInt<long int> > & lattice, int idx,
                                const linalg::math_rowvector<arithmetic::NInt<long int> > & vec, const LatticeReduction::EnumCallbackFunction & cf)
    {
        arithmetic::NIntContext<long int> ic;
        MatrixConversionConst2<arithmetic::NIntContext<long int> > m(ic, lattice);
        VectorConversionConst2<arithmetic::NIntContext<long int> > v(ic, vec);
        cf(m.matrix(), idx, v.vector());
    }
    
    const char * LatticeReduction::stop_reduction::what() const PLLL_INTERNAL_NOTHROW_POSTFIX_ENFORCE
    {
        return "Stopping lattice reduction.";
    }
    
    const char * LatticeReduction::stop_enumeration::what() const PLLL_INTERNAL_NOTHROW_POSTFIX_ENFORCE
    {
        return "Stopping lattice enumeration.";
    }
    
    const char * reduction_error::what() const PLLL_INTERNAL_NOTHROW_POSTFIX_ENFORCE
    {
        return d_message.c_str();
    }
    
    const char * change_interface_exception::what() const PLLL_INTERNAL_NOTHROW_POSTFIX_ENFORCE
    {
        return "Changing lattice reduction interface.";
    }
    
    const char * feature_not_implemented::what() const PLLL_INTERNAL_NOTHROW_POSTFIX_ENFORCE
    {
        return "Feature not implemented!";
    }
    
    LRIInterface * CreateLRIInterface(LatticeReduction::VerboseOutputLevel vl, LatticeReduction::VerboseFunction vf, linalg::math_matrix<arithmetic::Integer> & lattice,
                                      LatticeReduction::Arithmetic arith, LatticeReduction::Integers ints, LatticeReduction::GramSchmidt gs,
                                      bool gsr, LatticeReduction::SVPMode svp, unsigned max_cores); // see lll2-multiplexer.cpp
    
    void DefaultVerboseFunction(LatticeReduction::VerboseLevel, const std::string & s)
    {
        std::cerr << s << "\n";
    }
    
    void GaussianFactorComputer::computeHPF(unsigned dim)
    // Make sure HPF(dim) is available in cache.
    {
        if (d_hermite_prefactors.size() > dim)
            return;
        d_hermite_prefactors.reserve(dim + 1);
        while (d_hermite_prefactors.size() <= dim)
        {
            unsigned d = d_hermite_prefactors.size();
            long double fact = 1.0l / ::sqrt(3.1415926535897932384626433832795028841971693993751l);
            d_hermite_prefactors.push_back(std::exp(lgammal((long double)d * (long double)0.5 + (long double)1) / (long double)d) * fact);
        }
    }
    
    void addStatistics(LatticeReduction::Statistics & stats, const LatticeReduction::Statistics & add)
    {
        stats.swaps = add.swaps;
        stats.adds = add.adds;
        stats.adds_pm1 = add.adds_pm1;
        stats.adds_pm2 = add.adds_pm2;
        stats.flips = add.flips;
        stats.trans = add.trans;
        stats.sizereductions = add.sizereductions;
        stats.deepinsertions = add.deepinsertions;
        stats.enumcalls = add.enumcalls;
        stats.enumfails = add.enumfails;
        stats.vectorinsertions = add.vectorinsertions;
        stats.vectorinsertions_rearrange = add.vectorinsertions_rearrange;
    }
    
    void LatticeReduction::Statistics::reset()
    {
        swaps = 0;
        adds = 0;
        adds_pm1 = 0;
        adds_pm2 = 0;
        flips = 0;
        trans = 0;
        sizereductions = 0;
        deepinsertions = 0;
        enumcalls = 0;
        enumfails = 0;
        vectorinsertions = 0;
        vectorinsertions_rearrange = 0;
    }
    
    LatticeReduction::Statistics::Statistics()
    {
        reset();
    }
    
    void Ranger::insertVector(unsigned ofs)
    {
        for (std::deque<Range>::iterator i = d_ranges.begin(); i != d_ranges.end(); ++i)
        {
            if (i->first > ofs)
                ++i->first;
            if (i->second >= ofs)
                ++i->second;
        }
    }
    
    void Ranger::removeVector(unsigned ofs)
    {
        for (std::deque<Range>::iterator i = d_ranges.begin(); i != d_ranges.end(); ++i)
        {
            if (i->first > ofs)
                --i->first;
            if (i->second >= ofs)
                --i->second;
            assert(i->second >= i->first);
        }
    }
    
    class LatticeReductionImpl
    {
    private:
        mutable linalg::math_matrix<arithmetic::Integer> d_lattice;
        
        LatticeReduction::Arithmetic d_arith;
        LatticeReduction::Integers d_ints;
        LatticeReduction::GramSchmidt d_gs;
        bool d_gsr;
        
        LatticeReduction::SVPMode d_svp;
        
        unsigned d_max_cores;
        
        LatticeReduction::CallbackFunction d_cf;
        LatticeReduction::CallbackFunction_LI d_cf2;
        double d_cf_int;
        LatticeReduction::MinCallbackFunction d_mcf;
        LatticeReduction::MinCallbackFunction_LI d_mcf2;
        LatticeReduction::EnumCallbackFunction d_ecf;
        LatticeReduction::EnumCallbackFunction_LI d_ecf2;
        
        LatticeReduction::AnnealCallbackFunction d_acf;
        LatticeReduction::LLL_AnnealFunction d_af_lll;
        LatticeReduction::BKZ_AnnealFunction d_af_bkz;
        bool d_anneal_lll, d_anneal_bkz;
        
        LatticeReduction::DIMethod d_di;
        LatticeReduction::DIMode d_di_mode;
        LatticeReduction::DIChoice d_di_choice;
        unsigned d_di_bs;
        
        Transform d_trans;
        
        unsigned d_begin, d_end;
        unsigned long d_min_prec;
        
        LatticeReduction::VerboseOutputLevel d_verboseoutputlevel;
        LatticeReduction::VerboseFunction d_verbosefunction;
        
        mutable STD_AUTO_PTR<LRIInterface> d_if;
        
        inline void setupInterface() const
        {
            if (d_if.get() == NULL)
            {
                d_if.reset(CreateLRIInterface(d_verboseoutputlevel, d_verbosefunction, d_lattice, d_arith, d_ints, d_gs, d_gsr, d_svp, d_max_cores));
                d_if->ensureMinimumPrecision(d_min_prec);
            }
        }
        
        inline bool hasInterface() const
        {
            return d_if.get() != NULL;
        }
    
        inline void clearInterface()
        {
            d_if.reset();
        }
    
    public:
        inline LatticeReductionImpl()
            : d_arith(LatticeReduction::A_Default), d_ints(LatticeReduction::I_Default), d_gs(LatticeReduction::G_Default), d_gsr(false),
              d_svp(LatticeReduction::SVP_Default), d_max_cores(0), d_cf(NULL), d_cf2(NULL), d_cf_int(60.0 * 5.0), d_mcf(NULL), d_mcf2(NULL),
              d_ecf(NULL), d_ecf2(NULL), d_acf(NULL), d_af_lll(NULL), d_af_bkz(NULL), d_anneal_lll(false), d_anneal_bkz(false),
              d_di(LatticeReduction::DI_Default), d_di_mode(LatticeReduction::DIM_Default),
              d_di_choice(LatticeReduction::DIC_Default), d_di_bs(0), d_begin(0), d_end(0), d_min_prec(0),
              d_verboseoutputlevel(LatticeReduction::VOL_Warnings), d_verbosefunction(DefaultVerboseFunction), d_if()
        {
        }
        
        inline void setLattice(const linalg::math_matrix<arithmetic::Integer> & lattice)
        {
            clearInterface();
            d_lattice = lattice;
            d_begin = 0;
            d_end = d_lattice.rows() ? d_lattice.rows() - 1 : 0;
            d_trans.reset(d_lattice);
        }
        
        template<typename IType>
        inline void setLattice(const linalg::math_matrix<arithmetic::NInt<IType> > & lattice)
        {
            clearInterface();
            d_lattice.resize(lattice.rows(), lattice.cols());
            for (unsigned i = 0; i < lattice.rows(); ++i)
                for (unsigned j = 0; j < lattice.cols(); ++j)
                    d_lattice(i, j) = arithmetic::convert<arithmetic::Integer>(lattice(i, j));
            d_begin = 0;
            d_end = d_lattice.rows() ? d_lattice.rows() - 1 : 0;
            d_trans.reset(d_lattice);
        }
        
#if __cplusplus >= 201103L
        inline void setLattice(linalg::math_matrix<arithmetic::Integer> && lattice)
        {
            clearInterface();
            d_lattice = std::move(lattice);
            d_begin = 0;
            d_end = d_lattice.rows() ? d_lattice.rows() - 1 : 0;
            d_trans.reset(d_lattice);
        }
#endif
        
        inline const linalg::math_matrix<arithmetic::Integer> & getLattice() const
        {
            return d_lattice;
        }
        
        inline void setArithmetic(LatticeReduction::Arithmetic arith)
        {
            if (d_arith == arith)
                return;
            clearInterface();
            d_arith = arith;
        }
        
        inline LatticeReduction::Arithmetic getArithmetic() const
        {
            return d_arith;
        }
        
        inline void setIntegers(LatticeReduction::Integers ints)
        {
            if (d_ints == ints)
                return;
            clearInterface();
            d_ints = ints;
        }
        
        inline LatticeReduction::Integers getIntegers() const
        {
            return d_ints;
        }
        
        inline bool ensurePrecision(unsigned long p)
        // Tries to ensure a minimum precision. If requirement can be realized, returns true, otherwise false.
        {
            d_min_prec = p;
            if (d_if.get())
                d_if->ensureMinimumPrecision(d_min_prec);
            switch (d_arith)
            {
            case LatticeReduction::A_Rational:     return true;
            case LatticeReduction::A_Real:         return true;
            case LatticeReduction::A_LongDouble:   return (p <= (unsigned long)std::numeric_limits<long>::max()) && ((long)p <= std::numeric_limits<long double>::digits);
            case LatticeReduction::A_Double:       return (p <= (unsigned long)std::numeric_limits<long>::max()) && ((long)p <= std::numeric_limits<double>::digits);
#ifndef PLLL_CONFIG_NO_ARITHMETIC_DOUBLEDOUBLE
            case LatticeReduction::A_DoubleDouble: return (p <= (unsigned long)std::numeric_limits<long>::max()) && ((long)p <= std::numeric_limits<dd_real>::digits);
#else
            case LatticeReduction::A_DoubleDouble: break;
#endif
#ifndef PLLL_CONFIG_NO_ARITHMETIC_QUADDOUBLE
            case LatticeReduction::A_QuadDouble:   return (p <= (unsigned long)std::numeric_limits<long>::max()) && ((long)p <= std::numeric_limits<qd_real>::digits);
#else
            case LatticeReduction::A_QuadDouble:   break;
#endif
            }
            return false;
        }
        
        inline void setGramSchmidt(LatticeReduction::GramSchmidt gs)
        {
            if (d_gs == gs)
                return;
            clearInterface();
            d_gs = gs;
        }
        
        inline void setGramSchmidtRestart(bool gsr)
        {
            if (d_gsr == gsr)
                return;
            clearInterface();
            d_gsr = gsr;
        }
        
        inline LatticeReduction::GramSchmidt getGramSchmidt() const
        {
            return d_gs;
        }
        
        inline bool getGramSchmidtRestart() const
        {
            return d_gsr;
        }
        
        inline void forceGSRebuild(bool makeSureAllComputed)
        {
            setupInterface();
            d_if->forceGSRebuild(makeSureAllComputed);
        }
        
        inline double getGSCoefficientD(unsigned i, unsigned j) const
        {
            assert(i > j);
            setupInterface();
            return d_if->getGSCoefficientD(i, j);
        }
        
        inline long double getGSCoefficientLD(unsigned i, unsigned j) const
        {
            assert(i > j);
            setupInterface();
            return d_if->getGSCoefficientLD(i, j);
        }
        
        inline arithmetic::Real getGSCoefficientR(unsigned i, unsigned j, const arithmetic::RealContext & rc) const
        {
            assert(i > j);
            setupInterface();
            return d_if->getGSCoefficientR(i, j, rc);
        }
        
        inline double getGSSqNormD(unsigned i) const
        {
            setupInterface();
            return d_if->getGSSqNormD(i);
        }
        
        inline long double getGSSqNormLD(unsigned i) const
        {
            setupInterface();
            return d_if->getGSSqNormLD(i);
        }
        
        inline arithmetic::Real getGSSqNormR(unsigned i, const arithmetic::RealContext & rc) const
        {
            setupInterface();
            return d_if->getGSSqNormR(i, rc);
        }
        
        inline void modFlip(unsigned i)
        {
            if (i >= d_lattice.rows())
                return;
            if (hasInterface())
                d_if->modFlip(d_trans, i);
            else
                plll::linalg::neg(d_lattice.row(i), d_lattice.row(i));
        }
        
        inline void modSwap(unsigned i, unsigned j)
        {
            if ((i >= d_lattice.rows()) || (j >= d_lattice.rows()))
                return;
            if (i == j)
                return;
            if (hasInterface())
            {
                d_if->modSwap(d_trans, i, j);
            }
            else
                swap(d_lattice.row(i), d_lattice.row(j));
        }
        
        inline void modAdd(unsigned i, unsigned j, const arithmetic::Integer & lambda)
        {
            if ((i >= d_lattice.rows()) || (j >= d_lattice.rows()))
                return;
            if (hasInterface())
            {
                d_if->modAdd(d_trans, i, j, lambda);
            }
            else
                plll::linalg::addmul(d_lattice.row(i), d_lattice.row(j), lambda);
        }
        
        inline void enableTransform(LatticeReduction::Transform mode)
        {
            d_trans.enable(mode, d_lattice);
        }
        
        inline void disableTransform()
        {
            d_trans.disable();
        }
        
        inline bool isTransformationRecorded() const
        {
            return d_trans.isEnabled();
        }
        
        inline LatticeReduction::Transform getTransformationMode() const
        {
            return d_trans.mode();
        }
        
        inline const linalg::math_matrix<arithmetic::Integer> * getTransformation() const
        {
            return d_trans.transform();
        }
        
        inline void setSVPMode(LatticeReduction::SVPMode mode)
        {
            if (mode == d_svp)
                return;
            clearInterface();
            d_svp = mode;
        }
        
        inline LatticeReduction::SVPMode getSVPMode() const
        {
            return d_svp;
        }
        
        inline void setMaximalCoreUsage(unsigned cores)
        {
            if (cores == d_max_cores)
                return;
            clearInterface();
            d_max_cores = cores;
        }
        
        inline unsigned getMaximalCoreUsage()
        {
            return d_max_cores;
        }
        
        inline void setCallbackFunction(const LatticeReduction::CallbackFunction & cf, const LatticeReduction::CallbackFunction_LI & cf2)
        {
            d_cf = cf;
            d_cf2 = cf2;
        }
        
        inline void setCallbackInterval(double iv)
        {
            d_cf_int = iv;
        }
        
        inline std::pair<LatticeReduction::CallbackFunction, LatticeReduction::CallbackFunction_LI> getCallbackFunction() const
        {
            return std::make_pair(d_cf, d_cf2);
        }
        
        inline double getCallbackInterval() const
        {
            return d_cf_int;
        }
        
        inline void setMinCallbackFunction(const LatticeReduction::MinCallbackFunction & mcf, const LatticeReduction::MinCallbackFunction_LI & mcf2)
        {
            d_mcf = mcf;
            d_mcf2 = mcf2;
        }
        
        inline std::pair<LatticeReduction::MinCallbackFunction, LatticeReduction::MinCallbackFunction_LI> getMinCallbackFunction() const
        {
            return std::make_pair(d_mcf, d_mcf2);
        }
        
        inline void setEnumCallbackFunction(const LatticeReduction::EnumCallbackFunction & ecf, const LatticeReduction::EnumCallbackFunction_LI & ecf2)
        {
            d_ecf = ecf;
            d_ecf2 = ecf2;
        }
        
        inline std::pair<LatticeReduction::EnumCallbackFunction, LatticeReduction::EnumCallbackFunction_LI> getEnumCallbackFunction()
        {
            return std::make_pair(d_ecf, d_ecf2);
        }
        
        inline void setAnneal(bool lll, bool bkz)
        {
            d_anneal_lll = lll;
            d_anneal_bkz = bkz;
        }
        
        inline void setAnnealFun(const LatticeReduction::AnnealCallbackFunction & acf, const LatticeReduction::LLL_AnnealFunction & af_lll, const LatticeReduction::BKZ_AnnealFunction & af_bkz)
        {
            d_acf = acf;
            d_af_lll = af_lll;
            d_af_bkz = af_bkz;
        }
        
        inline const LatticeReduction::LLL_AnnealFunction & getAnnealLLL() const
        {
            return d_af_lll;
        }
        
        inline const LatticeReduction::BKZ_AnnealFunction & getAnnealBKZ() const
        {
            return d_af_bkz;
        }
        
        inline bool isAnnealLLL() const
        {
            return d_anneal_lll;
        }
        
        inline bool isAnnealBKZ() const
        {
            return d_anneal_bkz;
        }
        
        inline void setDI(LatticeReduction::DIMethod method)
        {
            d_di = method;
        }
        
        inline void setDIM(LatticeReduction::DIMode mode)
        {
            d_di_mode = mode;
        }
        
        inline void setDIC(LatticeReduction::DIChoice choice, unsigned bs)
        {
            d_di_choice = choice;
            d_di_bs = bs;
        }
        
        inline LatticeReduction::DIMethod getDI() const
        {
            return d_di;
        }
        
        inline LatticeReduction::DIMode getDIM() const
        {
            return d_di_mode;
        }
        
        inline LatticeReduction::DIChoice getDIC() const
        {
            return d_di_choice;
        }
        
        inline unsigned getDIBS() const
        {
            return d_di_bs;
        }
        
        void setRange(unsigned begin, unsigned end)
        {
            d_begin = begin;
            d_end = end;
            if (d_begin >= d_lattice.rows())
                d_begin = 0;
            if (d_begin > d_end)
                d_end = d_begin;
            if (d_end >= d_lattice.rows())
                d_end = d_lattice.rows() ? d_lattice.rows() - 1 : 0;
        }
        
        inline std::pair<unsigned, unsigned> getRange() const
        {
            return std::make_pair(d_begin, d_end);
        }
        
        inline void sort(bool projected)
        {
            if (d_lattice.rows() == 0)
                return;
            if (projected)
            {
                setupInterface();
                d_if->sortProjected(d_begin, d_end, d_trans);
            }
            else
            {
                // Implement a simple insertion sort
                linalg::base_rowvector<arithmetic::Integer> norms(d_lattice.rows());
                for (unsigned i = d_begin; i <= d_end; ++i)
                    normSq(norms[i], d_lattice.row(i));
                for (unsigned i = d_begin; i < d_end; ++i)
                {
                    unsigned mini = i;
                    for (unsigned j = i + 1; j <= d_end; ++j)
                        if (norms[j] < norms[mini])
                            mini = j;
                    for (unsigned j = mini; j > i; --j)
                    {
                        modSwap(j - 1, j);
                        swap(norms[j - 1], norms[j]);
                    }
                }
            }
        }
        
        inline void sizereduction()
        {
            if (d_lattice.rows() == 0)
                return;
            setupInterface();
            d_if->sizereduction(d_begin, d_end, d_trans);
        }
        
        inline void lll(double alpha, LatticeReduction::LLLMode mode)
        {
            if (d_lattice.rows() == 0)
                return;
            setupInterface();
            d_if->lll(d_begin, d_end, d_trans, alpha, mode, d_cf, d_cf2, d_cf_int, d_mcf, d_mcf2,
                      NULL, d_anneal_lll, d_acf, d_af_lll, d_di, d_di_mode, d_di_choice, d_di_bs);
        }
        
        inline void bkz(double alpha, unsigned blocksize, LatticeReduction::BKZMode mode)
        {
            if (d_lattice.rows() == 0)
                return;
            setupInterface();
            d_if->bkz(d_begin, d_end, d_trans, alpha, blocksize, mode, d_cf, d_cf2, d_cf_int, d_mcf, d_mcf2,
                      NULL, d_ecf, d_ecf2, d_anneal_bkz, d_acf, d_af_bkz, d_di, d_di_mode, d_di_choice, d_di_bs);
        }
        
        inline void hkz(bool dual)
        {
            if (d_lattice.rows() == 0)
                return;
            setupInterface();
            d_if->hkz(d_begin, d_end, d_trans, dual, d_cf, d_cf2, d_cf_int, d_mcf, d_mcf2, NULL, d_ecf, d_ecf2);
        }
        
        inline void svp(bool make_basis, bool extreme, bool dual)
        {
            if (d_lattice.rows() == 0)
                return;
            setupInterface();
            d_if->svp(d_begin, d_end, d_trans, make_basis, extreme, dual, d_cf, d_cf2, d_cf_int, d_mcf, d_mcf2, NULL, d_ecf, d_ecf2);
        }
        
        inline bool isSizeReduced() const
        {
            if (d_lattice.rows() == 0)
                return true;
            setupInterface();
            return d_if->isSizeReduced(d_begin, d_end);
        }
        
        inline bool isLLLBasis(double alpha, LatticeReduction::LLLMode mode) const
        {
            if (d_lattice.rows() == 0)
                return true;
            setupInterface();
            return d_if->isLLLBasis(d_begin, d_end, alpha, mode, d_di, d_di_choice, d_di_bs);
        }
        
        inline bool isBKZBasis(double alpha, unsigned blocksize, LatticeReduction::BKZMode mode) const
        {
            if (d_lattice.rows() == 0)
                return true;
            setupInterface();
            return d_if->isBKZBasis(d_begin, d_end, alpha, blocksize, mode, d_di, d_di_choice, d_di_bs);
        }
        
        inline bool isHKZBasis(bool dual) const
        {
            if (d_lattice.rows() == 0)
                return true;
            setupInterface();
            return d_if->isHKZBasis(d_begin, d_end, dual);
        }
        
        inline bool isSVPBasis(bool dual) const
        {
            if (d_lattice.rows() == 0)
                return true;
            setupInterface();
            return d_if->isSVPBasis(d_begin, d_end, dual);
        }
        
        const LatticeReduction::Statistics & getStatistics() const
        {
            setupInterface();
            return d_if->getStatistics();
        }
        
        void resetStatistics()
        {
            if (d_if.get())
                d_if->resetStatistics();
        }
        
        void setVerbose(LatticeReduction::VerboseOutputLevel l, const LatticeReduction::VerboseFunction & f)
        {
            d_verboseoutputlevel = l;
            if (f)
                d_verbosefunction = f;
            if (d_if.get())
                d_if->setupVerbose(d_verboseoutputlevel, d_verbosefunction);
        }
        
        LatticeReduction::VerboseOutputLevel getVerboseLevel() const
        {
            return d_verboseoutputlevel;
        }
        
        const LatticeReduction::VerboseFunction & getVerboseFunction() const
        {
            return d_verbosefunction;
        }
    };
    
    LatticeReduction::LatticeReduction()
        : d_impl(new LatticeReductionImpl())
    {
    }
    
    LatticeReduction::LatticeReduction(const linalg::math_matrix<arithmetic::Integer> & lattice)
        : d_impl(new LatticeReductionImpl())
    {
        d_impl->setLattice(lattice);
    }
    
    template<typename IType>
    LatticeReduction::LatticeReduction(const linalg::math_matrix<arithmetic::NInt<IType> > & lattice)
        : d_impl(new LatticeReductionImpl())
    {
        d_impl->setLattice(lattice);
    }
    
    template LatticeReduction::LatticeReduction(const linalg::math_matrix<arithmetic::NInt<int> > &);
    template LatticeReduction::LatticeReduction(const linalg::math_matrix<arithmetic::NInt<long int> > &);
    template LatticeReduction::LatticeReduction(const linalg::math_matrix<arithmetic::NInt<long long> > &);
    
#if __cplusplus >= 201103L
    LatticeReduction::LatticeReduction(linalg::math_matrix<arithmetic::Integer> && lattice)
        : d_impl(new LatticeReductionImpl())
    {
        d_impl->setLattice(std::move(lattice));
    }
    
    LatticeReduction::LatticeReduction(LatticeReduction && lr)
        : d_impl(lr.d_impl)
    {
    }
    
    LatticeReduction & LatticeReduction::operator = (LatticeReduction && lr)
    {
        if (this != &lr)
            d_impl = lr.d_impl;
        return *this;
    }
#endif
    
    LatticeReduction::~LatticeReduction()
    {
    }
    
    void LatticeReduction::setLattice(const linalg::math_matrix<arithmetic::Integer> & lattice)
    {
        d_impl->setLattice(lattice);
    }
    
    template<typename IType>
    void LatticeReduction::setLattice(const linalg::math_matrix<arithmetic::NInt<IType> > & lattice)
    {
        d_impl->setLattice(lattice);
    }
    
    template void LatticeReduction::setLattice(const linalg::math_matrix<arithmetic::NInt<int> > &);
    template void LatticeReduction::setLattice(const linalg::math_matrix<arithmetic::NInt<long int> > &);
    template void LatticeReduction::setLattice(const linalg::math_matrix<arithmetic::NInt<long long> > &);
    
#if __cplusplus >= 201103L
    void LatticeReduction::setLattice(linalg::math_matrix<arithmetic::Integer> && lattice)
    {
        d_impl->setLattice(std::move(lattice));
    }
#endif
    
    const linalg::math_matrix<arithmetic::Integer> & LatticeReduction::getLattice() const
    {
        return d_impl->getLattice();
    }
    
    unsigned LatticeReduction::rank() const
    {
        return d_impl->getLattice().rows();
    }
    
    unsigned LatticeReduction::dimension() const
    {
        return d_impl->getLattice().cols();
    }
    
    void LatticeReduction::setArithmetic(LatticeReduction::Arithmetic a)
    {
        d_impl->setArithmetic(a);
    }
    
    LatticeReduction::Arithmetic LatticeReduction::getArithmetic() const
    {
        return d_impl->getArithmetic();
    }
    
    bool LatticeReduction::ensurePrecision(unsigned long p)
    // Tries to ensure a minimum precision. If requirement can be realized, returns true, otherwise false.
    {
        return d_impl->ensurePrecision(p);
    }
    
    void LatticeReduction::setIntegers(LatticeReduction::Integers a)
    {
        d_impl->setIntegers(a);
    }
    
    LatticeReduction::Integers LatticeReduction::getIntegers() const
    {
        return d_impl->getIntegers();
    }
    
    void LatticeReduction::setGramSchmidt(LatticeReduction::GramSchmidt gs)
    {
        d_impl->setGramSchmidt(gs);
    }
    
    void LatticeReduction::setGramSchmidtRestart(bool r)
    {
        d_impl->setGramSchmidtRestart(r);
    }
    
    LatticeReduction::GramSchmidt LatticeReduction::getGramSchmidt() const
    {
        return d_impl->getGramSchmidt();
    }
    
    bool LatticeReduction::getGramSchmidtRestart() const
    {
        return d_impl->getGramSchmidtRestart();
    }
    
    void LatticeReduction::forceGSRebuild(bool makeSureAllComputed)
    {
        d_impl->forceGSRebuild(makeSureAllComputed);
    }
    
    double LatticeReduction::getGSCoefficientD(unsigned i, unsigned j) const
    {
        return d_impl->getGSCoefficientD(i, j);
    }
    
    long double LatticeReduction::getGSCoefficientLD(unsigned i, unsigned j) const
    {
        return d_impl->getGSCoefficientLD(i, j);
    }
    
    arithmetic::Real LatticeReduction::getGSCoefficientR(unsigned i, unsigned j) const
    {
        return d_impl->getGSCoefficientR(i, j, arithmetic::getThreadRealContext());
    }
    
    arithmetic::Real LatticeReduction::getGSCoefficientR(unsigned i, unsigned j, const arithmetic::RealContext & rc) const
    {
        return d_impl->getGSCoefficientR(i, j, rc);
    }
    
    double LatticeReduction::getGSSqNormD(unsigned i) const
    {
        return d_impl->getGSSqNormD(i);
    }
    
    long double LatticeReduction::getGSSqNormLD(unsigned i) const
    {
        return d_impl->getGSSqNormLD(i);
    }
    
    arithmetic::Real LatticeReduction::getGSSqNormR(unsigned i) const
    {
        return d_impl->getGSSqNormR(i, arithmetic::getThreadRealContext());
    }
    
    arithmetic::Real LatticeReduction::getGSSqNormR(unsigned i, const arithmetic::RealContext & rc) const
    {
        return d_impl->getGSSqNormR(i, rc);
    }
    
    void LatticeReduction::modFlip(unsigned i)
    {
        d_impl->modFlip(i);
    }
    
    void LatticeReduction::modSwap(unsigned i, unsigned j)
    {
        d_impl->modSwap(i, j);
    }
    
    void LatticeReduction::modAdd(unsigned i, unsigned j, const arithmetic::Integer & lambda)
    {
        d_impl->modAdd(i, j, lambda);
    }
    
    void LatticeReduction::enableTransform(LatticeReduction::Transform t)
    {
        d_impl->enableTransform(t);
    }
    
    void LatticeReduction::disableTransform()
    {
        d_impl->disableTransform();
    }
    
    bool LatticeReduction::isTransformationRecorded() const
    {
        return d_impl->isTransformationRecorded();
    }
    
    LatticeReduction::Transform LatticeReduction::getTransformationMode() const
    {
        return d_impl->getTransformationMode();
    }
    
    const linalg::math_matrix<arithmetic::Integer> * LatticeReduction::getTransformation() const
    {
        return d_impl->getTransformation();
    }
    
    void LatticeReduction::setSVPMode(LatticeReduction::SVPMode m)
    {
        d_impl->setSVPMode(m);
    }
    
    LatticeReduction::SVPMode LatticeReduction::getSVPMode() const
    {
        return d_impl->getSVPMode();
    }
    
    void LatticeReduction::setMaximalCoreUsage(unsigned m)
    {
        d_impl->setMaximalCoreUsage(m);
    }
    
    unsigned LatticeReduction::getMaximalCoreUsage()
    {
        return d_impl->getMaximalCoreUsage();
    }
    
    void LatticeReduction::setCallbackFunction(const LatticeReduction::CallbackFunction & f, const LatticeReduction::CallbackFunction_LI & f2)
    {
        d_impl->setCallbackFunction(f, f2);
    }
    
    void LatticeReduction::setCallbackFunction(const LatticeReduction::CallbackFunction_LI & f)
    {
        d_impl->setCallbackFunction(LatticeReduction::CallbackFunction(), f);
    }
    
    void LatticeReduction::setCallbackInterval(double iv)
    {
        d_impl->setCallbackInterval(iv);
    }
    
    std::pair<LatticeReduction::CallbackFunction, LatticeReduction::CallbackFunction_LI> LatticeReduction::getCallbackFunction() const
    {
        return d_impl->getCallbackFunction();
    }
    
    double LatticeReduction::getCallbackInterval() const
    {
        return d_impl->getCallbackInterval();
    }
    
    void LatticeReduction::setMinCallbackFunction(const LatticeReduction::MinCallbackFunction & f, const LatticeReduction::MinCallbackFunction_LI & f2)
    {
        d_impl->setMinCallbackFunction(f, f2);
    }
    
    void LatticeReduction::setMinCallbackFunction(const LatticeReduction::MinCallbackFunction_LI & f)
    {
        d_impl->setMinCallbackFunction(LatticeReduction::MinCallbackFunction(), f);
    }
    
    std::pair<LatticeReduction::MinCallbackFunction, LatticeReduction::MinCallbackFunction_LI> LatticeReduction::getMinCallbackFunction() const
    {
        return d_impl->getMinCallbackFunction();
    }
    
    void LatticeReduction::setEnumCallbackFunction(const LatticeReduction::EnumCallbackFunction & ecf, const LatticeReduction::EnumCallbackFunction_LI & ecf2)
    {
        d_impl->setEnumCallbackFunction(ecf, ecf2);
    }
    
    void LatticeReduction::setEnumCallbackFunction(const LatticeReduction::EnumCallbackFunction_LI & ecf)
    {
        d_impl->setEnumCallbackFunction(LatticeReduction::EnumCallbackFunction(), ecf);
    }
    
    std::pair<LatticeReduction::EnumCallbackFunction, LatticeReduction::EnumCallbackFunction_LI> LatticeReduction::getEnumCallbackFunction()
    {
        return d_impl->getEnumCallbackFunction();
    }
    
    void LatticeReduction::setDefaultAnnealing()
    {
        d_impl->setAnnealFun(NULL, NULL, NULL);
        d_impl->setAnneal(true, true);
    }
    
    void LatticeReduction::setDefaultAnnealingLLL()
    {
        d_impl->setAnnealFun(NULL, NULL, d_impl->getAnnealBKZ());
        d_impl->setAnneal(true, d_impl->isAnnealBKZ());
    }
    
    void LatticeReduction::setDefaultAnnealingBKZ()
    {
        d_impl->setAnnealFun(NULL, d_impl->getAnnealLLL(), NULL);
        d_impl->setAnneal(d_impl->isAnnealLLL(), true);
    }
    
    void LatticeReduction::setAnnealing(const LatticeReduction::AnnealCallbackFunction & f, const LatticeReduction::LLL_AnnealFunction & lf, const LatticeReduction::BKZ_AnnealFunction & bf)
    {
        d_impl->setAnnealFun(f, lf, bf);
        d_impl->setAnneal(true, true);
    }
    
    void LatticeReduction::setAnnealingLLL(const LatticeReduction::AnnealCallbackFunction & f, const LatticeReduction::LLL_AnnealFunction & lf)
    {
        d_impl->setAnnealFun(f, lf, d_impl->getAnnealBKZ());
        d_impl->setAnneal(true, d_impl->isAnnealBKZ());
    }
    
    void LatticeReduction::setAnnealingBKZ(const LatticeReduction::AnnealCallbackFunction & f, const LatticeReduction::BKZ_AnnealFunction & bf)
    {
        d_impl->setAnnealFun(f, d_impl->getAnnealLLL(), bf);
        d_impl->setAnneal(d_impl->isAnnealLLL(), true);
    }
    
    void LatticeReduction::disableAnnealing()
    {
        d_impl->setAnneal(false, false);
    }
    
    void LatticeReduction::disableAnnealingLLL()
    {
        d_impl->setAnneal(false, d_impl->isAnnealLLL());
    }
    
    void LatticeReduction::disableAnnealingBKZ()
    {
        d_impl->setAnneal(d_impl->isAnnealLLL(), false);
    }
    
    bool LatticeReduction::isAnnealingLLLEnabled() const
    {
        return d_impl->isAnnealLLL();
    }
    
    bool LatticeReduction::isAnnealingBKZEnabled() const
    {
        return d_impl->isAnnealBKZ();
    }
    
    void LatticeReduction::setDeepInsertionMethod(LatticeReduction::DIMethod method, LatticeReduction::DIMode mode)
    {
        d_impl->setDI(method);
        d_impl->setDIM(mode);
    }
    
    void LatticeReduction::setDeepInsertionMode(LatticeReduction::DIMode mode)
    {
        d_impl->setDIM(mode);
    }
    
    void LatticeReduction::setDeepInsertionChoice(LatticeReduction::DIChoice choice, unsigned blocksize)
    {
        d_impl->setDIC(choice, blocksize);
    }
    
    LatticeReduction::DIMethod LatticeReduction::getDeepInsertionMethod() const
    {
        return d_impl->getDI();
    }
    
    LatticeReduction::DIMode LatticeReduction::getDeepInsertionMode() const
    {
        return d_impl->getDIM();
    }
    
    LatticeReduction::DIChoice LatticeReduction::getDeepInsertionChoice() const
    {
        return d_impl->getDIC();
    }
    
    unsigned LatticeReduction::getDeepInsertionBlocksize() const
    {
        return d_impl->getDIBS();
    }
    
    void LatticeReduction::setRange(unsigned begin, unsigned end)
    {
        d_impl->setRange(begin, end);
    }
    
    std::pair<unsigned, unsigned> LatticeReduction::getRange() const
    {
        return d_impl->getRange();
    }
    
    void LatticeReduction::sort(bool projected)
    {
        d_impl->sort(projected);
    }
    
    void LatticeReduction::sizereduction()
    {
        d_impl->sizereduction();
    }
    
    void LatticeReduction::lll(double alpha, LatticeReduction::LLLMode mode)
    {
        d_impl->lll(alpha, mode);
    }
    
    void LatticeReduction::bkz(double alpha, unsigned blocksize, LatticeReduction::BKZMode mode)
    {
        d_impl->bkz(alpha, blocksize, mode);
    }
    
    void LatticeReduction::hkz(bool dual)
    {
        d_impl->hkz(dual);
    }
    
    void LatticeReduction::svp(bool make_basis, bool extreme, bool dual)
    {
        d_impl->svp(make_basis, extreme, dual);
    }
    
    bool LatticeReduction::isSizeReduced() const
    {
        return d_impl->isSizeReduced();
    }
    
    bool LatticeReduction::isLLLBasis(double alpha, LatticeReduction::LLLMode mode) const
    {
        return d_impl->isLLLBasis(alpha, mode);
    }
    
    bool LatticeReduction::isBKZBasis(double alpha, unsigned blocksize, LatticeReduction::BKZMode mode) const
    {
        return d_impl->isBKZBasis(alpha, blocksize, mode);
    }
    
    bool LatticeReduction::isHKZBasis(bool dual) const
    {
        return d_impl->isHKZBasis(dual);
    }
    
    bool LatticeReduction::isSVPBasis(bool dual) const
    {
        return d_impl->isSVPBasis(dual);
    }
    
    const LatticeReduction::Statistics & LatticeReduction::getStatistics() const
    {
        return d_impl->getStatistics();
    }
    
    void LatticeReduction::resetStatistics()
    {
        d_impl->resetStatistics();
    }
    
    void LatticeReduction::setVerbose(VerboseOutputLevel l, const VerboseFunction & f)
    {
        d_impl->setVerbose(l, f);
    }
    
    LatticeReduction::VerboseOutputLevel LatticeReduction::getVerboseOutputLevel()
    {
        return d_impl->getVerboseLevel();
    }
    
    const LatticeReduction::VerboseFunction & LatticeReduction::getVerboseFunction()
    {
        return d_impl->getVerboseFunction();
    }
    
    template<class IntTypeContext>
    bool DefaultLLLAnnealFunction(IntTypeContext & ic, linalg::math_matrix<typename IntTypeContext::Integer> & A, int k,
                                  arithmetic::RealContext & rc, arithmetic::RandomNumberGenerator & rng, arithmetic::Real & T)
    {
        arithmetic::Real prob = arithmetic::convert(normSq(A.row(k - 1)), rc) / arithmetic::convert(normSq(A.row(k)), rc);
//    prob *= arithmetic::convert(k, rc) / arithmetic::convert(A.rows(), rc);
        if (prob > arithmetic::convert(1, rc))
            setOne(prob);
        prob *= T;
        arithmetic::Real rnd;
        rng.randomUniform(rnd);
//            std::cout << prob << " " << rnd << " " << arithmetic::convert<arithmetic::Real>(normSq(A.row(k - 1))) / arithmetic::convert<arithmetic::Real>(normSq(A.row(k))) << "\n";
        return rnd <= prob;
    }
    
    bool DefaultAnnealCallbackFunction(const linalg::math_matrix<arithmetic::Integer> & A, arithmetic::RealContext & arc, arithmetic::Real & T)
    {
        if (sign(T) < 0)
        {
            T = arithmetic::convert(0.36, arc);
            return true;
        }
        else
        {
            if (T == arithmetic::convert(0.36, arc))
                T = arithmetic::convert(0.30, arc);
            else if (T == arithmetic::convert(0.30, arc))
                T = arithmetic::convert(0.20, arc);
            else if (T == arithmetic::convert(0.20, arc))
                T = arithmetic::convert(0.10, arc);
            else if (T == arithmetic::convert(0.10, arc))
                T = arithmetic::convert(0.05, arc);
            else if (T == arithmetic::convert(0.05, arc))
                return false;
            return true;
        }
    }
    
    template<class IntTypeContext>
    bool DefaultBKZAnnealFunction(IntTypeContext & ic, linalg::math_matrix<typename IntTypeContext::Integer> & A, int k, int windowsize,
                                  linalg::math_rowvector<typename IntTypeContext::Integer> & lincomb,
                                  arithmetic::RealContext & rc, arithmetic::RandomNumberGenerator & rng, arithmetic::Real & T)
    {
        if (windowsize < 2)
            return false;
    
        arithmetic::Real prob = arithmetic::convert(normSq(A.row(k)), rc) / arithmetic::convert(normSq(A.row(k + 1)), rc);
        if (prob > arithmetic::convert(1, rc))
            setOne(prob);
        prob *= T;
        arithmetic::Real rnd;
        rng.randomUniform(rnd);
//            std::cout << prob << " " << rnd << " " << arithmetic::convert<arithmetic::Real>(normSq(lattice.row(k - 1))) / arithmetic::convert<arithmetic::Real>(normSq(lattice.row(k))) << "\n";
    
        setZero(lincomb[0]);
        setOne(lincomb[1]);
        for (unsigned i = 2; i < lincomb.size(); ++i)
            setZero(lincomb[i]);
        return rnd <= prob;
    }
}

// Enforce template instantiations

#ifndef PLLL_CONFIG_NO_ARITHMETIC_BIGINT
namespace plll
{
    template bool DefaultLLLAnnealFunction<arithmetic::IntegerContext>(arithmetic::IntegerContext & ic,
                                                                       linalg::math_matrix<arithmetic::Integer> &,
                                                                       int,
                                                                       arithmetic::RealContext &,
                                                                       arithmetic::RandomNumberGenerator &,
                                                                       arithmetic::Real &);
    template bool DefaultBKZAnnealFunction<arithmetic::IntegerContext>(arithmetic::IntegerContext & ic,
                                                                       linalg::math_matrix<arithmetic::Integer> &,
                                                                       int, int,
                                                                       linalg::math_rowvector<arithmetic::Integer> &,
                                                                       arithmetic::RealContext &,
                                                                       arithmetic::RandomNumberGenerator &,
                                                                       arithmetic::Real &);
}
#endif

#ifndef PLLL_CONFIG_NO_ARITHMETIC_LONGINT
namespace plll
{
    template bool DefaultLLLAnnealFunction<arithmetic::NIntContext<long> >(arithmetic::NIntContext<long> & ic,
                                                                           linalg::math_matrix<arithmetic::NInt<long> > &,
                                                                           int,
                                                                           arithmetic::RealContext &,
                                                                           arithmetic::RandomNumberGenerator &,
                                                                           arithmetic::Real &);
    template bool DefaultBKZAnnealFunction<arithmetic::NIntContext<long> >(arithmetic::NIntContext<long> & ic,
                                                                           linalg::math_matrix<arithmetic::NInt<long> > &,
                                                                           int, int,
                                                                           linalg::math_rowvector<arithmetic::NInt<long> > &,
                                                                           arithmetic::RealContext &,
                                                                           arithmetic::RandomNumberGenerator &,
                                                                           arithmetic::Real &);
}
#endif

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

#ifndef PLLL_INCLUDE_GUARD__LLL2_GENERIC_CPP
#define PLLL_INCLUDE_GUARD__LLL2_GENERIC_CPP

#include "lll2-internal.hpp"

#include "lattice.cpp"
#include "lllimpl.cpp"
#include "enumimpl.hpp"
#include "enumimpl-preproc.cpp"
#include "bkzimpl.cpp"
#include "svpimpl.cpp"
#include "verifyimpl.cpp"
#include "matrixconversion.hpp"
#include "lll2-callback.cpp"
#include "lll2-deepinsertion.cpp"

#include <boost/bind.hpp>

namespace plll
{
    template<class RealTypeContext, class IntTypeContext>
    GSInterface<RealTypeContext, IntTypeContext> * createGSInterface(Verbose & verbose, RealTypeContext & rc, IntTypeContext & ic,
                                                                     LatticeReduction::GramSchmidt gs, bool restarts,
                                                                     linalg::math_matrix<typename IntTypeContext::Integer> & A);
    
    template<class RealTypeContext, class IntTypeContext>
    class LRIImplementation : public LRIInterface
    {
        class GSIImplementation;
        friend class GSIImplementation;
        
    private:
        mutable RealTypeContext d_rc;
        mutable IntTypeContext d_ic;
        MatrixConversion<IntTypeContext> d_intlattice;
        Verbose d_verbose;
        mutable STD_AUTO_PTR<GSInterface<RealTypeContext, IntTypeContext> > d_gs;
        mutable STD_AUTO_PTR<Lattice<RealTypeContext, IntTypeContext> > d_lattice_object;
        mutable Transform * d_lattice_object_trans;
        LatticeReduction::SVPMode d_svp;
        unsigned d_max_cores;
        
        void handle_exception(std::exception * e) const
        {
            if (reduction_error * ee = dynamic_cast<reduction_error*>(e))
            {
                d_verbose(LatticeReduction::VL_Error) << "FATAL ERROR during reduction: " << ee->what();
            }
            else if (change_interface_exception * ee = dynamic_cast<change_interface_exception*>(e))
            {
                d_verbose(LatticeReduction::VL_Error) << ee->what();
            }
            else if (feature_not_implemented * ee = dynamic_cast<feature_not_implemented*>(e))
            {
                d_verbose(LatticeReduction::VL_Error) << ee->what();
            }
            else if (dynamic_cast<LatticeReduction::stop_reduction*>(e) != NULL)
            {
                d_verbose(LatticeReduction::VL_Information) << "Stopping reduction requested";
            }
            else if (dynamic_cast<LatticeReduction::stop_enumeration*>(e) != NULL)
            {
                d_verbose(LatticeReduction::VL_Error) << "FATAL ERROR: Enumeration stopping request caught outside enumeration code.";
            }
            else if (dynamic_cast<std::bad_alloc*>(e) != NULL)
            {
                d_verbose(LatticeReduction::VL_Error) << "FATAL ERROR: Memory allocation failed!";
            }
            else
                d_verbose(LatticeReduction::VL_Error) << "FATAL ERROR: Unhandled standard exception: " << e->what();
        }
        
        void handle_unknown_exception() const
        {
            d_verbose(LatticeReduction::VL_Error) << "FATAL ERROR: Unknown exception! Re-throwing.";
        }
        
        void setupLattice(Transform * transform) const
        {
            if (d_lattice_object.get() == NULL)
                d_lattice_object.reset(new Lattice<RealTypeContext, IntTypeContext>(d_gs.get(), d_rc, d_ic, d_stats, 0, d_gs->getDimension() - 1));
            else
            {
                if (transform == d_lattice_object_trans)
                    return;
            }
            if (d_lattice_object_trans)
            {
                d_lattice_object->removeAllNotifiers();
                d_lattice_object_trans->releaseTransformObject();
                d_lattice_object_trans = NULL;
            }
            if (transform)
            {
                d_lattice_object_trans = transform;
                d_lattice_object->addNotifier(transform->createTransformObject<IntTypeContext>(d_ic));
            }
        }
        
        void freeLattice() const
        {
            if (d_lattice_object_trans)
            {
                if (d_lattice_object.get())
                    d_lattice_object->removeAllNotifiers();
                d_lattice_object_trans->releaseTransformObject();
                d_lattice_object_trans = NULL;
            }
            d_lattice_object.reset();
        }
        
        class GSIImplementation : public LatticeReduction::GramSchmidtInformer
        {
        private:
            LRIImplementation & d_i;
            
        public:
            GSIImplementation(LRIImplementation & i)
                : d_i(i)
            {
            }
            
            virtual ~GSIImplementation()
            {
            }
            
            virtual double getGSCoefficientD(unsigned i, unsigned j) const
            {
                d_i.setupLattice(NULL);
                d_i.d_lattice_object->update(i + 1);
                return arithmetic::convert<double>(d_i.d_lattice_object->getCoeff(i, j));
            }
            
            virtual long double getGSCoefficientLD(unsigned i, unsigned j) const
            {
                d_i.setupLattice(NULL);
                d_i.d_lattice_object->update(i + 1);
                return arithmetic::convert<long double>(d_i.d_lattice_object->getCoeff(i, j));
            }
            
            virtual arithmetic::Real getGSCoefficientR(unsigned i, unsigned j, const arithmetic::RealContext & rc) const
            {
                d_i.setupLattice(NULL);
                d_i.d_lattice_object->update(i + 1);
                return arithmetic::convert(d_i.d_lattice_object->getCoeff(i, j), rc);
            }
            
            virtual double getGSSqNormD(unsigned i) const
            {
                d_i.setupLattice(NULL);
                d_i.d_lattice_object->update(i + 1);
                return arithmetic::convert<double>(d_i.d_lattice_object->getNormSq(i));
            }
            
            virtual long double getGSSqNormLD(unsigned i) const
            {
                d_i.setupLattice(NULL);
                d_i.d_lattice_object->update(i + 1);
                return arithmetic::convert<long double>(d_i.d_lattice_object->getNormSq(i));
            }
            
            virtual arithmetic::Real getGSSqNormR(unsigned i, const arithmetic::RealContext & rc) const
            {
                d_i.setupLattice(NULL);
                d_i.d_lattice_object->update(i + 1);
                return arithmetic::convert(d_i.d_lattice_object->getNormSq(i), rc);
            }
            
            virtual double computeProjectionLengthD(unsigned k, unsigned b,
                                                    const linalg::math_rowvector<arithmetic::Integer> & vec) const
            {
                d_i.setupLattice(NULL);
                typename RealTypeContext::Real result(d_i.d_rc);
                VectorConversionConst<IntTypeContext> v(d_i.d_ic, vec);
                d_i.d_lattice_object->computeProjectionLength(result, k, b, v.vector());
                return arithmetic::convert<double>(result);
            }
            
            virtual long double computeProjectionLengthLD(unsigned k, unsigned b,
                                                          const linalg::math_rowvector<arithmetic::Integer> & vec) const
            {
                d_i.setupLattice(NULL);
                typename RealTypeContext::Real result(d_i.d_rc);
                VectorConversionConst<IntTypeContext> v(d_i.d_ic, vec);
                d_i.d_lattice_object->computeProjectionLength(result, k, b, v.vector());
                return arithmetic::convert<long double>(result);
            }
            
            virtual arithmetic::Real computeProjectionLengthR(unsigned k, unsigned b,
                                                              const linalg::math_rowvector<arithmetic::Integer> & vec,
                                                              const arithmetic::RealContext & rc) const
            {
                d_i.setupLattice(NULL);
                typename RealTypeContext::Real result(d_i.d_rc);
                VectorConversionConst<IntTypeContext> v(d_i.d_ic, vec);
                d_i.d_lattice_object->computeProjectionLength(result, k, b, v.vector());
                return arithmetic::convert(result, rc);
            }
        };
        
        GSIImplementation d_gsi;
        
    public:
        LRIImplementation(LatticeReduction::VerboseOutputLevel vol, LatticeReduction::VerboseFunction vf, linalg::math_matrix<arithmetic::Integer> & lattice,
                          LatticeReduction::GramSchmidt gs, bool gsr, LatticeReduction::SVPMode svp, unsigned max_cores,
                          const RealTypeContext & rc, const IntTypeContext & ic)
            : LRIInterface(lattice), d_rc(rc), d_ic(ic), d_intlattice(d_ic, d_lattice, true, true), d_verbose(vol, vf),
              d_gs(createGSInterface(d_verbose, d_rc, d_ic, gs, gsr, d_intlattice.matrix())),
              d_lattice_object(), d_lattice_object_trans(NULL), d_svp(svp), d_max_cores(max_cores), d_gsi(*this)
        {
        }
        
        virtual ~LRIImplementation()
        {
        }
        
        virtual void setupVerbose(LatticeReduction::VerboseOutputLevel vol, LatticeReduction::VerboseFunction vf)
        {
            d_verbose.setup(vol, vf);
        }
        
        virtual void ensureMinimumPrecision(unsigned long p)
        {
            if (RealTypeContext::is_variable_precision && (p > d_rc.getRealPrecision()))
                d_rc.setRealPrecision((p + 7) & ~7); // make multiple of 8
        }
        
        virtual const LatticeReduction::GramSchmidtInformer * getInformer() const
        {
            return &d_gsi;
        }
        
        virtual void forceGSRebuild(bool makeSureAllComputed)
        {
            if (d_lattice_object.get())
            {
                d_lattice_object->reset();
                if (makeSureAllComputed)
                    d_lattice_object->update(d_lattice_object->dimension());
            }
        }
        
        virtual double getGSCoefficientD(unsigned i, unsigned j) const
        {
            setupLattice(NULL);
            d_lattice_object->update(i + 1);
            return arithmetic::convert<double>(d_lattice_object->getCoeff(i, j));
        }
        
        virtual long double getGSCoefficientLD(unsigned i, unsigned j) const
        {
            setupLattice(NULL);
            d_lattice_object->update(i + 1);
            return arithmetic::convert<long double>(d_lattice_object->getCoeff(i, j));
        }
        
        virtual arithmetic::Real getGSCoefficientR(unsigned i, unsigned j, const arithmetic::RealContext & rc) const
        {
            setupLattice(NULL);
            d_lattice_object->update(i + 1);
            return arithmetic::convert(d_lattice_object->getCoeff(i, j), rc);
        }
        
        virtual double getGSSqNormD(unsigned i) const
        {
            setupLattice(NULL);
            d_lattice_object->update(i + 1);
            return arithmetic::convert<double>(d_lattice_object->getNormSq(i));
        }
        
        virtual long double getGSSqNormLD(unsigned i) const
        {
            setupLattice(NULL);
            d_lattice_object->update(i + 1);
            return arithmetic::convert<long double>(d_lattice_object->getNormSq(i));
        }
        
        virtual arithmetic::Real getGSSqNormR(unsigned i, const arithmetic::RealContext & rc) const
        {
            setupLattice(NULL);
            d_lattice_object->update(i + 1);
            return arithmetic::convert(d_lattice_object->getNormSq(i), rc);
        }
        
        virtual void modFlip(Transform & transform, unsigned i)
        {
            setupLattice(&transform);
            d_lattice_object->flip(i);
        }
        
        virtual void modSwap(Transform & transform, unsigned i, unsigned j)
        {
            setupLattice(&transform);
            d_lattice_object->swap(i, j);
        }
        
        virtual void modAdd(Transform & transform, unsigned i, unsigned j, const arithmetic::Integer & lambda)
        {
            setupLattice(&transform);
            d_lattice_object->add(i, arithmetic::convert(lambda, d_ic), j);
        }
        
        template<class Enumerator, class CallbackFunction>
        static void removeZeroVectors(Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction> & workspace,
                                      Lattice<RealTypeContext, IntTypeContext> & lattice)
        {
            unsigned r = 0, c = 0;
            while (r < lattice.dimension())
            {
                if (lattice.isRowZero(r))
                {
                    lattice.removeZeroVector(r);
                    ++c;
                }
                else
                    ++r;
            }
            if (c > 0)
                if (workspace.verbose().yieldsOutput(LatticeReduction::VL_Information))
                    workspace.verbose()(LatticeReduction::VL_Information) << "Removed " << c << " additional zero vectors!";
            lattice.compactify();
        }
        
        static void setupPruning(StandardEnumSettings & settings, bool enable, double p = 0.5, bool extreme = false)
        {
            if (enable)
            {
                if (extreme)
                    settings.d_pruning = StandardEnumSettings::PM_GNR110_Poly8;
                else
                {
                    settings.d_pruning = StandardEnumSettings::PM_Piecewise;
                    settings.d_pruning_parameter = p;
                    // ...
                }
            }
            else
                settings.d_pruning = StandardEnumSettings::PM_None;
        }
        
        template<class Enumerator>
        void prepareEnumerator(Enumerator & enumerator, bool LLLonly = false) const
        {
            typedef typename Enumerator::PreprocessorStackT::Entry EntryT;
            EntryT entry, lllentry;
            entry.d_minalpha.setContext(d_rc);
            arithmetic::convert(entry.d_minalpha, 0.999, d_rc);
            lllentry.d_method = EntryT::RM_LLL;
            lllentry.d_minalpha.setContext(d_rc);
            lllentry.d_minalpha = entry.d_minalpha;
            
            // Add simple LLL
            enumerator.getPreprocessorStack().begin_add(entry);
            enumerator.getPreprocessorStack().end_add();
            if (LLLonly)
                return;
            
            // For higher dimensions, add better preprocessing
            entry.d_min_dim = 21;
            entry.d_max_dim = 30;
            entry.d_method = EntryT::RM_BKZ;
            entry.d_bkz_mode = LatticeReduction::BKZ_HanrotPujolStehleSVP;
            entry.d_windowsize = 10;
            entry.d_force_tours = true;
            entry.d_bkz_tours = arithmetic::convert<arithmetic::Integer>(100);
            entry.d_enumboundselectionmethod = EBS_Standard;
            entry.d_enum_repetitions = 1;
            setupPruning(entry.d_enumsettings, false); // no pruning
            enumerator.getPreprocessorStack().begin_add(entry);
            {
                enumerator.getPreprocessorStack().begin_add(lllentry);
                enumerator.getPreprocessorStack().end_add();
            }
            enumerator.getPreprocessorStack().end_add();
            // For higher dimensions, add better preprocessing
            entry.d_min_dim = 31;
            entry.d_max_dim = 40;
            entry.d_method = EntryT::RM_BKZ;
            entry.d_bkz_mode = LatticeReduction::BKZ_HanrotPujolStehleSVP;
            entry.d_windowsize = 16;
            entry.d_force_tours = true;
            entry.d_bkz_tours = arithmetic::convert<arithmetic::Integer>(400);
            entry.d_enumboundselectionmethod = EBS_Standard;
            entry.d_enum_repetitions = 1;
            setupPruning(entry.d_enumsettings, false); // no pruning
            enumerator.getPreprocessorStack().begin_add(entry);
            {
                enumerator.getPreprocessorStack().begin_add(lllentry);
                enumerator.getPreprocessorStack().end_add();
            }
            enumerator.getPreprocessorStack().end_add();
            // For higher dimensions, add better preprocessing
            entry.d_min_dim = 41;
            entry.d_max_dim = std::numeric_limits<unsigned>::max();
            entry.d_method = EntryT::RM_BKZ;
            entry.d_bkz_mode = LatticeReduction::BKZ_HanrotPujolStehleSVP;
            entry.d_windowsize = 20;
            entry.d_force_tours = true;
            entry.d_bkz_tours = arithmetic::convert<arithmetic::Integer>(1200);
            entry.d_enumboundselectionmethod = EBS_Standard;
//        entry.d_enumboundselectionmethod = EBS_GaussianHeuristic105;
            entry.d_enum_repetitions = 1;
            setupPruning(entry.d_enumsettings, false); // no pruning
            enumerator.getPreprocessorStack().begin_add(entry);
            {
                enumerator.getPreprocessorStack().begin_add(lllentry);
                enumerator.getPreprocessorStack().end_add();
                entry.d_min_dim = 10;
                entry.d_max_dim = std::numeric_limits<unsigned>::max();
                entry.d_method = EntryT::RM_BKZ;
                entry.d_bkz_mode = LatticeReduction::BKZ_HanrotPujolStehleSVP;
                entry.d_windowsize = 10;
                entry.d_force_tours = true;
                entry.d_bkz_tours = arithmetic::convert<arithmetic::Integer>(100);
                entry.d_enumboundselectionmethod = EBS_Standard;
                entry.d_enum_repetitions = 1;
                setupPruning(entry.d_enumsettings, false); // no pruning
                enumerator.getPreprocessorStack().begin_add(entry);
                {
                    enumerator.getPreprocessorStack().begin_add(lllentry);
                    enumerator.getPreprocessorStack().end_add();
                }
                enumerator.getPreprocessorStack().end_add();
            }
            enumerator.getPreprocessorStack().end_add();
/*
            // For higher dimensions, add better preprocessing
            entry.d_min_dim = 50;
            entry.d_max_dim = std::numeric_limits<unsigned>::max();
            entry.d_method = EntryT::RM_BKZ;
            entry.d_bkz_mode = LatticeReduction::BKZ_HanrotPujolStehleSVP;
            entry.d_windowsize = 40;
            entry.d_force_tours = true;
            entry.d_bkz_tours = arithmetic::convert<arithmetic::Integer>(2400);
            entry.d_enumboundselectionmethod = EBS_Standard;
//        entry.d_enumboundselectionmethod = EBS_GaussianHeuristic105;
            entry.d_enum_repetitions = 1;
            setupPruning(entry.d_enumsettings, true, 0.1); // pruning
            enumerator.getPreprocessorStack().begin_add(entry);
            {
                enumerator.getPreprocessorStack().begin_add(lllentry);
                enumerator.getPreprocessorStack().end_add();
                entry.d_min_dim = 10;
                entry.d_max_dim = std::numeric_limits<unsigned>::max();
                entry.d_method = EntryT::RM_BKZ;
                entry.d_bkz_mode = LatticeReduction::BKZ_HanrotPujolStehleSVP;
                entry.d_windowsize = 16;
                entry.d_force_tours = true;
                entry.d_bkz_tours = arithmetic::convert<arithmetic::Integer>(250);
                entry.d_enumboundselectionmethod = EBS_Standard;
                entry.d_enum_repetitions = 1;
                setupPruning(entry.d_enumsettings, false); // no pruning
                enumerator.getPreprocessorStack().begin_add(entry);
                {
                    enumerator.getPreprocessorStack().begin_add(lllentry);
                    enumerator.getPreprocessorStack().end_add();
                }
                enumerator.getPreprocessorStack().end_add();
            }
            enumerator.getPreprocessorStack().end_add();
*/
        }
        
        /////////////// PROJECTED SORTING ///////////////
        
        void sortProjected(unsigned & begin, unsigned & end, TransformNotifier<IntTypeContext> * T)
        {
            freeLattice();
            try
            {
                Lattice<RealTypeContext, IntTypeContext> lattice(d_gs.get(), d_rc, d_ic, d_stats, begin, end);
                lattice.addNotifier(T);
                Workspace<RealTypeContext, IntTypeContext, Enumerator<RealTypeContext, IntTypeContext>, EmptyCallback<RealTypeContext, IntTypeContext> >
                    workspace(d_verbose, EmptyCallback<RealTypeContext, IntTypeContext>(), 0, d_max_cores, NULL);
                // Implement a simple insertion sort
                typename RealTypeContext::Real len(d_rc), minlen(d_rc);
                for (unsigned i = begin; i < end; ++i)
                {
                    lattice.update(end + 1);
                    minlen = lattice.getNormSq(i);
                    unsigned mini = i;
                    for (unsigned j = i + 1; j <= end; ++j)
                    {
                        lattice.computeProjectionLengthBV(len, i, j);
                        if (len < minlen)
                        {
                            mini = j;
                            minlen = len;
                        }
                    }
                    for (unsigned j = mini; j > i; --j)
                        lattice.swap(j, j - 1);
                }
                removeZeroVectors(workspace, lattice);
                begin = lattice.range().begin();
                end = lattice.range().end();
            }
            catch (std::exception & e)
            {
                handle_exception(&e);
            }
            catch (...)
            {
                handle_unknown_exception();
                throw;
            }
        }
        
        /////////////// SIZE REDUCTION ///////////////
        
        void sizereduction(unsigned & begin, unsigned & end, TransformNotifier<IntTypeContext> * T)
        {
            freeLattice();
            try
            {
                Lattice<RealTypeContext, IntTypeContext> lattice(d_gs.get(), d_rc, d_ic, d_stats, begin, end);
                lattice.addNotifier(T);
                Workspace<RealTypeContext, IntTypeContext, Enumerator<RealTypeContext, IntTypeContext>, EmptyCallback<RealTypeContext, IntTypeContext> >
                    workspace(d_verbose, EmptyCallback<RealTypeContext, IntTypeContext>(), 0, d_max_cores, NULL);
                workspace.sizereduction(lattice);
                removeZeroVectors(workspace, lattice);
                begin = lattice.range().begin();
                end = lattice.range().end();
            }
            catch (std::exception & e)
            {
                handle_exception(&e);
            }
            catch (...)
            {
                handle_unknown_exception();
                throw;
            }
        }
        
        template<typename CallbackObject, typename LatticeObject>
        void setupCallback_impl(CallbackObject & cb,
                                LatticeReduction::CallbackFunction cf, LatticeReduction::CallbackFunction_LI cf2, double cf_int,
                                LatticeReduction::MinCallbackFunction mcf, LatticeReduction::MinCallbackFunction_LI mcf2,
                                unsigned & begin, LatticeObject & lattice, arithmetic::IntegerContext *)
        {
            if ((cf != NULL) || (cf2 != NULL))
            {
                if (cf != NULL)
                    cb.setCallback(cf, cf_int, begin);
                else
                    cb.setCallback(boost::bind(CallbackAdaptor, _1, cf2), cf_int, begin);
            }
            if ((mcf != NULL) || (mcf2 != NULL))
            {
                if (mcf != NULL)
                    cb.setMinCallback(lattice, mcf);
                else
                    cb.setMinCallback(lattice, boost::bind(MinCallbackAdaptor, _1, _2, _3, mcf2));
            }
        }
        
        template<typename CallbackObject, typename LatticeObject>
        void setupCallback_impl(CallbackObject & cb,
                                LatticeReduction::CallbackFunction cf, LatticeReduction::CallbackFunction_LI cf2, double cf_int,
                                LatticeReduction::MinCallbackFunction mcf, LatticeReduction::MinCallbackFunction_LI mcf2,
                                unsigned & begin, LatticeObject & lattice, arithmetic::NIntContext<long int> *)
        {
            if ((cf != NULL) || (cf2 != NULL))
            {
                if (cf2 != NULL)
                    cb.setCallback(cf2, cf_int, begin);
                else
                    cb.setCallback(boost::bind(CallbackAdaptor_LI, _1, cf), cf_int, begin);
            }
            if ((mcf != NULL) || (mcf2 != NULL))
            {
                if (mcf2 != NULL)
                    cb.setMinCallback(lattice, mcf2);
                else
                    cb.setMinCallback(lattice, boost::bind(MinCallbackAdaptor_LI, _1, _2, _3, mcf));
            }
        }
        
        template<typename CallbackObject, typename LatticeObject>
        void setupCallback(CallbackObject & cb,
                           LatticeReduction::CallbackFunction cf, LatticeReduction::CallbackFunction_LI cf2, double cf_int,
                           LatticeReduction::MinCallbackFunction mcf, LatticeReduction::MinCallbackFunction_LI mcf2,
                           MaxBitsCallbackFunction mbcf,
                           unsigned & begin, LatticeObject & lattice)
        {
            setupCallback_impl(cb, cf, cf2, cf_int, mcf, mcf2, begin, lattice, static_cast<IntTypeContext*>(NULL));
            if (mbcf != NULL)
                cb.setMaxBitsCallback(lattice, mbcf);
        }
        
        /////////////// LLL ///////////////
        
        template<class Enumerator, class CallbackFunction, class Reorderer>
        void lll_internal(Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction> & workspace,
                          Lattice<RealTypeContext, IntTypeContext> & lattice, LatticeReduction::LLLMode mode,
                          bool anneal, LatticeReduction::AnnealCallbackFunction acf, LatticeReduction::LLL_AnnealFunction af, Reorderer & reorder)
        {
            if (anneal)
            {
                lattice.range().dupRange();
                typedef typename Enumerator::PreprocessorStackT::Entry EntryT;
                EntryT entry;
                entry.d_method = EntryT::RM_None;
                workspace.enumerator().getPreprocessorStack().begin_add(entry);
                workspace.enumerator().getPreprocessorStack().end_add();
                workspace.enumerator().getPreprocessorStack().finalize_adding();
                workspace.enumerator().getPreprocessorStack().begin_preprocess();
                switch (mode)
                {
                case LatticeReduction::LLL_Classic:
                {
                    if (af == NULL)
                        workspace.partialLLL_Anneal(lattice, acf, &DefaultLLLAnnealFunction<IntTypeContext>, reorder);
                    else
                        workspace.partialLLL_Anneal(lattice, acf, LLLAnnealWrapper<IntTypeContext>(af, &d_gsi), reorder);
                }
                break;
                case LatticeReduction::LLL_Unprojected:
                case LatticeReduction::LLL_Siegel:
                    lattice.range().popRange();
                    workspace.verbose()(LatticeReduction::VL_Warning) << "No annealing implemented for this LLL mode!";
                    throw feature_not_implemented();
                }
                workspace.enumerator().getPreprocessorStack().end_preprocess();
            }
            else
            {
                typedef typename Enumerator::PreprocessorStackT::Entry EntryT;
                EntryT entry;
                entry.d_method = EntryT::RM_LLL;
                entry.d_lll_mode = mode;
                entry.d_enumalg = d_svp;
                setupPruning(entry.d_enumsettings, false); // no pruning
                workspace.enumerator().getPreprocessorStack().begin_add(entry);
                workspace.enumerator().getPreprocessorStack().end_add();
                workspace.enumerator().getPreprocessorStack().finalize_adding();
                workspace.enumerator().preprocess(workspace, lattice, reorder);
            }
            return;
        }
        
        template<class Enumerator, class CallbackFunction>
        void lll_internal(Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction> & workspace,
                          Lattice<RealTypeContext, IntTypeContext> & lattice, LatticeReduction::LLLMode mode,
                          bool anneal, LatticeReduction::AnnealCallbackFunction acf, LatticeReduction::LLL_AnnealFunction af,
                          LatticeReduction::DIMethod di, LatticeReduction::DIMode di_mode, LatticeReduction::DIChoice di_choice, unsigned di_bs)
        {
            switch (di)
            {
            case LatticeReduction::DI_None:
            {
                NoReorder<RealTypeContext, IntTypeContext> r;
                lll_internal(workspace, lattice, mode, anneal, acf, af, r);
                break;
            }
            case LatticeReduction::DI_Classic:
            {
                DoClassicDeepInsertion<RealTypeContext, IntTypeContext> r(di_mode, di_choice, di_bs, d_stats);
                lll_internal(workspace, lattice, mode, anneal, acf, af, r);
                break;
            }
            case LatticeReduction::DI_MinimizePotential1:
            {
                DoMinPotDeepInsertion<RealTypeContext, IntTypeContext> r(d_intlattice.matrix().rows(), di_mode, di_choice, di_bs, d_stats, true);
                lll_internal(workspace, lattice, mode, anneal, acf, af, r);
                break;
            }
            case LatticeReduction::DI_MinimizePotential2:
            {
                DoMinPotDeepInsertion<RealTypeContext, IntTypeContext> r(d_intlattice.matrix().rows(), di_mode, di_choice, di_bs, d_stats, false);
                lll_internal(workspace, lattice, mode, anneal, acf, af, r);
                break;
            }
            }
        }
        
        void lll(unsigned & begin, unsigned & end, TransformNotifier<IntTypeContext> * T, double alpha, LatticeReduction::LLLMode mode,
                 LatticeReduction::CallbackFunction cf, LatticeReduction::CallbackFunction_LI cf2, double cf_int,
                 LatticeReduction::MinCallbackFunction mcf, LatticeReduction::MinCallbackFunction_LI mcf2, MaxBitsCallbackFunction mbcf,
                 bool anneal, LatticeReduction::AnnealCallbackFunction acf, LatticeReduction::LLL_AnnealFunction af,
                 LatticeReduction::DIMethod di, LatticeReduction::DIMode di_mode, LatticeReduction::DIChoice di_choice, unsigned di_bs)
        {
            freeLattice();
            try
            {
                Workspace<RealTypeContext, IntTypeContext, Enumerator<RealTypeContext, IntTypeContext>, Callback<RealTypeContext, IntTypeContext> >
                    workspace(d_verbose, Callback<RealTypeContext, IntTypeContext>(end - begin + 1), 0, d_max_cores, NULL, NULL);
                Lattice<RealTypeContext, IntTypeContext> lattice(d_gs.get(), d_rc, d_ic, d_stats, begin, end);
                lattice.addNotifier(T);
                setupCallback(workspace.getCallbackObject(), cf, cf2, cf_int, mcf, mcf2, mbcf, begin, lattice);
                d_gs->adjustAlpha(alpha);
                
                lll_internal(workspace, lattice, mode, anneal, acf, af, di, di_mode, di_choice, di_bs);
                removeZeroVectors(workspace, lattice);
                
                begin = lattice.range().begin();
                end = lattice.range().end();
            }
            catch (std::exception & e)
            {
                handle_exception(&e);
            }
            catch (...)
            {
                handle_unknown_exception();
                throw;
            }
        }
        
        /////////////// BKZ ///////////////
        
        template<class Enumerator, class CallbackFunction, class Reorderer>
        inline void bkz_internal(Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction> & workspace,
                                 Lattice<RealTypeContext, IntTypeContext> & lattice,
                                 unsigned blocksize, LatticeReduction::BKZMode mode,
                                 bool anneal, LatticeReduction::AnnealCallbackFunction acf, LatticeReduction::BKZ_AnnealFunction af,
                                 Reorderer & reorder)
        {
            if (anneal)
            {
                lattice.range().dupRange();
                typedef typename Enumerator::PreprocessorStackT::Entry EntryT;
                EntryT entry;
                entry.d_method = EntryT::RM_None;
                entry.d_enumalg = d_svp;
                workspace.enumerator().getPreprocessorStack().begin_add(entry);
                prepareEnumerator(workspace.enumerator());
                workspace.enumerator().getPreprocessorStack().end_add();
                workspace.enumerator().getPreprocessorStack().finalize_adding();
                workspace.enumerator().getPreprocessorStack().begin_preprocess();
                switch (mode)
                {
                case LatticeReduction::BKZ_SchnorrEuchner:
                {
                    if (af == NULL)
                        workspace.partialBKZ_Anneal(lattice, blocksize, acf, &DefaultBKZAnnealFunction<IntTypeContext>, reorder, false);
                    else
                        workspace.partialBKZ_Anneal(lattice, blocksize, acf, BKZAnnealWrapper<IntTypeContext>(af, &d_gsi), reorder, false);
                }
                break;
                case LatticeReduction::BKZ_Simplified:
                {
                    if (af == NULL)
                        workspace.partialBKZ_Anneal(lattice, blocksize, acf, &DefaultBKZAnnealFunction<IntTypeContext>, reorder, true);
                    else
                        workspace.partialBKZ_Anneal(lattice, blocksize, acf, BKZAnnealWrapper<IntTypeContext>(af, &d_gsi), reorder, true);
                }
                break;
                case LatticeReduction::BKZ_HanrotPujolStehleSVP:
                case LatticeReduction::BKZ_HanrotPujolStehleHKZ:
                case LatticeReduction::BKZ_SemiBlock2k: // Schnorr's Semi-block-2k-reduction (blocksize is made even by rounding up)
                case LatticeReduction::BKZ_PrimalDual: // Koy's primal-dual BKZ
                case LatticeReduction::BKZ_SlideReduction: // Gama-Nguyen Slide Reduction (Gama, Nguyen: "Finding Short Lattice Vectors within Mordell's Inequality")
                case LatticeReduction::BKZ_ImprovedSlideReduction: // Schnorr's improved and accelerated Slide Reduction
                case LatticeReduction::BKZ_ImprovedSlideReduction2: // Schnorr's improved and accelerated Slide Reduction (using a larger DSVP and a single SVP)
                case LatticeReduction::BKZ_ImprovedSlideReduction3: // Schnorr's improved and accelerated Slide Reduction (using a single larger SVP and a DSVP)
                case LatticeReduction::BKZ_SamplingReduction:
                case LatticeReduction::BKZ_Experimental:
                    lattice.range().popRange();
                    workspace.verbose()(LatticeReduction::VL_Warning) << "No annealing implemented for this BKZ mode!";
                    throw feature_not_implemented();
                }
                workspace.enumerator().getPreprocessorStack().end_preprocess();
            }
            else
            {
                typedef typename Enumerator::PreprocessorStackT::Entry EntryT;
                EntryT entry;
                entry.d_method = EntryT::RM_BKZ;
                entry.d_bkz_mode = mode;
                entry.d_windowsize = blocksize;
                entry.d_enumalg = d_svp;
                setupPruning(entry.d_enumsettings, false); // no pruning
                workspace.enumerator().getPreprocessorStack().begin_add(entry);
                prepareEnumerator(workspace.enumerator());
                workspace.enumerator().getPreprocessorStack().end_add();
                workspace.enumerator().getPreprocessorStack().finalize_adding();
                workspace.enumerator().preprocess(workspace, lattice, reorder);
            }
            return;
        }
        
        template<class Enumerator, class CallbackFunction>
        inline void bkz_internal(Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction> & workspace,
                                 Lattice<RealTypeContext, IntTypeContext> & lattice,
                                 unsigned blocksize, LatticeReduction::BKZMode mode,
                                 bool anneal, LatticeReduction::AnnealCallbackFunction acf, LatticeReduction::BKZ_AnnealFunction af,
                                 LatticeReduction::DIMethod di, LatticeReduction::DIMode di_mode, LatticeReduction::DIChoice di_choice, unsigned di_bs)
        {
            switch (di)
            {
            case LatticeReduction::DI_None:
            {
                NoReorder<RealTypeContext, IntTypeContext> r;
                bkz_internal(workspace, lattice, blocksize, mode, anneal, acf, af, r);
                break;
            }
            case LatticeReduction::DI_Classic:
            {
                DoClassicDeepInsertion<RealTypeContext, IntTypeContext> r(di_mode, di_choice, di_bs, d_stats);
                bkz_internal(workspace, lattice, blocksize, mode, anneal, acf, af, r);
                break;
            }
            case LatticeReduction::DI_MinimizePotential1:
            {
                DoMinPotDeepInsertion<RealTypeContext, IntTypeContext> r(lattice.dimension(), di_mode, di_choice, di_bs, d_stats, true);
                bkz_internal(workspace, lattice, blocksize, mode, anneal, acf, af, r);
                break;
            }
            case LatticeReduction::DI_MinimizePotential2:
            {
                DoMinPotDeepInsertion<RealTypeContext, IntTypeContext> r(lattice.dimension(), di_mode, di_choice, di_bs, d_stats, false);
                bkz_internal(workspace, lattice, blocksize, mode, anneal, acf, af, r);
                break;
            }
            }
        }
        
        void bkz(unsigned & begin, unsigned & end, TransformNotifier<IntTypeContext> * T, double alpha, unsigned blocksize,
                 LatticeReduction::BKZMode mode, LatticeReduction::CallbackFunction cf, LatticeReduction::CallbackFunction_LI cf2,
                 double cf_int, LatticeReduction::MinCallbackFunction mcf, LatticeReduction::MinCallbackFunction_LI mcf2,
                 MaxBitsCallbackFunction mbcf, LatticeReduction::EnumCallbackFunction ecf, LatticeReduction::EnumCallbackFunction_LI ecf2,
                 bool anneal, LatticeReduction::AnnealCallbackFunction acf, LatticeReduction::BKZ_AnnealFunction af,
                 LatticeReduction::DIMethod di, LatticeReduction::DIMode di_mode, LatticeReduction::DIChoice di_choice, unsigned di_bs)
        {
            freeLattice();
            try
            {
                Workspace<RealTypeContext, IntTypeContext, Enumerator<RealTypeContext, IntTypeContext>, Callback<RealTypeContext, IntTypeContext> >
                    workspace(d_verbose, Callback<RealTypeContext, IntTypeContext>(end - begin + 1), blocksize, d_max_cores, ecf, ecf2);
                Lattice<RealTypeContext, IntTypeContext> lattice(d_gs.get(), d_rc, d_ic, d_stats, begin, end);
                lattice.addNotifier(T);
                setupCallback(workspace.getCallbackObject(), cf, cf2, cf_int, mcf, mcf2, mbcf, begin, lattice);
                d_gs->adjustAlpha(alpha);
                
                bkz_internal(workspace, lattice, blocksize, mode, anneal, acf, af, di, di_mode, di_choice, di_bs);
                removeZeroVectors(workspace, lattice);
                
                begin = lattice.range().begin();
                end = lattice.range().end();
            }
            catch (std::exception & e)
            {
                handle_exception(&e);
            }
            catch (...)
            {
                handle_unknown_exception();
                throw;
            }
        }
        
        /////////////// HKZ ///////////////
        
        void hkz(unsigned & begin, unsigned & end, TransformNotifier<IntTypeContext> * T, bool dual,
                 LatticeReduction::CallbackFunction cf, LatticeReduction::CallbackFunction_LI cf2, double cf_int,
                 LatticeReduction::MinCallbackFunction mcf, LatticeReduction::MinCallbackFunction_LI mcf2,
                 MaxBitsCallbackFunction mbcf,
                 LatticeReduction::EnumCallbackFunction ecf, LatticeReduction::EnumCallbackFunction_LI ecf2)
        {
            freeLattice();
            try
            {
                Workspace<RealTypeContext, IntTypeContext, Enumerator<RealTypeContext, IntTypeContext>, Callback<RealTypeContext, IntTypeContext> >
                    workspace(d_verbose, Callback<RealTypeContext, IntTypeContext>(end - begin + 1), end - begin + 1, d_max_cores, ecf, ecf2);
                Lattice<RealTypeContext, IntTypeContext> lattice(d_gs.get(), d_rc, d_ic, d_stats, begin, end);
                lattice.addNotifier(T);
                setupCallback(workspace.getCallbackObject(), cf, cf2, cf_int, mcf, mcf2, mbcf, begin, lattice);
                d_gs->adjustAlpha(0.999);
                
                typedef typename Enumerator<RealTypeContext, IntTypeContext>::PreprocessorStackT::Entry EntryT;
                EntryT entry;
                entry.d_method = dual ? EntryT::RM_HKZDual : EntryT::RM_HKZ;
                entry.d_enumalg = d_svp;
                setupPruning(entry.d_enumsettings, false); // no pruning
                workspace.enumerator().getPreprocessorStack().begin_add(entry);
                prepareEnumerator(workspace.enumerator());
                workspace.enumerator().getPreprocessorStack().end_add();
                workspace.enumerator().getPreprocessorStack().finalize_adding();
                workspace.enumerator().preprocess(workspace, lattice);
                removeZeroVectors(workspace, lattice);
                
                begin = lattice.range().begin();
                end = lattice.range().end();
            }
            catch (std::exception & e)
            {
                handle_exception(&e);
            }
            catch (...)
            {
                handle_unknown_exception();
                throw;
            }
        }
        
        /////////////// SVP ///////////////
        
        inline void svp(unsigned & begin, unsigned & end, TransformNotifier<IntTypeContext> * T, bool make_basis, bool extreme, bool dual,
                        LatticeReduction::CallbackFunction cf, LatticeReduction::CallbackFunction_LI cf2, double cf_int,
                        LatticeReduction::MinCallbackFunction mcf, LatticeReduction::MinCallbackFunction_LI mcf2,
                        MaxBitsCallbackFunction mbcf,
                        LatticeReduction::EnumCallbackFunction ecf, LatticeReduction::EnumCallbackFunction_LI ecf2)
        {
            freeLattice();
            try
            {
                Workspace<RealTypeContext, IntTypeContext, Enumerator<RealTypeContext, IntTypeContext>, Callback<RealTypeContext, IntTypeContext> >
                    workspace(d_verbose, Callback<RealTypeContext, IntTypeContext>(end - begin + 1), end - begin + 1, d_max_cores, ecf, ecf2);
                Lattice<RealTypeContext, IntTypeContext> lattice(d_gs.get(), d_rc, d_ic, d_stats, begin, end);
                lattice.addNotifier(T);
                setupCallback(workspace.getCallbackObject(), cf, cf2, cf_int, mcf, mcf2, mbcf, begin, lattice);
                d_gs->adjustAlpha(0.999);
                
                lattice.range().dupRange();
                typedef typename Enumerator<RealTypeContext, IntTypeContext>::PreprocessorStackT::Entry EntryT;
                EntryT entry, entry2;
                entry.d_method = dual ? EntryT::RM_SVPDual : EntryT::RM_SVP;
                entry.d_enumalg = d_svp;
                setOne(entry.d_minalpha);
                if (extreme)
                {
                    entry.d_enumboundselectionmethod = EBS_GaussianHeuristic105;
                    entry.d_enum_repetitions = 256;
                    setupPruning(entry.d_enumsettings, true, 0, true); // extreme pruning
                    workspace.enumerator().getPreprocessorStack().begin_add(entry);
                    prepareEnumerator(workspace.enumerator(), true); // add LLL step
                    entry2.d_method = EntryT::RM_BKZ;
                    entry2.d_bkz_mode = LatticeReduction::BKZ_HanrotPujolStehleSVP;
                    entry2.d_windowsize = (end - begin + 1) / 5 + 10; // 110 => 30, 120 => 32
                    workspace.enumerator().getPreprocessorStack().begin_add(entry2);
                    prepareEnumerator(workspace.enumerator());
                    workspace.enumerator().getPreprocessorStack().end_add();
                    workspace.enumerator().getPreprocessorStack().end_add();
                }
                else
                {
                    setupPruning(entry.d_enumsettings, false); // no pruning
                    workspace.enumerator().getPreprocessorStack().begin_add(entry);
                    prepareEnumerator(workspace.enumerator());
                    workspace.enumerator().getPreprocessorStack().end_add();
                }
                workspace.enumerator().getPreprocessorStack().finalize_adding();
                workspace.enumerator().preprocess(workspace, lattice);
                removeZeroVectors(workspace, lattice);
                
                begin = lattice.range().begin();
                end = lattice.range().end();
            }
            catch (std::exception & e)
            {
                handle_exception(&e);
            }
            catch (...)
            {
                handle_unknown_exception();
                throw;
            }
        }
        
        /////////////// VERIFICATION ///////////////
        
        virtual bool isSizeReduced(unsigned begin, unsigned end) const
        {
            try
            {
                Lattice<RealTypeContext, IntTypeContext> lattice(d_gs.get(), d_rc, d_ic, d_stats, begin, end, true);
                Workspace<RealTypeContext, IntTypeContext, Enumerator<RealTypeContext, IntTypeContext>, EmptyCallback<RealTypeContext, IntTypeContext> >
                    workspace(d_verbose, EmptyCallback<RealTypeContext, IntTypeContext>(), 0, d_max_cores, NULL);
                return isSR(workspace, lattice);
            }
            catch (std::exception & e)
            {
                handle_exception(&e);
                return false;
            }
            catch (...)
            {
                handle_unknown_exception();
                throw;
            }
        }
        
        template<class Enumerator, class CallbackFunction>
        bool isDI_Classic(Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction> & workspace,
                          Lattice<RealTypeContext, IntTypeContext> & lattice,
                          LatticeReduction::DIChoice di_choice, unsigned di_bs) const
        {
            lattice.update(d_intlattice.matrix().rows());
            typename RealTypeContext::Real res(lattice.rc());
            
            for (unsigned j = lattice.range().begin() + 1; j <= lattice.range().end(); ++j)
            {
                unsigned i = lattice.range().begin();
                for (; i < j; ++i)
                {
                    if (i >= lattice.range().begin() + di_bs)
                    {
                        if (di_choice == LatticeReduction::DIC_First)
                        {
                            i = j;
                            break;
                        }
                        if ((di_choice != LatticeReduction::DIC_All) && (i + di_bs < j))
                            i = j - di_bs;
                    }
                    
                    lattice.computeProjectionLengthBV(res, i, j);
                    if (res < lattice.getLLLalpha() * lattice.getNormSq(i))
                    {
                        if (workspace.verbose().yieldsOutput(LatticeReduction::VL_Information))
                            workspace.verbose()(LatticeReduction::VL_Information) << "Deep Insertion condition failed for (" << i << ", " << j << ")\n";
                        return false;
                    }
                }
            }
            return true;
        }
        
        template<class Enumerator, class CallbackFunction>
        bool isDI_MinPot(Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction> & workspace,
                         Lattice<RealTypeContext, IntTypeContext> & lattice,
                         LatticeReduction::DIChoice di_choice, unsigned di_bs, int mode = 0) const
        {
            lattice.update(d_intlattice.matrix().rows());
            typedef typename arithmetic::HugeExponent<RealTypeContext>::Type Huge;
            typename RealTypeContext::Real t(lattice.rc());
            Huge minpot(lattice.rc()), pot(lattice.rc());
            for (unsigned j = lattice.range().begin() + 1; j <= lattice.range().end(); ++j)
            {
                unsigned i = j - 1;
                minpot = Huge(lattice.getLLLalpha());
                setOne(pot);
                
                unsigned blockbegin = lattice.range().begin();
                if ((di_choice == LatticeReduction::DIC_Block) && (blockbegin + di_bs < i))
                    blockbegin = i - di_bs;
                if (blockbegin == i)
                    continue;
                
                for (; i >= blockbegin; --i)
                {
                    // Update potential
                    pot *= square(Huge(lattice.getNormSq(i - 1)));
                    lattice.computeProjectionLengthBV(t, i - 1, i);
                    pot /= square(Huge(t));
                    
                    // Check whether we want to actually compare the potential for this position
                    bool doComp = true;
                    if (i > lattice.range().begin() + di_bs)
                    {
                        if (di_choice == LatticeReduction::DIC_First)
                            doComp = false;
                        if ((di_choice != LatticeReduction::DIC_All) && (i + di_bs <= j))
                            doComp = false;
                    }
                    else
                        if ((di_choice == LatticeReduction::DIC_Block) && (i + di_bs <= j))
                            doComp = false;
                    
                    // Compare potential
                    if (doComp)
                        if (pot < minpot)
                        {
                            if (workspace.verbose().yieldsOutput(LatticeReduction::VL_Information))
                                workspace.verbose()(LatticeReduction::VL_Information) << "MinPot Deep Insertion condition failed for (" << i << ", " << j << ")\n";
                            return false;
                        }
                }
            }
            return true;
        }
        
        template<class Enumerator, class CallbackFunction>
        bool isDI(Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction> & workspace,
                  Lattice<RealTypeContext, IntTypeContext> & lattice,
                  LatticeReduction::DIMethod di, LatticeReduction::DIChoice di_choice, unsigned di_bs) const
        {
            switch (di)
            {
            case LatticeReduction::DI_None: return true;
            case LatticeReduction::DI_Classic: return isDI_Classic(workspace, lattice, di_choice, di_bs);
            case LatticeReduction::DI_MinimizePotential1: return isDI_MinPot(workspace, lattice, di_choice, di_bs, 0);
            case LatticeReduction::DI_MinimizePotential2: return isDI_MinPot(workspace, lattice, di_choice, di_bs, 1);
            }
            return true; // unnecessary, but g++ complains if this is missing...
        }
        
        virtual bool isLLLBasis(unsigned begin, unsigned end, double alpha, LatticeReduction::LLLMode mode,
                                LatticeReduction::DIMethod di, LatticeReduction::DIChoice di_choice, unsigned di_bs) const
        {
            try
            {
                Lattice<RealTypeContext, IntTypeContext> lattice(d_gs.get(), d_rc, d_ic, d_stats, begin, end, true);
                Workspace<RealTypeContext, IntTypeContext, Enumerator<RealTypeContext, IntTypeContext>, EmptyCallback<RealTypeContext, IntTypeContext> >
                    workspace(d_verbose, EmptyCallback<RealTypeContext, IntTypeContext>(), 0, d_max_cores, NULL);
                d_gs->adjustAlpha(alpha);
                if (!(isSR(workspace, lattice) && isDI(workspace, lattice, di, di_choice, di_bs)))
                    return false;
                switch (mode)
                {
                case LatticeReduction::LLL_Classic:
                    return isLLL(workspace, lattice,
                                 Workspace<RealTypeContext, IntTypeContext, Enumerator<RealTypeContext, IntTypeContext>, EmptyCallback<RealTypeContext, IntTypeContext> >::LovaszCondition);
                case LatticeReduction::LLL_Unprojected:
                    return isLLL(workspace, lattice,
                                 Workspace<RealTypeContext, IntTypeContext, Enumerator<RealTypeContext, IntTypeContext>, EmptyCallback<RealTypeContext, IntTypeContext> >::UnprojectedLovaszCondition);
                case LatticeReduction::LLL_Siegel:
                    return isLLL(workspace, lattice,
                                 Workspace<RealTypeContext, IntTypeContext, Enumerator<RealTypeContext, IntTypeContext>, EmptyCallback<RealTypeContext, IntTypeContext> >::SiegelCondition);
                }
                return false;
            }
            catch (std::exception & e)
            {
                handle_exception(&e);
                return false;
            }
            catch (...)
            {
                handle_unknown_exception();
                throw;
            }
        }
        
        virtual bool isBKZBasis(unsigned begin, unsigned end, double alpha, unsigned blocksize, LatticeReduction::BKZMode mode,
                                LatticeReduction::DIMethod di, LatticeReduction::DIChoice di_choice, unsigned di_bs) const
        {
            try
            {
                Lattice<RealTypeContext, IntTypeContext> lattice(d_gs.get(), d_rc, d_ic, d_stats, begin, end, true);
                Workspace<RealTypeContext, IntTypeContext, Enumerator<RealTypeContext, IntTypeContext>, EmptyCallback<RealTypeContext, IntTypeContext> >
                    workspace(d_verbose, EmptyCallback<RealTypeContext, IntTypeContext>(), blocksize, d_max_cores);
                d_gs->adjustAlpha(alpha);
                typedef typename Enumerator<RealTypeContext, IntTypeContext>::PreprocessorStackT::Entry EntryT;
                EntryT entry;
                entry.d_method = EntryT::RM_None;
                entry.d_enumalg = d_svp;
                workspace.enumerator().getPreprocessorStack().begin_add(entry);
                prepareEnumerator(workspace.enumerator());
                workspace.enumerator().getPreprocessorStack().end_add();
                workspace.enumerator().getPreprocessorStack().finalize_adding();
                workspace.enumerator().getPreprocessorStack().begin_preprocess();
                bool result = true;
                switch (mode)
                {
                case LatticeReduction::BKZ_SchnorrEuchner:
                case LatticeReduction::BKZ_Simplified:
                case LatticeReduction::BKZ_SamplingReduction:
                    result = isSR(workspace, lattice) && isDI(workspace, lattice, di, di_choice, di_bs) && isBKZ(workspace, lattice, blocksize);
                    break;
                case LatticeReduction::BKZ_HanrotPujolStehleSVP:
                case LatticeReduction::BKZ_HanrotPujolStehleHKZ:
                case LatticeReduction::BKZ_Experimental:
                    // Don't know how to check this!!! ??? ...
                    result = isSR(workspace, lattice) && isDI(workspace, lattice, di, di_choice, di_bs)
                        && isLLL(workspace, lattice,
                                 Workspace<RealTypeContext, IntTypeContext, Enumerator<RealTypeContext, IntTypeContext>, EmptyCallback<RealTypeContext, IntTypeContext> >::LovaszCondition);
                    break;
                case LatticeReduction::BKZ_SemiBlock2k: // Schnorr's Semi-block-2k-reduction (blocksize is made even by rounding up)
                    result = isSR(workspace, lattice) && isDI(workspace, lattice, di, di_choice, di_bs) && isBKZ_SemiBlock2k(workspace, lattice, blocksize);
                    break;
                case LatticeReduction::BKZ_PrimalDual: // Koy's primal-dual BKZ
                    result = isSR(workspace, lattice) && isDI(workspace, lattice, di, di_choice, di_bs) && isBKZ_PrimalDual(workspace, lattice, blocksize);
                    break;
                case LatticeReduction::BKZ_SlideReduction: // Gama-Nguyen Slide Reduction (Gama, Nguyen: "Finding Short Lattice Vectors within Mordell's Inequality")
                    result = isSR(workspace, lattice) && isDI(workspace, lattice, di, di_choice, di_bs) && isBKZ_Slide(workspace, lattice, blocksize);
                    break;
                case LatticeReduction::BKZ_ImprovedSlideReduction: // Schnorr's improved and accelerated Slide Reduction
                case LatticeReduction::BKZ_ImprovedSlideReduction2: // Schnorr's improved and accelerated Slide Reduction (using a larger DSVP and a single SVP)
                case LatticeReduction::BKZ_ImprovedSlideReduction3: // Schnorr's improved and accelerated Slide Reduction (using a single larger SVP and a DSVP)
                    result = isSR(workspace, lattice) && isDI(workspace, lattice, di, di_choice, di_bs) && isBKZ_ImprovedSlide(workspace, lattice, blocksize);
                    break;
                }
                workspace.enumerator().getPreprocessorStack().end_preprocess();
                // Default:
                return result;
            }
            catch (std::exception & e)
            {
                handle_exception(&e);
                return false;
            }
            catch (...)
            {
                handle_unknown_exception();
                throw;
            }
        }
        
        virtual bool isHKZBasis(unsigned begin, unsigned end, bool dual) const
        {
            try
            {
                Lattice<RealTypeContext, IntTypeContext> lattice(d_gs.get(), d_rc, d_ic, d_stats, begin, end, true);
                Workspace<RealTypeContext, IntTypeContext, Enumerator<RealTypeContext, IntTypeContext>, EmptyCallback<RealTypeContext, IntTypeContext> >
                    workspace(d_verbose, EmptyCallback<RealTypeContext, IntTypeContext>(), lattice.range().dimension(), d_max_cores);
                d_gs->adjustAlpha(1);
                if (!isSR(workspace, lattice))
                    return false;
                typedef typename Enumerator<RealTypeContext, IntTypeContext>::PreprocessorStackT::Entry EntryT;
                EntryT entry;
                entry.d_method = EntryT::RM_None;
                entry.d_enumalg = d_svp;
                workspace.enumerator().getPreprocessorStack().begin_add(entry);
                prepareEnumerator(workspace.enumerator());
                workspace.enumerator().getPreprocessorStack().end_add();
                workspace.enumerator().getPreprocessorStack().finalize_adding();
                workspace.enumerator().getPreprocessorStack().begin_preprocess();
                bool result = isHKZ(workspace, lattice, dual);
                workspace.enumerator().getPreprocessorStack().end_preprocess();
                return result;
            }
            catch (std::exception & e)
            {
                handle_exception(&e);
                return false;
            }
            catch (...)
            {
                handle_unknown_exception();
                throw;
            }
        }
        
        virtual bool isSVPBasis(unsigned begin, unsigned end, bool dual) const
        {
            try
            {
                Lattice<RealTypeContext, IntTypeContext> lattice(d_gs.get(), d_rc, d_ic, d_stats, begin, end, true);
                Workspace<RealTypeContext, IntTypeContext, Enumerator<RealTypeContext, IntTypeContext>, EmptyCallback<RealTypeContext, IntTypeContext> >
                    workspace(d_verbose, EmptyCallback<RealTypeContext, IntTypeContext>(), lattice.range().dimension(), d_max_cores);
                d_gs->adjustAlpha(1);
                typedef typename Enumerator<RealTypeContext, IntTypeContext>::PreprocessorStackT::Entry EntryT;
                EntryT entry;
                entry.d_method = EntryT::RM_None;
                entry.d_enumalg = d_svp;
                workspace.enumerator().getPreprocessorStack().begin_add(entry);
                prepareEnumerator(workspace.enumerator());
                workspace.enumerator().getPreprocessorStack().end_add();
                workspace.enumerator().getPreprocessorStack().finalize_adding();
                workspace.enumerator().getPreprocessorStack().begin_preprocess();
                bool result = isSVP(workspace, lattice, dual);
                workspace.enumerator().getPreprocessorStack().end_preprocess();
                return result;
            }
            catch (std::exception & e)
            {
                handle_exception(&e);
                return false;
            }
            catch (...)
            {
                handle_unknown_exception();
                throw;
            }
        }
        
        /////////////// WRAPPERS ///////////////
        
        virtual void sortProjected(unsigned & begin, unsigned & end, Transform & transform)
        {
            sortProjected(begin, end, transform.createTransformObject<IntTypeContext>(d_ic));
            d_intlattice.update_original();
            transform.releaseTransformObject();
        }
        
        virtual void sizereduction(unsigned & begin, unsigned & end, Transform & transform)
        {
            sizereduction(begin, end, transform.createTransformObject<IntTypeContext>(d_ic));
            d_intlattice.update_original();
            transform.releaseTransformObject();
        }
        
        virtual void lll(unsigned & begin, unsigned & end, Transform & transform, double alpha, LatticeReduction::LLLMode mode,
                         LatticeReduction::CallbackFunction cf, LatticeReduction::CallbackFunction_LI cf2, double cf_int,
                         LatticeReduction::MinCallbackFunction mcf, LatticeReduction::MinCallbackFunction_LI mcf2,
                         MaxBitsCallbackFunction mbcf, bool anneal, LatticeReduction::AnnealCallbackFunction acf,
                         LatticeReduction::LLL_AnnealFunction af, LatticeReduction::DIMethod di, LatticeReduction::DIMode di_mode,
                         LatticeReduction::DIChoice di_choice, unsigned di_bs)
        {
            lll(begin, end, transform.createTransformObject<IntTypeContext>(d_ic), alpha, mode, cf, cf2, cf_int,
                mcf, mcf2, mbcf, anneal, acf, af, di, di_mode, di_choice, di_bs);
            d_intlattice.update_original();
            transform.releaseTransformObject();
        }
        
        virtual void bkz(unsigned & begin, unsigned & end, Transform & transform, double alpha, unsigned blocksize, LatticeReduction::BKZMode mode,
                         LatticeReduction::CallbackFunction cf, LatticeReduction::CallbackFunction_LI cf2, double cf_int,
                         LatticeReduction::MinCallbackFunction mcf, LatticeReduction::MinCallbackFunction_LI mcf2, MaxBitsCallbackFunction mbcf,
                         LatticeReduction::EnumCallbackFunction ecf, LatticeReduction::EnumCallbackFunction_LI ecf2,
                         bool anneal, LatticeReduction::AnnealCallbackFunction acf, LatticeReduction::BKZ_AnnealFunction af,
                         LatticeReduction::DIMethod di, LatticeReduction::DIMode di_mode, LatticeReduction::DIChoice di_choice, unsigned di_bs)
        {
            bkz(begin, end, transform.createTransformObject<IntTypeContext>(d_ic), alpha, blocksize, mode, cf, cf2, cf_int,
                mcf, mcf2, mbcf, ecf, ecf2, anneal, acf, af, di, di_mode, di_choice, di_bs);
            d_intlattice.update_original();
            transform.releaseTransformObject();
        }
        
        virtual void hkz(unsigned & begin, unsigned & end, Transform & transform, bool dual,
                         LatticeReduction::CallbackFunction cf, LatticeReduction::CallbackFunction_LI cf2, double cf_int,
                         LatticeReduction::MinCallbackFunction mcf, LatticeReduction::MinCallbackFunction_LI mcf2,
                         MaxBitsCallbackFunction mbcf,
                         LatticeReduction::EnumCallbackFunction ecf, LatticeReduction::EnumCallbackFunction_LI ecf2)
        {
            hkz(begin, end, transform.createTransformObject<IntTypeContext>(d_ic), dual,
                cf, cf2, cf_int, mcf, mcf2, mbcf, ecf, ecf2);
            d_intlattice.update_original();
            transform.releaseTransformObject();
        }
        
        virtual void svp(unsigned & begin, unsigned & end, Transform & transform, bool make_basis, bool extreme, bool dual,
                         LatticeReduction::CallbackFunction cf, LatticeReduction::CallbackFunction_LI cf2, double cf_int,
                         LatticeReduction::MinCallbackFunction mcf, LatticeReduction::MinCallbackFunction_LI mcf2,
                         MaxBitsCallbackFunction mbcf,
                         LatticeReduction::EnumCallbackFunction ecf, LatticeReduction::EnumCallbackFunction_LI ecf2)
        {
            svp(begin, end, transform.createTransformObject<IntTypeContext>(d_ic), make_basis, extreme, dual,
                cf, cf2, cf_int, mcf, mcf2, mbcf, ecf, ecf2);
            d_intlattice.update_original();
            transform.releaseTransformObject();
        }
    };
    
    template<class RealTypeContext, class IntTypeContext>
    LRIInterface * CreateLRIInterfaceWithContexts(LatticeReduction::VerboseOutputLevel vol, LatticeReduction::VerboseFunction vf,
                                                  linalg::math_matrix<arithmetic::Integer> & lattice, LatticeReduction::GramSchmidt gs, bool gsr,
                                                  LatticeReduction::SVPMode svp, unsigned max_cores,
                                                  const RealTypeContext & rc, const IntTypeContext & ic)
    {
        return new LRIImplementation<RealTypeContext, IntTypeContext>(vol, vf, lattice, gs, gsr, svp, max_cores, rc, ic);
    }
}

#endif

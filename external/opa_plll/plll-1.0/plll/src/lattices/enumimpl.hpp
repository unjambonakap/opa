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

#ifndef PLLL_INCLUDE_GUARD__ENUM_IMPL_HPP
#define PLLL_INCLUDE_GUARD__ENUM_IMPL_HPP

#include "myalloc.hpp"
#include <stack>
#include <limits>

namespace plll
{
    class StandardEnumSettings
    {
    public:
        enum PruningMethod {
            PM_None, // no pruning
            PM_Linear, PM_Step, PM_Piecewise, // pruning functions from Gama, Nguyen, Regev: "Lattice
            // Enumeration using Extreme Pruning", EUROCRYPT 2010
            PM_SchnorrHoerner1, PM_SchnorrHoerner2, // pruning functions from Schnorr, Hoerner:
            // "Attacking the Chor-Rivest Cryptosystem by
            // improved lattice reduction"
            PM_GNR110_Poly8  // "optimal" XP curve for dimension 110 by Gama, Regev and Nguyen,
            // interpolated by a degree 8 polynomial as in Kuo, Schneider, Dagdelen,
            // Reichelt, Buchmann, Cheng, Yang: "Extreme Enumeration on GPU and in
            // Clouds -- How Many Dollars You Need to Break SVP Challenges".
        };
    
        PruningMethod d_pruning;
        double d_pruning_parameter; // needed for PM_Step, PM_Piecewise and PM_SchnorrHoerner2
    
        StandardEnumSettings()
            : d_pruning(PM_None), d_pruning_parameter(1.0)
        {
        }
    };

    enum EnumBoundSelectionMethod { EBS_Standard, EBS_GaussianHeuristic105 };

    template<class RealTypeContext, class IntTypeContext,
             class EnumSettings> // EnumSettings contain infos about pruning etc. which are enumerator specific
    class Preprocessor
    {
    public:
        enum ReductionMethod { RM_None, RM_LLL, RM_BKZ, RM_SVP, RM_HKZ, RM_SVPDual, RM_HKZDual };
        typedef void (*ToursComputer)(arithmetic::Integer &, const Lattice<RealTypeContext, IntTypeContext> &);
        typedef void (*EnumSettingsComputer)(EnumSettings &, const Lattice<RealTypeContext, IntTypeContext> &);
    
        unsigned d_min_dim, d_max_dim; // apply this preprocessing step only if enumeration dimension is in range
        typename RealTypeContext::Real d_minalpha; // minimal value for reduction parameter alpha !!! ??? ... NOT IMPLEMENTED CURRENTLY!!!
        ReductionMethod d_method;
        LatticeReduction::LLLMode d_lll_mode;
        LatticeReduction::BKZMode d_bkz_mode;
        unsigned d_windowsize; // for BKZ-type methods (RM_BKZ*)
        bool d_force_tours; // for RM_BKZ with modes BKZ_HanrotPujolStehle*
        arithmetic::Integer d_bkz_tours; // for RM_BKZ with modes BKZ_HanrotPujolStehle*
        ToursComputer d_bkz_tours_computer; // for RM_BKZ with modes BKZ_HanrotPujolStehle*
        EnumBoundSelectionMethod d_enumboundselectionmethod; // for RM_BKZ, RM_SVP*, RM_HKZ*
        LatticeReduction::SVPMode d_enumalg;
        mutable EnumSettings d_enumsettings; // for RM_BKZ, RM_SVP*, RM_HKZ*
        EnumSettingsComputer d_enum_settings_computer; // for RM_BKZ, RM_SVP*, RM_HKZ*
        unsigned d_enum_repetitions; // for extreme pruning: repeat enumeration n times and take shortest result
    
        Preprocessor()
            : d_min_dim(0), d_max_dim(std::numeric_limits<unsigned>::max()),
              d_method(RM_None), d_lll_mode(LatticeReduction::LLL_Default), d_bkz_mode(LatticeReduction::BKZ_Default),
              d_windowsize(10), d_force_tours(false), d_bkz_tours_computer(NULL),
              d_enumboundselectionmethod(EBS_Standard), d_enumalg(LatticeReduction::SVP_ParallelKannanSchnorrEuchner),
              d_enum_settings_computer(NULL), d_enum_repetitions(1)
        {
            setZero(d_minalpha);
        }
    };

    template<class RealTypeContext, class IntTypeContext,
             class EnumSettings> // EnumSettings contain infos about pruning, ...
    class PreprocessorStack
    {
    public:
        typedef Preprocessor<RealTypeContext, IntTypeContext, EnumSettings> Entry;
    
    private:
        class PaddedEntry;

        typedef std::list<PaddedEntry> EntryCollector;
    
        typedef std::pair<EntryCollector *, typename std::list<PaddedEntry>::const_iterator> EntryPosition;
    
        EntryCollector d_root;
        bool d_can_add;
    
        mutable std::stack<EntryPosition> d_stack;
    
    public:
        PreprocessorStack();
    
        // Add preprocessors
    
        void begin_add(const Entry & pp);
        void end_add();
    
        // Stop adding phase; begin reduction phase
    
        void finalize_adding();
    
        // Reduction phase: begin/end preprocessing
    
        inline void begin_preprocess();
        inline void end_preprocess();
        inline bool is_toplevel() const;
    
        // Only allowed during begin_preprocess() and end_preprocess() calls:
    
        inline bool has_more() const;
        inline const Entry & current() const;
        inline void go_to_next();
    
        template<class Enumerator, // a enumerator object (ParallelEnumerator, SerialEnumerator, ...)
                 class CallbackFunction,
                 class Reorder>
        void preprocess(Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction> &,
                        Lattice<RealTypeContext, IntTypeContext> &, Reorder &);
    
        template<class Enumerator, // a enumerator object (ParallelEnumerator, SerialEnumerator, ...)
                 class CallbackFunction>
        inline void preprocess(Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction> &,
                               Lattice<RealTypeContext, IntTypeContext> &);
    };

    template<class RealTypeContext, class IntTypeContext>
    class ParallelEnumerator;
    template<class RealTypeContext, class IntTypeContext>
    class SerialEnumerator;
    template<class RealTypeContext, class IntTypeContext>
    class Sieve;
    template<class RealTypeContext, class IntTypeContext>
    class SimpleEnumerator;
    template<class RealTypeContext, class IntTypeContext>
    class SchnorrEnumerator;
    template<class RealTypeContext, class IntTypeContext>
    class VoronoiCellComputer;

    template<class RealTypeContext, class IntTypeContext>
    void StandardEnumBoundSelection(GaussianFactorComputer & gfc, Lattice<RealTypeContext, IntTypeContext> & lattice,
                                    typename RealTypeContext::Real & bound);

    template<class RealTypeContext, class IntTypeContext>
    class GaussianHeuristicEnumBoundSelection
    {
    private:
        typename RealTypeContext::Real d_factor;
    
    public:
        GaussianHeuristicEnumBoundSelection(RealTypeContext & rc, double factor = 1.05)
            : d_factor(factor, rc)
        {
        }
    
        void operator() (GaussianFactorComputer & gfc, Lattice<RealTypeContext, IntTypeContext> & lattice,
                         typename RealTypeContext::Real & bound) const;
    };

    template<class RealTypeContext, class IntTypeContext>
    class Enumerator
    {
    public:
        typedef PreprocessorStack<RealTypeContext, IntTypeContext, StandardEnumSettings> PreprocessorStackT;
        
    private:
        unsigned d_enumdim, d_maxthreads;
        LatticeReduction::EnumCallbackFunction d_ecf;
        LatticeReduction::EnumCallbackFunction_LI d_ecf2;
        
        ParallelEnumerator<RealTypeContext, IntTypeContext> * d_parallel_enum;
        SerialEnumerator<RealTypeContext, IntTypeContext> * d_serial_enum;
        Sieve<RealTypeContext, IntTypeContext> * d_sieve;
        SimpleEnumerator<RealTypeContext, IntTypeContext> * d_simple_enum;
        SchnorrEnumerator<RealTypeContext, IntTypeContext> * d_schnorr_enum;
        VoronoiCellComputer<RealTypeContext, IntTypeContext> * d_voronoi;
        
        enum { c_enum_threshold = 30 }; // all dimensions < 30 are enumerated by SimpleEnumerator<>
        
        void fullInit();
        
        Verbose & d_verbose;
        GaussianFactorComputer & d_gaussian_factors;
        
        PreprocessorStackT d_preprocessor_stack;
        
    public:
        
        template<class Enum, class CallbackF>
        Enumerator(Workspace<RealTypeContext, IntTypeContext, Enum, CallbackF> & workspace,
                   unsigned enumdimension, unsigned max_threads,
                   LatticeReduction::EnumCallbackFunction ecf = NULL, LatticeReduction::EnumCallbackFunction_LI ecf2 = NULL);
        
        bool enumerate(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned begin, unsigned end,
                       linalg::math_rowvector<typename IntTypeContext::Integer> & result,
                       typename RealTypeContext::Real & bound, LatticeReduction::SVPMode algorithm,
                       const StandardEnumSettings & settings, bool DontFallback = false);
        // Finds a shortest vector in the lattice generated by the orthogonal projections of the vectors
        // A.row(begin) to A.row(end) into the orthogonal complement of the vectors A.row(0) to
        // A.row(begin-1). Uses the Kannan-Schnorr-Euchner enumeration method.
        
        void setCallback(LatticeReduction::EnumCallbackFunction cf, LatticeReduction::EnumCallbackFunction_LI cf2);
        
        template<class CallbackFunction, class Reorder>
        inline void preprocess(Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction> & workspace,
                               Lattice<RealTypeContext, IntTypeContext> & lattice, Reorder & reorder)
        {
            d_preprocessor_stack.preprocess(workspace, lattice, reorder);
        }
        
        template<class CallbackFunction>
        inline void preprocess(Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction> & workspace,
                               Lattice<RealTypeContext, IntTypeContext> & lattice)
        {
            d_preprocessor_stack.preprocess(workspace, lattice);
        }
        
        template<class CallbackFunction>
        inline void computeEnumBound(Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction> & workspace,
                                     Lattice<RealTypeContext, IntTypeContext> & lattice, typename RealTypeContext::Real & bound) const
        {
            switch (d_preprocessor_stack.current().d_enumboundselectionmethod)
            {
            default:
            case EBS_Standard:
                StandardEnumBoundSelection(workspace.gaussianFactors(), lattice, bound);
                break;
            case EBS_GaussianHeuristic105:
                GaussianHeuristicEnumBoundSelection<RealTypeContext, IntTypeContext>(lattice.rc(), 1.05)
                    (workspace.gaussianFactors(), lattice, bound);
                break;
            }
        }
    
        inline unsigned getCurrentNumberOfEnumRepetitions() const
        {
            return d_preprocessor_stack.current().d_enum_repetitions;
        }
    
        inline LatticeReduction::SVPMode getCurrentAlgorithm() const
        {
            return d_preprocessor_stack.current().d_enumalg;
        }
    
        inline const StandardEnumSettings & getCurrentEnumSettings() const
        {
            return d_preprocessor_stack.current().d_enumsettings;
        }
    
        inline bool isToplevelEnum() const
        {
            return d_preprocessor_stack.is_toplevel();
        }
    
        inline void updateSettings(Lattice<RealTypeContext, IntTypeContext> & lattice)
        {
            if (d_preprocessor_stack.current().d_enum_settings_computer)
                d_preprocessor_stack.current().d_enum_settings_computer(d_preprocessor_stack.current().d_enumsettings, lattice);
        }

        inline PreprocessorStackT & getPreprocessorStack()
        {
            return d_preprocessor_stack;
        }
    };

    template<class RealTypeContext, class IntTypeContext>
    class MoveShortestFirst
    {
    private:
        LatticeReduction::Statistics & d_stats;
    
    public:
        MoveShortestFirst(LatticeReduction::Statistics & stats)
            : d_stats(stats)
        {
        }
    
        template<class Enumerator, class CallbackFunction>
        bool operator() (Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction> &,
                         Lattice<RealTypeContext, IntTypeContext> &, unsigned &, signed, bool, bool);
    };
}

#endif

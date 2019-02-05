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

#ifndef PLLL_INCLUDE_GUARD__ENUMIMPL_PREPROC_CPP
#define PLLL_INCLUDE_GUARD__ENUMIMPL_PREPROC_CPP

namespace plll
{
    template<class RealTypeContext, class IntTypeContext>
    template<class Enumerator, class CallbackFunction>
    bool MoveShortestFirst<RealTypeContext, IntTypeContext>::operator()
        (Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction> & workspace,
         Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned & k, signed modified_up_to,
         bool before_sizereduce, bool size_reduce_changed)
    {
        // Check if something should be done?
        if (!before_sizereduce && !size_reduce_changed)
            return false;
    
        // Do deep insertions
        if (modified_up_to < (signed)k)
            modified_up_to = (signed)k;
        unsigned new_k = k;
        bool modified = false, recompute_gs = false;
        typename RealTypeContext::Real res(lattice.rc());
        for (unsigned begin_k = k; begin_k <= (unsigned)modified_up_to; ++begin_k)
        {
            if (recompute_gs)
                lattice.update(begin_k + 1);
            lattice.computeProjectionLengthBV(res, lattice.range().begin(), begin_k);
            if (isZero(res))
                continue;
            if ((lattice.getNormSq(lattice.range().begin()) <= res) || isZero(lattice.getNormSq(lattice.range().begin())))
                continue;
        
            k = begin_k;
            if (k - lattice.range().begin() > 1)
                ++d_stats.deepinsertions;
            while (k > lattice.range().begin())
            {
                lattice.swap(k, k - 1);
                --k;
            }
            modified = true;
            if (k < new_k)
                new_k = k;
            recompute_gs = true;
        }
        k = new_k;
        return modified;
    }

    template<class RealTypeContext, class IntTypeContext, class EnumSettings>
    class PreprocessorStack<RealTypeContext, IntTypeContext, EnumSettings>::PaddedEntry
    {
    public:
        Entry d_entry;
        STD_AUTO_PTR<EntryCollector> d_subset;
        
        PaddedEntry()
            : d_subset(new EntryCollector())
        {
        }
    
        PaddedEntry(const Entry & e)
            : d_entry(e), d_subset(new EntryCollector())
        {
        }
    
        PaddedEntry(const PaddedEntry & pe)
            : d_entry(pe.d_entry), d_subset(new EntryCollector(*pe.d_subset))
        {
        }
    
        ~PaddedEntry()
        {
        }
    
        PaddedEntry & operator= (const PaddedEntry & pe)
        {
            if (&pe != this)
                d_entry.reset(new EntryCollector(*pe.d_subset));
            return *this;
        }
    };

    template<class RealTypeContext, class IntTypeContext, class EnumSettings>
    PreprocessorStack<RealTypeContext, IntTypeContext, EnumSettings>::PreprocessorStack()
        : d_can_add(true)
    {
        d_stack.push(std::make_pair(&d_root, d_root.end()));
    }

    template<class RealTypeContext, class IntTypeContext, class EnumSettings>
    void PreprocessorStack<RealTypeContext, IntTypeContext, EnumSettings>::begin_add(const Entry & pp)
    {
        assert(d_can_add);
        assert(!d_stack.empty());
        d_stack.top().first->push_back(PaddedEntry(pp));
        d_stack.push(std::make_pair(d_stack.top().first->back().d_subset.get(), d_stack.top().first->back().d_subset->end()));
    }

    template<class RealTypeContext, class IntTypeContext, class EnumSettings>
    void PreprocessorStack<RealTypeContext, IntTypeContext, EnumSettings>::end_add()
    {
        assert(d_can_add);
        assert(!d_stack.empty());
        d_stack.pop();
        assert(!d_stack.empty());
    }

    template<class RealTypeContext, class IntTypeContext, class EnumSettings>
    void PreprocessorStack<RealTypeContext, IntTypeContext, EnumSettings>::finalize_adding()
    {
        assert(d_can_add);
        assert(!d_stack.empty());
        d_stack.pop();
        assert(d_stack.empty());
        d_can_add = false;
    }

    template<class RealTypeContext, class IntTypeContext, class EnumSettings>
    inline void PreprocessorStack<RealTypeContext, IntTypeContext, EnumSettings>::begin_preprocess()
    {
        assert(!d_can_add);
        if (d_stack.empty())
            d_stack.push(EntryPosition(&d_root, d_root.begin()));
        else
            d_stack.push(EntryPosition(d_stack.top().second->d_subset.get(), d_stack.top().second->d_subset->begin()));
    }

    template<class RealTypeContext, class IntTypeContext, class EnumSettings>
    inline void PreprocessorStack<RealTypeContext, IntTypeContext, EnumSettings>::end_preprocess()
    {
        assert(!d_can_add);
        assert(!d_stack.empty());
        d_stack.pop();
    }

    template<class RealTypeContext, class IntTypeContext, class EnumSettings>
    inline bool PreprocessorStack<RealTypeContext, IntTypeContext, EnumSettings>::is_toplevel() const
    {
        return d_stack.size() == 1;
    }

    template<class RealTypeContext, class IntTypeContext, class EnumSettings>
    inline bool PreprocessorStack<RealTypeContext, IntTypeContext, EnumSettings>::has_more() const
    {
        return d_stack.top().second != d_stack.top().first->end();
    }

    template<class RealTypeContext, class IntTypeContext, class EnumSettings>
    inline const typename PreprocessorStack<RealTypeContext, IntTypeContext, EnumSettings>::Entry &
    PreprocessorStack<RealTypeContext, IntTypeContext, EnumSettings>::current() const
    {
        return d_stack.top().second->d_entry;
    }

    template<class RealTypeContext, class IntTypeContext, class EnumSettings>
    inline void PreprocessorStack<RealTypeContext, IntTypeContext, EnumSettings>::go_to_next()
    {
        ++d_stack.top().second;
    }

    template<class RealTypeContext, class IntTypeContext, class EnumSettings>
    template<class Enumerator, class CallbackFunction, class Reorder>
    void PreprocessorStack<RealTypeContext, IntTypeContext, EnumSettings>::preprocess
    (Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction> & workspace,
     Lattice<RealTypeContext, IntTypeContext> & lattice, Reorder & reorder)
    {
        if (workspace.verbose().yieldsOutput(LatticeReduction::VL_Chatter))
            workspace.verbose()(LatticeReduction::VL_Chatter) << "BEGIN PREPROCESS (dim=" << lattice.range().dimension() << ")";
        begin_preprocess();
        for (; has_more(); go_to_next())
        {
            const Entry & entry = current();
            // Check dimension
            if ((lattice.range().dimension() < entry.d_min_dim) || (lattice.range().dimension() > entry.d_max_dim))
                continue;
//            // Set minimal alpha
//            internalAlpha.setContext(lattice.rc());
//            internalAlpha = entry.d_minalpha;
//            if (internalAlpha < ralpha)
//                internalAlpha = ralpha;
            // Apply method
            if (workspace.verbose().yieldsOutput(LatticeReduction::VL_Chatter))
                workspace.verbose()(LatticeReduction::VL_Chatter) << "  DO PREPROCESS(method = " << entry.d_method << ", windowsize = " << entry.d_windowsize << ")";
            lattice.range().dupRange();
            switch (entry.d_method)
            {
            default:
            case Entry::RM_None:
                lattice.range().popRange();
                break;
            
            case Entry::RM_LLL:
            {
                NoAnnealLLL<RealTypeContext, IntTypeContext> noanneal;
                switch (entry.d_lll_mode)
                {
                case LatticeReduction::LLL_Classic:
                    workspace.partialLLL(lattice, lattice.range().begin() + 1, reorder, noanneal,
                                         Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::LovaszCondition);
                    break;
                case LatticeReduction::LLL_Unprojected:
                    workspace.partialLLL(lattice, lattice.range().begin() + 1, reorder, noanneal,
                                         Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::UnprojectedLovaszCondition);
                    break;
                case LatticeReduction::LLL_Siegel:
                    workspace.partialLLL(lattice, lattice.range().begin() + 1, reorder, noanneal,
                                         Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::SiegelCondition);
                    break;
                }
            }
            break;
        
            case Entry::RM_SVP:
                workspace.computeSVPbasis(lattice, true);
                break;
            case Entry::RM_SVPDual:
                workspace.computeDSVPbasis(lattice);
                break;
            case Entry::RM_HKZ:
                workspace.computeHKZbasis(lattice);
                break;
            case Entry::RM_HKZDual:
                workspace.computeDHKZbasis(lattice);
                break;
        
            case Entry::RM_BKZ:
            {
                // Start BKZ-style reduction (enumeration settings will be set on further solveSVP() calls)
                switch (entry.d_bkz_mode)
                {
                case LatticeReduction::BKZ_SchnorrEuchner:
                    workspace.partialBKZ_Classical(lattice, entry.d_windowsize, reorder);
                    break;
                case LatticeReduction::BKZ_Simplified:
                    workspace.partialBKZ_Simplified(lattice, entry.d_windowsize, reorder);
                    break;
                case LatticeReduction::BKZ_HanrotPujolStehleHKZ:
                case LatticeReduction::BKZ_HanrotPujolStehleSVP:
                    if (entry.d_force_tours)
                    {
                        arithmetic::Integer tours = entry.d_bkz_tours;
                        if (entry.d_bkz_tours_computer)
                            entry.d_bkz_tours_computer(tours, lattice);
                        workspace.partialBKZ_Terminating(lattice, entry.d_windowsize,
                                                         entry.d_bkz_mode == LatticeReduction::BKZ_HanrotPujolStehleHKZ, tours, reorder);
                    }
                    else
                    {
                        workspace.partialBKZ_Terminating(lattice, entry.d_windowsize,
                                                         entry.d_bkz_mode == LatticeReduction::BKZ_HanrotPujolStehleHKZ, reorder);
                    }
                    break;
                case LatticeReduction::BKZ_PrimalDual:
                    workspace.partialBKZ_PrimalDual(lattice, entry.d_windowsize, reorder);
                    break;
                case LatticeReduction::BKZ_SlideReduction:
                    workspace.partialBKZ_Slide(lattice, entry.d_windowsize, reorder);
                    break;
                case LatticeReduction::BKZ_ImprovedSlideReduction:
                    workspace.partialBKZ_ImprovedSlide(lattice, entry.d_windowsize,
                                                       Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::partialBKZ_ImprovedSlide_Original, reorder);
                    break;
                case LatticeReduction::BKZ_ImprovedSlideReduction2:
                    workspace.partialBKZ_ImprovedSlide(lattice, entry.d_windowsize,
                                                       Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::partialBKZ_ImprovedSlide_LargerDSVP, reorder);
                    break;
                case LatticeReduction::BKZ_ImprovedSlideReduction3:
                    workspace.partialBKZ_ImprovedSlide(lattice, entry.d_windowsize,
                                                       Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::partialBKZ_ImprovedSlide_LargerSVP, reorder);
                    break;
                case LatticeReduction::BKZ_SemiBlock2k:
                    workspace.partialBKZ_SemiBlock2k(lattice, entry.d_windowsize, reorder);
                    break;
                case LatticeReduction::BKZ_SamplingReduction:
                    workspace.partialBKZ_SamplingReduction(lattice, entry.d_windowsize, reorder);
                    break;
                case LatticeReduction::BKZ_Experimental:
                    workspace.partialBKZ_Experimental(lattice, entry.d_windowsize, reorder, 0 /* FIX TO 0 !!! ??? ... */);
                    break;
                }
            }
            break;
            }
        }
        end_preprocess();
        if (workspace.verbose().yieldsOutput(LatticeReduction::VL_Chatter))
            workspace.verbose()(LatticeReduction::VL_Chatter) << "END PREPROCESS";
    }
    
    template<class RealTypeContext, class IntTypeContext, class EnumSettings>
    template<class Enumerator, class CallbackFunction>
    inline void PreprocessorStack<RealTypeContext, IntTypeContext, EnumSettings>::
    preprocess(Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction> & workspace, Lattice<RealTypeContext, IntTypeContext> & lattice)
    {
        MoveShortestFirst<RealTypeContext, IntTypeContext> reorder(lattice.getStatistics());
        preprocess(workspace, lattice, reorder);
    }
}

#endif

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

#ifndef PLLL_INCLUDE_GUARD__LLL2_INTERNAL_HPP
#define PLLL_INCLUDE_GUARD__LLL2_INTERNAL_HPP

#include <plll/config.hpp>
#include <plll.hpp>

#include <cassert>
#include <vector>
#include <deque>
#include <sstream>

#include <plll/rational.hpp>

#if __cplusplus >= 201103L
  #define STD_AUTO_PTR std::unique_ptr
#else
  #define STD_AUTO_PTR std::auto_ptr
#endif

namespace plll
{
    void CallbackAdaptor(const linalg::math_matrix<arithmetic::Integer> & lattice, const LatticeReduction::CallbackFunction_LI & cf);
    void CallbackAdaptor_LI(const linalg::math_matrix<arithmetic::NInt<long int> > & lattice, const LatticeReduction::CallbackFunction & cf);
    void MinCallbackAdaptor(const linalg::math_matrix<arithmetic::Integer> & lattice, unsigned idx,
                            const arithmetic::Integer & len, const LatticeReduction::MinCallbackFunction_LI & cf);
    void MinCallbackAdaptor_LI(const linalg::math_matrix<arithmetic::NInt<long int> > & lattice, unsigned idx,
                               const arithmetic::NInt<long int> & len, const LatticeReduction::MinCallbackFunction & cf);
    void EnumCallbackAdaptor(const linalg::math_matrix<arithmetic::Integer> & lattice, int idx,
                             const linalg::math_rowvector<arithmetic::Integer> & vec, const LatticeReduction::EnumCallbackFunction_LI & cf);
    void EnumCallbackAdaptor_LI(const linalg::math_matrix<arithmetic::NInt<long int> > & lattice, int idx,
                                const linalg::math_rowvector<arithmetic::NInt<long int> > & vec, const LatticeReduction::EnumCallbackFunction & cf);
    
    class reduction_error : public LatticeReduction::stop_reduction
    {
    private:
        std::string d_message;
        
    public:
        reduction_error(const char * message) PLLL_INTERNAL_NOTHROW_POSTFIX_ENFORCE
            : d_message(message)
        {
        }
        
        reduction_error(const std::string & message) PLLL_INTERNAL_NOTHROW_POSTFIX_ENFORCE
            : d_message(message)
        {
        }
        
        virtual ~reduction_error() PLLL_INTERNAL_NOTHROW_POSTFIX_ENFORCE
        {
        }
        
        virtual const char * what() const PLLL_INTERNAL_NOTHROW_POSTFIX_ENFORCE;
    };
    
    class change_interface_exception : public LatticeReduction::stop_reduction
    {
    public:
        virtual ~change_interface_exception() PLLL_INTERNAL_NOTHROW_POSTFIX_ENFORCE
        {
        }
        
        virtual const char * what() const PLLL_INTERNAL_NOTHROW_POSTFIX_ENFORCE;
    };
    
    class feature_not_implemented : public LatticeReduction::stop_reduction
    {
    public:
        virtual ~feature_not_implemented() PLLL_INTERNAL_NOTHROW_POSTFIX_ENFORCE
        {
        }
        
        virtual const char * what() const PLLL_INTERNAL_NOTHROW_POSTFIX_ENFORCE;
    };
    
    template<class RealTypeContext, class IntTypeContext>
    class Lattice;
    
    template<class RealTypeContext, class IntTypeContext>
    class GSInterface
    {
    public:
        virtual ~GSInterface() { }
        
        virtual std::pair<linalg::math_matrix<typename RealTypeContext::Real> *, linalg::math_rowvector<typename RealTypeContext::Real> *> getStorage() = 0;
        virtual void update(unsigned i) = 0; // make sure all indices < i are available
        virtual void reset() = 0;
        virtual void adjustAlpha(double LLLalpha) = 0;
        
        virtual unsigned getDimension() const = 0;
        virtual bool isRowZero(unsigned i) const = 0;
        virtual const typename RealTypeContext::Real & getZDF() const = 0;
        virtual const typename RealTypeContext::Real & getLLLalpha() const = 0;
        
        virtual void swap(unsigned i, unsigned j) = 0;
        virtual void add(unsigned i, const typename IntTypeContext::Integer & m, unsigned j) = 0;
        virtual void flip(unsigned i) = 0;
        virtual void trans(unsigned i, unsigned j,
                           const typename IntTypeContext::Integer & B00, const typename IntTypeContext::Integer & B01,
                           const typename IntTypeContext::Integer & B10, const typename IntTypeContext::Integer & B11) = 0;
        virtual void insertVectorLC(unsigned ofs, const linalg::math_rowvector<typename IntTypeContext::Integer> & result) = 0;
        virtual void insertVector(unsigned ofs, const linalg::math_rowvector<typename IntTypeContext::Integer> & result) = 0;
        virtual void removeZeroVector(unsigned ofs) = 0;
        virtual void compactify() = 0;
        
        virtual bool canGetIntegerVector() const = 0;
        virtual void getLinearCombination(linalg::math_rowvector<typename IntTypeContext::Integer> & result, unsigned ofs, const linalg::math_rowvector<typename IntTypeContext::Integer> & lincomb) const = 0;
        // Will only do something if canGetIntegerVector() returns true.
        virtual void getVector(linalg::math_rowvector<typename IntTypeContext::Integer> & result, unsigned ofs) const = 0;
        // Will only do something if canGetIntegerVector() returns true.
        virtual std::pair<linalg::math_matrix<typename IntTypeContext::Integer> *, unsigned> getMatrixIndex(unsigned ofs) const = 0;
        // Returns matrix (if available) and index of specified vector <ofs> (might be different from
        // <ofs>). Will only do something if canGetIntegerVector() returns true.
        virtual const typename IntTypeContext::Integer & getNormSqUP(unsigned i) const = 0;
        // Will only do something if canGetIntegerVector() returns true.
        virtual void getNormSqUP(typename IntTypeContext::Integer & r, unsigned i) const = 0;
        // Will only do something if canGetIntegerVector() returns true.
        virtual void getDotProduct(typename IntTypeContext::Integer & r, unsigned i, unsigned j) const = 0;
        // Will only do something if canGetIntegerVector() returns true.
        
        virtual void computeDotProductProjected(typename RealTypeContext::Real & r, unsigned k, unsigned i, unsigned j) const = 0;
        // Computes the dot product of the projections of B_i and B_j onto the
        // orthogonal complement of B_0, ..., B_{k-1}. (Here, B_0, ..., B_{n-1} is the basis of the
        // lattice.)
        virtual void computeProjectionLength(typename RealTypeContext::Real & result, unsigned k, unsigned b, const linalg::math_rowvector<typename IntTypeContext::Integer> & vec) const = 0;
        // Computes the squared length of the projection of \sum_{i=0}^{m-1} vec[i] B_{b+i} onto the
        // orthogonal complement of B_0, ..., B_{k-1}. (Here, B_0, ..., B_{n-1} is the basis of the
        // lattice.)
        virtual void computeProjectionLengthBV(typename RealTypeContext::Real & result, unsigned k, unsigned b) const = 0;
        // Computes the squared length of the projection of B_b onto the orthogonal complement of B_0, ...,
        // B_{k-1}. (Here, B_0, ..., B_{n-1} is the basis of the lattice.)
        virtual void computeProjectionLength(typename RealTypeContext::Real & result, unsigned k, const linalg::math_rowvector<typename IntTypeContext::Integer> & vec) const = 0;
        
        virtual bool sizereduce(Lattice<RealTypeContext, IntTypeContext> &, unsigned i) = 0;
        
        virtual GSInterface * clone(RealTypeContext &, IntTypeContext &) = 0;
        virtual void copyFrom(GSInterface * source) = 0;
        
        virtual void changeOfPrecision() = 0; // inform that RealTypeContext's precision changed; has to update the matrix/vector set by setStorage()
        virtual bool avoidInsertions() = 0; // inform that insertions/removals should be avoided if possible
        
        virtual void print(std::ostream &) = 0; // can print more internal details if wanted 
    };
    
    template<class IntTypeContext>
    class TransformNotifier
    {
    public:
        virtual ~TransformNotifier() { }
        
        virtual void swap(unsigned i, unsigned j) = 0;
        virtual void add(unsigned i, const typename IntTypeContext::Integer & m, unsigned j) = 0;
        virtual void flip(unsigned i) = 0;
        virtual void trans(unsigned i, unsigned j,
                           const typename IntTypeContext::Integer & B00, const typename IntTypeContext::Integer & B01,
                           const typename IntTypeContext::Integer & B10, const typename IntTypeContext::Integer & B11) = 0;
        virtual void insertVectorLC(unsigned ofs, const linalg::math_rowvector<typename IntTypeContext::Integer> & result) = 0;
        virtual void insertVector(unsigned ofs, const linalg::math_rowvector<typename IntTypeContext::Integer> & result) = 0;
        virtual bool canInsertVector() = 0;
        virtual void removeZeroVector(unsigned ofs) = 0;
        virtual void compactify() = 0;
    };
    
    class Ranger
    // Implemented in lattice.cpp
    {
    public:
        typedef std::pair<unsigned, unsigned> Range;
        
    private:
        std::deque<Range> d_ranges;
        
    public:
        inline bool empty() const;
        inline unsigned begin() const;
        inline unsigned end() const;
        inline unsigned dimension() const;
        
        inline void popRange();
        inline Ranger & addRange(unsigned b, unsigned e);
        inline Ranger & addRange(const std::pair<unsigned, unsigned> & b);
        inline Ranger & dupRange();
        
        void insertVector(unsigned ofs);
        void removeVector(unsigned ofs);
        
        inline friend std::ostream & operator << (std::ostream &, const Ranger &);
    };
    
    template<class RealTypeContext, class IntTypeContext>
    class Lattice
    // Implemented in lattice.cpp
    {
    private:
        RealTypeContext & d_rc;
        IntTypeContext & d_ic;
        Ranger d_range;
        linalg::math_matrix<typename RealTypeContext::Real> * d_coeffs;
        linalg::math_rowvector<typename RealTypeContext::Real> * d_sqnorms;
        GSInterface<RealTypeContext, IntTypeContext> * d_gs;
        LatticeReduction::Statistics & d_stats;
        std::deque<TransformNotifier<IntTypeContext>*> d_notifier;
        bool d_do_not_modify_flag;
        
        Lattice(GSInterface<RealTypeContext, IntTypeContext> * gs, RealTypeContext & rc, IntTypeContext & ic, const Lattice & lattice, bool own_stats, LatticeReduction::Statistics & stats);
        
    public:
        typedef RealTypeContext RealContext;
        
        class DuplicateStorage;
        friend class DuplicateStorage;
        
        class DuplicateStorage
        {
            friend class Lattice<RealTypeContext, IntTypeContext>;
            
        private:
            RealTypeContext d_rc;
            IntTypeContext d_ic;
            STD_AUTO_PTR<GSInterface<RealTypeContext, IntTypeContext> > d_gs;
            LatticeReduction::Statistics d_stats;
            
            // Prohibit copying
            DuplicateStorage(const DuplicateStorage &);
            DuplicateStorage & operator = (const DuplicateStorage &);
            
        public:
            DuplicateStorage(const Lattice<RealTypeContext, IntTypeContext> & source, bool has_own_stats = false)
                : d_rc(source.rc()), d_ic(source.ic()),
                  d_gs(source.gs()->clone(d_rc, d_ic)),
                  d_lattice(d_gs.get(), d_rc, d_ic, source, has_own_stats, d_stats)
            {
            }
            
            Lattice<RealTypeContext, IntTypeContext> d_lattice;
        };
        
        void copyFrom(const DuplicateStorage & storage, bool add_own_stats = false);
        
        Lattice(GSInterface<RealTypeContext, IntTypeContext> * gs, RealTypeContext & rc, IntTypeContext & ic, LatticeReduction::Statistics & stats, unsigned begin, unsigned end, bool do_not_modify_flag = false);
        ~Lattice();
        
        /////////////////////////////////////////////
        // Get information
        
        inline LatticeReduction::Statistics & getStatistics();
        inline Ranger & range();
        inline const Ranger & range() const;
        inline unsigned dimension() const;
        inline GSInterface<RealTypeContext, IntTypeContext> * gs() const;
        inline RealTypeContext & rc();
        inline const RealTypeContext & rc() const;
        inline IntTypeContext & ic();
        inline const IntTypeContext & ic() const;
        inline const typename RealTypeContext::Real & getZDF() const;
        inline const typename RealTypeContext::Real & getLLLalpha() const;
        inline LatticeReduction::Statistics & getStats();
        
        inline bool allowedToModify() const;
        inline bool avoidInsertions() const; // Returns true if insertions/removals of vectors are
        // supported, but should be avoided since they cause
        // numerical problems
        
        /////////////////////////////////////////////
        // Manage transformation notifiers
        
        void addNotifier(TransformNotifier<IntTypeContext> * t);
        void removeNotifier(TransformNotifier<IntTypeContext> * tt);
        void removeAllNotifiers();
        
        /////////////////////////////////////////////
        // Get information on Gram-Schmidt coefficients
        
        inline void update(unsigned i);
        inline void reset();
        
        inline const typename RealTypeContext::Real & getCoeff(unsigned i, unsigned j) const;
        inline const typename RealTypeContext::Real & getNormSq(unsigned i) const;
        
        /////////////////////////////////////////////
        // Get information on vectors
        
        inline bool isRowZero(unsigned i) const;
        
        /////////////////////////////////////////////
        // Modify vectors
        
        void swap(unsigned i, unsigned j);
        
        void add(unsigned i, const typename IntTypeContext::Integer & m, unsigned j);
        
        void flip(unsigned i);
        
        void trans(unsigned i, unsigned j,
                   const typename IntTypeContext::Integer & B00, const typename IntTypeContext::Integer & B01,
                   const typename IntTypeContext::Integer & B10, const typename IntTypeContext::Integer & B11);
        // Assuming that the determinant is 1:
        // [ B00 B01 ]   [ ..... < row i > ..... ]
        // [ B10 B11 ] * [ ..... < row j > ..... ]
        
        void insertVectorLC(unsigned ofs, const linalg::math_rowvector<typename IntTypeContext::Integer> & result);
        // Inserts a new vector at index ofs, which will equal \sum_{i=0}^{result.size()-1} result[i] * b[ofs+i]
        
        void insertVector(unsigned ofs, const linalg::math_rowvector<typename IntTypeContext::Integer> & result);
        // Inserts a new vector <result> at index ofs
        
        void removeZeroVector(unsigned ofs);
        // Removes the zero vector at position ofs
        
        void compactify();
        
        inline bool sizereduce(unsigned i);
        
        /////////////////////////////////////////////
        // Compute lengths of projections
        
        inline void computeDotProductProjected(typename RealTypeContext::Real & r, unsigned k, unsigned i, unsigned j) const;
        // Computes the dot product of the projections of B_i and B_j onto the
        // orthogonal complement of B_0, ..., B_{k-1}. (Here, B_0, ..., B_{n-1} is the basis of the
        // lattice.)
        
        inline void computeProjectionLength(typename RealTypeContext::Real & result, unsigned k, unsigned b, const linalg::math_rowvector<typename IntTypeContext::Integer> & vec) const;
        // Computes the squared length of the projection of \sum_{i=0}^{m-1} vec[i] B_{b+i} onto the
        // orthogonal complement of B_0, ..., B_{k-1}. (Here, B_0, ..., B_{n-1} is the basis of the
        // lattice.)
        
        inline void computeProjectionLengthBV(typename RealTypeContext::Real & result, unsigned k, unsigned b) const;
        // Computes the squared length of the projection of B_b onto the orthogonal complement of B_0, ...,
        // B_{k-1}. (Here, B_0, ..., B_{n-1} is the basis of the lattice.)
        
        inline void computeProjectionLength(typename RealTypeContext::Real & result, unsigned k, const linalg::math_rowvector<typename IntTypeContext::Integer> & vec) const;
        
        /////////////////////////////////////////////
        // Get vectors (if possible). These functions only work if canGetIntegerVector() returns true.
        
        inline bool canGetIntegerVector() const;
        inline void getLinearCombination(linalg::math_rowvector<typename IntTypeContext::Integer> & result, unsigned ofs, const linalg::math_rowvector<typename IntTypeContext::Integer> & lincomb) const;
        inline void getVector(linalg::math_rowvector<typename IntTypeContext::Integer> & result, unsigned ofs) const;
        inline linalg::math_matrix<typename IntTypeContext::Integer> * getMatrix() const;
        inline std::pair<linalg::math_matrix<typename IntTypeContext::Integer> *, unsigned> getMatrixIndex(unsigned ofs) const;
        // Returns matrix (if available) and index of specified vector <ofs> (might be different from
        // <ofs>).
        inline const typename IntTypeContext::Integer & getNormSqUP(unsigned i) const;
        inline void getNormSqUP(typename IntTypeContext::Integer & r, unsigned i) const;
        // Computes the squared norm (unprojected)
        inline void getDotProduct(typename IntTypeContext::Integer & r, unsigned i, unsigned j) const;
        // Computes arbitrary dot products (unprojected)
        
        /////////////////////////////////////////////
        // Misc
        
        inline void changeOfPrecision(); // inform that RealTypeContext's precision changed
    };
    
    template<class RealTypeContext, class IntTypeContext>
    struct LatticeHelper
    {
    public:
        class RangePopper
        {
        private:
            Lattice<RealTypeContext, IntTypeContext> & d_lattice;
            
        public:
            inline RangePopper(Lattice<RealTypeContext, IntTypeContext> & lattice) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_lattice(lattice)
            {
            }
            
            inline ~RangePopper() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                if (!d_lattice.range().empty()) // avoid exceptions!
                    d_lattice.range().popRange();
            }
        };
        
        class RangePopper2
        {
        private:
            Lattice<RealTypeContext, IntTypeContext> & d_lattice;
            bool d_popped;
            
        public:
            inline RangePopper2(Lattice<RealTypeContext, IntTypeContext> & lattice) PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
                : d_lattice(lattice), d_popped(false)
            {
            }
            
            void pop()
            {
                if (!d_popped)
                {
                    d_lattice.range().popRange();
                    d_popped = true;
                }
            }
            
            inline ~RangePopper2() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                if (!d_popped && !d_lattice.range().empty()) // avoid exceptions!
                    d_lattice.range().popRange();
            }
        };
        
        class AddNotifier
        {
        private:
            Lattice<RealTypeContext, IntTypeContext> & d_lattice;
            TransformNotifier<IntTypeContext> * d_notifier;
            
        public:
            inline AddNotifier(Lattice<RealTypeContext, IntTypeContext> & lattice, TransformNotifier<IntTypeContext> * notifier)
                : d_lattice(lattice), d_notifier(notifier)
            {
                d_lattice.addNotifier(d_notifier);
            }
            
            inline ~AddNotifier() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
            {
                try
                {
                    d_lattice.removeNotifier(d_notifier);
                }
                catch (...)
                {
                }
            }
        };
    };
    
    template<class RealTypeContext, class IntTypeContext>
    std::ostream & operator << (std::ostream &, const Lattice<RealTypeContext, IntTypeContext> &); // (implemented in lattice.cpp)
    
    class GaussianFactorComputer
    {
    private:
        mutable std::vector<long double> d_hermite_prefactors;
    
    public:
        void computeHPF(unsigned dim); // (implemented in lll2.cpp)
        // Make sure HPF(dim) is available in cache.
    
        inline long double HPF(unsigned dim)
        {
            if (d_hermite_prefactors.size() <= dim)
                computeHPF(dim);
            return d_hermite_prefactors[dim];
        }
    };

    class Verbose
    {
    private:
        LatticeReduction::VerboseOutputLevel d_outputlevel;
        LatticeReduction::VerboseFunction d_function;
    
        mutable std::ostringstream d_s;
    
    public:
        class VerboseStream
        {
            friend class Verbose;
        
        private:
            LatticeReduction::VerboseFunction d_f;
            LatticeReduction::VerboseLevel d_level;
            std::ostringstream & d_s;
            
            inline VerboseStream(LatticeReduction::VerboseFunction f, LatticeReduction::VerboseLevel level, std::ostringstream & s)
                : d_f(f), d_level(level), d_s(s)
            {
            }
        
        public:
            inline ~VerboseStream()
            {
                if (d_f)
                    d_f(d_level, d_s.str());
            }
        
            template<class T>
            VerboseStream & operator << (const T & v)
            {
                d_s << v;
                return *this;
            }
        };
        
        inline Verbose(LatticeReduction::VerboseOutputLevel outputlevel, LatticeReduction::VerboseFunction function)
            : d_outputlevel(outputlevel), d_function(function)
        {
        }
        
        inline Verbose(const Verbose & v)
            : d_outputlevel(v.d_outputlevel), d_function(v.d_function)
        {
        }
        
        void setup(LatticeReduction::VerboseOutputLevel outputlevel, LatticeReduction::VerboseFunction function)
        {
            d_outputlevel = outputlevel;
            d_function = function;
        }
    
        static inline bool yieldsOutput(LatticeReduction::VerboseOutputLevel outputlevel, LatticeReduction::VerboseLevel level)
        {
//    enum VerboseOutputLevel { VOL_None, VOL_Warnings, VOL_Informative, VOL_Full };
//    enum VerboseLevel { VL_Error, VL_Warning, VL_Information, VL_Chatter };
            static bool VOB_VL[4][4] = { { false, false, false, false },
                                         { true, true, false, false },
                                         { true, true, true, false },
                                         { true, true, true, true } };
            return VOB_VL[outputlevel][level];
        }
    
        inline bool yieldsOutput(LatticeReduction::VerboseLevel level) const
        {
            return yieldsOutput(d_outputlevel, level);
        }
    
        inline VerboseStream operator() (LatticeReduction::VerboseLevel level) const
        {
            d_s.str(std::string());
            return VerboseStream(yieldsOutput(level) ? d_function : NULL, level, d_s);
        }
    
        inline VerboseStream operator() () const
        {
            d_s.str(std::string());
            return VerboseStream(d_function, LatticeReduction::VL_Chatter, d_s);
        }
    };
    
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    class Workspace
    {
    private:
        GaussianFactorComputer d_gaussian_factors; // initialized with lattice dimension, i.e. should be thread-safe
        Verbose d_verbose; // verbose output
        CallbackFunction d_callback; // callback
        Enumerator d_enum; // enumerator
        
    public:
        Workspace(const Verbose & verbose, const CallbackFunction & callback,
                  unsigned enumdimension, unsigned max_threads,
                  LatticeReduction::EnumCallbackFunction ecf = NULL,
                  LatticeReduction::EnumCallbackFunction_LI ecf2 = NULL)
            : d_verbose(verbose), d_callback(callback),
              d_enum(*this, enumdimension, max_threads, ecf, ecf2)
        {
            d_gaussian_factors.computeHPF(enumdimension);
        }
        
        Workspace(LatticeReduction::VerboseOutputLevel verbose_outputlevel, LatticeReduction::VerboseFunction verbose_function,
                  const CallbackFunction & callback,
                  unsigned enumdimension, unsigned max_threads,
                  LatticeReduction::EnumCallbackFunction ecf = NULL, LatticeReduction::EnumCallbackFunction_LI ecf2 = NULL)
            : d_verbose(verbose_outputlevel, verbose_function), d_callback(callback),
              d_enum(*this, enumdimension, max_threads, ecf, ecf2)
        {
            d_gaussian_factors.computeHPF(enumdimension);
        }
        
        inline Verbose & verbose()
        {
            return d_verbose;
        }
        
        inline Enumerator & enumerator()
        {
            return d_enum;
        }
        
        inline void callback(Lattice<RealTypeContext, IntTypeContext> & l, unsigned i)
        {
            d_callback(l, i);
        }
        
        inline CallbackFunction & getCallbackObject()
        {
            return d_callback;
        }
        
        inline long double HPF(unsigned dim)
        {
            return d_gaussian_factors.HPF(dim);
        }
        
        inline GaussianFactorComputer & gaussianFactors()
        {
            return d_gaussian_factors;
        }
        
        ////////////////////////////////////////////////////////////////////////
        // Lattice helpers
        
        GSInterface<RealTypeContext, IntTypeContext> * getDualGSI(RealTypeContext & rc, IntTypeContext & ic,
                                                                  const Lattice<RealTypeContext, IntTypeContext> & lattice); // (implemented in lattice.cpp)
        
        void applyRandomUnimodularTransformation(Lattice<RealTypeContext, IntTypeContext> & lattice);
        
        void rearrange(Lattice<RealTypeContext, IntTypeContext> & lattice, linalg::math_rowvector<typename IntTypeContext::Integer> & result,
                       bool minimizeCoeffs = true);
        // Rearranges A.row(begin) to A.row(end) such that the vector described by the linear combination in
        // result (with end-begin+1 entries) is at position A.row(begin). Note that during the process,
        // result might be arbitrarily changed. Also note that the GCD of the entries of result must be 1.
        //
        // Uses trans(), swap() and flip().
        
        void rearrangeDual(Lattice<RealTypeContext, IntTypeContext> & lattice, linalg::math_rowvector<typename IntTypeContext::Integer> & result,
                           bool minimizeCoeffs = true);
        // Rearranges A.row(begin) to A.row(end) such that in the dual of the projected sublattice of these
        // vectors, the first vector equals the given linear combination of the dual lattice's "canonical"
        // reversed basis. Note that during the process, result might be arbitrarily changed. Also note that
        // the GCD of the entries of result must be 1.
        // 
        // Uses trans(), swap() and flip().
        
        template<class Reorderer>
        void solveSVPOnce(Lattice<RealTypeContext, IntTypeContext> & lattice,
                          linalg::math_rowvector<typename IntTypeContext::Integer> & result, Reorderer & reorder);
        
        template<class Reorderer>
        void solveSVPRepeat(Lattice<RealTypeContext, IntTypeContext> & lattice,
                            linalg::math_rowvector<typename IntTypeContext::Integer> & result, Reorderer & reorder, unsigned number_of_iterations);
        
        ////////////////////////////////////////////////////////////////////////
        // Helper functions
        
        static void computeProjectedDeterminant(typename RealTypeContext::Real & result, const Lattice<RealTypeContext, IntTypeContext> & lattice,
                                                const std::pair<unsigned, unsigned> & block);
        
        // Arguments for partialLLL:
        static bool LovaszCondition(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned k);
        static bool UnprojectedLovaszCondition(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned k);
        static bool SiegelCondition(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned k);
        
        enum ImprovedSlideReturn { ISR_Continue = true, ISR_Stop };
        
        // Arguments for partialBKZ_ImprovedSlide:
        static ImprovedSlideReturn partialBKZ_ImprovedSlide_Original(Workspace & workspace, unsigned ell,
                                                                     Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize);
        static ImprovedSlideReturn partialBKZ_ImprovedSlide_LargerDSVP(Workspace & workspace, unsigned ell,
                                                                       Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize);
        static ImprovedSlideReturn partialBKZ_ImprovedSlide_LargerSVP(Workspace & workspace, unsigned ell,
                                                                      Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize);
        
        ////////////////////////////////////////////////////////////////////////
        // Algorithms
        
        void sizereduction(Lattice<RealTypeContext, IntTypeContext> & lattice);
        // Applies size reduction: A.row(begin) to A.row(end) will be size-reduced.
        
        template<class Reorderer, class Annealer, class LCondition>
        unsigned partialLLL(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned beginstage,
                            Reorderer & reorder, Annealer & anneal, LCondition & lcondition);
        // Applies LLL to a subset of the vectors.
        
        template<class Reorderer>
        unsigned partialLLL(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned beginstage, Reorderer & reorder);
        
        template<class Reorderer>
        unsigned partialLLL2(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned beginstage, Reorderer & reorder);
        
        template<class Reorderer>
        unsigned partialLLL3(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned beginstage, Reorderer & reorder);
        
        template<class AnnealFunction, class Reorderer>
        void partialLLL_Anneal(Lattice<RealTypeContext, IntTypeContext> & lattice, LatticeReduction::AnnealCallbackFunction ACF,
                               AnnealFunction AF, Reorderer & reorder);
        
        template<class Reorderer, class Annealer>
        void partialBKZ_Classical(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize, Reorderer & reorder, Annealer & anneal);
        // Applies Classical (Schnorr-Euchner) BKZ to a subset of the vectors.
        
        template<class Reorderer>
        void partialBKZ_Classical(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize, Reorderer & reorder);
        
        template<class Reorderer, class Annealer>
        void partialBKZ_Simplified(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize, Reorderer & reorder, Annealer & anneal);
        // Applies Simplified (Hanrot-Pujol-Stehle) BKZ to a subset of the vectors.
        
        template<class Reorderer>
        void partialBKZ_Simplified(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize, Reorderer & reorder);
        
        template<class Reorderer>
        void partialBKZ_Terminating(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize, bool doHKZ,
                                    const arithmetic::Integer & tours, Reorderer & reorder);
        // doHKZ: true = apply HKZ to local bases, false = apply SVP to local bases
        // tours: number of tours (i.e. runs)
        //
        // Applies Terminating BKZ (Hanrot-Pujol-Stehle) to a subset of the vectors.
        
        template<class Reorderer>
        void partialBKZ_Terminating(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize, bool doHKZ, Reorderer & reorder);
        // doHKZ: true = apply HKZ to local bases, false = apply SVP to local bases
        //
        // Applies Terminating BKZ (Hanrot-Pujol-Stehle) to a subset of the vectors.
        
        template<class Reorderer, class Annealer>
        void partialBKZ_SemiBlock2k(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize, Reorderer & reorder, Annealer & anneal);
        
        template<class Reorderer>
        void partialBKZ_SemiBlock2k(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize, Reorderer & reorder);
        
        template<class Reorderer, class Annealer>
        void partialBKZ_PrimalDual(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize, Reorderer & reorder, Annealer & anneal);
        
        template<class Reorderer>
        void partialBKZ_PrimalDual(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize, Reorderer & reorder);
        
        template<class Reorderer>
        void partialBKZ_PrimalDualBHO(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize, Reorderer & reorder);
        
        template<class Reorderer, class Annealer>
        void partialBKZ_Slide(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize, Reorderer & reorder, Annealer & anneal);
        // Applies Gama-Nguyen Slide Reduction to a subset of the vectors.
        
        template<class Reorderer>
        void partialBKZ_Slide(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize, Reorderer & reorder);
        
        template<class Reduction, class Reorderer, class Annealer>
        void partialBKZ_ImprovedSlide(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize,
                                      Reduction reduction, Reorderer & reorder, Annealer & anneal);
        // Applies Schnorr's Improved Slide Reduction to a subset of the vectors.
        
        template<class Reduction, class Reorderer>
        void partialBKZ_ImprovedSlide(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize,
                                      Reduction reduction, Reorderer & reorder);
        
        template<class Reorderer, class Annealer>
        void partialBKZ_SamplingReduction(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize,
                                          Reorderer & reorder, Annealer & anneal);
        
        template<class Reorderer>
        void partialBKZ_SamplingReduction(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize, Reorderer & reorder);
        
        template<class Reorderer, class Annealer>
        void partialBKZ_Experimental(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize,
                                     Reorderer & reorder, Annealer & anneal, unsigned id);
        
        template<class Reorderer>
        void partialBKZ_Experimental(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize,
                                     Reorderer & reorder, unsigned id = 0);
        
        template<class AnnealFunction, class Reorderer>
        void partialBKZ_Anneal(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize,
                               LatticeReduction::AnnealCallbackFunction ACF, AnnealFunction AF, Reorderer & reorder, bool simplified);
        
        template<class Reorderer>
        void solveSVP(Lattice<RealTypeContext, IntTypeContext> & lattice, linalg::math_rowvector<typename IntTypeContext::Integer> & result,
                      Reorderer & reorder, bool dontModify = false);
        // Finds a shortest vector in the lattice generated by the orthogonal projections of the vectors
        // A.row(begin) to A.row(end) into the orthogonal complement of the vectors A.row(0) to
        // A.row(begin-1).
        //
        // In case no shortest vector was found, result is a vector of length 0.
        
        template<class Reorderer>
        unsigned rearrangeLLL(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned reducebegin,
                              linalg::math_rowvector<typename IntTypeContext::Integer> & result, Reorderer & reorder);
        // Inserts a linear combination of the vectors A.row(begin), ..., A.row(begin+result.size()-1)
        // before A.row(begin) and applies LLL on the range A.row(reducebegin), ..., A.row(end+1) to remove
        // the induced dependency.
        
        template<class Reorderer>
        inline void solveSVPRearrange(Lattice<RealTypeContext, IntTypeContext> & lattice, Reorderer & reorder, bool minimizeCoeffs = true);
        // Finds a shortest vector in the lattice generated by the orthogonal projections of the vectors
        // A.row(begin) to A.row(end) into the orthogonal complement of the vectors A.row(0) to
        // A.row(begin-1). Modifies the basis so that A.row(begin) is said vector. Only A.row(begin) to
        // A.row(end) are modified.
        
        template<class Reorderer>
        inline unsigned solveSVPRearrangeLLL(Lattice<RealTypeContext, IntTypeContext> & lattice,
                                             unsigned reducebegin, Reorderer & reorder);
        // Finds a shortest vector in the lattice generated by the orthogonal projections of the vectors
        // A.row(begin) to A.row(end) into the orthogonal complement of the vectors A.row(0) to
        // A.row(begin-1). Modifies the basis so that A.row(begin) is said vector. Only A.row(begin) to
        // A.row(end) are modified.
        
        inline void computeSVPbasis(Lattice<RealTypeContext, IntTypeContext> & lattice, bool make_basis);
        // Computes a "SVP basis" (i.e. first vector is shortest) for the lattice generated by the orthogonal
        // projections of the vectors A.row(begin) to A.row(end) into the orthogonal complement of the
        // vectors A.row(0) to A.row(begin-1). Modifies the basis so that A.row(begin) to A.row(end) is said
        // basis (or more precisely, its preimage).
        
        inline bool computeDSVPbasis(Lattice<RealTypeContext, IntTypeContext> & lattice);
        // Computes a "Dual-SVP basis" (i.e. first vector of canonical reversed projected dual is shortest)
        // for the lattice generated by the orthogonal projections of the vectors A.row(begin) to A.row(end)
        // into the orthogonal complement of the vectors A.row(0) to A.row(begin-1). Modifies the basis so
        // that A.row(begin) to A.row(end) is said basis (or more precisely, its preimage).
        //
        // The returned bool indicates whether a new shortest vector was inserted (true) or not
        // (false).
        
        inline void computeHKZbasis(Lattice<RealTypeContext, IntTypeContext> & lattice);
        // Computes a Hermite-Korkine-Zolotarev basis for the lattice generated by the orthogonal
        // projections of the vectors A.row(begin) to A.row(end) into the orthogonal complement of the
        // vectors A.row(0) to A.row(begin-1). Modifies the basis so that A.row(begin) to A.row(end) is said
        // basis (or more precisely, its preimage).
        
        inline void computeDHKZbasis(Lattice<RealTypeContext, IntTypeContext> & lattice);
        // Computes a "Dual-HKZ basis" (i.e. canonical reversed projected dual is HKZ reduced) for the
        // lattice generated by the orthogonal projections of the vectors A.row(begin) to A.row(end) into
        // the orthogonal complement of the vectors A.row(0) to A.row(begin-1). Modifies the basis so that
        // A.row(begin) to A.row(end) is said basis (or more precisely, its preimage).
    };
    
    template<class RealTypeContext, class IntTypeContext>
    class EmptyCallback
    {
    public:
        inline void operator() (Lattice<RealTypeContext, IntTypeContext> &, unsigned) const
        {
        }
    };
    
    template<class RealTypeContext, class IntTypeContext>
    class NoReorder
    {
    public:
        template<class Enumerator, class CallbackFunction>
        bool operator() (Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction> &,
                         Lattice<RealTypeContext, IntTypeContext> &, unsigned &, signed,
                         bool /*before_sizereduce*/, bool /*size_reduce_changed*/) const
        {
            return false;
        }
    };
    
    template<class RealTypeContext, class IntTypeContext>
    class NoAnnealLLL
    {
    public:
        bool operator() (Lattice<RealTypeContext, IntTypeContext> &, int) const
        {
            return false;
        }
    };
    
    template<class RealTypeContext, class IntTypeContext>
    class NoAnnealBKZ
    {
    public:
        bool operator() (Lattice<RealTypeContext, IntTypeContext> & lattice, int k, int windowsize,
                         linalg::math_rowvector<typename IntTypeContext::Integer> & lincomb) const
        {
            return false;
        }
    };
    
    typedef boost::function<void(unsigned dimension, unsigned maxbits)> MaxBitsCallbackFunction;
    /* Such callback functions will be called when the maximal bit count of all entries of the
       matrix changed. */
}

#include "matrixconversion.hpp"

namespace plll
{
    template<class IntTypeContext>
    class LLLAnnealWrapper
    {
    private:
        LatticeReduction::LLL_AnnealFunction d_af;
        LatticeReduction::GramSchmidtInformer * d_gsi;
        
    public:
        LLLAnnealWrapper(LatticeReduction::LLL_AnnealFunction af, LatticeReduction::GramSchmidtInformer * gsi)
            : d_af(af), d_gsi(gsi)
        {
        }
        
        bool operator() (IntTypeContext & ic, linalg::math_matrix<typename IntTypeContext::Integer> & A, int k,
                         arithmetic::RealContext & rc, arithmetic::RandomNumberGenerator & rng, arithmetic::Real & T)
        {
            MatrixConversionConst2<IntTypeContext> m(ic, A);
            return d_af(m.matrix(), k, rc, rng, T, d_gsi);
        }
    };
    
    template<class IntTypeContext>
    class BKZAnnealWrapper
    {
    private:
        LatticeReduction::BKZ_AnnealFunction d_af;
        LatticeReduction::GramSchmidtInformer * d_gsi;
        
    public:
        BKZAnnealWrapper(LatticeReduction::BKZ_AnnealFunction af, LatticeReduction::GramSchmidtInformer * gsi)
            : d_af(af), d_gsi(gsi)
        {
        }
        
        bool operator() (IntTypeContext & ic, linalg::math_matrix<typename IntTypeContext::Integer> & A, int k, int windowsize,
                         linalg::math_rowvector<typename IntTypeContext::Integer> & lincomb,
                         arithmetic::RealContext & rc, arithmetic::RandomNumberGenerator & rng, arithmetic::Real & T)
        {
            MatrixConversionConst2<IntTypeContext> m(ic, A);
            VectorConversion2<IntTypeContext> v(ic, lincomb, true, true);
            return d_af(m.matrix(), k, windowsize, v.vector(), rc, rng, T, d_gsi);
        }
    };
    
    bool DefaultAnnealCallbackFunction(const linalg::math_matrix<arithmetic::Integer> & A, arithmetic::RealContext & arc, arithmetic::Real & T);
    
    template<class IntTypeContext>
    bool DefaultLLLAnnealFunction(IntTypeContext & ic, linalg::math_matrix<typename IntTypeContext::Integer> & A, int k,
                                  arithmetic::RealContext & rc, arithmetic::RandomNumberGenerator & rng, arithmetic::Real & T);
    
    template<class IntTypeContext>
    bool DefaultBKZAnnealFunction(IntTypeContext & ic, linalg::math_matrix<typename IntTypeContext::Integer> & A, int k, int windowsize,
                                  linalg::math_rowvector<typename IntTypeContext::Integer> & lincomb,
                                  arithmetic::RealContext & rc, arithmetic::RandomNumberGenerator & rng, arithmetic::Real & T);
    
    void DefaultVerboseFunction(LatticeReduction::VerboseLevel, const std::string &);
    
    template<class IntTypeContext>
    void setUnit(linalg::math_matrix<typename IntTypeContext::Integer> & T);
// Sets T to the unit matrix. Assumes it is a square matrix.
    
    void addStatistics(LatticeReduction::Statistics &, const LatticeReduction::Statistics &);
    
////////////////////////////////////////////////////////////////////////
/// TRANSFORMATION NOTIFIERS ///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
    
    template<class IntTypeContext>
    class TransformData : public TransformNotifier<IntTypeContext>
    {
    private:
        linalg::math_matrix<typename IntTypeContext::Integer> * d_T;
        bool d_inv;
        
    public:
        TransformData(linalg::math_matrix<typename IntTypeContext::Integer> * T = NULL, bool inv = false);
        virtual ~TransformData() { }
        void setMatrix(linalg::math_matrix<typename IntTypeContext::Integer> * T, bool inv);
        
        virtual void swap(unsigned i, unsigned j);
        // Swap rows i and j.
        virtual void add(unsigned i, const typename IntTypeContext::Integer & m, unsigned j);
        // Add m-times the j-th row to the i-th row
        virtual void flip(unsigned i);
        // Flip the signs of the i-th row
        virtual void trans(unsigned i, unsigned j,
                           const typename IntTypeContext::Integer & B00, const typename IntTypeContext::Integer & B01,
                           const typename IntTypeContext::Integer & B10, const typename IntTypeContext::Integer & B11);
        // [ B00 B01 ]   [ ..... < row i > ..... ]
        // [ B10 B11 ] * [ ..... < row j > ..... ]
        virtual void insertVectorLC(unsigned ofs, const linalg::math_rowvector<typename IntTypeContext::Integer> & result);
        // Inserts a new vector at index ofs, which will equal \sum_{i=0}^{result.size()-1} result[i] * b[ofs+i]
        virtual void insertVector(unsigned ofs, const linalg::math_rowvector<typename IntTypeContext::Integer> & result);
        // Inserts a new vector <result> at index ofs
        virtual bool canInsertVector() { return d_inv; }
        virtual void removeZeroVector(unsigned ofs);
        // Removes the zero vector at position ofs
        virtual void compactify();
    };
    
    template<class RealTypeContext, class IntTypeContext>
    class Dumper : public TransformNotifier<IntTypeContext>
    {
    private:
        Verbose & d_verbose;
        Lattice<RealTypeContext, IntTypeContext> & d_lattice;
        
    public:
        Dumper(Verbose & v, Lattice<RealTypeContext, IntTypeContext> & lattice)
            : d_verbose(v), d_lattice(lattice)
        {
            d_lattice.addNotifier(this);
        }
        
        virtual ~Dumper()
        {
            d_lattice.removeNotifier(this);
        }
        
        virtual void swap(unsigned i, unsigned j)
        {
            d_verbose(LatticeReduction::VL_Information) << "swap(r" << i << ", r" << j << ")";
        }
        
        virtual void add(unsigned i, const typename IntTypeContext::Integer & m, unsigned j)
        {
            d_verbose(LatticeReduction::VL_Information) << "add(r" << i << " += " << m << " * r" << j << ")";
        }
        
        virtual void flip(unsigned i)
        {
            d_verbose(LatticeReduction::VL_Information) << "flip(r" << i << ")";
        }
        
        virtual void trans(unsigned i, unsigned j,
                           const typename IntTypeContext::Integer & B00, const typename IntTypeContext::Integer & B01,
                           const typename IntTypeContext::Integer & B10, const typename IntTypeContext::Integer & B11)
        {
            d_verbose(LatticeReduction::VL_Information) << "add([r" << i << ", r" << j << "] = [" << B00 << " " << B01 << " ; " << B10 << " " << B11 << "] * [r" << i << ", r" << j << "])";
        }
        
        virtual void insertVectorLC(unsigned ofs, const linalg::math_rowvector<typename IntTypeContext::Integer> & result)
        {
            d_verbose(LatticeReduction::VL_Information) << "insert(r" << ofs << ", " << result << ")";
        }
        
        virtual void insertVector(unsigned ofs, const linalg::math_rowvector<typename IntTypeContext::Integer> & result)
        {
            d_verbose(LatticeReduction::VL_Information) << "insert_rowvector(r" << ofs << ", " << result << ")";
        }
        
        virtual void removeZeroVector(unsigned ofs)
        {
            d_verbose(LatticeReduction::VL_Information) << "remove(r" << ofs << ")";
        }
        
        virtual void compactify()
        {
        }
    };
    
////////////////////////////////////////////////////////////////////////
/// LATTICE REDUCTION INTERFACE ////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
    
    class Transform
    {
    private:
        LatticeReduction::Transform d_mode;
        STD_AUTO_PTR<linalg::math_matrix<arithmetic::Integer> > d_transform;
        
        class TransformImpl;
        STD_AUTO_PTR<TransformImpl> d_impl;
        friend class TransformImpl;
        
        Transform(const Transform &);
        Transform & operator = (const Transform &);
        
    public:
        Transform();
        ~Transform();
        void reset(const linalg::math_matrix<arithmetic::Integer> & A);
        void enable(LatticeReduction::Transform mode, const linalg::math_matrix<arithmetic::Integer> & A);
        void disable();
        bool isEnabled() const;
        LatticeReduction::Transform mode() const;
        const linalg::math_matrix<arithmetic::Integer> * transform() const;
        
        template<class IntTypeContext>
        TransformNotifier<IntTypeContext> * createTransformObject(IntTypeContext &);
        void releaseTransformObject();
    };
    
    class LRIInterface
    {
    protected:
        linalg::math_matrix<arithmetic::Integer> & d_lattice;
        mutable LatticeReduction::Statistics d_stats;
        
    public:
        LRIInterface(linalg::math_matrix<arithmetic::Integer> & lattice)
            : d_lattice(lattice)
        {
        }
        
        virtual ~LRIInterface() { }
        
        virtual void setupVerbose(LatticeReduction::VerboseOutputLevel, LatticeReduction::VerboseFunction) = 0;
        virtual const LatticeReduction::GramSchmidtInformer * getInformer() const = 0;
        virtual void ensureMinimumPrecision(unsigned long) = 0;
        virtual void forceGSRebuild(bool) = 0;
        virtual double getGSCoefficientD(unsigned i, unsigned j) const = 0;
        virtual long double getGSCoefficientLD(unsigned i, unsigned j) const = 0;
        virtual arithmetic::Real getGSCoefficientR(unsigned i, unsigned j, const arithmetic::RealContext & rc) const = 0;
        virtual double getGSSqNormD(unsigned i) const = 0;
        virtual long double getGSSqNormLD(unsigned i) const = 0;
        virtual arithmetic::Real getGSSqNormR(unsigned i, const arithmetic::RealContext & rc) const = 0;
        virtual void modFlip(Transform & transform, unsigned i) = 0;
        virtual void modSwap(Transform & transform, unsigned i, unsigned j) = 0;
        virtual void modAdd(Transform & transform, unsigned, unsigned, const arithmetic::Integer &) = 0;
        
        virtual void sortProjected(unsigned & begin, unsigned & end, Transform & transform) = 0;
        virtual void sizereduction(unsigned & begin, unsigned & end, Transform & transform) = 0;
        
        virtual void lll(unsigned & begin, unsigned & end, Transform & transform, double alpha, LatticeReduction::LLLMode mode,
                         LatticeReduction::CallbackFunction cf, LatticeReduction::CallbackFunction_LI cf2, double cf_int,
                         LatticeReduction::MinCallbackFunction mcf, LatticeReduction::MinCallbackFunction_LI mcf2,
                         MaxBitsCallbackFunction mbcf, bool anneal, LatticeReduction::AnnealCallbackFunction acf,
                         LatticeReduction::LLL_AnnealFunction af, LatticeReduction::DIMethod di,
                         LatticeReduction::DIMode di_mode, LatticeReduction::DIChoice di_choice, unsigned di_bs) = 0;
        virtual void bkz(unsigned & begin, unsigned & end, Transform & transform, double alpha, unsigned blocksize, LatticeReduction::BKZMode mode,
                         LatticeReduction::CallbackFunction cf, LatticeReduction::CallbackFunction_LI cf2, double cf_int,
                         LatticeReduction::MinCallbackFunction mcf, LatticeReduction::MinCallbackFunction_LI mcf2,
                         MaxBitsCallbackFunction mbcf, LatticeReduction::EnumCallbackFunction ecf,
                         LatticeReduction::EnumCallbackFunction_LI ecf2, bool anneal, LatticeReduction::AnnealCallbackFunction acf,
                         LatticeReduction::BKZ_AnnealFunction af, LatticeReduction::DIMethod di,
                         LatticeReduction::DIMode di_mode, LatticeReduction::DIChoice di_choice, unsigned di_bs) = 0;
        virtual void hkz(unsigned & begin, unsigned & end, Transform & transform, bool dual,
                         LatticeReduction::CallbackFunction cf, LatticeReduction::CallbackFunction_LI cf2,
                         double cf_int, LatticeReduction::MinCallbackFunction mcf, LatticeReduction::MinCallbackFunction_LI mcf2,
                         MaxBitsCallbackFunction mbcf, LatticeReduction::EnumCallbackFunction ecf,
                         LatticeReduction::EnumCallbackFunction_LI ecf2) = 0;
        virtual void svp(unsigned & begin, unsigned & end, Transform & transform, bool make_basis, bool extreme, bool dual,
                         LatticeReduction::CallbackFunction cf, LatticeReduction::CallbackFunction_LI cf2, double cf_int,
                         LatticeReduction::MinCallbackFunction mcf, LatticeReduction::MinCallbackFunction_LI mcf2,
                         MaxBitsCallbackFunction mbcf, LatticeReduction::EnumCallbackFunction ecf,
                         LatticeReduction::EnumCallbackFunction_LI ecf2) = 0;
        virtual bool isSizeReduced(unsigned begin, unsigned end) const = 0;
        virtual bool isLLLBasis(unsigned begin, unsigned end, double alpha, LatticeReduction::LLLMode mode,
                                LatticeReduction::DIMethod di, LatticeReduction::DIChoice di_choice, unsigned di_bs) const = 0;
        virtual bool isBKZBasis(unsigned begin, unsigned end, double alpha, unsigned blocksize, LatticeReduction::BKZMode mode,
                                LatticeReduction::DIMethod di, LatticeReduction::DIChoice di_choice, unsigned di_bs) const = 0;
        virtual bool isHKZBasis(unsigned begin, unsigned end, bool dual) const = 0;
        virtual bool isSVPBasis(unsigned begin, unsigned end, bool dual) const = 0;
        
        const LatticeReduction::Statistics & getStatistics() const
        {
            return d_stats;
        }
        
        void resetStatistics()
        {
            d_stats.reset();
        }
    };
}

#endif

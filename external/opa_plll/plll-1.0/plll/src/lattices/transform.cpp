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

#ifndef PLLL_INCLUDE_GUARD__TRANSFORM_IMPL_CPP
#define PLLL_INCLUDE_GUARD__TRANSFORM_IMPL_CPP

#include <plll/matrix.hpp>
#include <plll/arithmetic.hpp>
#include "lll2-internal.hpp"
#include <list>

namespace plll
{
//////////////////////////////////////////////////////////////////////////////////////////////////////////
// Transformation Data Collectors

    template<class IntTypeContext>
    class TransformDataChangeNotifyer : public TransformNotifier<IntTypeContext>
    {
    private:
        bool d_changed;

    public:
        inline TransformDataChangeNotifyer()
            : d_changed(false)
        {
        }
    
        inline void resetChangeNotifyer()
        {
            d_changed = false;
        }
    
        inline bool wasChanged() const
        {
            return d_changed;
        }
    
        virtual void swap(unsigned i, unsigned j)
        // Swap rows i and j.
        {
            d_changed = true;
        }
    
        virtual void add(unsigned i, const typename IntTypeContext::Integer & m, unsigned j)
        // Add m-times the j-th row to the i-th row
        {
            d_changed = true;
        }
    
        virtual void flip(unsigned i)
        // Flip the signs of the i-th row
        {
            d_changed = true;
        }
    
        virtual void trans(unsigned i, unsigned j,
                           const typename IntTypeContext::Integer & B00, const typename IntTypeContext::Integer & B01,
                           const typename IntTypeContext::Integer & B10, const typename IntTypeContext::Integer & B11)
        // [ B00 B01 ]   [ ..... < row i > ..... ]
        // [ B10 B11 ] * [ ..... < row j > ..... ]
        {
            d_changed = true;
        }
    
        virtual void insertVectorLC(unsigned ofs, const linalg::math_rowvector<typename IntTypeContext::Integer> & result)
        // Inserts a new vector at index ofs, which will equal \sum_{i=0}^{result.size()-1} result[i] * b[ofs+i]
        {
            d_changed = true;
        }
    
        virtual void insertVector(unsigned ofs, const linalg::math_rowvector<typename IntTypeContext::Integer> & result)
        // Inserts a new vector <result> at index ofs
        {
            d_changed = true;
        }
    
        virtual bool canInsertVector()
        {
            return true;
        }
    
        virtual void removeZeroVector(unsigned ofs)
        // Removes the zero vector at position ofs
        {
            d_changed = true;
        }
    
        virtual void compactify()
        {
        }
    };

    template<class IntTypeContext>
    class TransformDataChangeLogger : public TransformNotifier<IntTypeContext>
    // Keep track of which vectors were changed.
    {
    private:
        unsigned d_changed_size;
        linalg::base_rowvector<bool> d_changed;
        unsigned d_min_change, d_max_change;
        bool d_something_changed;
    
    public:
        TransformDataChangeLogger(unsigned dim)
            : d_min_change(0), d_max_change(0), d_something_changed(false)
        {
            d_changed_size = dim;
            d_changed.resize(dim);
            for (unsigned i = 0; i < d_changed.size(); ++i)
                d_changed[i] = false;
        }
    
        // Keep track of exact changes
    
        void resetChangeLogger()
        {
            d_something_changed = false;
            for (unsigned i = 0; i < d_changed_size; ++i)
                d_changed[i] = false;
        }
    
        inline bool wasSomethingChanged() const
        {
            return d_something_changed;
        }
    
        inline unsigned getMinChanged() const
        {
            return d_min_change;
        }
    
        inline unsigned getMaxChanged() const
        {
            return d_max_change;
        }
    
        inline bool wasChanged(unsigned idx) const
        {
            return d_changed[idx];
        }
    
        // Usual transformation interface
    
        virtual void swap(unsigned i, unsigned j)
        // Swap rows i and j.
        {
            if (i > j)
                std::swap(i, j);
            d_changed[i] = true;
            d_changed[j] = true;
            if (d_something_changed)
            {
                if (d_min_change > i)
                    d_min_change = i;
                if (d_max_change < j)
                    d_max_change = j;
            }
            else
            {
                d_something_changed = true;
                d_min_change = i;
                d_max_change = j;
            }
        }
    
        virtual void add(unsigned i, const typename IntTypeContext::Integer & m, unsigned j)
        // Add m-times the j-th row to the i-th row
        {
            d_changed[i] = true;
            if (d_something_changed)
            {
                if (d_min_change > i)
                    d_min_change = i;
                if (d_max_change < i)
                    d_max_change = i;
            }
            else
            {
                d_something_changed = true;
                d_min_change = i;
                d_max_change = i;
            }
        }
    
        virtual void flip(unsigned i)
        // Flip the signs of the i-th row
        {
            d_changed[i] = true;
            if (d_something_changed)
            {
                if (d_min_change > i)
                    d_min_change = i;
                if (d_max_change < i)
                    d_max_change = i;
            }
            else
            {
                d_something_changed = true;
                d_min_change = i;
                d_max_change = i;
            }
        }
    
        virtual void trans(unsigned i, unsigned j,
                           const typename IntTypeContext::Integer & B00, const typename IntTypeContext::Integer & B01,
                           const typename IntTypeContext::Integer & B10, const typename IntTypeContext::Integer & B11)
        // [ B00 B01 ]   [ ..... < row i > ..... ]
        // [ B10 B11 ] * [ ..... < row j > ..... ]
        {
            if (i > j)
                std::swap(i, j);
            d_changed[i] = true;
            d_changed[j] = true;
            if (d_something_changed)
            {
                if (d_min_change > i)
                    d_min_change = i;
                if (d_max_change < j)
                    d_max_change = j;
            }
            else
            {
                d_something_changed = true;
                d_min_change = i;
                d_max_change = j;
            }
        }
    
        virtual void insertVectorLC(unsigned ofs, const linalg::math_rowvector<typename IntTypeContext::Integer> & result)
        // Inserts a new vector at index ofs, which will equal \sum_{i=0}^{result.size()-1} result[i] * b[ofs+i]
        {
            if (d_changed_size == d_changed.size())
                d_changed.resize(d_changed.size() + 1);
            for (unsigned i = d_changed_size; i > ofs; --i)
                d_changed[i] = d_changed[i - 1];
            d_changed[ofs] = true;
            ++d_changed_size;
            if (d_something_changed)
            {
                if (d_min_change > ofs)
                    d_min_change = ofs;
                if (d_max_change < ofs)
                    d_max_change = ofs;
                else
                    ++d_max_change;
            }
            else
            {
                d_something_changed = true;
                d_min_change = ofs;
                d_max_change = ofs;
            }
        }
    
        virtual void insertVector(unsigned ofs, const linalg::math_rowvector<typename IntTypeContext::Integer> & result)
        // Inserts a new vector <result> at index ofs
        {
            if (d_changed_size == d_changed.size())
                d_changed.resize(d_changed.size() + 1);
            for (unsigned i = d_changed_size; i > ofs; --i)
                d_changed[i] = d_changed[i - 1];
            d_changed[ofs] = true;
            ++d_changed_size;
            if (d_something_changed)
            {
                if (d_min_change > ofs)
                    d_min_change = ofs;
                if (d_max_change < ofs)
                    d_max_change = ofs;
                else
                    ++d_max_change;
            }
            else
            {
                d_something_changed = true;
                d_min_change = ofs;
                d_max_change = ofs;
            }
        }
    
        virtual bool canInsertVector()
        {
            return true;
        }
    
        virtual void removeZeroVector(unsigned ofs)
        // Removes the zero vector at position ofs
        {
            for (unsigned i = ofs + 1; i < d_changed_size; ++i)
                d_changed[i - 1] = d_changed[i];
            --d_changed_size;
            if (d_min_change > ofs)
                --d_min_change;
            if (d_max_change > ofs)
                --d_max_change;
        }
    
        virtual void compactify()
        {
        }
    };

    template<class IntTypeContext>
    class LLLBridgeNotifier : public TransformNotifier<IntTypeContext>
    // Tracks whether a swap between the given index i and i+1 appeared
    {
    private:
        unsigned d_index;
        bool d_bridge_appeared;
    
    public:
        inline LLLBridgeNotifier(unsigned idx)
            : d_index(idx), d_bridge_appeared(false)
        {
        }
    
        inline bool bridgeSwapAppeared() const
        {
            return d_bridge_appeared;
        }
    
        virtual ~LLLBridgeNotifier()
        {
        }
    
        virtual void swap(unsigned i, unsigned j)
        {
            if (std::min(i, j) == d_index)
                d_bridge_appeared = true;
        }
    
        virtual void add(unsigned, const typename IntTypeContext::Integer &, unsigned)
        {
        }
    
        virtual void flip(unsigned)
        {
        }
    
        virtual void trans(unsigned, unsigned,
                           const typename IntTypeContext::Integer &, const typename IntTypeContext::Integer &,
                           const typename IntTypeContext::Integer &, const typename IntTypeContext::Integer &)
        {
        }
    
        virtual void insertVectorLC(unsigned, const linalg::math_rowvector<typename IntTypeContext::Integer> &)
        {
        }
    
        virtual void insertVector(unsigned, const linalg::math_rowvector<typename IntTypeContext::Integer> &)
        {
        }
    
        virtual bool canInsertVector()
        {
            return true;
        }
    
        virtual void removeZeroVector(unsigned)
        {
        }
    
        virtual void compactify()
        {
        }
    };

    template<class IntTypeContext>
    class SimpleTransformationCollector : public TransformNotifier<IntTypeContext>
    {
    private:
        unsigned d_dimension;
        mutable linalg::math_matrix<typename IntTypeContext::Integer> d_T;
    
        inline void compactifyImpl() const
        {
            d_T.resize(d_dimension, d_T.cols());
        }
    
    public:
        inline SimpleTransformationCollector(unsigned dimension)
            : d_dimension(dimension)
        {
            d_T.resize(dimension, dimension);
            setUnit<IntTypeContext>(d_T);
        }
    
        virtual ~SimpleTransformationCollector()
        {
        }
    
        inline const linalg::math_matrix<typename IntTypeContext::Integer> & transform() const
        {
            if (d_dimension < d_T.rows())
                compactifyImpl();
            return d_T;
        }
    
        virtual void swap(unsigned i, unsigned j)
        // Swap rows i and j.
        {
            linalg::swap(d_T.row(i), d_T.row(j));
        }
    
        virtual void add(unsigned i, const typename IntTypeContext::Integer & m, unsigned j)
        // Add m-times the j-th row to the i-th row
        {
            if (isPMOne(m))
            {
                if (sign(m) > 0)
                    d_T.row(i) += d_T.row(j);
                else
                    d_T.row(i) -= d_T.row(j);
            }
            else
                plll::linalg::addmul(d_T.row(i), d_T.row(j), m); // A.row(i) += A.row(j) * m;
        }
    
        virtual void flip(unsigned i)
        // Flip the signs of the i-th row
        {
            plll::linalg::neg(d_T.row(i), d_T.row(i));
        }
    
        virtual void trans(unsigned i, unsigned j,
                           const typename IntTypeContext::Integer & B00, const typename IntTypeContext::Integer & B01,
                           const typename IntTypeContext::Integer & B10, const typename IntTypeContext::Integer & B11)
        // [ B00 B01 ]   [ ..... < row i > ..... ]
        // [ B10 B11 ] * [ ..... < row j > ..... ]
        {
            typename IntTypeContext::Integer t0, t1;
            for (unsigned k = 0; k < d_T.cols(); ++k)
            {
                t0 = B00 * d_T(i, k);
                t0 += B01 * d_T(j, k);
                t1 = B10 * d_T(i, k);
                t1 += B11 * d_T(j, k);
                d_T(i, k) = t0;
                d_T(j, k) = t1;
            }
        }
    
        virtual void insertVectorLC(unsigned ofs, const linalg::math_rowvector<typename IntTypeContext::Integer> & result)
        // Inserts a new vector at index ofs, which will equal \sum_{i=0}^{result.size()-1} result[i] * b[ofs+i]
        {
            bool add = d_dimension == d_T.rows();
            if (add)
                d_T.resize(d_T.rows() + 1, d_T.cols());
            ++d_dimension;
            for (unsigned i = d_T.rows() - 1; i > ofs; --i)
                linalg::swap(d_T.row(i), d_T.row(i - 1));
            for (unsigned i = 0; i < result.size(); ++i)
                if (!isZero(result[i]))
                    plll::linalg::addmul(d_T.row(ofs), d_T.row(ofs + i + 1), result[i]); // v += result[i] * A.row(ofs + i);
        }
    
        virtual void insertVector(unsigned ofs, const linalg::math_rowvector<typename IntTypeContext::Integer> & result)
        // Inserts a new vector <result> at index ofs
        {
            // This should not be called!
            assert(false);
        }
    
        virtual bool canInsertVector()
        {
            return false;
        }
    
        virtual void removeZeroVector(unsigned ofs)
        // Removes the zero vector at position ofs
        {
            --d_dimension;
        }
    
        virtual void compactify()
        {
            if (d_dimension < d_T.rows())
                compactifyImpl();
        }
    };

    template<class RealTypeContext, class IntTypeContext>
    class DualTransformationApplier : public TransformNotifier<IntTypeContext>
    {
    private:
        std::pair<unsigned, unsigned> d_range; // range for d_lattice
        unsigned d_dimension;
        Lattice<RealTypeContext, IntTypeContext> & d_lattice;
        Lattice<RealTypeContext, IntTypeContext> & d_dual;
    
        /*
          If R_n is the "reverse identity matrix", satisfying R_n^2 = I_n and R_n \neq I_n (for n > 1),
          then if B = R * Q is the RQ-decomposition of the basis B, then
      
          R^{-s} Q^{-s},    where A^{-s} := R_n (A^T)^{-1} R_n,
      
          is the canonical reversed basis of the dual of the lattice generated by B.
      
          If we transform B^{-s} by some transform T, i.e. if we have T B^{-s}, then
      
          T B^{-s} = ( T^{-s} B )^{-s}.
      
          Therefore, we apply T^{-s} to B to achieve the transformation on the non-dualed lattice.
        */
        
        static inline bool notequal(const typename RealTypeContext::Real & a, const typename RealTypeContext::Real & b,
                                    const RealTypeContext & rc, helper::BoolToType<false>)
        {
            return a != b;
        }
        
        static inline bool notequal(const typename RealTypeContext::Real & a, const typename RealTypeContext::Real & b,
                                    const RealTypeContext & rc, helper::BoolToType<true>)
        {
            return compareAbsValues(a - b, rc.getEpsilon() << 10) > 0;
        }
        
        static inline bool notequal(const typename RealTypeContext::Real & a, const typename RealTypeContext::Real & b,
                                    const RealTypeContext & rc)
        {
            return notequal(a, b, rc, helper::BoolToType<RealTypeContext::has_constants>());
        }
        
        void verify()
        // Debug function
        {
            d_lattice.range().addRange(d_range);
            d_lattice.update(d_range.second + 1);
            STD_AUTO_PTR<GSInterface<RealTypeContext, IntTypeContext> > gs = getDualGSI(d_lattice.rc(), d_lattice);
            Lattice<RealTypeContext, IntTypeContext> dl(gs.get(), d_lattice.rc(), d_lattice.getStats(), 0, gs->getDimension() - 1);
            bool ok = true;
            for (unsigned i = 0; i < d_dimension; ++i)
            {
                for (unsigned j = 0; j < i; ++j)
                    if (notequal(d_dual.getCoeff(i, j), dl.getCoeff(i, j), d_lattice.rc()))
                    {
                        std::cout << i << " " << j << ": " << d_dual.getCoeff(i, j) << " vs " << dl.getCoeff(i, j) << "\n";
                        ok = false;
                    }
                if (notequal(d_dual.getNormSq(i), dl.getNormSq(i), d_lattice.rc()))
                {
                    std::cout << i << ": " << d_dual.getNormSq(i) << " vs " << dl.getNormSq(i) << "\n";
                    ok = false;
                }
            }
            if (!ok)
            {
                std::cout << "Should be:\n" << d_dual << "\n";
                std::cout << "Is:\n" << dl << "\n";
            }
            assert(ok);
            d_lattice.range().popRange();
        }
    
        inline unsigned c(unsigned i)
        {
            return d_range.first + d_dimension - 1 - i;
        }
    
    public:
        inline DualTransformationApplier(Lattice<RealTypeContext, IntTypeContext> & lattice,
                                         Lattice<RealTypeContext, IntTypeContext> & dual)
            // Assumes that <dual> is the dual of the projected sublattice corresponding to the ranbge
            // <lattice.range()>.
            : d_range(lattice.range().begin(), lattice.range().end()),
              d_dimension(lattice.range().dimension()), d_lattice(lattice), d_dual(dual)
        {
            d_dual.addNotifier(this);
//        verify();
        }
        
        virtual ~DualTransformationApplier() PLLL_INTERNAL_NOTHROW_POSTFIX_INLINE
        {
            try
            {
                d_dual.removeNotifier(this);
            }
            catch (...)
            {
            }
        }
        
        virtual void swap(unsigned i, unsigned j)
        /*
          Note that transpose and inverse of a swapping matrix are the swapping matrix itself.
      
          [ 1     i       j       ]^{-s}   [       i       j     1 ]         [ 1   n-j-1   n-i-1     ]
          [   \   |       |       ]        [       |       |   /   ]         [   \   |       |       ]
          [     1 |       |       ]        [       |       | 1     ]         [     1 |       |       ]
          [       0.......1------i]        [       1.......0- n-j-1]         [       0.......1- n-j-1]
          [       . 1     .       ]        [       .     1 .       ]         [       . 1     .       ]
          [       .   \   .       ]      = [       .   /   .       ] * R_n = [       .   \   .       ]
          [       .     1 .       ]        [       . 1     .       ]         [       .     1 .       ]
          [       1.......0------j]        [       0.......1- n-i-1]         [       1.......0- n-i-1]
          [                 1     ]        [     1                 ]         [                 1     ]
          [                   \   ]        [   /                   ]         [                   \   ]
          [                     1 ]        [ 1                     ]         [                     1 ]
        */
        {
            d_lattice.swap(c(j), c(i));
//        verify();
        }
    
        virtual void add(unsigned i, const typename IntTypeContext::Integer & m, unsigned j)
        /*
          [ 1     i       j       ]^{-s}         [ 1     i               ]         [ 1   n-j-1   n-i-1     ]
          [   \   |       |       ]              [   \   |               ]         [   \   |       |       ]
          [     1 |       |       ]              [     1 |               ]         [     1 |       |       ]
          [       1.......m------i]              [       1               ]         [       1......-m  n-j-1]
          [         1     .       ]              [       | 1             ]         [         1     .       ]
          [           \   .       ]      = R_n * [       |   \           ] * R_n = [           \   .       ]
          [             1 .       ]              [       |     1         ]         [             1 .       ]
          [               1------j]              [      -m.......1------j]         [               1- n-i-1]
          [                 1     ]              [                 1     ]         [                 1     ]
          [                   \   ]              [                   \   ]         [                   \   ]
          [                     1 ]              [                     1 ]         [                     1 ]
        */
        {
            d_lattice.add(c(j), -m, c(i));
//        verify();
        }
    
        virtual void flip(unsigned i)
        {
            d_lattice.flip(c(i));
//        verify();
        }
    
        virtual void trans(unsigned i, unsigned j,
                           const typename IntTypeContext::Integer & B00, const typename IntTypeContext::Integer & B01,
                           const typename IntTypeContext::Integer & B10, const typename IntTypeContext::Integer & B11)
        /*
          [ 1     i       j       ]^{-s}         [ 1     i       j       ]         [ 1   n-j-1   n-i-1     ]
          [   \   |       |       ]              [   \   |       |       ]         [   \   |       |       ]
          [     1 |       |       ]              [     1 |       |       ]         [     1 |       |       ]
          [      B00.....B01-----i]              [      B11....-B10-----i]         [      B00....-B01 n-j-1]
          [       . 1     .       ]              [       . 1     .       ]         [       . 1     .       ]
          [       .   \   .       ]      = R_n * [       .   \   .       ] * R_n = [       .   \   .       ]
          [       .     1 .       ]              [       .     1 .       ]         [       .     1 .       ]
          [      B10.....B11-----j]              [     -B01.....B00-----j]         [     -B10.....B11 n-i-1]
          [                 1     ]              [                 1     ]         [                 1     ]
          [                   \   ]              [                   \   ]         [                   \   ]
          [                     1 ]              [                     1 ]         [                     1 ]
        */
        {
            d_lattice.trans(c(j), c(i), B00, -B01, -B10, B11);
//        verify();
        }
    
        // Vector insertions and removals are not supported at the moment!
    
        virtual void insertVectorLC(unsigned ofs, const linalg::math_rowvector<typename IntTypeContext::Integer> & result)
        {
            assert(false);
        }
    
        virtual void insertVector(unsigned ofs, const linalg::math_rowvector<typename IntTypeContext::Integer> & result)
        {
            assert(false);
        }
    
        virtual bool canInsertVector()
        {
            return false;
        }
    
        virtual void removeZeroVector(unsigned ofs)
        {
            assert(false);
        }
    
        virtual void compactify()
        {
        }
    };

    template<class IntTypeContext>
    class TransformationStoreAndReplay : public TransformNotifier<IntTypeContext>
    {
    private:
        enum CommandType { CT_Swap, CT_Add, CT_Flip, CT_Trans, CT_Insert, CT_InsertLC, CT_Remove };
    
        struct Command
        {
            CommandType cmd;
            unsigned i, j;
            linalg::math_rowvector<typename IntTypeContext::Integer> args;
        
            inline Command(CommandType cmd_, unsigned i_)
                : cmd(cmd_), i(i_), j(0)
            {
            }
        
            inline Command(CommandType cmd_, unsigned i_, unsigned j_)
                : cmd(cmd_), i(i_), j(j_)
            {
            }
        
            inline Command(CommandType cmd_, unsigned i_, unsigned j_, const typename IntTypeContext::Integer & a)
                : cmd(cmd_), i(i_), j(j_)
            {
                args.resize(1);
                args[0] = a;
            }
        
            inline Command(CommandType cmd_, unsigned i_, unsigned j_,
                           const typename IntTypeContext::Integer & a, const typename IntTypeContext::Integer & b,
                           const typename IntTypeContext::Integer & c, const typename IntTypeContext::Integer & d)
                : cmd(cmd_), i(i_), j(j_)
            {
                args.resize(4);
                args[0] = a;
                args[1] = b;
                args[2] = c;
                args[3] = d;
            }
        
            inline Command(CommandType cmd_, unsigned i_, const linalg::math_rowvector<typename IntTypeContext::Integer> & args_)
                : cmd(cmd_), i(i_), j(0), args(args_)
            {
            }
        };
    
        std::list<Command> d_commands;
    
    public:
        virtual ~TransformationStoreAndReplay()
        {
        }
    
        virtual void swap(unsigned i, unsigned j)
        {
            d_commands.push_back(Command(CT_Swap, i, j));
        }
    
        virtual void add(unsigned i, const typename IntTypeContext::Integer & m, unsigned j)
        {
            d_commands.push_back(Command(CT_Add, i, j, m));
        }
    
        virtual void flip(unsigned i)
        {
            d_commands.push_back(Command(CT_Flip, i));
        }
    
        virtual void trans(unsigned i, unsigned j,
                           const typename IntTypeContext::Integer & B00, const typename IntTypeContext::Integer & B01,
                           const typename IntTypeContext::Integer & B10, const typename IntTypeContext::Integer & B11)
        {
            d_commands.push_back(Command(CT_Trans, i, j, B00, B01, B10, B11));
        }
    
        virtual void insertVectorLC(unsigned ofs, const linalg::math_rowvector<typename IntTypeContext::Integer> & result)
        {
            d_commands.push_back(Command(CT_InsertLC, ofs, result));
        }
    
        virtual void insertVector(unsigned ofs, const linalg::math_rowvector<typename IntTypeContext::Integer> & result)
        {
            d_commands.push_back(Command(CT_Insert, ofs, result));
        }
    
        virtual bool canInsertVector()
        {
            return true;
        }
    
        virtual void removeZeroVector(unsigned ofs)
        {
            d_commands.push_back(Command(CT_Remove, ofs));
        }
    
        virtual void compactify()
        {
        }
    
        template<class RealTypeContext>
        void replay(Lattice<RealTypeContext, IntTypeContext> & lattice) const
        {
            for (typename std::list<Command>::const_iterator i = d_commands.begin(); i != d_commands.end(); ++i)
                switch (i->cmd)
                {
                case CT_Swap:
                    lattice.swap(i->i, i->j);
                    break;
                case CT_Add:
                    lattice.add(i->i, i->args[0], i->j);
                    break;
                case CT_Flip:
                    lattice.flip(i->i);
                    break;
                case CT_Trans:
                    lattice.trans(i->i, i->j, i->args[0], i->args[1], i->args[2], i->args[3]);
                    break;
                case CT_Insert:
                    lattice.insertVector(i->i, i->args);
                    break;
                case CT_InsertLC:
                    lattice.insertVectorLC(i->i, i->args);
                    break;
                case CT_Remove:
                    lattice.removeZeroVector(i->i);
                    break;
                }
        }
    };

    template<class IntTypeContext>
    TransformData<IntTypeContext>::TransformData(linalg::math_matrix<typename IntTypeContext::Integer> * T, bool inv)
        : d_T(T), d_inv(inv)
    {
    }

    template<class IntTypeContext>
    void TransformData<IntTypeContext>::setMatrix(linalg::math_matrix<typename IntTypeContext::Integer> * T, bool inv)
    {
        d_T = T;
        d_inv = inv;
    }

    template<class IntTypeContext>
    void TransformData<IntTypeContext>::swap(unsigned i, unsigned j)
    // Swap rows i and j.
    {
        if (d_inv)
            linalg::swap(d_T->col(i), d_T->col(j));
        else
            linalg::swap(d_T->row(i), d_T->row(j));
    }

    template<class IntTypeContext>
    void TransformData<IntTypeContext>::add(unsigned i, const typename IntTypeContext::Integer & m, unsigned j)
    // Add m-times the j-th row to the i-th row
    {
        if (d_inv)
        {
            if (isPMOne(m))
            {
                if (sign(m) > 0)
                    d_T->col(j) -= d_T->col(i);
                else
                    d_T->col(j) += d_T->col(i);
            }
            else
                plll::linalg::addmul(d_T->col(j), d_T->col(i), -m);
        }
        else
        {
            if (isPMOne(m))
            {
                if (sign(m) > 0)
                    d_T->row(i) += d_T->row(j);
                else
                    d_T->row(i) -= d_T->row(j);
            }
            else
                plll::linalg::addmul(d_T->row(i), d_T->row(j), m); // A.row(i) += A.row(j) * m;
        }
    }

    template<class IntTypeContext>
    void TransformData<IntTypeContext>::flip(unsigned i)
    // Flip the signs of the i-th row
    {
        if (d_inv)
            plll::linalg::neg(d_T->col(i), d_T->col(i));
        else
            plll::linalg::neg(d_T->row(i), d_T->row(i));
    }

    template<class IntTypeContext>
    void TransformData<IntTypeContext>::trans(unsigned i, unsigned j,
                                              const typename IntTypeContext::Integer & B00, const typename IntTypeContext::Integer & B01,
                                              const typename IntTypeContext::Integer & B10, const typename IntTypeContext::Integer & B11)
    // [ B00 B01 ]   [ ..... < row i > ..... ]
    // [ B10 B11 ] * [ ..... < row j > ..... ]
    {
        typename IntTypeContext::Integer t0, t1;
        if (d_inv)
            // The inverse of [ B00 B01 ] is [  B11  -B01 ] since the determinant is assumed to be +1.
            //                [ B10 B11 ]    [ -B10   B00 ]
            // Hence, we have to compute
            //          [   :        :   ]
            //          [ col i    col j ] * [  B11  -B01 ]
            //          [   :        :   ]   [ -B10   B00 ]
            for (unsigned k = 0; k < d_T->rows(); ++k)
            {
                t0 = B11 * (*d_T)(k, i);
                t0 -= B10 * (*d_T)(k, j);
                t1 = B01 * (*d_T)(k, i);
                t1 -= B00 * (*d_T)(k, j);
                (*d_T)(k, i) = t0;
                (*d_T)(k, j) = -t1;
            }
        else
            for (unsigned k = 0; k < d_T->cols(); ++k)
            {
                t0 = B00 * (*d_T)(i, k);
                t0 += B01 * (*d_T)(j, k);
                t1 = B10 * (*d_T)(i, k);
                t1 += B11 * (*d_T)(j, k);
                (*d_T)(i, k) = t0;
                (*d_T)(j, k) = t1;
            }
    }

    template<class IntTypeContext>
    void TransformData<IntTypeContext>::insertVectorLC(unsigned ofs, const linalg::math_rowvector<typename IntTypeContext::Integer> & result)
    // Inserts a new vector at index ofs, which will equal \sum_{i=0}^{result.size()-1} result[i] * b[ofs+i]
    {
        if (d_inv)
        {
            d_T->resize(d_T->rows(), d_T->cols() + 1);
            for (unsigned i = d_T->cols() - 1; i > ofs; --i)
                linalg::swap(d_T->col(i), d_T->col(i - 1));
            for (unsigned i = 0; i < d_T->rows(); ++i)
                setZero((*d_T)(i, ofs));
        }
        else
        {
            linalg::math_rowvector<typename IntTypeContext::Integer> v(d_T->cols());
            for (unsigned i = 0; i < result.size(); ++i)
                if (!isZero(result[i]))
                    plll::linalg::addmul(v, d_T->row(ofs + i), result[i]); // v += result[i] * A.row(ofs + i);
            d_T->resize(d_T->rows() + 1, d_T->cols());
            for (unsigned i = d_T->rows() - 1; i > ofs; --i)
                linalg::swap(d_T->row(i), d_T->row(i - 1));
            d_T->row(ofs) = v;
        }
    }

    template<class IntTypeContext>
    void TransformData<IntTypeContext>::insertVector(unsigned ofs, const linalg::math_rowvector<typename IntTypeContext::Integer> & result)
    // Inserts a new vector <result> at index ofs
    {
        if (d_inv)
        {
            d_T->resize(d_T->rows(), d_T->cols() + 1);
            for (unsigned i = d_T->cols() - 1; i > ofs; --i)
                linalg::swap(d_T->col(i), d_T->col(i - 1));
            for (unsigned i = 0; i < d_T->rows(); ++i)
                setZero((*d_T)(i, ofs));
        }
        else
            assert(false); // This should never be called!
    }

    template<class IntTypeContext>
    void TransformData<IntTypeContext>::removeZeroVector(unsigned ofs)
    // Removes the zero vector at position ofs
    {
        if (d_inv)
        {
            for (unsigned i = ofs + 1; i < d_T->cols(); ++i)
                linalg::swap(d_T->col(i), d_T->col(i - 1));
            d_T->resize(d_T->rows(), d_T->cols() - 1);
        }
        else
        {
            for (unsigned i = ofs + 1; i < d_T->rows(); ++i)
                linalg::swap(d_T->row(i), d_T->row(i - 1));
            d_T->resize(d_T->rows() - 1, d_T->cols());
        }
    }

    template<class IntTypeContext>
    void TransformData<IntTypeContext>::compactify()
    {
    }

    template<class IntTypeContext>
    void setUnit(linalg::math_matrix<typename IntTypeContext::Integer> & A)
    // Sets T to the unit matrix. Assumes it is a square matrix.
    {
        assert(A.rows() == A.cols());
        unsigned n = A.rows();
        for (unsigned i = 0; i < n; ++i)
            for (unsigned j = 0; j < n; ++j)
                if (i == j)
                    setOne(A(i, j));
                else
                    setZero(A(i, j));
    }

    template<class IntTypeContext, bool is_cputype>
    class TransformDataMaxBitsCounterImpl;

    template<class IntTypeContext>
    class TransformDataMaxBitsCounterImpl<IntTypeContext, false> : private TransformNotifier<IntTypeContext>
    {
    private:
        linalg::math_matrix<typename IntTypeContext::Integer> * d_basis;
        unsigned d_max_bits_size;
        linalg::base_rowvector<unsigned> d_max_bits;
        mutable unsigned d_total_max_bits, d_max_index, d_max_since_changed;
        mutable bool d_recompute, d_changed;
    
        void updateMax() const
        {
            unsigned old_max = d_total_max_bits;
            d_total_max_bits = d_max_bits[0];
            d_max_index = 0;
            for (unsigned ii = 1; ii < d_max_bits_size; ++ii)
                if (d_max_bits[ii] > d_total_max_bits)
                {
                    d_total_max_bits = d_max_bits[ii];
                    d_max_index = ii;
                }
            if (old_max != d_total_max_bits)
            {
                d_max_since_changed = std::max(d_max_since_changed, d_total_max_bits);
                d_changed = true;
            }
        }
    
        void update(unsigned i)
        {
            d_max_bits[i] = 0;
            for (unsigned k = 0; k < d_basis->cols(); ++k)
            {
                long v = approxLog2((*d_basis)(i, k));
                if (d_max_bits[i] < v)
                    d_max_bits[i] = v;
            }
            if (d_max_bits[i] > d_total_max_bits)
            {
                d_total_max_bits = d_max_bits[i];
                d_max_index = i;
                d_max_since_changed = std::max(d_max_since_changed, d_total_max_bits);
                d_changed = true;
            }
            else if ((d_max_index == i) && (d_max_bits[i] < d_total_max_bits))
                d_recompute = true;
        }
    
        // Usual transformation interface
    
        virtual void swap(unsigned i, unsigned j)
        // Swap rows i and j.
        {
            std::swap(d_max_bits[i], d_max_bits[j]);
            if (d_max_index == i)
                d_max_index = j;
            else if (d_max_index == j)
                d_max_index = i;
        }
    
        virtual void add(unsigned i, const typename IntTypeContext::Integer &, unsigned)
        // Add m-times the j-th row to the i-th row
        {
            update(i);
        }
    
        virtual void flip(unsigned i)
        // Flip the signs of the i-th row
        {
        }
    
        virtual void trans(unsigned i, unsigned j,
                           const typename IntTypeContext::Integer &, const typename IntTypeContext::Integer &,
                           const typename IntTypeContext::Integer &, const typename IntTypeContext::Integer &)
        // [ B00 B01 ]   [ ..... < row i > ..... ]
        // [ B10 B11 ] * [ ..... < row j > ..... ]
        {
            update(i);
            update(j);
        }
    
        virtual void insertVectorLC(unsigned ofs, const linalg::math_rowvector<typename IntTypeContext::Integer> &)
        {
            if (d_max_bits_size == d_max_bits.size())
                d_max_bits.resize(d_max_bits.size() + 1);
            for (unsigned i = d_max_bits_size; i > ofs; --i)
                d_max_bits[i] = d_max_bits[i - 1];
            ++d_max_bits_size;
            if (ofs <= d_max_index)
                ++d_max_index;
            d_max_bits[ofs] = 0;
            update(ofs);
        }
    
        virtual void insertVector(unsigned ofs, const linalg::math_rowvector<typename IntTypeContext::Integer> &)
        {
            if (d_max_bits_size == d_max_bits.size())
                d_max_bits.resize(d_max_bits.size() + 1);
            for (unsigned i = d_max_bits_size; i > ofs; --i)
                d_max_bits[i] = d_max_bits[i - 1];
            ++d_max_bits_size;
            if (ofs <= d_max_index)
                ++d_max_index;
            d_max_bits[ofs] = 0;
            update(ofs);
        }
    
        virtual bool canInsertVector()
        {
            return true;
        }
    
        virtual void removeZeroVector(unsigned ofs)
        // Removes the zero vector at position ofs
        {
            for (unsigned i = ofs + 1; i < d_max_bits_size; ++i)
                d_max_bits[i - 1] = d_max_bits[i];
            --d_max_bits_size;
            if (d_max_index > ofs)
                --d_max_index;
        }
    
        virtual void compactify()
        {
        }
    
    public:
        inline TransformDataMaxBitsCounterImpl()
            : d_basis(NULL), d_max_bits_size(0), d_total_max_bits(0), d_max_index(0), d_max_since_changed(0), d_recompute(false), d_changed(false)
        {
        }

        template<class RealTypeContext>
        void addTo(Lattice<RealTypeContext, IntTypeContext> & lattice)
        {
            d_basis = lattice.getMatrixIndex(0).first;
            if (d_basis)
            {
                lattice.addNotifier(this);
                d_changed = true;
                d_max_since_changed = 0;
                d_max_bits_size = lattice.dimension();
                d_max_bits.resize(d_max_bits_size);
                for (unsigned i = 0; i < d_max_bits_size; ++i)
                    update(i);
            }
        }
    
        inline bool changed() const
        // Query whether maximum number of bits changed since last reset.
        {
            if (d_recompute)
            {
                updateMax();
                d_recompute = false;
            }
            return d_changed;
        }
    
        inline void resetChange()
        {
            d_changed = false;
            d_max_since_changed = 0;
        }
    
        inline unsigned maxBits() const
        // Returns maximal number of bits needed since last change.
        {
            if (d_recompute)
            {
                updateMax();
                d_recompute = false;
            }
            return std::max(d_max_since_changed, d_total_max_bits);
        }
    };

    template<class IntTypeContext>
    class TransformDataMaxBitsCounterImpl<IntTypeContext, true> : private TransformNotifier<IntTypeContext>
    {
    private:
        linalg::math_matrix<typename IntTypeContext::Integer> * d_basis;
        unsigned d_maxs_size;
        linalg::base_rowvector<typename IntTypeContext::Integer> d_maxs;
        mutable typename IntTypeContext::Integer d_total_max, d_max_since_changed;
        mutable unsigned d_max_index;
        mutable bool d_recompute, d_changed;
    
        void updateMax() const
        {
            typename IntTypeContext::Integer old_max = d_total_max;
            d_total_max = d_maxs[0];
            d_max_index = 0;
            for (unsigned ii = 1; ii < d_maxs_size; ++ii)
                if (d_maxs[ii] > d_total_max)
                {
                    d_total_max = d_maxs[ii];
                    d_max_index = ii;
                }
            if (old_max != d_total_max)
            {
                if (d_max_since_changed < d_total_max)
                    d_max_since_changed = d_total_max;
                d_changed = true;
            }
        }
    
        void update(unsigned i)
        {
            setZero(d_maxs[i]);
            for (unsigned k = 0; k < d_basis->cols(); ++k)
            {
                if (d_maxs[i] < (*d_basis)(i, k))
                    d_maxs[i] = (*d_basis)(i, k);
                else if (d_maxs[i] < -(*d_basis)(i, k))
                    d_maxs[i] = -(*d_basis)(i, k);
            }
            if (d_maxs[i] > d_total_max)
            {
                d_total_max = d_maxs[i];
                d_max_index = i;
                if (d_max_since_changed < d_total_max)
                    d_max_since_changed = d_total_max;
                d_changed = true;
            }
            else if ((d_max_index == i) && (d_maxs[i] < d_total_max))
                d_recompute = true;
        }
    
        // Usual transformation interface
    
        virtual void swap(unsigned i, unsigned j)
        // Swap rows i and j.
        {
            std::swap(d_maxs[i], d_maxs[j]);
            if (d_max_index == i)
                d_max_index = j;
            else if (d_max_index == j)
                d_max_index = i;
        }
    
        virtual void add(unsigned i, const typename IntTypeContext::Integer &, unsigned)
        // Add m-times the j-th row to the i-th row
        {
            update(i);
        }
    
        virtual void flip(unsigned i)
        // Flip the signs of the i-th row
        {
        }
    
        virtual void trans(unsigned i, unsigned j,
                           const typename IntTypeContext::Integer &, const typename IntTypeContext::Integer &,
                           const typename IntTypeContext::Integer &, const typename IntTypeContext::Integer &)
        // [ B00 B01 ]   [ ..... < row i > ..... ]
        // [ B10 B11 ] * [ ..... < row j > ..... ]
        {
            update(i);
            update(j);
        }
    
        virtual void insertVectorLC(unsigned ofs, const linalg::math_rowvector<typename IntTypeContext::Integer> &)
        {
            if (d_maxs_size == d_maxs.size())
                d_maxs.resize(d_maxs.size() + 1);
            for (unsigned i = d_maxs_size; i > ofs; --i)
                d_maxs[i] = d_maxs[i - 1];
            ++d_maxs_size;
            if (ofs <= d_max_index)
                ++d_max_index;
            setZero(d_maxs[ofs]);
            update(ofs);
        }
    
        virtual void insertVector(unsigned ofs, const linalg::math_rowvector<typename IntTypeContext::Integer> &)
        {
            if (d_maxs_size == d_maxs.size())
                d_maxs.resize(d_maxs.size() + 1);
            for (unsigned i = d_maxs_size; i > ofs; --i)
                d_maxs[i] = d_maxs[i - 1];
            ++d_maxs_size;
            if (ofs <= d_max_index)
                ++d_max_index;
            setZero(d_maxs[ofs]);
            update(ofs);
        }
    
        virtual bool canInsertVector()
        {
            return true;
        }
    
        virtual void removeZeroVector(unsigned ofs)
        // Removes the zero vector at position ofs
        {
            for (unsigned i = ofs + 1; i < d_maxs_size; ++i)
                d_maxs[i - 1] = d_maxs[i];
            --d_maxs_size;
            if (d_max_index > ofs)
                --d_max_index;
        }
    
        virtual void compactify()
        {
        }
    
    public:
        inline TransformDataMaxBitsCounterImpl()
            : d_basis(NULL), d_maxs_size(0), d_max_index(0), d_recompute(false), d_changed(false)
        {
        }
    
        template<class RealTypeContext>
        void addTo(Lattice<RealTypeContext, IntTypeContext> & lattice)
        {
            d_basis = lattice.getMatrixIndex(0).first;
            if (d_basis)
            {
                lattice.addNotifier(this);
                d_changed = true;
                setZero(d_total_max);
                setZero(d_max_since_changed);
                d_maxs_size = lattice.dimension();
                d_maxs.resize(d_maxs_size);
                for (unsigned i = 0; i < d_maxs_size; ++i)
                    update(i);
            }
        }
    
        inline bool changed() const
        // Query whether maximum number of bits changed since last reset.
        {
            if (d_recompute)
            {
                updateMax();
                d_recompute = false;
            }
            return d_changed;
        }
    
        inline void resetChange()
        {
            d_changed = false;
            setZero(d_max_since_changed);
        }
    
        inline unsigned maxBits() const
        // Returns maximal number of bits needed since last change.
        {
            if (d_recompute)
            {
                updateMax();
                d_recompute = false;
            }
            return approxLog2(std::max(d_max_since_changed, d_total_max));
        }
    };

    template<class IntTypeContext>
    class TransformDataMaxBitsCounter
    {
    private:
        TransformDataMaxBitsCounterImpl<IntTypeContext, IntTypeContext::is_cputype> d_impl;
    
    public:
        inline TransformDataMaxBitsCounter()
        {
        }
    
        template<class RealTypeContext>
        inline void addTo(Lattice<RealTypeContext, IntTypeContext> & lattice)
        {
            d_impl.addTo(lattice);
        }
    
        inline bool changed() const
        // Query whether maximum number of bits changed since last reset.
        {
            return d_impl.changed();
        }
    
        inline void resetChange()
        {
            d_impl.resetChange();
        }
    
        inline unsigned maxBits() const
        // Returns maximal number of bits needed since last change.
        {
            return d_impl.maxBits();
        }
    };
}

#endif

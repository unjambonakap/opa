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

#ifndef PLLL_INCLUDE_GUARD__SVP_IMPL_CPP
#define PLLL_INCLUDE_GUARD__SVP_IMPL_CPP

namespace plll
{
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    void Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    applyRandomUnimodularTransformation(Lattice<RealTypeContext, IntTypeContext> & lattice)
    {
        typename LatticeHelper<RealTypeContext, IntTypeContext>::RangePopper rangepopper(lattice);
        unsigned dim = lattice.range().dimension();
        if (dim < 2)
            return;
        // Initialize random generator
        arithmetic::RandomNumberGenerator rng;
        rng.randomizeSeed();
        typename IntTypeContext::UniformRNG urng(rng);
        // Select number of transformations
        unsigned iterations = dim * 2 + rng.random(dim);
        // Apply random transformations
        typename IntTypeContext::Integer m, b, bb;
        arithmetic::convert(b, 3, lattice.ic());
        bb = b << 1;
        ++bb;
        for (unsigned i = 0; i < iterations; ++i)
        {
            // Select two vectors
            unsigned j = rng.random(dim), k = j;
            while (j == k)
                k = rng.random(dim);
            // Select modifier
            urng.random(m, bb);
            m -= b;
            // Transform
            lattice.add(lattice.range().begin() + j, m, lattice.range().begin() + k);
        }
        // Now apply random permutation
        for (unsigned i = 1; i < dim; ++i)
        {
            unsigned j = rng.random(dim - i + 1);
            while (j)
            {
                lattice.swap(lattice.range().begin() + i + j - 1, lattice.range().begin() + i + j - 2);
                --j;
            }
        }
    }

    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    template<class Reorderer>
    void Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    solveSVPOnce(Lattice<RealTypeContext, IntTypeContext> & lattice, linalg::math_rowvector<typename IntTypeContext::Integer> & result, Reorderer & reorder)
    // Finds a shortest vector in the lattice generated by the orthogonal projections of the vectors
    // A.row(begin) to A.row(end) into the orthogonal complement of the vectors A.row(0) to
    // A.row(begin-1).
    //
    // In case no shortest vector was found. result is a vector of length 0.
    {
        typename LatticeHelper<RealTypeContext, IntTypeContext>::RangePopper rangepopper(lattice);
        // First step: preprocessing
        
        // Perform preprocessing
        enumerator().preprocess(*this, lattice);
        
        // Resize result storage
        result.resize(lattice.range().dimension());
        
        // Compute shortest vector
        if (lattice.range().dimension() <= 1)
        {
            // The trivial cases: dimension 0 [the empty vector] and 1 [any basis vector]
            if (lattice.range().dimension() == 1)
                setOne(result[0]);
        }
        else
        {
            // Recompute settings if necessary
            enumerator().updateSettings(lattice);
        
            // Make sure GS is really up to date
            lattice.update(lattice.range().end() + 1);
        
            // Compute enumeration bound
            typename RealTypeContext::Real bound(lattice.rc());
            enumerator().computeEnumBound(*this, lattice, bound);
        
            // Do enumeration
            if (!enumerator().enumerate(lattice, lattice.range().begin(), lattice.range().end(), result, bound,
                                        enumerator().getCurrentAlgorithm(), enumerator().getCurrentEnumSettings(), enumerator().isToplevelEnum()))
                // In case nothing was found, rescale result array to size 0
                result.resize(0);
        }
    }
    
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    template<class Reorderer>
    void Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    solveSVPRepeat(Lattice<RealTypeContext, IntTypeContext> & lattice, linalg::math_rowvector<typename IntTypeContext::Integer> & result,
                   Reorderer & reorder, unsigned number_of_iterations)
    // Finds a shortest vector in the lattice generated by the orthogonal projections of the vectors
    // A.row(begin) to A.row(end) into the orthogonal complement of the vectors A.row(0) to
    // A.row(begin-1).
    //
    // In case no shortest vector was found. result is a vector of length 0.
    {
        typename LatticeHelper<RealTypeContext, IntTypeContext>::RangePopper rangepopper(lattice);
        assert(number_of_iterations > 1);
    
        // Perform iterations
        linalg::math_rowvector<typename IntTypeContext::Integer> result_buffer;
        bool restore_state = false;
        unsigned iteration_stored = 0;
        STD_AUTO_PTR<typename Lattice<RealTypeContext, IntTypeContext>::DuplicateStorage> state;
        typename RealTypeContext::Real current_len(lattice.rc());
        current_len = lattice.getNormSq(lattice.range().begin());
        for (unsigned iteration = 0; iteration < number_of_iterations; ++iteration)
        {
            // If this is not the first iteration, perform random transformation
            if (iteration)
            {
                // If this is the first iteration, make sure we have something to restore. If it is not,
                // restore.
                if (restore_state == false)
                {
                    // Store state
                    restore_state = true;
                    state.reset(new typename Lattice<RealTypeContext, IntTypeContext>::DuplicateStorage(lattice));
                }
                else
                {
                    if (iteration_stored + 1 < iteration)
                    {
                        // Restore state
                        lattice.copyFrom(*state);
                    }
                }
            
                // Reduce with LLL
                MoveShortestFirst<RealTypeContext, IntTypeContext> reorder2(lattice.getStatistics());
                lattice.range().dupRange();
                partialLLL(lattice, lattice.range().begin() + 1, reorder2);
                
                // Insert shortest vector so far
                lattice.insertVectorLC(lattice.range().begin(), result_buffer);
                
                // Apply random transformation to vectors after the newly inserted one
                lattice.range().addRange(lattice.range().begin() + 1, lattice.range().end());
                applyRandomUnimodularTransformation(lattice);
                
                // Reduce with LLL
                lattice.range().dupRange();
                partialLLL(lattice, lattice.range().begin() + 1, reorder2);
            }
            
            lattice.range().dupRange();
            solveSVPOnce(lattice, result, reorder);
            
            // If nothing was found, check first basis vector: maybe it is shorter than before
            if (result.size() == 0)
            {
                if (lattice.getNormSq(lattice.range().begin()) < current_len)
                {
                    result.resize(lattice.range().dimension());
                    setOne(result[0]);
                    for (unsigned i = 1; i < result.size(); ++i)
                        setOne(result[i]);
                }
            }
            // Found something?
            if (result.size() > 0)
            {
                // Compute length
                typename RealTypeContext::Real len(lattice.rc());
                lattice.computeProjectionLength(len, lattice.range().begin(), lattice.range().begin(), result);
                verbose()(LatticeReduction::VL_Information) << iteration << ": found " << len;
                // Take new vector if it is shorter
                if (len < current_len)
                {
                    linalg::swap(result_buffer, result);
                    current_len = len;
                    verbose()(LatticeReduction::VL_Information) << "    new shortest vector has length " << current_len;
                    state.reset();
                    if (iteration + 1 < number_of_iterations)
                    {
                        restore_state = true;
                        state.reset(new typename Lattice<RealTypeContext, IntTypeContext>::DuplicateStorage(lattice));
                        iteration_stored = iteration;
                    }
                    else
                    {
                        restore_state = false;
                    }
                }
            }
            else
                verbose()(LatticeReduction::VL_Information) << iteration << ": found nothing";
        }
        
        // Restore state
        if (restore_state && (state.get() != NULL))
            lattice.copyFrom(*state);
        // Return minimum
        linalg::swap(result_buffer, result);
    }
    
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    template<class Reorderer>
    void Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    solveSVP(Lattice<RealTypeContext, IntTypeContext> & lattice, linalg::math_rowvector<typename IntTypeContext::Integer> & result,
             Reorderer & reorder, bool dontModify)
    // Finds a shortest vector in the lattice generated by the orthogonal projections of the vectors
    // A.row(begin) to A.row(end) into the orthogonal complement of the vectors A.row(0) to
    // A.row(begin-1).
    //
    // In case no shortest vector was found. result is a vector of length 0.
    {
        typename LatticeHelper<RealTypeContext, IntTypeContext>::RangePopper rangepopper(lattice);
        if (!lattice.allowedToModify() || dontModify)
        {
            // Create and work with copy
            typename Lattice<RealTypeContext, IntTypeContext>::DuplicateStorage work_lattice(lattice);
            SimpleTransformationCollector<IntTypeContext> dt(lattice.dimension());
            work_lattice.d_lattice.addNotifier(&dt);
            work_lattice.d_lattice.range().dupRange();
            solveSVP(work_lattice.d_lattice, result, reorder);
            if (result.size() > 0)
            {
                // Transform result back
                linalg::math_rowvector<typename IntTypeContext::Integer> v(result.size());
                const linalg::math_matrix<typename IntTypeContext::Integer> & T = dt.transform();
                for (unsigned i = 0; i < result.size(); ++i)
                    if (!isZero(result[i]))
                    {
                        for (unsigned j = 0; j < result.size(); ++j)
                            v[j] += result[i] * T(work_lattice.d_lattice.range().begin() + i, work_lattice.d_lattice.range().begin() + j);
                    }
                linalg::swap(result, v);
            }
            return;
        }
        
        // Find out number of iterations
        unsigned number_of_iterations = enumerator().getCurrentNumberOfEnumRepetitions();
        if (number_of_iterations == 0)
        {
            result.resize(0);
            return;
        }
        else if (number_of_iterations == 1)
        {
            lattice.range().dupRange();
            solveSVPOnce(lattice, result, reorder);
            return;
        }
        else
        {
            lattice.range().dupRange();
            solveSVPRepeat(lattice, result, reorder, number_of_iterations);
            return;
        }
    }
    
    template<class RealTypeContext, class IntTypeContext>
    void rearrangeHelper(Lattice<RealTypeContext, IntTypeContext> & lattice, linalg::math_rowvector<typename IntTypeContext::Integer> & result,
                         unsigned i, unsigned j,
                         typename IntTypeContext::Integer & d, typename IntTypeContext::Integer & x, typename IntTypeContext::Integer & y)
    {
        if (!isZero(result[i]))
        {
            if (isZero(result[j]))
            {
                // Swap A.row(begin + i) with A.row(begin + i)
                lattice.swap(lattice.range().begin() + j, lattice.range().begin() + i);
                arithmetic::swap(result[i], result[j]);
            }
            else
            {
                // Both are non-zero
                XGCD(d, x, y, result[j], result[i]);
                result[j] /= d;
                result[i] /= d;
                y = -y;
                lattice.trans(lattice.range().begin() + j, lattice.range().begin() + i, result[j], result[i], y, x);
                result[j] = d;
                setZero(result[i]);
            }
        }
    }

    template<class RealTypeContext, class IntTypeContext>
    void rearrangeDualHelper(Lattice<RealTypeContext, IntTypeContext> & lattice, linalg::math_rowvector<typename IntTypeContext::Integer> & result,
                             unsigned i, unsigned j,
                             typename IntTypeContext::Integer & d, typename IntTypeContext::Integer & x, typename IntTypeContext::Integer & y)
    {
        // Move (result[i], result[j]) to (0, d), where d = gcd(result[i],result[j])
        if (!isZero(result[i]))
        {
            if (isZero(result[j]))
            {
                // Do a simple swap
                lattice.swap(lattice.range().end() - i, lattice.range().end() - j);
                // Update v
                arithmetic::swap(result[i], result[j]);
            }
            else
            {
                // First, compute d and Bezout coefficients x,y such that x * result[i] + y * result[j] == d
                XGCD(d, x, y, result[i], result[j]);
            
                // Then, apply transformation [result[j]/d  -result[i]/d] to vectors end-i, end-i+1
                //                            [x             y          ]
                result[j] /= d;
                result[i] /= d;
                result[i] = -result[i];
                lattice.trans(lattice.range().end() - i, lattice.range().end() - j, result[j], result[i], x, y); 
            
                // Update v
                setZero(result[i]);
                result[j] = d;
            }
        }
    }
    
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    void Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    rearrange(Lattice<RealTypeContext, IntTypeContext> & lattice, linalg::math_rowvector<typename IntTypeContext::Integer> & result, bool minimizeCoeffs)
    // Rearranges A.row(begin) to A.row(end) such that the vector described by the linear combination in
    // result (with end-begin+1 entries) is at position A.row(begin). Note that during the process,
    // result might be arbitrarily changed. Also note that the GCD of the entries of result must be 1.
    //
    // Uses trans(), swap() and flip().
    {
        typename LatticeHelper<RealTypeContext, IntTypeContext>::RangePopper rangepopper(lattice);
        const unsigned dim = lattice.range().dimension();
        assert(result.size() == dim);
        if (dim == 0)
            return;
        ++lattice.getStats().vectorinsertions_rearrange;
    
        // Work back from the end
        typename IntTypeContext::Integer d, x, y;
        if (minimizeCoeffs)
        {
            // This is the approach used by Gama and Nguyen in Algorithm 3 of "Finding Short Lattice
            // Vectors within Mordell's Inequality".
            for (unsigned a = 1; a < result.size(); a <<= 1)
                for (unsigned b = 0; b < result.size() - a; b += 2*a)
                    rearrangeHelper(lattice, result, b + a, b, d, x, y);
        }
        else
        {
            for (unsigned i = dim - 1; i > 0; --i)
                rearrangeHelper(lattice, result, i, i - 1, d, x, y);
        }
        
        // This is only needed when dim == 1, or maybe when no GCD was computed
        if (sign(result[0]) < 0)
        {
            // Flip sign of A.row(begin)
            lattice.flip(lattice.range().begin());
            result[0] = -result[0];
        }
    }
    
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    void Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    rearrangeDual(Lattice<RealTypeContext, IntTypeContext> & lattice, linalg::math_rowvector<typename IntTypeContext::Integer> & result, bool minimizeCoeffs)
    // Rearranges A.row(begin) to A.row(end) such that in the dual of the projected sublattice of these
    // vectors, the first vector equals the given linear combination of the dual lattice's "canonical"
    // reversed basis. Note that during the process, result might be arbitrarily changed. Also note that
    // the GCD of the entries of result must be 1.
    // 
    // Uses trans(), swap() and flip().
    {
        typename LatticeHelper<RealTypeContext, IntTypeContext>::RangePopper rangepopper(lattice);
        const unsigned dim = lattice.range().dimension();
        assert(result.size() == dim);
        if (dim == 0)
            return;
        ++lattice.getStats().vectorinsertions_rearrange;
        
        /*
          Let R be the d_coeff (sub-)matrix of the original lattice. We want
              R^{-1} [result[i]  result[i-1]] = R^{-1} T^{-1} T [result[i]  result[i-1]]^t
                                              = R^{-1} T^{-1} [0  d]^t = (T R)^{-1} [0  d]^t,
          where d = gcd(result[i], result[i-1]), and x, y are integers such that x*result[i] +
          y*result[i-1] = d. Then we must have
               T = [ result[i-1]/d  -result[i]/d ]
                   [ x               y           ].
          
          In case result[i] = 0, there is nothing to do. In case result[i-1] = 0 and result[i] != 0, we
          have result[i] = d = gcd(result[i], 0) = 1*result[i] + 0*0, whence x = 1, y = 0 can be
          chosen. Therefore,
               T = [ 0  -1 ]
                   [ 1   0 ].
          Since we don't care that the determinant is +1, a simple swap suffices.
        */
    
        typename IntTypeContext::Integer d, x, y;
        if (minimizeCoeffs)
        {
            // This is the approach used by Gama and Nguyen in Algorithm 3 of "Finding Short Lattice
            // Vectors within Mordell's Inequality".
            for (unsigned a = 1; a < result.size(); a <<= 1)
                for (unsigned b = 0; b < result.size() - a; b += 2*a)
                    rearrangeDualHelper(lattice, result, b + a, b, d, x, y);
        }
        else
        {
            for (unsigned i = result.size() - 1; i > 0; --i)
                rearrangeDualHelper(lattice, result, i, i - 1, d, x, y);
        }
        
        // This is only needed when dim == 1, or maybe when no GCD was computed
        if (sign(result[0]) < 0)
        {
            // Flip sign of A.row(begin)
            lattice.flip(lattice.range().end());
            result[0] = -result[0];
        }
    }

    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    template<class Reorderer>
    unsigned Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    rearrangeLLL(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned reducebegin,
                 linalg::math_rowvector<typename IntTypeContext::Integer> & result, Reorderer & reorder)
    // Inserts a linear combination of the vectors A.row(begin), ..., A.row(begin+result.size()-1)
    // before A.row(begin) and applies LLL on the range A.row(reducebegin), ..., A.row(end+1) to remove
    // the induced dependency.
    {
        typename LatticeHelper<RealTypeContext, IntTypeContext>::RangePopper rangepopper(lattice);
        unsigned origdim = lattice.range().dimension();
        if (result.size() > 0)
        {
            // Is result = +-e_i?
            bool is_e_i = true;
            int idx = -1;
            for (unsigned i = 0; i < result.size(); ++i)
                if (!isZero(result[i]))
                {
                    if (idx != -1)
                    {
                        is_e_i = false;
                        break;
                    }
                    idx = i;
                }
            assert(idx != -1);
            if (idx == -1)
                // result is the zero vector!
                return lattice.range().begin();
            
            if (is_e_i)
            {
                // Just move the (begin+idx)-th vector to the begin-th position and apply LLL to "clean up"
                assert(isPMOne(result[idx]));
                while (idx > 0)
                {
                    lattice.swap(lattice.range().begin() + idx - 1, lattice.range().begin() + idx);
                    --idx;
                }
                ++lattice.getStats().vectorinsertions_rearrange;
            }
            else
                lattice.insertVectorLC(lattice.range().begin(), result);
        }
        unsigned beg = lattice.range().begin() + 1;
        lattice.range().addRange(reducebegin, lattice.range().end());
        unsigned ms = partialLLL(lattice, beg, reorder);
        // Early abort?
        if (lattice.range().dimension() > origdim)
            throw reduction_error("Removing dependency using LLL failed!");
        return ms;
    }
    
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    template<class Reorderer>
    inline void Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    solveSVPRearrange(Lattice<RealTypeContext, IntTypeContext> & lattice, Reorderer & reorder, bool minimizeCoeffs)
    // Finds a shortest vector in the lattice generated by the orthogonal projections of the vectors
    // A.row(begin) to A.row(end) into the orthogonal complement of the vectors A.row(0) to
    // A.row(begin-1). Modifies the basis so that A.row(begin) is said vector. Only A.row(begin) to
    // A.row(end) are modified.
    {
        typename LatticeHelper<RealTypeContext, IntTypeContext>::RangePopper rangepopper(lattice);
        linalg::math_rowvector<typename IntTypeContext::Integer> result;
        lattice.range().dupRange();
        solveSVP(lattice, result, reorder);
        lattice.range().dupRange();
        rearrange(lattice, result, minimizeCoeffs);
    }
    
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    template<class Reorderer>
    inline unsigned Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    solveSVPRearrangeLLL(Lattice<RealTypeContext, IntTypeContext> & lattice,
                         unsigned reducebegin, Reorderer & reorder)
    // Finds a shortest vector in the lattice generated by the orthogonal projections of the vectors
    // A.row(begin) to A.row(end) into the orthogonal complement of the vectors A.row(0) to
    // A.row(begin-1). Modifies the basis so that A.row(begin) is said vector. Only A.row(begin) to
    // A.row(end) are modified.
    {
        typename LatticeHelper<RealTypeContext, IntTypeContext>::RangePopper rangepopper(lattice);
        linalg::math_rowvector<typename IntTypeContext::Integer> result;
        lattice.range().dupRange();
        solveSVP(lattice, result, reorder);
        lattice.range().dupRange();
        return rearrangeLLL(lattice, reducebegin, result, reorder);
    }
    
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    inline void Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    computeSVPbasis(Lattice<RealTypeContext, IntTypeContext> & lattice, bool make_basis)
    // Computes a "SVP basis" (i.e. first vector is shortest) for the lattice generated by the orthogonal
    // projections of the vectors A.row(begin) to A.row(end) into the orthogonal complement of the
    // vectors A.row(0) to A.row(begin-1). Modifies the basis so that A.row(begin) to A.row(end) is said
    // basis (or more precisely, its preimage).
    {
        typename LatticeHelper<RealTypeContext, IntTypeContext>::RangePopper rangepopper(lattice);
        // Prepare GS
        lattice.update(lattice.range().end() + 1);
        
        // Find shortest vector and put to beginning.
        NoReorder<RealTypeContext, IntTypeContext> reorder;
        linalg::math_rowvector<typename IntTypeContext::Integer> result;
        lattice.range().dupRange();
        solveSVP(lattice, result, reorder);
        if (make_basis)
        {
            if (lattice.avoidInsertions())
            {
                lattice.range().dupRange();
                rearrange(lattice, result);
            }
            else
            {
                // Use LLL to turn this into a basis
                lattice.range().dupRange();
                rearrangeLLL(lattice, lattice.range().begin(), result, reorder);
            }
        }
        else
        {
            // Insert vector
            lattice.insertVectorLC(lattice.range().begin(), result);
            // Do size reduction
            lattice.sizereduce(lattice.range().begin());
            // Do not remove dependencies
        }
    }
    
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    inline bool Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    computeDSVPbasis(Lattice<RealTypeContext, IntTypeContext> & lattice)
    // Computes a "Dual-SVP basis" (i.e. first vector of canonical reversed projected dual is shortest)
    // for the lattice generated by the orthogonal projections of the vectors A.row(begin) to A.row(end)
    // into the orthogonal complement of the vectors A.row(0) to A.row(begin-1). Modifies the basis so
    // that A.row(begin) to A.row(end) is said basis (or more precisely, its preimage).
    //
    // The returned bool indicates whether a new shortest vector was inserted (true) or not (false).
    {
        typename LatticeHelper<RealTypeContext, IntTypeContext>::RangePopper rangepopper(lattice);
        // Prepare GS
        lattice.update(lattice.range().end() + 1);
        
        // Create interface for canonical reverse dual basis
        STD_AUTO_PTR<GSInterface<RealTypeContext, IntTypeContext> > gs(getDualGSI(lattice.rc(), lattice.ic(), lattice));
        Lattice<RealTypeContext, IntTypeContext> dl(gs.get(), lattice.rc(), lattice.ic(), lattice.getStats(), 0, gs->getDimension() - 1);
        SimpleTransformationCollector<IntTypeContext> dt(dl.dimension());
        dl.addNotifier(&dt);
        // Find shortest vector
        typename RealTypeContext::Real ol(lattice.rc()), l(lattice.rc());
        ol = dl.getNormSq(0);
        linalg::math_rowvector<typename IntTypeContext::Integer> result;
        NoReorder<RealTypeContext, IntTypeContext> reorder;
        dl.range().dupRange();
        solveSVP(dl, result, reorder);
        bool ins = false;
        if (result.size() > 0)
        {
            dl.computeProjectionLength(l, 0, 0, result);
            if (l < lattice.getLLLalpha() * ol) // use LLL alpha here!
            {
                result = result * dt.transform();
                lattice.range().dupRange();
                rearrangeDual(lattice, result);
                ins = true;
            }
        }
        return ins;
    }
    
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    inline void Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    computeHKZbasis(Lattice<RealTypeContext, IntTypeContext> & lattice)
    // Computes a Hermite-Korkine-Zolotarev basis for the lattice generated by the orthogonal
    // projections of the vectors A.row(begin) to A.row(end) into the orthogonal complement of the
    // vectors A.row(0) to A.row(begin-1). Modifies the basis so that A.row(begin) to A.row(end) is said
    // basis (or more precisely, its preimage).
    {
        typename LatticeHelper<RealTypeContext, IntTypeContext>::RangePopper rangepopper(lattice);
        // Prepare GS
        lattice.update(lattice.range().end() + 1);
        bool avoidInsertions = lattice.avoidInsertions();
        
        for (unsigned i = lattice.range().begin(); i < lattice.range().end(); ++i)
        // No need to do this for i = end: a one-dimensional lattice's basis is always a HKZ basis.
        {
            // Callback
            callback(lattice, i);
            
            // Find and insert shortest vector using LLL
            NoReorder<RealTypeContext, IntTypeContext> reorder;
            lattice.range().addRange(i, lattice.range().end());
            if (avoidInsertions)
                solveSVPRearrange(lattice, reorder);
            else
                solveSVPRearrangeLLL(lattice, i, reorder);
        }
        
        // Only size reduction has to be done for i = end as well
        lattice.sizereduce(lattice.range().end());
    }
    
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    inline void Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    computeDHKZbasis(Lattice<RealTypeContext, IntTypeContext> & lattice)
    // Computes a "Dual-HKZ basis" (i.e. canonical reversed projected dual is HKZ reduced) for the
    // lattice generated by the orthogonal projections of the vectors A.row(begin) to A.row(end) into
    // the orthogonal complement of the vectors A.row(0) to A.row(begin-1). Modifies the basis so that
    // A.row(begin) to A.row(end) is said basis (or more precisely, its preimage).
    {
        typename LatticeHelper<RealTypeContext, IntTypeContext>::RangePopper rangepopper(lattice);
        // Prepare GS
        lattice.update(lattice.range().end() + 1);
        
        // Create interface for canonical reverse dual basis
        STD_AUTO_PTR<GSInterface<RealTypeContext, IntTypeContext> > gs(getDualGSI(lattice.rc(), lattice.ic(), lattice));
        Lattice<RealTypeContext, IntTypeContext> dl(gs.get(), lattice.rc(), lattice.ic(), lattice.getStats(), 0, gs->getDimension() - 1);
        DualTransformationApplier<RealTypeContext, IntTypeContext> dt(lattice, dl);
        
        // Compute HKZ basis of dual
        dl.range().dupRange();
        computeHKZbasis(dl);
    }
}

#endif
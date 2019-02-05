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

#ifndef PLLL_INCLUDE_GUARD__VERIFY_IMPL_CPP
#define PLLL_INCLUDE_GUARD__VERIFY_IMPL_CPP

namespace plll
{
    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    bool isSR(Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction> & workspace,
              Lattice<RealTypeContext, IntTypeContext> & lattice)
    {
        // Compute GS coefficients
        lattice.update(lattice.range().end() + 1);
    
        // Check if size reduced
        for (unsigned i = lattice.range().begin(); i <= lattice.range().end(); ++i)
            for (unsigned j = lattice.range().begin(); j < i; ++j)
                if (arithmetic::abs(lattice.getCoeff(i, j)) > lattice.getZDF())
                {
                    if (workspace.verbose().yieldsOutput(LatticeReduction::VL_Information))
                        workspace.verbose()(LatticeReduction::VL_Information) << "Not size reduced (" << i << ", " << j << ") [" << lattice.getCoeff(i, j) << "]\n";
                    return false;
                }
        return true;
    }

    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction, class LCondition>
    bool isLLL(Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction> & workspace,
               Lattice<RealTypeContext, IntTypeContext> & lattice, LCondition & condition)
    {
        // Compute GS coefficients
        lattice.update(lattice.range().end() + 1);
    
        // Check for Lovász condition
        for (unsigned k = lattice.range().begin() + 1; k <= lattice.range().end(); ++k)
            if (condition(lattice, k))
            {
                if (workspace.verbose().yieldsOutput(LatticeReduction::VL_Information))
                    workspace.verbose()(LatticeReduction::VL_Information) << "LLL condition failed for [" << k - 1 << "," << k << "]\n";
                return false;
            }
        return true;
    }

    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    bool isShortestVector(Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction> & workspace,
                          Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned begin, unsigned end)
    {
        // Find shortest vector
        linalg::math_rowvector<typename IntTypeContext::Integer> result;
        lattice.range().addRange(begin, end);
        NoReorder<RealTypeContext, IntTypeContext> noreorder;
        workspace.solveSVP(lattice, result, noreorder);
        if (result.size() == 0)
            return true;
    
        // Compute length of projection
        typename RealTypeContext::Real l(lattice.rc());
        lattice.computeProjectionLength(l, begin, begin, result);
    
        // Compare
        if (lattice.canGetIntegerVector() && (begin == 0))
        {
            // Compare (integral) norms
            linalg::math_rowvector<typename IntTypeContext::Integer> res;
            lattice.getLinearCombination(res, lattice.range().begin(), result);
            return normSq(res) >= lattice.getNormSqUP(0);
        }
        else
            return l >= lattice.getNormSq(begin);
    }

    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    bool isShortestVector(Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction> & workspace,
                          Lattice<RealTypeContext, IntTypeContext> & lattice,
                          unsigned begin, unsigned end, const typename RealTypeContext::Real & alpha)
    {
        // Find shortest vector
        linalg::math_rowvector<typename IntTypeContext::Integer> result;
        lattice.range().addRange(begin, end);
        NoReorder<RealTypeContext, IntTypeContext> noreorder;
        workspace.solveSVP(lattice, result, noreorder);
    
        // Compute length of projection
        typename RealTypeContext::Real l(lattice.rc());
        lattice.computeProjectionLength(l, begin, begin, result);
    
        // Compare
        return l >= lattice.getNormSq(begin) * alpha;
    }

    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    bool isBKZ(Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction> & workspace,
               Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize)
    {
        // Compute GS coefficients
        lattice.update(lattice.range().end() + 1);
    
        // Check for BKZ condition
        for (unsigned k = lattice.range().begin(); k < lattice.range().end(); ++k)
        {
            unsigned bend = std::min(k + blocksize - 1, lattice.range().end());
        
            if (!isShortestVector(workspace, lattice, k, bend, lattice.getLLLalpha()))
            {
                if (workspace.verbose().yieldsOutput(LatticeReduction::VL_Information))
                    workspace.verbose()(LatticeReduction::VL_Information) << "BKZ condition failed for [" << k << "," << bend << "]\n";
                return false;
            }
        }
        return true;
    }

    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    bool isBKZ_SemiBlock2k(Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction> & workspace,
                           Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize)
    {
        unsigned blockno = (lattice.range().dimension() + blocksize - 1) / blocksize;
        // Check HKZ conditions
        for (unsigned b = 0; b < blockno; ++b)
        {
            lattice.range().addRange(getBlock(b, lattice, blocksize));
            typename LatticeHelper<RealTypeContext, IntTypeContext>::RangePopper rangepopper(lattice);
            if (!isHKZ(workspace, lattice))
                return false;
        }
        // Check block interconnect conditions
        typename RealTypeContext::Real factor(lattice.rc()), t(lattice.rc());
        // Set t to alpha - 1/4
        setOne(t);
        t >>= 2;
        t = -t;
        t += lattice.getLLLalpha();
        // Now set factor to 1/t
        setOne(factor);
        factor /= t;
        for (unsigned b = 1; b < blockno; ++b)
        {
            unsigned i = lattice.range().begin() + b * blocksize;
            if (lattice.getNormSq(i - 1) > factor * lattice.getNormSq(i))
            {
                if (workspace.verbose().yieldsOutput(LatticeReduction::VL_Information))
                    workspace.verbose()(LatticeReduction::VL_Information) << "Block interconnect condition failed for blocks " << b-1 << " and " << b << "!\n";
                return false;
            }
        }
        // Check for determinant condition
        typename RealTypeContext::Real d1(lattice.rc()), d2(lattice.rc()), betak(lattice.rc());
        workspace.computeProjectedDeterminant(d2, lattice, getBlock(0, lattice, blocksize));
        // We use the upper bound ((1 + blocksize/2)^{2 \log 2 + 1/blocksize})^k for beta^k. (See:
        // Schnorr: "Blockwise Lattice Reduction Revisited", before Definition 1 in Section 2.)
        convert(betak, std::pow((double)(blocksize + 2) * 0.5, blocksize * std::log(4.0) + 1.0), lattice.rc());
        for (unsigned b = 1; b < blockno; ++b)
        {
            arithmetic::swap(d1, d2);
            workspace.computeProjectedDeterminant(d2, lattice, getBlock(b, lattice, blocksize));
            if (lattice.getLLLalpha() * d1 > betak * d2)
            {
                if (workspace.verbose().yieldsOutput(LatticeReduction::VL_Information))
                    workspace.verbose()(LatticeReduction::VL_Information) << "Determinant condition failed for blocks " << b-1 << " and " << b << "!\n";
                return false;
            }
        }
        return true;
    }

    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    bool isBKZ_PrimalDual(Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction> & workspace,
                          Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize)
    {
        unsigned blockno = (lattice.range().dimension() + blocksize - 1) / blocksize;
        // Check primal conditions
        for (unsigned b = 0; b < blockno; ++b)
        {
            lattice.range().addRange(getBlock(b, lattice, blocksize));
            typename LatticeHelper<RealTypeContext, IntTypeContext>::RangePopper rangepopper(lattice);
            if (!isHKZ(workspace, lattice))
                return false;
        }
        // Check dual conditions
        typename RealTypeContext::Real factor(lattice.rc()), t(lattice.rc());
        // Set t to alpha - 1/4
        setOne(t);
        t >>= 2;
        t = -t;
        t += lattice.getLLLalpha();
        // Now set factor to 1/t
        setOne(factor);
        factor /= t;
        for (unsigned b = 0; b < blockno; ++b)
        {
            std::pair<unsigned, unsigned> block = getBlock(b, lattice, blocksize);
            if ((block.second <= block.first) || (block.second + 1 > lattice.range().end()))
                continue;
            // Create interface for canonical reverse dual basis
            lattice.range().addRange(block);
            typename LatticeHelper<RealTypeContext, IntTypeContext>::RangePopper rangepopper(lattice);
            STD_AUTO_PTR<GSInterface<RealTypeContext, IntTypeContext> > gs(workspace.getDualGSI(lattice.rc(), lattice.ic(), lattice));
            Lattice<RealTypeContext, IntTypeContext> dl(gs.get(), lattice.rc(), lattice.ic(), lattice.getStats(), 0, gs->getDimension() - 1);
            // Find shortest vector
            linalg::math_rowvector<typename IntTypeContext::Integer> result;
            dl.range().dupRange();
            NoReorder<RealTypeContext, IntTypeContext> noreorder;
            workspace.solveSVP(dl, result, noreorder);
            if (result.size() == 0)
                continue;
            // Compute length of projection multiplied by alpha-1/4
            typename RealTypeContext::Real l(dl.rc()), t(dl.rc());
            dl.computeProjectionLength(t, 0, 0, result);
            setOne(l);
            l >>= 2;
            l = -l;
            l += lattice.getLLLalpha();
            l /= t;
            // Check for Siegel condition
            if (l > lattice.getNormSq(block.second))
            {
                if (workspace.verbose().yieldsOutput(LatticeReduction::VL_Information))
                    workspace.verbose()(LatticeReduction::VL_Information) << "Siegel condition connecting blocks " << b << " and " << b + 1 << " failed!\n";
                return false;
            }
        }
        return true;
    }

    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    bool isBKZ_Slide(Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction> & workspace,
                     Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize)
    {
        unsigned blockno = (lattice.range().dimension() + blocksize - 1) / blocksize;
        // Check primal conditions
        for (unsigned b = 0; b < blockno; ++b)
        {
            lattice.range().addRange(getBlock(b, lattice, blocksize));
            typename LatticeHelper<RealTypeContext, IntTypeContext>::RangePopper rangepopper(lattice);
            if (!isHKZ(workspace, lattice))
                return false;
        }
        // Check dual conditions
        for (unsigned b = 0; b < blockno; ++b)
        {
            std::pair<unsigned, unsigned> block = getBlock(b, lattice, blocksize, 1);
            if (block.second <= block.first)
                continue;
            // Create interface for canonical reverse dual basis
            lattice.range().addRange(block);
            typename LatticeHelper<RealTypeContext, IntTypeContext>::RangePopper rangepopper(lattice);
            STD_AUTO_PTR<GSInterface<RealTypeContext, IntTypeContext> > gs(workspace.getDualGSI(lattice.rc(), lattice.ic(), lattice));
            Lattice<RealTypeContext, IntTypeContext> dl(gs.get(), lattice.rc(), lattice.ic(), lattice.getStats(), 0, gs->getDimension() - 1);
            // Check
            if (!isShortestVector(workspace, dl, 0, gs->getDimension() - 1, lattice.getLLLalpha()))
            {
                if (workspace.verbose().yieldsOutput(LatticeReduction::VL_Information))
                    workspace.verbose()(LatticeReduction::VL_Information) << "Dual condition for block [" << block.first << "," << block.second << "] failed!\n";
                return false;
            }
        }
        return true;
    }

    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    bool isBKZ_ImprovedSlide(Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction> & workspace,
                             Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned blocksize)
    {
        unsigned blockno = (lattice.range().dimension() + blocksize - 1) / blocksize;
        if (blockno == 0)
            return true;
    
        // Check primal condition for block 0
        lattice.range().addRange(getBlock(0, lattice, blocksize));
        {
            typename LatticeHelper<RealTypeContext, IntTypeContext>::RangePopper rangepopper(lattice);
            if (!isSVP(workspace, lattice, false))
            {
                if (workspace.verbose().yieldsOutput(LatticeReduction::VL_Information))
                    workspace.verbose()(LatticeReduction::VL_Information) << "Block [" << lattice.range().begin() << ", " << getBlock(0, lattice, blocksize).second << "] is not SVP reduced!\n";
                return false;
            }
        }
        
        if (blockno == 1)
            return true;
        
        // Compute determinants
        linalg::math_rowvector<typename RealTypeContext::Real> D(blockno);
        for (unsigned i = 0; i < blockno; ++i)
        {
            D[i].setContext(lattice.rc());
            workspace.computeProjectedDeterminant(D[i], lattice, getBlock(i, lattice, blocksize));
        }
        
        // Determine maximal quotient
        typename RealTypeContext::Real maximal(lattice.rc()), t(lattice.rc());
        maximal = D[0] / D[1];
        for (unsigned i = 2; i < blockno; ++i)
        {
            t = D[i - 1] / D[i];
            if (t > maximal)
                maximal = t;
        }
        if (workspace.verbose().yieldsOutput(LatticeReduction::VL_Chatter))
        {
            for (unsigned i = 0; i < blockno; ++i)
            {
                Verbose::VerboseStream s = workspace.verbose()(LatticeReduction::VL_Chatter);
                s << "det[" << i << "] = " << D[i];
                if (i + 1 < blockno)
                    s << ", q = " << D[i] / D[i+1];
            }
            workspace.verbose()(LatticeReduction::VL_Chatter) << "Maximal value: " << maximal << "\n";
        }
    
        // Determine all ell's for which quotient is maximal
        for (unsigned ell = 0; ell + 1 < blockno; ++ell)
        {
            bool is_ok = false;
            t = D[ell] / D[ell + 1];
            if (t >= maximal)
            {
                if (workspace.verbose().yieldsOutput(LatticeReduction::VL_Chatter))
                    workspace.verbose()(LatticeReduction::VL_Chatter) << "Maximal value at " << ell << "\n";
                // Check primal condition for block ell
                lattice.range().addRange(getBlock(ell, lattice, blocksize));
                typename LatticeHelper<RealTypeContext, IntTypeContext>::RangePopper2 rangepopper(lattice);
                if (!isSVP(workspace, lattice, false))
                {
                    if (workspace.verbose().yieldsOutput(LatticeReduction::VL_Information))
                        workspace.verbose()(LatticeReduction::VL_Information) << "Block " << ell << " is not SVP reduced!\n";
                }
                else
                {
                    rangepopper.pop();
                    // Check primal condition for block ell+1
                    lattice.range().addRange(getBlock(ell + 1, lattice, blocksize));
                    typename LatticeHelper<RealTypeContext, IntTypeContext>::RangePopper2 rangepopper2(lattice);
                    if (!isSVP(workspace, lattice, false))
                    {
                        if (workspace.verbose().yieldsOutput(LatticeReduction::VL_Information))
                            workspace.verbose()(LatticeReduction::VL_Information) << "Block " << ell + 1 << " is not SVP reduced!\n";
                    }
                    else
                    {
                        rangepopper2.pop();
                        // Check dual condition for slided block ell
                        std::pair<unsigned, unsigned> block = getBlock(ell, lattice, blocksize, 1);
                        if (block.second <= block.first)
                        {
                        }
                        else
                        {
                            // Create interface for canonical reverse dual basis
                            lattice.range().addRange(block);
                            typename LatticeHelper<RealTypeContext, IntTypeContext>::RangePopper rangepopper3(lattice);
                            STD_AUTO_PTR<GSInterface<RealTypeContext, IntTypeContext> > gs(workspace.getDualGSI(lattice.rc(), lattice.ic(), lattice));
                            Lattice<RealTypeContext, IntTypeContext> dl(gs.get(), lattice.rc(), lattice.ic(), lattice.getStats(), 0, gs->getDimension() - 1);
                            // Check
                            if (!isShortestVector(workspace, dl, 0, gs->getDimension() - 1, lattice.getLLLalpha()))
                            {
                                if (workspace.verbose().yieldsOutput(LatticeReduction::VL_Information))
                                    workspace.verbose()(LatticeReduction::VL_Information) << "Dual condition for block [" << block.first << "," << block.second << "] failed!\n";
                            }
                            else
                                is_ok = true;
                        }
                    }
                }
            }
            if (is_ok)
                // We found one ell maximizing the quotient which didn't fail
                return true;
        }
        // Only if all checks for maximal determinant quotient fails, fail total test.
        return false;
    }

    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    bool isHKZ(Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction> & workspace,
               Lattice<RealTypeContext, IntTypeContext> & lattice, bool dual = false)
    {
        // Compute GS coefficients
        lattice.update(lattice.range().end() + 1);
    
        if (dual)
        {
            // Create interface for canonical reverse dual basis
            STD_AUTO_PTR<GSInterface<RealTypeContext, IntTypeContext> > gs(workspace.getDualGSI(lattice.rc(), lattice.ic(), lattice));
            Lattice<RealTypeContext, IntTypeContext> dl(gs.get(), lattice.rc(), lattice.ic(), lattice.getStats(), 0, gs->getDimension() - 1);
            // Check
            return isHKZ(workspace, dl);
        }
        else
        {
            // Check for HKZ condition
            for (unsigned k = lattice.range().begin(); k <= lattice.range().end(); ++k)
            {
                if (!isShortestVector(workspace, lattice, k, lattice.range().end()))
                {
                    if (workspace.verbose().yieldsOutput(LatticeReduction::VL_Information))
                        workspace.verbose()(LatticeReduction::VL_Information) << "HKZ condition failed for [" << k << "," << lattice.range().end() << "]\n";
                    return false;
                }
            }
            return true;
        }
    }

    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    bool isSVP(Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction> & workspace,
               Lattice<RealTypeContext, IntTypeContext> & lattice, bool dual = false)
    {
        // Compute GS coefficients
        lattice.update(lattice.range().end() + 1);
    
        if (dual)
        {
            // Create interface for canonical reverse dual basis
            STD_AUTO_PTR<GSInterface<RealTypeContext, IntTypeContext> > gs(workspace.getDualGSI(lattice.rc(), lattice.ic(), lattice));
            Lattice<RealTypeContext, IntTypeContext> dl(gs.get(), lattice.rc(), lattice.ic(), lattice.getStats(), 0, gs->getDimension() - 1);
            // Check
            return isShortestVector(workspace, dl, dl.range().begin(), dl.range().end());
        }
        else
            // Check
            return isShortestVector(workspace, lattice, lattice.range().begin(), lattice.range().end());
    }

/*
  b_0,   ...    basis
  b_0^*, ...    GS-basis
  
  b_i^* = b_i - \sum_{j=0}^{i-1} \mu_{ij} b_j^*    with    \mu_{ij} = <b_i, b_j^*> / <b_j^*, b_j^*>
  
  pi_i = projection onto orthogonal complement of b_0, ..., b_{i-1}

  pi_i(v) = v - \sum_{j=0}^{i-1} <v, b_j^*> / <b_j^*, b_j^*> b_j^*

  hence   pi_i(b_k) = v - \sum_{j=0}^{i-1} \mu_{kj} b_j^*
                    = b_k^* + \sum_{j=0}^{k-1} \mu_{kj} b_j^* - \sum_{j=0}^{i-1} \mu_{kj} b_j^*
                    = b_k^* + \sum_{j=i}^{k-1} \mu_{kj} b_j^*
  i.e. w.r.t. the basis b_i^*, ... of the orthogonal complement of b_0, ..., b_{i-1}, the vector pi_i(b_k) is represented by
                    (mu_{k,i}, mu_{k,i+1}, mu_{k,i+2}, ..., mu_{k,k-1}, 1, 0, 0, ...),
                     with a 1 at position k-i (indexed from 0)

  therefore, the image of b = \sum_{j=i}^... lambda_j b_j under pi_i w.r.t. the basis b_i^*, ... is given by the vector
  whose k-th entry is
                     \sum_{j=i}^{i+k-1} \lambda_j \mu_{j,i+k} + \lambda_k,
  and the norm of pi_i(b) is therefore
                     \sum_{k=0}^... \left( \sum_{j=i}^{i+k-1} \lambda_j \mu_{j,i+k} + \lambda_k \right)^2 \|b_{i+k}\|^2
                   = \sum_{k=i}^... \left( \sum_{j=i}^{k-1} \lambda_{j-i} \mu_{j,k} + \lambda_{k-i} \right)^2 \|b_k\|^2
 */
}

#endif

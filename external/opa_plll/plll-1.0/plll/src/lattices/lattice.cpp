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

#ifndef PLLL_INCLUDE_GUARD__LATTICE_CPP
#define PLLL_INCLUDE_GUARD__LATTICE_CPP

#include <deque>
#include "transform.cpp"
#include "gramschmidt-generic.cpp"
#include <plll/linalg.hpp>

namespace plll
{
    /////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////
    //////  IMPLEMENTATION  /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////

    inline bool Ranger::empty() const
    {
        return d_ranges.empty();
    }

    inline unsigned Ranger::begin() const
    {
        return d_ranges.back().first;
    }

    inline unsigned Ranger::end() const
    {
        return d_ranges.back().second;
    }

    inline unsigned Ranger::dimension() const
    {
        return d_ranges.back().second - d_ranges.back().first + 1;
    }

    inline void Ranger::popRange()
    {
        d_ranges.pop_back();
    }

    inline Ranger & Ranger::addRange(unsigned b, unsigned e)
    {
        assert(b <= e);
        d_ranges.push_back(Range(b, e));
        return *this;
    }

    inline Ranger & Ranger::addRange(const std::pair<unsigned, unsigned> & b)
    {
        assert(b.first <= b.second);
        d_ranges.push_back(Range(b.first, b.second));
        return *this;
    }

    inline Ranger & Ranger::dupRange()
    {
        assert(!d_ranges.empty());
        d_ranges.push_back(d_ranges.back());
        return *this;
    }

    inline std::ostream & operator << (std::ostream & s, const Ranger & range)
    {
        return s << "[" << range.begin() << "," << range.end() << "]";
    }

    /////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////
    //////  IMPLEMENTATION  /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////
    
    template<class RealTypeContext, class IntTypeContext, class GSComputer>
    class DefaultGSI : public GSInterface<RealTypeContext, IntTypeContext>
    {
    private:
        struct Storage;
        
        STD_AUTO_PTR<Storage> d_storage;
        RealTypeContext & d_rc;
        IntTypeContext & d_ic;
        linalg::math_matrix<typename IntTypeContext::Integer> & d_A;
        unsigned d_used_vectors;
        GSComputer d_gs;
        typename RealTypeContext::Real d_zerodotfive, d_LLLalpha;
        
        struct Storage
        {
            linalg::math_matrix<typename IntTypeContext::Integer> A;
            
            Storage(const DefaultGSI<RealTypeContext, IntTypeContext, GSComputer> & src)
                : A(src.d_A)
            {
            }
        };
        
        void initZDFAlpha(double LLLalpha)
        {
            d_zerodotfive.setContext(d_rc);
            setOne(d_zerodotfive);
            d_zerodotfive >>= 1;
            d_LLLalpha.setContext(d_rc);
            arithmetic::convert(d_LLLalpha, LLLalpha, d_rc);
            d_gs.registerZDFAlpha(d_zerodotfive, d_LLLalpha);
            d_gs.adjustZerodotfive(d_zerodotfive, d_LLLalpha);
        }
        
        DefaultGSI(const DefaultGSI & source, RealTypeContext & rc, IntTypeContext & ic)
            : d_storage(new Storage(source)), d_rc(rc), d_ic(ic), d_A(d_storage->A),
              d_used_vectors(source.d_used_vectors),
              d_gs(source.d_gs, d_rc, d_ic, d_A),
              d_zerodotfive(d_rc), d_LLLalpha(d_rc)
        {
            d_zerodotfive = source.d_zerodotfive;
            d_LLLalpha = source.d_LLLalpha;
            d_gs.registerZDFAlpha(d_zerodotfive, d_LLLalpha);
        }
        
    public:
        DefaultGSI(Verbose & v, RealTypeContext & rc, IntTypeContext & ic, linalg::math_matrix<typename IntTypeContext::Integer> & A, double LLLalpha = 0.999)
            : d_storage(), d_rc(rc), d_ic(ic), d_A(A), d_used_vectors(d_A.rows()), d_gs(v, d_rc, d_ic, d_A)
        {
            initZDFAlpha(LLLalpha);
        }
        
        DefaultGSI(Verbose & v, RealTypeContext & rc, IntTypeContext & ic, linalg::math_matrix<typename IntTypeContext::Integer> & A, bool recompute, double LLLalpha = 0.999)
            : d_storage(), d_rc(rc), d_ic(ic), d_A(A), d_used_vectors(d_A.rows()), d_gs(v, d_rc, d_ic, d_A, recompute)
        {
            initZDFAlpha(LLLalpha);
        }
        
        virtual void adjustAlpha(double LLLalpha)
        // Sets a new value of LLLalpha and resets zerodotfive.
        {
            initZDFAlpha(LLLalpha);
        }
        
        virtual ~DefaultGSI()
        {
        }
        
        virtual std::pair<linalg::math_matrix<typename RealTypeContext::Real> *, linalg::math_rowvector<typename RealTypeContext::Real> *> getStorage()
        {
            return std::make_pair(d_gs.getCoeffsStorage(), d_gs.getSqNormsStorage());
        }
        
        virtual unsigned getDimension() const
        {
            return d_used_vectors;
        }
        
        virtual const typename RealTypeContext::Real & getZDF() const
        {
            return d_zerodotfive;
        }
        
        virtual const typename RealTypeContext::Real & getLLLalpha() const
        {
            return d_LLLalpha;
        }
        
        virtual void update(unsigned i)
        {
            d_gs.update(i);
        }
        
        virtual void reset()
        {
            d_gs.reset();
        }
        
        virtual bool isRowZero(unsigned i) const
        {
            for (unsigned j = 0; j < d_A.cols(); ++j)
                if (!isZero(d_A(i, j)))
                    return false;
            return true;
        }
    
        virtual void swap(unsigned i, unsigned j)
        {
            assert((i < d_used_vectors) && (j < d_used_vectors));
            if (i == j)
                // Avoid "zero" operations
                return;
            linalg::swap(d_A.row(i), d_A.row(j));
            d_gs.swap(i, j);
        }
        
        virtual void add(unsigned i, const typename IntTypeContext::Integer & m, unsigned j)
        {
            assert(i != j);
            assert((i < d_used_vectors) && (j < d_used_vectors));
            if (isZero(m))
                // Avoid "zero" operations
                return;
            if (isPMOne(m))
            {
                if (sign(m) > 0)
                    d_A.row(i) += d_A.row(j);
                else
                    d_A.row(i) -= d_A.row(j);
            }
            else
                plll::linalg::addmul(d_A.row(i), d_A.row(j), m); // d_A.row(i) += d_A.row(j) * m;
            d_gs.add(i, m, j);
        }
        
        virtual void flip(unsigned i)
        {
            assert(i < d_used_vectors);
            plll::linalg::neg(d_A.row(i), d_A.row(i));
            d_gs.flip(i);
        }
        
        virtual void trans(unsigned i, unsigned j,
                           const typename IntTypeContext::Integer & B00, const typename IntTypeContext::Integer & B01,
                           const typename IntTypeContext::Integer & B10, const typename IntTypeContext::Integer & B11)
        // Assuming that the determinant is 1:
        // [ B00 B01 ]   [ ..... < row i > ..... ]
        // [ B10 B11 ] * [ ..... < row j > ..... ]
        {
            assert(i != j);
            assert((i < d_used_vectors) && (j < d_used_vectors));
            assert(isOne(B00 * B11 - B01 * B10));
            typename IntTypeContext::Integer t0, t1;
            for (unsigned k = 0; k < d_A.cols(); ++k)
            {
                t0 = B00 * d_A(i, k);
                t0 += B01 * d_A(j, k);
                t1 = B10 * d_A(i, k);
                t1 += B11 * d_A(j, k);
                d_A(i, k) = t0;
                d_A(j, k) = t1;
            }
            d_gs.trans(i, j, B00, B01, B10, B11);
        }
        
        virtual void insertVectorLC(unsigned ofs, const linalg::math_rowvector<typename IntTypeContext::Integer> & result)
        // Inserts a new vector at index ofs, which will equal \sum_{i=0}^{result.size()-1} result[i] * b[ofs+i]
        {
            assert(ofs <= d_used_vectors);
            bool add = d_used_vectors == d_A.rows();
            assert(ofs + result.size() <= d_used_vectors);
            linalg::math_rowvector<typename IntTypeContext::Integer> v(d_A.cols());
            for (unsigned i = 0; i < result.size(); ++i)
                if (!isZero(result[i]))
                    plll::linalg::addmul(v, d_A.row(ofs + i), result[i]); // v += result[i] * d_A.row(ofs + i);
            if (add)
                d_A.resize(d_A.rows() + 1, d_A.cols());
            for (unsigned i = d_used_vectors; i > ofs; --i)
                linalg::swap(d_A.row(i), d_A.row(i - 1));
            d_A.row(ofs) = v;
            ++d_used_vectors;
            d_gs.insertVector(ofs, result);
        }
        
        virtual void insertVector(unsigned ofs, const linalg::math_rowvector<typename IntTypeContext::Integer> & result)
        // Inserts a new vector <result> at index ofs
        {
            assert(ofs <= d_used_vectors);
            bool add = d_used_vectors == d_A.rows();
            assert(result.size() == d_A.cols());
            if (add)
                d_A.resize(d_A.rows() + 1, d_A.cols());
            for (unsigned i = d_used_vectors; i > ofs; --i)
                linalg::swap(d_A.row(i), d_A.row(i - 1));
            d_A.row(ofs) = result;
            ++d_used_vectors;
            d_gs.insertVector(ofs, result);
        }
        
        virtual void removeZeroVector(unsigned ofs)
        // Removes the zero vector at position ofs
        {
            assert(ofs < d_used_vectors);
            for (unsigned i = ofs + 1; i < d_used_vectors; ++i)
                linalg::swap(d_A.row(i), d_A.row(i - 1));
            --d_used_vectors;
            d_gs.removeZeroVector(ofs);
        }
        
        virtual void compactify()
        {
            if (d_used_vectors < d_A.rows())
                d_A.resize(d_used_vectors, d_A.cols());
            d_gs.compactify();
        }
        
        virtual bool canGetIntegerVector() const
        {
            return true;
        }
        
        virtual void getLinearCombination(linalg::math_rowvector<typename IntTypeContext::Integer> & result, unsigned ofs, const linalg::math_rowvector<typename IntTypeContext::Integer> & lincomb) const
        {
            result.resize(d_A.cols());
            for (unsigned i = 0; i < result.size(); ++i)
                setZero(result[i]);
            for (unsigned i = 0; i < lincomb.size(); ++i)
                if (!isZero(lincomb[i]))
                {
                    if (isPMOne(lincomb[i]))
                    {
                        if (sign(lincomb[i]) > 0)
                            result += d_A.row(ofs + i);
                        else
                            result -= d_A.row(ofs + i);
                    }
                    else
                        result += lincomb[i] * d_A.row(ofs + i);
                }
        }
        
        virtual void getVector(linalg::math_rowvector<typename IntTypeContext::Integer> & result, unsigned ofs) const
        {
            result = d_A.row(ofs);
        }
        
        virtual std::pair<linalg::math_matrix<typename IntTypeContext::Integer> *, unsigned> getMatrixIndex(unsigned ofs) const
        // Returns matrix (if available) and index of specified vector <ofs> (might be different from
        // <ofs>).
        {
            return std::make_pair(&d_A, ofs);
        }
        
        virtual const typename IntTypeContext::Integer & getNormSqUP(unsigned i) const
        {
            return d_gs.getNormSqUP(i);
        }
        
        virtual void getNormSqUP(typename IntTypeContext::Integer & r, unsigned i) const
        {
            d_gs.getNormSqUP(r, i);
        }
        
        virtual void getDotProduct(typename IntTypeContext::Integer & r, unsigned i, unsigned j) const
        {
            d_gs.getDotProduct(r, i, j);
        }
        
        virtual void computeDotProductProjected(typename RealTypeContext::Real & r, unsigned k, unsigned i, unsigned j) const
        {
            d_gs.computeDotProductProjected(r, k, i, j);
        }
        
        virtual void computeProjectionLength(typename RealTypeContext::Real & result, unsigned k, unsigned b, const linalg::math_rowvector<typename IntTypeContext::Integer> & vec) const
        // Computes the squared length of the projection of \sum_{i=0}^{m-1} vec[i] B_{b+i} onto the
        // orthogonal complement of B_0, ..., B_{k-1}. (Here, B_0, ..., B_{n-1} is the basis of the
        // lattice.)
        {
            d_gs.computeProjectionLength(result, k, b, vec);
        }
        
        virtual void computeProjectionLengthBV(typename RealTypeContext::Real & result, unsigned k, unsigned b) const
        // Computes the squared length of the projection of B_b onto the orthogonal complement of B_0, ...,
        // B_{k-1}. (Here, B_0, ..., B_{n-1} is the basis of the lattice.)
        {
            d_gs.computeProjectionLengthBV(result, k, b);
        }
        
        virtual void computeProjectionLength(typename RealTypeContext::Real & result, unsigned k, const linalg::math_rowvector<typename IntTypeContext::Integer> & vec) const
        {
            d_gs.computeProjectionLength(result, k, vec);
        }
        
        virtual bool sizereduce(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned i)
        {
            if (i == 0)
                return false;
            return d_gs.sizereduce(lattice, d_zerodotfive, d_LLLalpha, 0, i - 1, i);
        }
        
        virtual GSInterface<RealTypeContext, IntTypeContext> * clone(RealTypeContext & rc, IntTypeContext & ic)
        {
            return new DefaultGSI<RealTypeContext, IntTypeContext, GSComputer>(*this, rc, ic);
        }
        
        virtual void copyFrom(GSInterface<RealTypeContext, IntTypeContext> * source_)
        {
            DefaultGSI * source = dynamic_cast<DefaultGSI<RealTypeContext, IntTypeContext, GSComputer>*>(source_);
            assert(source != NULL);
            d_rc = source->d_rc;
            d_ic = source->d_ic;
            d_A = source->d_A;
            d_used_vectors = source->d_used_vectors;
            d_gs.copyFrom(source->d_gs);
            d_zerodotfive.setContext(d_rc);
            d_zerodotfive = source->d_zerodotfive;
            d_LLLalpha.setContext(d_rc);
            d_LLLalpha = source->d_LLLalpha;
        }
    
        virtual void changeOfPrecision()
        {
            d_gs.changeOfPrecision();
        }
    
        virtual bool avoidInsertions()
        {
            return false;
        }
    
        virtual void print(std::ostream & s)
        {
            s << "Basis:\n";
            s << d_A << "\n";
        }
    };

    /////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////
    //////  IMPLEMENTATION  /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////

    template<class RealTypeContext, class IntTypeContext>
    class GSContainer : public GSInterface<RealTypeContext, IntTypeContext>
    /*
      Keeps track of a GS decomposition
    */
    {
    private:
        linalg::math_matrix<typename RealTypeContext::Real> d_coeffs;
        linalg::math_rowvector<typename RealTypeContext::Real> d_sqnorms;
        linalg::base_rowvector<bool> d_dependent;
        unsigned d_dimension;
        RealTypeContext & d_rc;
        IntTypeContext & d_ic;
        typename RealTypeContext::Real d_zerodotfive, d_LLLalpha;
        typename IntTypeContext::Integer d_tmp;
    
        GSContainer(const GSContainer & source, RealTypeContext & rc, IntTypeContext & ic)
            : d_coeffs(source.d_coeffs), d_sqnorms(source.d_sqnorms), d_dependent(source.d_dependent),
              d_dimension(source.d_dimension), d_rc(rc), d_ic(ic), d_zerodotfive(d_rc), d_LLLalpha(d_rc)
        {
            d_zerodotfive = source.d_zerodotfive;
            d_LLLalpha = source.d_LLLalpha;
            setZero(d_tmp);
        }
    
    public:
        GSContainer(RealTypeContext & rc, IntTypeContext & ic, unsigned dimension = 0, double alpha = 0.999)
            : d_dimension(dimension), d_rc(rc), d_ic(ic), d_zerodotfive(d_rc), d_LLLalpha(d_rc)
        {
            setOne(d_zerodotfive);
            d_zerodotfive >>= 1;
            arithmetic::convert(d_LLLalpha, alpha, d_rc);
            d_coeffs.resize(d_dimension, d_dimension, linalg::Initialize(d_rc));
            d_sqnorms.resize(d_dimension, linalg::Initialize(d_rc));
            d_dependent.resize(d_dimension);
            for (unsigned i = 0; i < d_dimension; ++i)
                d_dependent[i] = false;
        }
    
        virtual ~GSContainer()
        {
        }
    
        virtual std::pair<linalg::math_matrix<typename RealTypeContext::Real> *, linalg::math_rowvector<typename RealTypeContext::Real> *> getStorage()
        {
            return std::make_pair(&d_coeffs, &d_sqnorms);
        }
    
        void updateDependent()
        {
            for (unsigned i = 0; i < d_dimension; ++i)
                d_dependent[i] = isZero(d_sqnorms[i]);
        }
    
        virtual void update(unsigned)
        {
        }
    
        virtual void reset()
        {
        }
    
        virtual void adjustAlpha(double LLLalpha)
        {
            arithmetic::convert(d_LLLalpha, LLLalpha, d_rc);
        }
    
        virtual unsigned getDimension() const
        {
            return d_dimension;
        }
    
        virtual bool isRowZero(unsigned i) const
        {
            for (unsigned j = 0; j < i; ++j)
                if (!isZero(d_coeffs(i, j)))
                    return false;
            return isZero(d_sqnorms[i]);
        }
    
        virtual const typename RealTypeContext::Real & getZDF() const
        {
            return d_zerodotfive;
        }
    
        virtual const typename RealTypeContext::Real & getLLLalpha() const
        {
            return d_LLLalpha;
        }
    
        virtual void swap(unsigned i, unsigned j)
        {
            if (i < j)
                // Make sure that j <= i.
                std::swap(i, j);
            if (i == j)
                // Do nothing.
                return;
            if (i - j == 1)
            {
                // Swap two adjacent vectors
                GSGeneric::swap(d_rc, d_ic, d_coeffs, d_sqnorms, d_dependent, d_dimension, i, j);
            }
            else
            {
                // First swap i to j+1
                for (unsigned k = i - 1; k > j; --k)
                    GSGeneric::swap(d_rc, d_ic, d_coeffs, d_sqnorms, d_dependent, d_dimension, k + 1, k);
            
                // Then swap j and j+1
                GSGeneric::swap(d_rc, d_ic, d_coeffs, d_sqnorms, d_dependent, d_dimension, j + 1, j);

                // Finally swap j+1 back to i
                for (unsigned k = j + 1; k < i; ++k)
                    GSGeneric::swap(d_rc, d_ic, d_coeffs, d_sqnorms, d_dependent, d_dimension, k + 1, k);
            }
        }
    
        virtual void add(unsigned i, const typename IntTypeContext::Integer & m, unsigned j)
        {
            GSGeneric::add(d_rc, d_ic, d_coeffs, d_sqnorms, d_dependent, d_dimension, i, m, j);
        }
    
        virtual void flip(unsigned i)
        {
            GSGeneric::flip(d_rc, d_ic, d_coeffs, d_sqnorms, d_dependent, d_dimension, i);
        }
    
        virtual void trans(unsigned i, unsigned j,
                           const typename IntTypeContext::Integer & B00, const typename IntTypeContext::Integer & B01,
                           const typename IntTypeContext::Integer & B10, const typename IntTypeContext::Integer & B11)
        {
            if (i + 1 == j)
            {
                // Use the generic method
                GSGeneric::trans(d_rc, d_ic, d_coeffs, d_sqnorms, d_dependent, d_dimension, i, j, B00, B01, B10, B11);
            }
            else if (i > j)
            {
                // Wrong order: change it
                trans(j, i, B01, B00, B11, B10);
            }
            else
            {
                // Correct order (i < j), but wrong distance! Solve by swapping j to i+1, apply
                // transformation, and swap back.
            
                // Swap j to i+1
                for (unsigned k = j - 1; k > i; --k)
                    GSGeneric::swap(d_rc, d_ic, d_coeffs, d_sqnorms, d_dependent, d_dimension, k + 1, k);
            
                // Apply transformation to i and i+1
                GSGeneric::trans(d_rc, d_ic, d_coeffs, d_sqnorms, d_dependent, d_dimension, i, i + 1, B00, B01, B10, B11);
            
                // Swap i+1 to j
                for (unsigned k = i + 1; k < j; ++k)
                    GSGeneric::swap(d_rc, d_ic, d_coeffs, d_sqnorms, d_dependent, d_dimension, k + 1, k);
            }
        }
    
        virtual void insertVectorLC(unsigned ofs, const linalg::math_rowvector<typename IntTypeContext::Integer> & result)
        {
            // We "cheat": first, we insert the vector *behind* the last vector of the linear
            // combination. Then, we add the previous vectors (as usual). Finally, we move the vector to
            // the right position using swaps.
            GSGeneric::insertVector(d_rc, d_ic, d_coeffs, d_sqnorms, d_dependent, d_dimension, d_dimension, ofs + result.size());
            for (unsigned i = 0; i < result.size(); ++i)
                if (!isZero(result[i]))
                    GSGeneric::add(d_rc, d_ic, d_coeffs, d_sqnorms, d_dependent, d_dimension, ofs + result.size(), result[i], ofs + i);
            for (unsigned i = ofs + result.size(); i > ofs; --i)
                GSGeneric::swap(d_rc, d_ic, d_coeffs, d_sqnorms, d_dependent, d_dimension, i - 1, i);
        }
    
        virtual void insertVector(unsigned ofs, const linalg::math_rowvector<typename IntTypeContext::Integer> & result)
        {
            // This operation is not supported. (We don't support integer vectors.)
        }
    
        virtual void removeZeroVector(unsigned ofs)
        {
            GSGeneric::removeVector(d_rc, d_ic, d_coeffs, d_sqnorms, d_dependent, d_dimension, d_dimension, ofs);
        }
    
        virtual void compactify()
        {
            // GSGeneric::compactify(d_rc, d_ic, d_coeffs, d_sqnorms, d_dependent, d_dimension);
        }
    
        virtual bool canGetIntegerVector() const
        {
            return false;
        }
    
        virtual void getLinearCombination(linalg::math_rowvector<typename IntTypeContext::Integer> & result, unsigned ofs, const linalg::math_rowvector<typename IntTypeContext::Integer> & lincomb) const
        {
            // Do nothing. (We don't have integer vectors.)
        }
    
        virtual void getVector(linalg::math_rowvector<typename IntTypeContext::Integer> & result, unsigned ofs) const
        {
            // Do nothing. (We don't have integer vectors.)
        }
    
        virtual std::pair<linalg::math_matrix<typename IntTypeContext::Integer> *, unsigned> getMatrixIndex(unsigned ofs) const
        {
            // Do nothing. (We don't have integer vectors.)
            return std::make_pair<linalg::math_matrix<typename IntTypeContext::Integer> *, unsigned>(NULL, 0);
        }
    
        virtual const typename IntTypeContext::Integer & getNormSqUP(unsigned i) const
        {
            // Do nothing. (We don't have integer vectors.)
            return d_tmp;
        }
    
        virtual void getNormSqUP(typename IntTypeContext::Integer & r, unsigned i) const
        {
            // Do nothing. (We don't have integer vectors.)
        }
    
        virtual void getDotProduct(typename IntTypeContext::Integer & r, unsigned i, unsigned j) const
        {
            // Do nothing. (We don't have integer vectors.)
        }
    
        virtual void computeDotProductProjected(typename RealTypeContext::Real & result, unsigned k, unsigned i, unsigned j) const
        {
            GSGeneric::computeDotProductProjected(d_rc, d_ic, d_coeffs, d_sqnorms, d_dependent, d_dimension, result, k, i, j);
        }
    
        virtual void computeProjectionLength(typename RealTypeContext::Real & result, unsigned k, unsigned b, const linalg::math_rowvector<typename IntTypeContext::Integer> & vec) const
        // Computes the squared length of the projection of \sum_{i=0}^{m-1} vec[i] B_{b+i} onto the
        // orthogonal complement of B_0, ..., B_{k-1}. (Here, B_0, ..., B_{n-1} is the basis of the
        // lattice.)
        {
            GSGeneric::computeProjectionLength(d_rc, d_ic, d_coeffs, d_sqnorms, d_dependent, d_dimension, result, k, b, vec);
        }
    
        virtual void computeProjectionLengthBV(typename RealTypeContext::Real & result, unsigned k, unsigned b) const
        // Computes the squared length of the projection of B_b onto the orthogonal complement of B_0, ...,
        // B_{k-1}. (Here, B_0, ..., B_{n-1} is the basis of the lattice.)
        {
            GSGeneric::computeProjectionLengthBV(d_rc, d_ic, d_coeffs, d_sqnorms, d_dependent, d_dimension, result, k, b);
        }
    
        virtual void computeProjectionLength(typename RealTypeContext::Real & result, unsigned k, const linalg::math_rowvector<typename IntTypeContext::Integer> & vec) const
        {
            // This operation is not supported. (We don't support integer vectors.)
        }
    
        virtual bool sizereduce(Lattice<RealTypeContext, IntTypeContext> & lattice, unsigned i)
        {
            if (i == 0)
                return false;
            return GSGeneric::sizereduce(d_rc, d_ic, d_coeffs, d_sqnorms, d_dependent, d_dimension, lattice, d_zerodotfive, 0, i - 1, i);
        }
    
        virtual GSInterface<RealTypeContext, IntTypeContext> * clone(RealTypeContext & rc, IntTypeContext & ic)
        {
            return new GSContainer(*this, rc, ic);
        }
    
        virtual void copyFrom(GSInterface<RealTypeContext, IntTypeContext> * source_)
        {
            GSContainer * source = dynamic_cast<GSContainer<RealTypeContext, IntTypeContext>*>(source_);
            assert(source != NULL);
            d_coeffs = source->d_coeffs;
            d_sqnorms = source->d_sqnorms;
            d_dependent = source->d_dependent;
            d_dimension = source->d_dimension;
            d_zerodotfive.setContext(d_rc);
            d_zerodotfive = source->d_zerodotfive;
            d_LLLalpha.setContext(d_rc);
            d_LLLalpha = source->d_LLLalpha;
        }
    
        virtual void changeOfPrecision()
        // inform that RealTypeContext's precision changed; has to update the matrix/vector set by setStor}
        {
            for (unsigned i = 0; i < d_dimension; ++i)
            {
                for (unsigned j = 0; j < d_dimension; ++j)
                    d_coeffs(i, j).setContext(d_rc);
                d_sqnorms[i].setContext(d_rc);
            }
        }
    
        virtual bool avoidInsertions()
        {
            return true;
        }
    
        virtual void print(std::ostream &)
        {
        }
    };

    /////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////
    //////  IMPLEMENTATION  /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////

    template<class RealTypeContext, class IntTypeContext, class Enumerator, class CallbackFunction>
    GSInterface<RealTypeContext, IntTypeContext> * Workspace<RealTypeContext, IntTypeContext, Enumerator, CallbackFunction>::
    getDualGSI(RealTypeContext & rc, IntTypeContext & ic, const Lattice<RealTypeContext, IntTypeContext> & lattice)
    {
        unsigned begin = lattice.range().begin(), dim = lattice.range().dimension();
        GSContainer<RealTypeContext, IntTypeContext> * gs =
            new GSContainer<RealTypeContext, IntTypeContext>(rc, ic, dim, arithmetic::convert<double>(lattice.getLLLalpha()));
        linalg::math_matrix<typename RealTypeContext::Real> * coeffs = gs->getStorage().first;
        linalg::math_rowvector<typename RealTypeContext::Real> * sqnorms = gs->getStorage().second;
//    // Make sure GS decomposition is available.  -- not possible since the lattice is const
//    lattice.update(lattice.range().end() + 1);
        // Compute dual basis. First, invert diagonal.
        typename RealTypeContext::Real t(rc);
        setOne(t);
        for (unsigned i = 0; i < dim; ++i)
            (*sqnorms)[i] = t / lattice.getNormSq(begin + dim - 1 - i);
        // Then invert rest of the matrix.
        for (unsigned i = 0; i < dim; ++i)
        {
            for (unsigned j = 0; j < i; ++j)
            {
                // Let t be the (i,j) entry of the inverse of *coeffs (with added 1-diagonal). Then t
                // should be stored at position (dim-j-1,dim-i-1).
                // 
                // If the inverse is stored in B with B = (b_{ij})_{ij}, and *coeffs is A with A =
                // (a_{ij})_{ij}, then
                //           b_{ij} = -a_{ij} - \sum_{k=j+1}^{i-1} a_{ik} b_{kj}.
                // 
                // This has to be computed for increasing i so that the needed values for the inverse
                // are already computed.
                t = lattice.getCoeff(begin + i, begin + j); // We add the negation later; first just add values.
                for (unsigned k = j + 1; k < i; ++k)
                    t += lattice.getCoeff(begin + i, begin + k) * (*coeffs)(dim - j - 1, dim - k - 1);
                // Store at correct position.
                (*coeffs)(dim - j - 1, dim - i - 1) = -t;
            }
        }
        // Now update dependent flags (all should be 'false') and we're done!
        gs->updateDependent();
        return gs;
    }

    /////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////
    //////  IMPLEMENTATION  /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////

    template<class RealTypeContext, class IntTypeContext>
    Lattice<RealTypeContext, IntTypeContext>::Lattice(GSInterface<RealTypeContext, IntTypeContext> * gs, RealTypeContext & rc, IntTypeContext & ic,
                                                      const Lattice & lattice, bool own_stats, LatticeReduction::Statistics & stats)
        : d_rc(rc), d_ic(ic), d_gs(gs), d_stats(own_stats ? stats : lattice.d_stats), d_do_not_modify_flag(false)
    {
        d_range = lattice.d_range;
        std::pair<linalg::math_matrix<typename RealTypeContext::Real> *, linalg::math_rowvector<typename RealTypeContext::Real> *> r = d_gs->getStorage();
        d_coeffs = r.first;
        d_sqnorms = r.second;
        // Do *NOT* copy notifiers! !!! ??? ...
    }

    template<class RealTypeContext, class IntTypeContext>
    void Lattice<RealTypeContext, IntTypeContext>::copyFrom(const DuplicateStorage & storage, bool add_own_stats)
    {
        if (add_own_stats)
            addStatistics(d_stats, storage.d_stats);
        d_rc = storage.d_lattice.d_rc;
        d_ic = storage.d_lattice.d_ic;
        d_gs->copyFrom(storage.d_lattice.d_gs);
        d_range = storage.d_lattice.d_range;
        // Do *NOT* copy notifiers! !!! ??? ...
    }

    template<class RealTypeContext, class IntTypeContext>
    Lattice<RealTypeContext, IntTypeContext>::Lattice(GSInterface<RealTypeContext, IntTypeContext> * gs, RealTypeContext & rc, IntTypeContext & ic,
                                                      LatticeReduction::Statistics & stats, unsigned begin, unsigned end, bool do_not_modify_flag)
        : d_rc(rc), d_ic(ic), d_coeffs(NULL), d_sqnorms(NULL), d_gs(gs), d_stats(stats), d_do_not_modify_flag(do_not_modify_flag)
    {
        assert(begin <= end);
        assert(end < gs->getDimension());
        d_range.addRange(begin, end);
        std::pair<linalg::math_matrix<typename RealTypeContext::Real> *, linalg::math_rowvector<typename RealTypeContext::Real> *> r = d_gs->getStorage();
        d_coeffs = r.first;
        d_sqnorms = r.second;
    }

    template<class RealTypeContext, class IntTypeContext>
    Lattice<RealTypeContext, IntTypeContext>::~Lattice()
    {
        if (!d_range.empty())
            d_range.popRange();
//    assert(d_range.empty()); -- In certain situations, such as when Lattice<> is cloned, it can
//    happen that the range stack is not empty! So don't check it via an assertion.
    }

    template<class RealTypeContext, class IntTypeContext>
    inline LatticeReduction::Statistics & Lattice<RealTypeContext, IntTypeContext>::getStatistics()
    {
        return d_stats;
    }

    template<class RealTypeContext, class IntTypeContext>
    inline Ranger & Lattice<RealTypeContext, IntTypeContext>::range()
    {
        return d_range;
    }

    template<class RealTypeContext, class IntTypeContext>
    inline const Ranger & Lattice<RealTypeContext, IntTypeContext>::range() const
    {
        return d_range;
    }

    template<class RealTypeContext, class IntTypeContext>
    inline unsigned Lattice<RealTypeContext, IntTypeContext>::dimension() const
    {
        return d_gs->getDimension();
    }

    template<class RealTypeContext, class IntTypeContext>
    inline GSInterface<RealTypeContext, IntTypeContext> * Lattice<RealTypeContext, IntTypeContext>::gs() const
    {
        return d_gs;
    }

    template<class RealTypeContext, class IntTypeContext>
    inline RealTypeContext & Lattice<RealTypeContext, IntTypeContext>::rc()
    {
        return d_rc;
    }

    template<class RealTypeContext, class IntTypeContext>
    inline const RealTypeContext & Lattice<RealTypeContext, IntTypeContext>::rc() const
    {
        return d_rc;
    }

    template<class RealTypeContext, class IntTypeContext>
    inline IntTypeContext & Lattice<RealTypeContext, IntTypeContext>::ic()
    {
        return d_ic;
    }

    template<class RealTypeContext, class IntTypeContext>
    inline const IntTypeContext & Lattice<RealTypeContext, IntTypeContext>::ic() const
    {
        return d_ic;
    }

    template<class RealTypeContext, class IntTypeContext>
    inline const typename RealTypeContext::Real & Lattice<RealTypeContext, IntTypeContext>::getZDF() const
    {
        return d_gs->getZDF();
    }

    template<class RealTypeContext, class IntTypeContext>
    inline const typename RealTypeContext::Real & Lattice<RealTypeContext, IntTypeContext>::getLLLalpha() const
    {
        return d_gs->getLLLalpha();
    }

    template<class RealTypeContext, class IntTypeContext>
    inline LatticeReduction::Statistics & Lattice<RealTypeContext, IntTypeContext>::getStats()
    {
        return d_stats;
    }

    template<class RealTypeContext, class IntTypeContext>
    inline bool Lattice<RealTypeContext, IntTypeContext>::avoidInsertions() const
    // Returns true if insertions/removals of vectors are
    {
        return d_gs->avoidInsertions();
    }

    template<class RealTypeContext, class IntTypeContext>
    inline bool Lattice<RealTypeContext, IntTypeContext>::allowedToModify() const
    {
        return !d_do_not_modify_flag;
    }

    template<class RealTypeContext, class IntTypeContext>
    inline void Lattice<RealTypeContext, IntTypeContext>::update(unsigned i)
    {
        d_gs->update(i);
    }

    template<class RealTypeContext, class IntTypeContext>
    inline void Lattice<RealTypeContext, IntTypeContext>::reset()
    {
        d_gs->reset();
    }

    template<class RealTypeContext, class IntTypeContext>
    inline const typename RealTypeContext::Real & Lattice<RealTypeContext, IntTypeContext>::getCoeff(unsigned i, unsigned j) const
    {
        assert(j < i);
        return (*d_coeffs)(i, j);
    }

    template<class RealTypeContext, class IntTypeContext>
    inline const typename RealTypeContext::Real & Lattice<RealTypeContext, IntTypeContext>::getNormSq(unsigned i) const
    {
        return (*d_sqnorms)[i];
    }

    template<class RealTypeContext, class IntTypeContext>
    void Lattice<RealTypeContext, IntTypeContext>::addNotifier(TransformNotifier<IntTypeContext> * t)
    {
        if (t)
            d_notifier.push_back(t);
    }

    template<class RealTypeContext, class IntTypeContext>
    void Lattice<RealTypeContext, IntTypeContext>::removeNotifier(TransformNotifier<IntTypeContext> * tt)
    {
        typename std::deque<TransformNotifier<IntTypeContext>*>::reverse_iterator t = d_notifier.rbegin(); // start searching at the end
        while (t != d_notifier.rend())
        {
            if (*t == tt)
            {
                d_notifier.erase(t.base());
                return;
            }
            ++t;
        }
        assert(!"Notifier not found!");
    }
    
    template<class RealTypeContext, class IntTypeContext>
    void Lattice<RealTypeContext, IntTypeContext>::removeAllNotifiers()
    {
        d_notifier.clear();
    }

    template<class RealTypeContext, class IntTypeContext>
    inline bool Lattice<RealTypeContext, IntTypeContext>::isRowZero(unsigned i) const
    {
        return d_gs->isRowZero(i);
    }

    template<class RealTypeContext, class IntTypeContext>
    void Lattice<RealTypeContext, IntTypeContext>::swap(unsigned i, unsigned j)
    {
        ++d_stats.swaps;
        d_gs->swap(i, j);
        for (typename std::deque<TransformNotifier<IntTypeContext>*>::iterator t = d_notifier.begin(); t != d_notifier.end(); ++t)
            (*t)->swap(i, j);
    }

    template<class RealTypeContext, class IntTypeContext>
    void Lattice<RealTypeContext, IntTypeContext>::add(unsigned i, const typename IntTypeContext::Integer & m, unsigned j)
    {
        ++d_stats.adds;
        if (isPMOne(m))
            ++d_stats.adds_pm1;
        if (isPMTwo(m))
            ++d_stats.adds_pm2;
        d_gs->add(i, m, j);
        for (typename std::deque<TransformNotifier<IntTypeContext>*>::iterator t = d_notifier.begin(); t != d_notifier.end(); ++t)
            (*t)->add(i, m, j);
    }

    template<class RealTypeContext, class IntTypeContext>
    void Lattice<RealTypeContext, IntTypeContext>::flip(unsigned i)
    {
        ++d_stats.flips;
        d_gs->flip(i);
        for (typename std::deque<TransformNotifier<IntTypeContext>*>::iterator t = d_notifier.begin(); t != d_notifier.end(); ++t)
            (*t)->flip(i);
    }

    template<class RealTypeContext, class IntTypeContext>
    void Lattice<RealTypeContext, IntTypeContext>::trans(unsigned i, unsigned j,
                                                         const typename IntTypeContext::Integer & B00, const typename IntTypeContext::Integer & B01,
                                                         const typename IntTypeContext::Integer & B10, const typename IntTypeContext::Integer & B11)
    // Assuming that the determinant is 1:
    // [ B00 B01 ]   [ ..... < row i > ..... ]
    // [ B10 B11 ] * [ ..... < row j > ..... ]
    {
        ++d_stats.trans;
        d_gs->trans(i, j, B00, B01, B10, B11);
        for (typename std::deque<TransformNotifier<IntTypeContext>*>::iterator t = d_notifier.begin(); t != d_notifier.end(); ++t)
            (*t)->trans(i, j, B00, B01, B10, B11);
    }

    template<class RealTypeContext, class IntTypeContext>
    void Lattice<RealTypeContext, IntTypeContext>::insertVectorLC(unsigned ofs, const linalg::math_rowvector<typename IntTypeContext::Integer> & result)
    // Inserts a new vector at index ofs, which will equal \sum_{i=0}^{result.size()-1} result[i] * b[ofs+i]
    {
        /*
          if (canGetIntegerVector()) // for debug reasons: use vector insertion to see if it works
          {
          linalg::math_rowvector<Integer> v;
          std::cout << "LC:     " << ofs << " " << result << "\n";
          getLinearCombination(v, ofs, result);
          insertVector(ofs, v);
          return;
          }
        */
    
        ++d_stats.vectorinsertions;
        d_gs->insertVectorLC(ofs, result);
        for (typename std::deque<TransformNotifier<IntTypeContext>*>::iterator t = d_notifier.begin(); t != d_notifier.end(); ++t)
            (*t)->insertVectorLC(ofs, result);
        d_range.insertVector(ofs);
    }

    template<class RealTypeContext, class IntTypeContext>
    void Lattice<RealTypeContext, IntTypeContext>::insertVector(unsigned ofs, const linalg::math_rowvector<typename IntTypeContext::Integer> & result)
    // Inserts a new vector <result> at index ofs
    {
        ++d_stats.vectorinsertions;
        // Do we need linear combination?
        bool need_combination = false;
        for (typename std::deque<TransformNotifier<IntTypeContext>*>::iterator t = d_notifier.begin(); t != d_notifier.end(); ++t)
            if (!(*t)->canInsertVector())
            {
                need_combination = true;
                break;
            }
        linalg::math_rowvector<typename IntTypeContext::Integer> resultLC;
        if (need_combination)
        {
            // Set up linear system
            typename std::pair<linalg::math_matrix<typename IntTypeContext::Integer> *, unsigned> mat = getMatrixIndex(ofs);
            assert(mat.first != NULL);
            unsigned dim = d_gs->getDimension();
            linalg::math_matrix<typename IntTypeContext::Integer> A;
            A.resize(mat.first->cols(), dim - ofs);
            for (unsigned r = 0; r < dim - ofs; ++r)
                A.col(r) = mat.first->row(mat.second + r);
            // Solve linear system to recover linear combination
            resultLC = solveInt(A, result.transpose()).transpose(); // we have A*resultLC==result in case resultLC.size() > 0
            // See if we can shorten the linear combination
            unsigned len = resultLC.size();
            while (len > 0)
                if (isZero(resultLC[len - 1]))
                    --len;
                else
                    break;
            resultLC.resize(len);
            // It should be non-trivial!
            assert(resultLC.size() > 0);
        }
    
        // Insert vector
        d_gs->insertVector(ofs, result);
        for (typename std::deque<TransformNotifier<IntTypeContext>*>::iterator t = d_notifier.begin(); t != d_notifier.end(); ++t)
            if ((*t)->canInsertVector())
                (*t)->insertVector(ofs, result);
            else
                (*t)->insertVectorLC(ofs, resultLC); // need linear combination!
        d_range.insertVector(ofs);
    }

    template<class RealTypeContext, class IntTypeContext>
    void Lattice<RealTypeContext, IntTypeContext>::removeZeroVector(unsigned ofs)
    // Removes the zero vector at position ofs
    {
        d_gs->removeZeroVector(ofs);
        for (typename std::deque<TransformNotifier<IntTypeContext>*>::iterator t = d_notifier.begin(); t != d_notifier.end(); ++t)
            (*t)->removeZeroVector(ofs);
        d_range.removeVector(ofs);
    }

    template<class RealTypeContext, class IntTypeContext>
    void Lattice<RealTypeContext, IntTypeContext>::compactify()
    {
        d_gs->compactify();
        for (typename std::deque<TransformNotifier<IntTypeContext>*>::iterator t = d_notifier.begin(); t != d_notifier.end(); ++t)
            (*t)->compactify();
    }

    template<class RealTypeContext, class IntTypeContext>
    inline bool Lattice<RealTypeContext, IntTypeContext>::sizereduce(unsigned i)
    {
        ++d_stats.sizereductions;
        return d_gs->sizereduce(*this, i);
    }

    template<class RealTypeContext, class IntTypeContext>
    inline void Lattice<RealTypeContext, IntTypeContext>::computeDotProductProjected(typename RealTypeContext::Real & result, unsigned k, unsigned i, unsigned j) const
    {
        d_gs->computeDotProductProjected(result, k, i, j);
    }

    template<class RealTypeContext, class IntTypeContext>
    inline void Lattice<RealTypeContext, IntTypeContext>::computeProjectionLength(typename RealTypeContext::Real & result, unsigned k, unsigned b,
                                                                                  const linalg::math_rowvector<typename IntTypeContext::Integer> & vec) const
    // Computes the squared length of the projection of \sum_{i=0}^{m-1} vec[i] B_{b+i} onto the
    // orthogonal complement of B_0, ..., B_{k-1}. (Here, B_0, ..., B_{n-1} is the basis of the
    // lattice.)
    {
        d_gs->computeProjectionLength(result, k, b, vec);
    }

    template<class RealTypeContext, class IntTypeContext>
    inline void Lattice<RealTypeContext, IntTypeContext>::computeProjectionLengthBV(typename RealTypeContext::Real & result, unsigned k, unsigned b) const
    // Computes the squared length of the projection of B_b onto the orthogonal complement of B_0, ...,
    // B_{k-1}. (Here, B_0, ..., B_{n-1} is the basis of the lattice.)
    {
        d_gs->computeProjectionLengthBV(result, k, b);
    }

    template<class RealTypeContext, class IntTypeContext>
    inline void Lattice<RealTypeContext, IntTypeContext>::computeProjectionLength(typename RealTypeContext::Real & result, unsigned k,
                                                                                  const linalg::math_rowvector<typename IntTypeContext::Integer> & vec) const
    {
        d_gs->computeProjectionLength(result, k, vec);
    }

    template<class RealTypeContext, class IntTypeContext>
    inline bool Lattice<RealTypeContext, IntTypeContext>::canGetIntegerVector() const
    {
        return d_gs->canGetIntegerVector();
    }

    template<class RealTypeContext, class IntTypeContext>
    inline void Lattice<RealTypeContext, IntTypeContext>::getLinearCombination(linalg::math_rowvector<typename IntTypeContext::Integer> & result, unsigned ofs,
                                                                               const linalg::math_rowvector<typename IntTypeContext::Integer> & lincomb) const
    {
        d_gs->getLinearCombination(result, ofs, lincomb);
    }

    template<class RealTypeContext, class IntTypeContext>
    inline void Lattice<RealTypeContext, IntTypeContext>::getVector(linalg::math_rowvector<typename IntTypeContext::Integer> & result, unsigned ofs) const
    {
        d_gs->getVector(result, ofs);
    }

    template<class RealTypeContext, class IntTypeContext>
    inline linalg::math_matrix<typename IntTypeContext::Integer> * Lattice<RealTypeContext, IntTypeContext>::getMatrix() const
    // Returns matrix (if available).
    {
        return d_gs->getMatrixIndex(0).first;
    }

    template<class RealTypeContext, class IntTypeContext>
    inline std::pair<linalg::math_matrix<typename IntTypeContext::Integer> *, unsigned> Lattice<RealTypeContext, IntTypeContext>::getMatrixIndex(unsigned ofs) const
    // Returns matrix (if available) and index of specified vector <ofs> (might be different from
    // <ofs>).
    {
        return d_gs->getMatrixIndex(ofs);
    }

    template<class RealTypeContext, class IntTypeContext>
    inline const typename IntTypeContext::Integer & Lattice<RealTypeContext, IntTypeContext>::getNormSqUP(unsigned i) const
    {
        return d_gs->getNormSqUP(i);
    }

    template<class RealTypeContext, class IntTypeContext>
    inline void Lattice<RealTypeContext, IntTypeContext>::getNormSqUP(typename IntTypeContext::Integer & r, unsigned i) const
    {
        d_gs->getNormSqUP(r, i);
    }

    template<class RealTypeContext, class IntTypeContext>
    inline void Lattice<RealTypeContext, IntTypeContext>::getDotProduct(typename IntTypeContext::Integer & r, unsigned i, unsigned j) const
    {
        d_gs->getDotProduct(r, i, j);
    }

    template<class RealTypeContext, class IntTypeContext>
    inline void Lattice<RealTypeContext, IntTypeContext>::changeOfPrecision()
    {
        d_gs->changeOfPrecision();
    }

    template<class RealTypeContext, class IntTypeContext>
    std::ostream & operator << (std::ostream & s, const Lattice<RealTypeContext, IntTypeContext> & l)
    {
        s << "Lattice: dimension " << l.dimension() << ", range [" << l.range() << "]\n";
        for (unsigned i = 0; i < l.dimension(); ++i)
        {
            s << i << ":";
            for (unsigned j = 0; j < i; ++j)
                s << " " << l.getCoeff(i, j);
            s << " => " << l.getNormSq(i) << "\n";
        }
        l.gs()->print(s);
        return s;
    }
}

#endif

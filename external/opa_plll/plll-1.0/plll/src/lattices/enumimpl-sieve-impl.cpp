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

#ifndef PLLL_INCLUDE_GUARD__ENUMIMPL_SIEVE_IMPL_CPP
#define PLLL_INCLUDE_GUARD__ENUMIMPL_SIEVE_IMPL_CPP

namespace plll
{
    template<class RealIntTypeContext, class IntTypeContext>
    class DPLinForm // Dot product linear form with basis given (explicitly or by dot products)
    {
    private:
        linalg::math_rowvector<typename RealIntTypeContext::Type> d_coeffs;
        mutable typename RealIntTypeContext::Type d_temp;
        RealIntTypeContext & d_ric;
    
    public:
        DPLinForm(RealIntTypeContext & ric)
            : d_temp(ric), d_ric(ric)
        {
        }
    
        void reset_basis(const linalg::math_rowvector<typename IntTypeContext::Type> & u, const linalg::math_matrix<typename RealIntTypeContext::Type> & basis)
        {
            d_coeffs.resize(u.size(), linalg::Initialize(d_ric));
        
            // First, prepare the vector u*basis
            linalg::math_rowvector<typename RealIntTypeContext::Type> uu;
            uu.resize(u.size(), linalg::Initialize(d_ric));
            arithmetic::convert(d_temp, u[0], d_ric);
            uu = d_temp * basis.row(0);
            for (unsigned i = 1; i < u.size(); ++i)
            {
                arithmetic::convert(d_temp, u[i], d_ric);
                uu += d_temp * basis.row(i);
            }
        
            // Now compute dot product of u*basis with all rows of basis
            for (unsigned i = 0; i < u.size(); ++i)
                dot(d_coeffs[i], uu, basis.row(i));
        }
    
        void reset_basis(const linalg::math_rowvector<typename RealIntTypeContext::Type> & u, const linalg::math_matrix<typename RealIntTypeContext::Type> & basis)
        {
            d_coeffs.resize(u.size(), linalg::Initialize(d_ric));
        
            // First, prepare the vector u*basis
            linalg::math_rowvector<typename RealIntTypeContext::Type> uu;
            uu.resize(u.size(), linalg::Initialize(d_ric));
            uu = u[0] * basis.row(0);
            for (unsigned i = 1; i < u.size(); ++i)
                uu += u[i] * basis.row(i);
        
            // Now compute dot product of u*basis with all rows of basis
            for (unsigned i = 0; i < u.size(); ++i)
                dot(d_coeffs[i], uu, basis.row(i));
        }
    
        void reset_dotprods(const linalg::math_rowvector<typename IntTypeContext::Type> & u, const linalg::math_matrix<typename RealIntTypeContext::Type> & dot_products)
        {
            d_coeffs.resize(u.size(), linalg::Initialize(d_ric));
            for (unsigned i = 0; i < u.size(); ++i)
            {
                arithmetic::convert(d_temp, u[0], d_ric);
                d_coeffs[i] = d_temp * dot_products(i, 0);
                for (unsigned j = 1; j < u.size(); ++j)
                {
                    arithmetic::convert(d_temp, u[j], d_ric);
                    d_coeffs[i] += d_temp * dot_products(i, j);
                }
            }
        }
    
        void reset_dotprods(const linalg::math_rowvector<typename RealIntTypeContext::Type> & u, const linalg::math_matrix<typename RealIntTypeContext::Type> & dot_products)
        {
            d_coeffs.resize(u.size(), linalg::Initialize(d_ric));
            for (unsigned i = 0; i < u.size(); ++i)
            {
                d_coeffs[i] = u[0] * dot_products(i, 0);
                for (unsigned j = 1; j < u.size(); ++j)
                    d_coeffs[i] += u[j] * dot_products(i, j);
            }
        }
    
        void operator() (typename RealIntTypeContext::Type & result, const linalg::math_rowvector<typename IntTypeContext::Type> & v) const
        // Computes dot product of <u> and <v> w.r.t <basis>.
        {
            arithmetic::convert(d_temp, v[0], d_ric);
            result = d_temp * d_coeffs[0];
            for (unsigned i = 1; i < d_coeffs.size(); ++i)
            {
                arithmetic::convert(d_temp, v[i], d_ric);
                result += d_temp * d_coeffs[i];
            }
        }
    
        void operator() (typename RealIntTypeContext::Type & result, const linalg::math_rowvector<typename RealIntTypeContext::Type> & v) const
        // Computes dot product of <u> and <v> w.r.t <basis>.
        {
            result = v[0] * d_coeffs[0];
            for (unsigned i = 1; i < d_coeffs.size(); ++i)
                result += v[i] * d_coeffs[i];
        }
    
        typename RealIntTypeContext::Type operator() (const linalg::math_rowvector<typename IntTypeContext::Type> & v) const
        // Computes dot product of <u> and <v> w.r.t <basis>.
        {
            typename RealIntTypeContext::Type r(d_ric);
            operator()(r, v);
            return r;
        }
    
        typename RealIntTypeContext::Type operator() (const linalg::math_rowvector<typename RealIntTypeContext::Type> & v) const
        // Computes dot product of <u> and <v> w.r.t <basis>.
        {
            typename RealIntTypeContext::Type r(d_ric);
            operator()(r, v);
            return r;
        }
    };

    template<class IntTypeContext>
    class DPLinForm<IntTypeContext, IntTypeContext> // Dot product linear form with basis given (explicitly or by dot products)
    {
        // This specialization is needed because in the general template, having
        // RealIntTypeContext::Type == IntTypeContext::Type leads to problems as the same function
        // prototypes are defined with two different implementations. (The difference is in fact that
        // one version first converts to a temporary, while the second does the operation directly.)
    
    private:
        linalg::math_rowvector<typename IntTypeContext::Type> d_coeffs;
        IntTypeContext & d_ic;
    
    public:
        DPLinForm(IntTypeContext & ric)
            : d_ic(ric)
        {
        }
    
        void reset_basis(const linalg::math_rowvector<typename IntTypeContext::Type> & u, const linalg::math_matrix<typename IntTypeContext::Type> & basis)
        {
            d_coeffs.resize(u.size(), linalg::Initialize(d_ic));
        
            // First, prepare the vector u*basis
            linalg::math_rowvector<typename IntTypeContext::Type> uu;
            uu.resize(u.size(), linalg::Initialize(d_ic));
            uu = u[0] * basis.row(0);
            for (unsigned i = 1; i < u.size(); ++i)
                uu += u[i] * basis.row(i);
        
            // Now compute dot product of u*basis with all rows of basis
            for (unsigned i = 0; i < u.size(); ++i)
                dot(d_coeffs[i], uu, basis.row(i));
        }
    
        void reset_dotprods(const linalg::math_rowvector<typename IntTypeContext::Type> & u, const linalg::math_matrix<typename IntTypeContext::Type> & dot_products)
        {
            d_coeffs.resize(u.size(), linalg::Initialize(d_ic));
            for (unsigned i = 0; i < u.size(); ++i)
            {
                d_coeffs[i] = u[0] * dot_products(i, 0);
                for (unsigned j = 1; j < u.size(); ++j)
                    d_coeffs[i] += u[j] * dot_products(i, j);
            }
        }
    
        void operator() (typename IntTypeContext::Type & result, const linalg::math_rowvector<typename IntTypeContext::Type> & v) const
        // Computes dot product of <u> and <v> w.r.t <basis>.
        {
            result = v[0] * d_coeffs[0];
            for (unsigned i = 1; i < d_coeffs.size(); ++i)
                result += v[i] * d_coeffs[i];
        }
    
        typename IntTypeContext::Type operator() (const linalg::math_rowvector<typename IntTypeContext::Type> & v) const
        // Computes dot product of <u> and <v> w.r.t <basis>.
        {
            typename IntTypeContext::Type r(d_ic);
            operator()(r, v);
            return r;
        }
    };

    template<class T, int Rows, class ST, bool MO>
    bool isZero(const linalg::base_rowvector<T, Rows, ST, MO> & v)
    {
        for (linalg::size_type i = 0; i < v.size(); ++i)
            if (!isZero(v[i]))
                return false;
        return true;
    }

    template<class RealTypeContext>
    class GaussianGenerator
    /*
      Given a RNG which delivers a stream of uniformly distributed values in [0, 1), returns a stream of
      standard normal distributed values.
    */
    {
    private:
        RealTypeContext & d_rc;
        typename RealTypeContext::UniformRNG & d_rng;
    
        typename RealTypeContext::Real d_x, d_y, d_s, d_one, d_m;
        bool d_has_one;
    
    public:
        GaussianGenerator(RealTypeContext & rc, typename RealTypeContext::UniformRNG & rng)
            : d_rc(rc), d_rng(rng), d_x(rc), d_y(rc), d_s(rc), d_one(rc), d_m(rc), d_has_one(false)
        {
            setOne(d_one);
        }
    
        void generate(typename RealTypeContext::Real & result)
        /*
          Use Marsaglia's Polar Method to generate two independent standard normal distributed
          values. Returns one now and the second on the next call. The result is stored in <result>.
        */
        {
            // First, check cache
            if (d_has_one)
            {
                // Return the other
                result = d_y * d_m;
                d_has_one = false;
            }
            else
            {
                // Generate a pair
                do
                {
                    d_rng.randomUniform(d_x);
                    d_x <<= 1;
                    d_x -= d_one;
                    d_rng.randomUniform(d_y);
                    d_y <<= 1;
                    d_y -= d_one;
                    square(d_s, d_x);
                    square(d_m, d_y);
                    d_s += d_m;
                }
                while (isZero(d_s) || (d_s >= d_one));
                log(d_m, d_s);
                d_m /= d_s;
                d_m <<= 1;
                abs(d_m, d_m);
                sqrt(d_m, d_m);
                // Return one
                result = d_x * d_m;
                d_has_one = true;
            }
        }
    };

    template<class RealTypeContext>
    void GenSphere(linalg::math_rowvector<typename RealTypeContext::Real> & u, GaussianGenerator<RealTypeContext> & grng,
                   typename RealTypeContext::UniformRNG & rng, RealTypeContext & rc, const typename RealTypeContext::Real & radius)
    /*
      Generates a uniformly distributed point on the u.size()-dimensional sphere of radius <radius> and
      stores it in u. For this, first generate a point on the sphere $\partial B_{u.size()}(0, radius)$
      using Muller's and Marsaglia's method. Then multiply with $r^{1/u.size()}$, where r is uniformly
      randomly sampled from [0, 1]. Finally, multiply by <radius>.
    */
    {
        typename RealTypeContext::Real t(rc);
        const unsigned dim = u.size();
        do
        {
            setZero(t);
            for (unsigned i = 0; i < dim; ++i)
            {
                grng.generate(u[i]);
                t += square(u[i]);
            }
        } while (isZero(t));
        t = sqrt(t);
        typename RealTypeContext::Real rad(rc);
        rng.randomUniform(rad);
        rad = power(rad, arithmetic::convert(1, rc) / arithmetic::convert(dim, rc));
        rad /= t;
        rad *= radius;
        u *= rad;
    }

    class Percentage
    // A quick and dirty percentage displayer
    {
    private:
        Verbose & d_verbose;
        long d_total, d_last, d_last_pc;
    
    public:
        Percentage(Verbose & v, long total)
            : d_verbose(v), d_total(total), d_last(-1), d_last_pc(-1)
        {
        }
    
        void update(long curr)
        {
            long pc = (curr * 100) / d_total;
            if (pc != d_last_pc)
            {
                d_verbose(LatticeReduction::VL_Information) << pc << "%";
                d_last_pc = pc;
            }
            d_last = curr;
        }
    };

    template<class DestRealTypeContext, class Lattice>
    void computeBasisApproximation(linalg::math_matrix<typename DestRealTypeContext::Real> & basis,
                                   DestRealTypeContext & drc, unsigned begin, unsigned end, 
                                   Lattice & lattice)
    /*
      Ensures that the basis is lower triangular, i.e. basis.row(i) has all non-zero entries in range
      coefficients 0 up to i, and basis(i, i) > 0. This ensures that the associated Gram-Schmidt basis
      is the standard basis.
    */
    {
        const unsigned dim = end - begin + 1;
        basis.resize(dim, dim, linalg::Initialize(drc));
        // Create basis
        for (unsigned i = 0; i < dim; ++i)
        {
            for (unsigned j = 0; j < i; ++j)
            {
                arithmetic::convert(basis(i, j), lattice.getCoeff(begin + i, begin + j), drc);
                setZero(basis(j, i));
            }
            arithmetic::convert(basis(i, i), lattice.getNormSq(begin + i), drc);
            sqrt(basis(i, i), basis(i, i));
        }
        for (unsigned i = 0; i < dim; ++i)
            for (unsigned j = 0; j < i; ++j)
                basis(i, j) *= basis(j, j);
    }

    template<class DestRealTypeContext, class Lattice>
    void computeBasisApproximation(linalg::math_matrix<typename DestRealTypeContext::Real> & basis,
                                   linalg::math_matrix<typename DestRealTypeContext::Real> & inverse_basis,
                                   DestRealTypeContext & drc, unsigned begin, unsigned end, 
                                   Lattice & lattice)
    /*
      Ensures that the basis is lower triangular, i.e. basis.row(i) has all non-zero entries in range
      coefficients 0 up to i, and basis(i, i) > 0. This ensures that the associated Gram-Schmidt basis
      is the standard basis. The matrix inverse_matrix satisfies matrix * inverse_matrix = I_dim (up to
      rounding errors).
    */
    {
        // For the inverse basis, we proceed by first computing the inverse of the unscaled version, and
        // then by scaling by the inverse of the diagonal matrix. For inversion of the unscaled basis,
        // we use backwards substitution (similarly as in getDualGSI() in lattice.cpp), which requires
        // no division since all diagonal elements are 1.
        const unsigned dim = end - begin + 1;
        basis.resize(dim, dim, linalg::Initialize(drc));
        inverse_basis.resize(dim, dim, linalg::Initialize(drc));
        // Create basis, with correct diagonal elements, but with unscaled off-diagonal elements
        for (unsigned i = 0; i < dim; ++i)
        {
            for (unsigned j = 0; j < i; ++j)
            {
                arithmetic::convert(basis(i, j), lattice.getCoeff(begin + i, begin + j), drc);
                setZero(basis(j, i));
            }
            arithmetic::convert(basis(i, i), lattice.getNormSq(begin + i), drc);
            sqrt(basis(i, i), basis(i, i));
        }
        // Create unscaled inverse basis
        for (unsigned j = 0; j < dim; ++j)
            for (unsigned i = j + 1; i < dim; ++i)
            {
                // Compare discussion in lattice.cpp in getDualGSI() for the computation
                inverse_basis(i, j) = basis(i, j);
                for (unsigned k = j + 1; k < i; ++k)
                    inverse_basis(i, j) += basis(i, k) * inverse_basis(k, j);
                inverse_basis(i, j) = -inverse_basis(i, j);
                setZero(inverse_basis(j, i));
            }
        // Scale with (inverse of) diagonal
        for (unsigned i = 0; i < dim; ++i)
            for (unsigned j = 0; j < i; ++j)
                basis(i, j) *= basis(j, j);
        for (unsigned i = 0; i < dim; ++i)
        {
            setOne(inverse_basis(i, i));
            inverse_basis(i, i) /= basis(i, i);
            for (unsigned j = 0; j < i; ++j)
                inverse_basis(i, j) *= inverse_basis(i, i);
        }
    }

    template<class TypeContext>
    class TypeVecPComparator
    {
    public:
        bool operator() (const linalg::math_rowvector<typename TypeContext::Type> & v,
                         const linalg::math_rowvector<typename TypeContext::Type> & w) const // tests whether v < w (lexicographically)
        {
            // assert(v.size() == w.size());
            for (unsigned i = 0; i < v.size(); ++i)
            {
                int c = compare(v[i], w[i]);
                if (c != 0)
                    return c < 0;
            }
            return false; // since we test for strictly less
        }
    };
}

#endif

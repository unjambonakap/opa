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

#include <plll.hpp>
#include <plll/helper.hpp>
#include "profiling.hpp"
#include <list>
#include <set>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <limits>
#include <cmath>
#include <queue>

using namespace plll;

namespace plll
{
    typedef arithmetic::Integer IntType;
    
    namespace arithmetic
    {
        const linalg::math_rowvector<Integer> & toInteger(const linalg::math_rowvector<Integer> & x)
        {
            return x;
        }
    }
}

TimedDataCollector * d_cvpp;
TimedDataCollector * d_venum;
TimedDataCollector * d_vrelevant;


linalg::math_matrix<IntType> dot_products;

class DPLinForm // Dot product linear form with basisx
{
private:
    linalg::math_rowvector<IntType> d_coeffs;
    
public:
    void reset(const linalg::math_rowvector<IntType> & u, const linalg::math_matrix<IntType> & dot_products)
    {
        d_coeffs.resize(u.size());
        for (unsigned i = 0; i < u.size(); ++i)
        {
            setZero(d_coeffs[i]);
            for (unsigned j = 0; j < u.size(); ++j)
                d_coeffs[i] += u[j] * dot_products(i, j);
        }
    }
    
    void operator() (IntType & result, const linalg::math_rowvector<IntType> & v)
    // Computes dot product of <u> and <v> w.r.t <basis>.
    {
        setZero(result);
        for (unsigned i = 0; i < d_coeffs.size(); ++i)
            result += v[i] * d_coeffs[i];
    }
    
    IntType operator() (const linalg::math_rowvector<IntType> & v)
    // Computes dot product of <u> and <v> w.r.t <basis>.
    {
        IntType r;
        operator()(r, v);
        return r;
    }
};

/*
  Note that we only store *half* the Voronoi cell, i.e. we leave away -v if v is included in it.
 */

void CVPP(unsigned dim, linalg::math_rowvector<IntType> & v,
          const std::list<std::pair<linalg::math_rowvector<IntType>, IntType> > & V)
{
    Profiler<TimedDataCollector> prof(d_cvpp[dim - 1]);
    DPLinForm dp;
    bool min_sign_val = false, sign_val;
    IntType min_val, min_val_d, val;
    IntType alpha, tmp;
    while (true)
    {
        setZero(min_val_d);
        setOne(min_val);
        dp.reset(v, dot_products);
        std::list<std::pair<linalg::math_rowvector<IntType>, IntType> >::const_iterator min_i = V.end();
        for (std::list<std::pair<linalg::math_rowvector<IntType>, IntType> >::const_iterator i = V.begin(); i != V.end(); ++i)
        {
            dp(val, i->first);
            if (sign(val) < 0)
            {
                neg(val, val);
                sign_val = true;
            }
            else
                sign_val = false;
            if (val * min_val_d < min_val * i->second)
            {
                min_val = val;
                min_val_d = i->second;
                min_sign_val = sign_val;
                min_i = i;
            }
        }
        if (min_val <= min_val_d)
            break;
        euclideanDivision(alpha, tmp, min_val, min_val_d);
        if (!isZero(tmp))
            ++alpha;
        // Note that alpha is always positive
        if (isOne(alpha))
        {
            if (min_sign_val)
                v += min_i->first;
            else
                v -= min_i->first;
        }
        else
        {
            if (min_sign_val)
                v += alpha * min_i->first;
            else
                v -= alpha * min_i->first;
        }
    }
}

unsigned getIndex(unsigned dim, const linalg::math_rowvector<IntType> & v)
{
    unsigned res = 0;
    for (unsigned i = 0; i < dim; ++i)
    {
        res <<= 1;
        if (bit(v[i], 0))
            res |= 1;
    }
    return res;
}

class SortIteratorSecond
{
public:
    template<class It>
    bool operator() (const It & x, const It & y) const
    {
        return x->second < y->second;
    }
};

class VecIntPairComparatorNeg
{
private:
    const std::vector<std::pair<linalg::math_rowvector<IntType>, IntType> > & d_A;
    
public:
    VecIntPairComparatorNeg(const std::vector<std::pair<linalg::math_rowvector<IntType>, IntType> > & A)
        : d_A(A)
    {
    }
    
    bool operator() (unsigned x, unsigned y) const
    {
        int c = compare(d_A[x].second, d_A[y].second);
        return (c == 0) ? x < y : c > 0;
    }
};

void VEnum(unsigned dim, std::vector<std::pair<linalg::math_rowvector<IntType>, IntType> > & A, std::vector<bool> & Aused,
           const std::list<std::pair<linalg::math_rowvector<IntType>, IntType> > & V, const std::vector<unsigned> & Vidx,
           const std::pair<linalg::math_rowvector<IntType>, IntType> & t)
{
    Profiler<TimedDataCollector> prof(d_venum[dim - 1]);
    for (unsigned i = 0; i < Aused.size(); ++i)
        Aused[i] = false;
    
    // Solve CVPP
    std::pair<linalg::math_rowvector<IntType>, IntType> y0 = t;
    CVPP(dim, y0.first, V);
    DPLinForm dp;
    dp.reset(y0.first, dot_products);
    dp(y0.second, y0.first);
    
    // Store into A and initialize priority queue
    unsigned idx = getIndex(dim, y0.first);
    A[idx] = y0;
    Aused[idx] = true;
    std::set<unsigned, VecIntPairComparatorNeg>
        Q = std::set<unsigned, VecIntPairComparatorNeg>(VecIntPairComparatorNeg(A));
    Q.insert(idx);
    
    // Do while loop
    IntType yn;
    IntType n;
    while (!Q.empty())
    {
        // Get shortest element from queue
        std::set<unsigned, VecIntPairComparatorNeg>::iterator i = Q.end();
        --i;
        unsigned ui = *i;
        // std::pair<linalg::math_rowvector<IntType>, IntType> u = A[ui];
        // std::cout << "u = " << A[ui].first << " (" << A[ui].second << ")\n";
        Q.erase(i);
        
        // Iterate through all elements of u+V and u-V
        dp.reset(A[ui].first, dot_products);
        unsigned yii = 0;
        for (std::list<std::pair<linalg::math_rowvector<IntType>, IntType> >::const_iterator yi = V.begin(); yi != V.end(); ++yi, ++yii)
        {
            // Compute norm of y = u +- *yi
            yn = A[ui].second + yi->second;
            dp(n, yi->first);
            n <<= 1;
            bool wasneg = sign(n) < 0;
            if (!wasneg)
                // We start with smaller of the two
                neg(n, n);
            // Compare norms
            if (A[ui].second <= yn + n)
            {
                yn += n;
                // Trick 17: compute index for y.first = A[ui].first + yi->first using indices of
                // A[ui].first (which is ui) and yi->first (which is Vidx[yii])
                unsigned idx = ui ^ Vidx[yii];
                // Compare to table
//                std::cout << "y = " << A[ui].first + yi->first << " (" << yn << ")\n";
                if (!Aused[idx])
                {
                    // Update table and priority queue
                    Aused[idx] = true;
                    if (wasneg)
                        A[idx].first = A[ui].first + yi->first;
                    else
                        A[idx].first = A[ui].first - yi->first;
                    A[idx].second = yn;
                    Q.insert(idx);
                }
                else if (yn < A[idx].second)
                {
                    // Update table and priority queue
                    i = Q.find(idx);
                    if (i != Q.end())
                        Q.erase(i);
                    if (wasneg)
                        A[idx].first = A[ui].first + yi->first;
                    else
                        A[idx].first = A[ui].first - yi->first;
                    A[idx].second = yn;
                    Q.insert(idx);
                }
            }
            else
            {
                neg(n, n);
                yn += n;
                if (A[ui].second <= yn)
                {
                    // Trick 17: compute index for y.first = A[ui].first + yi->first using indices of
                    // A[ui].first (which is ui) and yi->first (which is Vidx[yii])
                    unsigned idx = ui ^ Vidx[yii];
                    // Compare to table
//                std::cout << "y = " << A[ui].first + yi->first << " (" << yn << ")\n";
                    if (!Aused[idx])
                    {
                        // Update table and priority queue
                        Aused[idx] = true;
                        if (wasneg)
                            A[idx].first = A[ui].first - yi->first;
                        else
                            A[idx].first = A[ui].first + yi->first;
                        A[idx].second = yn;
                        Q.insert(idx);
                    }
                    else if (yn < A[idx].second)
                    {
                        // Update table and priority queue
                        i = Q.find(idx);
                        if (i != Q.end())
                            Q.erase(i);
                        if (wasneg)
                            A[idx].first = A[ui].first - yi->first;
                        else
                            A[idx].first = A[ui].first + yi->first;
                        A[idx].second = yn;
                        Q.insert(idx);
                    }
                }
            }
        }
    }
    /*
    for (unsigned i = 0; i < A.size(); ++i)
    {
        std::cout << i << ": ";
        if (A[i].first.size())
            std::cout << A[i].first << " " << getIndex(dim, A[i].first) << "\n";
        else
            std::cout << "none\n";
    }
    */
}

void VRelevant(unsigned dim, std::list<std::pair<linalg::math_rowvector<IntType>, IntType> > & dest,
               const std::list<std::pair<linalg::math_rowvector<IntType>, IntType> > & src,
               const std::pair<linalg::math_rowvector<IntType>, IntType> & b, const arithmetic::Real & gs_norm_sq)
{
    Profiler<TimedDataCollector> prof(d_vrelevant[dim - 1]);
    std::vector<std::pair<linalg::math_rowvector<IntType>, IntType> > T, A; // "T[i].first.size() == 0" means "infinity"
    std::vector<bool> Tused, Aused;
    T.resize(1u << dim);
    A.resize(1u << dim);
    Tused.resize(1u << dim);
    Aused.resize(1u << dim);
    for (unsigned i = 0; i < T.size(); ++i)
    {
        T[i].first.resize(dim);
        A[i].first.resize(dim);
        Tused[i] = false;
    }
    unsigned infinity_count = 1u << dim;
    bool recompute_max = false;
    IntType max, mm;
    setZero(max);
    std::vector<unsigned> srcidx;
    srcidx.resize(src.size());
    unsigned ii = 0;
    for (std::list<std::pair<linalg::math_rowvector<IntType>, IntType> >::const_iterator i = src.begin(); i != src.end(); ++i, ++ii)
        srcidx[ii] = getIndex(dim, i->first);
    for (unsigned m = 0; ; ++m)
    {
        if (infinity_count == 0)
        {
            arithmetic::RealContext rc;
            if (arithmetic::convert(m, rc) * gs_norm_sq >= arithmetic::convert(max, rc))
                break;
        }
        std::pair<linalg::math_rowvector<IntType>, IntType> t = b;
        mm = arithmetic::convert<IntType>(m);
        t.first *= mm;
        mm = square(mm);
        t.second *= mm;
        VEnum(dim, A, Aused, src, srcidx, t);
        for (unsigned i = 0; i < A.size(); ++i)
            if (Aused[i])
            {
                unsigned idx = getIndex(dim, A[i].first);
                if (!Tused[idx])
                {
                    Tused[idx] = true;
                    if (infinity_count == (1u << dim))
                        max = A[i].second;
                    else
                        if (max < A[i].second)
                            max = A[i].second;
                    --infinity_count;
                    T[idx] = A[i];
                    Tused[idx] = true;
                }
                else if (T[idx].second > A[i].second)
                {
                    if (T[idx].second == max)
                        recompute_max = true;
                    T[idx] = A[i];
                }
            }
        if (recompute_max)
        {
            recompute_max = false;
            setZero(max);
            for (unsigned i = 0; i < T.size(); ++i)
                if (Tused[i])
                    if (max < T[i].second)
                        max = T[i].second;
        }
    }
    dest.clear();
    for (unsigned i = 0; i < T.size(); ++i)
        if (Tused[i])
            dest.push_back(T[i]);
}

void VFilter(std::list<std::pair<linalg::math_rowvector<IntType>, IntType> > & dest,
             const std::list<std::pair<linalg::math_rowvector<IntType>, IntType> > & src)
{
    dest.clear();
    // Create sorted version of src by norm
    std::vector<std::list<std::pair<linalg::math_rowvector<IntType>, IntType> >::const_iterator> sorted;
    sorted.reserve(src.size());
    for (std::list<std::pair<linalg::math_rowvector<IntType>, IntType> >::const_iterator i = src.begin(); i != src.end(); ++i)
        if (!isZero(i->second))
            sorted.push_back(i);
    std::sort(sorted.begin(), sorted.end(), SortIteratorSecond());
    // Go through sorted list
    IntType n1, n2;
    DPLinForm dp;
    for (unsigned i = 0; i < sorted.size(); ++i)
    {
        dp.reset(sorted[i]->first, dot_products);
        bool fail1 = false, fail2 = false;
        for (std::list<std::pair<linalg::math_rowvector<IntType>, IntType> >::iterator j = dest.begin(); j != dest.end(); ++j)
        {
            // Note that ||2*sorted[i]->first +- j->first||^2 = 4 * sorted[i]->second + j->second +-
            // 4 * <sorted[i]->first, j->first>. Therefore, comparing this norm to ||j->first||^2
            // only involves comparing sorted[i]->second +- <sorted[i]->first, j->first> to 0.
            n1 = sorted[i]->second;
            dp(n2, j->first);
            // Check whether ||2*sorted[i]->first - j->first||^2 > ||j->first||^2 is violated
            if (n1 <= n2)
            {
                fail1 = true;
                if (fail2)
                    break;
            }
            // Check whether ||2*sorted[i]->first + j->first||^2 > ||j->first||^2 is violated
            if (n1 <= -n2)
            {
                fail2 = true;
                if (fail1)
                    break;
            }
        }
        if (!fail1 || !fail2)
            dest.push_back(*sorted[i]);
    }
}

int main(int argc, char ** argv)
{
    arithmetic::initArithmeticThreadAllocators();

    if (argc != 2)
    {
        std::cout << "Syntax: " << argv[0] << " <filename>\n";
        return -1;
    }
    std::string fn = argv[1];
    
    linalg::math_matrix<arithmetic::Integer> A;
    std::ifstream f(fn.c_str());
    f >> A;
    if (A.rows() == 0)
    {
        std::cout << "Cannot load '" << fn << "'!\n";
        return -1;
    }
    
    LatticeReduction lr(A);
    lr.setArithmetic(LatticeReduction::A_LongDouble);
    lr.setGramSchmidt(LatticeReduction::G_NumStable);
    lr.bkz(0.999, A.rows() / 2);
    lr.bkz(0.999, A.rows() / 2, LatticeReduction::BKZ_SlideReduction);
    A = lr.getLattice();
    
    dot_products.resize(A.rows(), A.rows());
    for (unsigned i = 0; i < A.rows(); ++i)
        for (unsigned j = 0; j <= i; ++j)
        {
            dot_products(i, j) = arithmetic::convert<IntType>(dot(A.row(i), A.row(j)));
            if (j < i)
                dot_products(j, i) = dot_products(i, j);
        }
    
    std::list<std::pair<linalg::math_rowvector<IntType>, IntType> > V, U;
    linalg::math_rowvector<IntType> e;
    e.resize(A.rows());
    setOne(e[0]);
    V.push_back(std::make_pair(e, dot_products(0, 0)));
    std::cout << "\nSTARTING VORONOI CELL COMPUTATION\n";
    
    linalg::math_matrix<std::string> names;
    d_cvpp = new TimedDataCollector[A.rows()];
    d_venum = new TimedDataCollector[A.rows()];
    d_vrelevant = new TimedDataCollector[A.rows()];
    names.resize(A.rows(), 3);
    for (unsigned i = 0; i < A.rows(); ++i)
    {
        std::ostringstream s;
        s << "CVPP d=" << i+1;
        names(i, 0) = s.str();
        s.str("");
        s << "VEnum d=" << i+1;
        names(i, 1) = s.str();
        s.str("");
        s << "VRelevant d=" << i+1;
        names(i, 2) = s.str();
        d_cvpp[i].setName(names(i, 0).c_str());
        d_venum[i].setName(names(i, 1).c_str());
        d_vrelevant[i].setName(names(i, 2).c_str());
    }
    
    std::cout << 1 << ": " << V.size() << "\n";
    /*
    for (std::list<std::pair<linalg::math_rowvector<IntType>, IntType> >::iterator i = V.begin(); i != V.end(); ++i)
        std::cout << "    " << i->first << " (" << i->second << ")\n";
    */
    for (unsigned k = 1; k < A.rows(); ++k)
    {
        arithmetic::Real gs_norm_sq = lr.getGSSqNormR(k);
        for (unsigned i = 0; i < e.size(); ++i)
            setZero(e[i]);
        setOne(e[k]);
        VRelevant(k + 1, U, V, std::make_pair(e, dot_products(k, k)), gs_norm_sq);
// We can uncomment this for "early abort" if we just want to solve SVP without having the Voronoi
// cell. In that case we have to search for the result in U though.
//        if (k + 1 == A.rows())
//            break;
        VFilter(V, U);
        std::cout << k+1 << ": " << V.size() << "\n";
        /*
        for (std::list<std::pair<linalg::math_rowvector<IntType>, IntType> >::iterator i = V.begin(); i != V.end(); ++i)
            std::cout << "    " << i->first << " (" << i->second << ")\n";
        */
    }
    
    std::list<std::pair<linalg::math_rowvector<IntType>, IntType> >::iterator dest_i = V.end();
    IntType dest_len;
    setZero(dest_len);
    for (std::list<std::pair<linalg::math_rowvector<IntType>, IntType> >::iterator i = V.begin(); i != V.end(); ++i)
        if ((dest_i == V.end()) || (i->second < dest_len))
        {
            dest_len = i->second;
            dest_i = i;
        }
    if (dest_i != V.end())
    {
        std::cout << "Shortest vector has squared norm " << dest_len << "\n";
        std::cout << arithmetic::toInteger(dest_i->first) * A << "\n";
    }
    else
        std::cout << "No shortest vector found!\n";
    
    delete[] d_cvpp;
    delete[] d_venum;
    delete[] d_vrelevant;
}

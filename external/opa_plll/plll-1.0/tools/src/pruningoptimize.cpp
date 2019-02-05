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
#include "profiling.hpp"
#include <plll/arguments.hpp>
#include <plll/linalg.hpp>
#include "polynomial.hpp"
#include "polynomial-real.hpp"
#include "taskmanager.hpp"
#include <list>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cctype>
#include <cmath>

using namespace plll;

void help(char * name)
{
    std::cout << name << " <options>\n";
    std::cout << "\n";
    std::cout << "Options:\n";
    std::cout << "  -dimension=<dim> (should be even if -montecarlo is not specified)\n";
    std::cout << "  -maxthreads=<count> (number of threads to use; default: number of cores)\n";
    std::cout << "  -randomize: use /dev/(u)random for RNG initialization and not just time\n";
    std::cout << "  -montecarlo: use Monte Carlo simulation without trying to minimize encapsulating body\n";
    std::cout << "\n";
    std::cout << "Curve selection:\n";
    std::cout << "  -read=<filename>: read curve as vector from file\n";
    std::cout << "  -linar: use linar curve\n";
    std::cout << "  -gnrpoly: use degree-8 polynomial approximating curve from Gama-Nguyen-Regev\n";
    std::cout << "  -parse=<vector>: give curve directly\n";
    std::cout << "  -simplify: \"simplify\" the given curve so that pruning[2*i] == pruning[2*i+1] for i=0,1,2,...\n";
    std::cout << "\n";
    std::cout << "Output:\n";
    std::cout << "  -shownodes: show number of nodes as text output and plot\n";
    std::cout << "  -dump: show pruning curve as vector\n";
    std::cout << "  -plot: show pruning curve as plot\n";
    std::cout << "  -write=<filename>: writes curve as vector into file\n";
    std::cout << "\n";
    std::cout << "Commands w/ options:\n";
    std::cout << "  -readqs <filename1> [<filename2> [...]]: read q values from lattice file(s)\n";
    std::cout << "  -onlyprob: only shows final probability of curve, and not node count\n";
    std::cout << "  -onlynodes: only shows final node count of curve, and not probability\n";
    std::cout << "  -optimize: tries to optimize curve\n";
    std::cout << "\n";
}

bool sanitize(linalg::math_rowvector<arithmetic::Real> & pruning)
// Ensures that pruning[pruning.size() - 1] == 1, and that pruning is non-decreasing. Returns true
// in case it already was sane.
{
    bool sane = true;
    if (pruning.size())
    {
        if (!isOne(pruning[pruning.size() - 1]))
        {
            setOne(pruning[pruning.size() - 1]);
            sane = false;
        }
        for (unsigned j = pruning.size() - 1; j > 0; --j)
        {
            if (pruning[j - 1] > pruning[j])
            {
                pruning[j - 1] = pruning[j];
                sane = false;
            }
            if (sign(pruning[j - 1]) < 0)
            {
                sane = false;
                setZero(pruning[j - 1]);
            }
        }
    }
    return sane;
}

void solve_regulafalsi(arithmetic::Real & lower, arithmetic::Real & upper, const realpoly<1, arithmetic::Real> & f, const arithmetic::Real & dest, unsigned iterations, arithmetic::RealContext & rc)
// Assumes that f(lower) <= dest <= f(upper). Does <iterations> Regula Falsi iterations,
// i.e. decreases upper-lower by a factor of 2^iterations, such that this inequality still holds.
{
    arithmetic::Real tmp1(rc), tmp2(rc);
    for (; iterations > 0; --iterations)
    {
        tmp1 = lower + upper;
        tmp1 >>= 1;
        f(tmp1).evaluate_to(tmp2);
        tmp2 -= dest;
        if (sign(tmp2) < 0)
            lower = tmp1;
        else if (sign(tmp2) > 0)
            upper = tmp1;
        else
        {
            lower = upper = tmp1;
            break;
        }
    }
}

void solve_newton(arithmetic::Real & sol, const realpoly<1, arithmetic::Real> & f, const realpoly<1, arithmetic::Real> & f_deriv, const arithmetic::Real & dest,
                  const arithmetic::Real & stopbound, arithmetic::RealContext & rc)
// Perform Newton iterations to solve f(sol) = dest. Assumes sol contains a suitable starting
// value. The polynomial f_deriv is assumed to be the derivative of f. Stops whenever |poly(sol) -
// dest| < stopbound.
{
    arithmetic::Real tmp1(rc), tmp2(rc);
    long loops = 0;
    while (true)
    {
        if (loops > 20)
            // Second abortion condition: too many iterations (probably due to floating point precision)
            return;
        // Evaluate poly(sol) - dest to tmp
        f(sol).evaluate_to(tmp1);
        tmp1 -= dest;
        // Check abortion condition
        if (compareAbsValues(tmp1, stopbound) < 0)
            return;
        // Evaluate deriv(sol) to tmp2
        f_deriv(sol).evaluate_to(tmp2);
        // Compute Newton iteration sol_new = sol_old - f(sol_old) / f'(sol_old)
        tmp1 /= tmp2;
        sol -= tmp1;
        ++loops;
    };
}

void computeCircularProbability(arithmetic::Real & prob, const arithmetic::Real & r, const arithmetic::Real & s, const arithmetic::Real & pi, arithmetic::RealContext & rc)
// Computes the probability that x, y with x^2 + y^2 = r satisfy x^2 <= s. Stores the result in
// prob. Expects pi to be an approximation of pi.
{
    // We have to compute the measure of t in [0, 1] with r * sin(2*pi*t)^2 <= s.
    prob = s / r;
    assert(!isNegative(prob));
    prob = sqrt(prob);
    prob = asin(prob);
    prob /= pi;
    prob <<= 1;
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

template<class A, class B>
std::ostream & operator << (std::ostream & s, const std::pair<A, B> & p)
{
    return s << "(" << p.first << ", " << p.second << ")";
}

arithmetic::Real computeProbability(const linalg::math_rowvector<arithmetic::Real> & pruning, arithmetic::RandomNumberGenerator & rng, arithmetic::RealContext & rc, bool montecarlo)
/*
  Compute the probability that a vector x of norm 1 satisfies \sum_{j=0}^i x[j]^2 \le pruning[j] for
  i=0..dimension-1.
 */
{
    arithmetic::Real res(rc);
    if (pruning.size() == 0)
    {
        setOne(res);
        return res;
    }
    setOne(res);
    const unsigned d = pruning.size();
    if (pruning[d - 1] < res)
    {
        // The last coefficient should be at least one!
        setZero(res);
        return res;
    }
    if (d == 1)
    {
        // Special case dimension 1: since the pruning coefficient is >= 1, the probability is 1
        return res;
    }
    arithmetic::Real pi(rc);
    rc.getPi(pi);
    if (montecarlo)
    {
        // Do a complete Monte Carlo approach
        arithmetic::RealContext::UniformRNG urng(rng);
        GaussianGenerator<arithmetic::RealContext> gg(rc, urng);
        setZero(res);
        unsigned total = 100000, inside = 0;
        arithmetic::Real tmp(rc), n(rc);
        linalg::math_rowvector<arithmetic::Real> v;
        v.resize(d);
        for (unsigned i = 0; i < v.size(); ++i)
            v[i].setContext(rc);
        for (unsigned i = 0; i < total; ++i)
        {
            // Generate point in unit ball
            while (true)
            {
                // Compute Gaussian distributed point and its norm
                setZero(n);
                for (unsigned j = 0; j < v.size(); ++j)
                {
                    gg.generate(v[j]);
                    square(tmp, v[j]);
                    n += tmp;
                }
                if (!isZero(n))
                    break;
            }
            // Divide by norm
            setOne(tmp);
            sqrt(n, n);
            tmp /= n;
            v *= tmp;
            
            // Test wether v satisfies conditions
            setZero(n);
            bool is_inside = true;
            for (unsigned j = 0; j + 1 < v.size(); ++j) // don't check v[v.size()-1], since the
                                                        // condition should be satisfied here (as
                                                        // pruning[v.size()-1] >= 1), but might not
                                                        // be due to numerical inaccuracies. If j
                                                        // runs until v.size() - 1, the error can be
                                                        // seen in practice, already for linear
                                                        // pruning curves in dimension 2.
            {
                square(tmp, v[j]);
                n += tmp;
                if (n > pruning[j])
                {
                    is_inside = false;
                    break;
                }
            }
            if (is_inside)
                ++inside;
        }
        // Return result
        std::cout << "Monte Carlo: " << inside << "/" << total << " => " << 100.0 * inside / (double)total << "%\n";
        arithmetic::convert(res, inside, rc);
        arithmetic::convert(tmp, total, rc);
        res /= tmp;
        return res;
    }
    if (d == 2)
    {
        // Special case dimension 2: here we can solve it exactly and directly.
        computeCircularProbability(res, pruning[1], pruning[0], pi, rc);
        return res;
    }
    
    assert((d & 1) == 0); // ONLY WORKS FOR DIMENSION EVEN!!!
    
    /*
      We use a hybrid Monte Carlo apporach.
      
      As a bounding volume, we increase pruning[2*i] or pruning[2*i+1] such that both are
      equal. Then, according to Gama, Nguyen, Regev: "Lattice Enumeration Using Extreme Pruning",
      the probability that a vector from the unit ball lies in the volume can be computed
      exactly. For this, note that the distribution of the vector (x[0]^2+x[1]^2, x[2]^2+x[3]^2,
      x[4]^2+x[5]^2, ...) is given by a Dirichlet distribution with parameters (1, ..., 1) (a vector
      of dimension d/2+1), which is a uniform distribution over the set of all vectors whose
      coordinates are non-negative and sum to at most 1.
      
      Thus we compute the probability that an (d/2)-dimensional vector x of non-negative
      coordinates which sum to at most 1 satisfies x[i] <= max(pruning[2*i], pruning[2*i+1]). The
      result is multiplied by the volume of a ball of dimension d: this yields the volume of
      the bounding set.
     */
    
    // Compute bounding in dimension d/2
    bool simple_function = true; // indicates that pruning[2*i] == pruning[2*i+1] for all i
    linalg::math_rowvector<arithmetic::Real> bounding;
    bounding.resize(d / 2);
    for (unsigned i = 0; i < bounding.size(); ++i)
        bounding[i].setContext(rc);
    for (unsigned i = bounding.size(); i > 0; --i)
    {
        if (pruning[2 * i - 2] != pruning[2 * i - 1])
        {
            simple_function = false;
            bounding[i - 1] = std::max(pruning[2 * i - 2], pruning[2 * i - 1]);
        }
        else
            bounding[i - 1] = pruning[2 * i - 1];
        // Make sure that bounding[] is non-decreasing
        if (i < bounding.size())
            if (bounding[i - 1] > bounding[i])
                bounding[i - 1] = bounding[i];
    }
    
    // Compute probability that if x is sampled from a sphere of radius 1, that x[0]^2+x[1]^2 <=
    // bounding[, x[2]^2+x[3]^2, ...).
    
    // Compute integrals
    linalg::base_rowvector<realpoly<1, arithmetic::Real> > integrals;
    arithmetic::Real tmp(rc), tmp2(rc);
    integrals.resize(bounding.size());
    integrals[0] = RealMonomial<arithmetic::Real>(0);
    for (unsigned i = 1; i < integrals.size(); ++i)
    {
        // Integrate
        integral<0>(integrals[i], integrals[i - 1]);
        // Evaluate the integral in bounding[bounding.size()-1]-bounding[bounding.size()-i-1] using
        // Horner's method. The negative of the result is stored in integrals[i][0].
        tmp = bounding[bounding.size() - 1] - bounding[bounding.size() - i - 1];
        integrals[i](tmp).evaluate_to(tmp2);
        integrals[i][0] = -tmp2;
    }
    // Evaluate last integral at bounding[bounding.size()-1] using Horner's method. Evaluate this into res.
    integrals[integrals.size() - 1](bounding[bounding.size() - 1]).evaluate_to(res);
    // Multiply res by (bounding.size()-1)!
    for (unsigned i = 2; i < bounding.size(); ++i)
    {
        arithmetic::convert(tmp, i, rc);
        res *= tmp;
    }
    // Now res can be multiplied by the probability obtained from the Monte Carlo algorithm.
    if (simple_function)
        // In case of simple functions, there is no need to run the Monte Carlo part...
        return res;
    
    // Apply Monte Carlo algorithm
    unsigned total = 1000;
    arithmetic::Real inside_p(rc);
    setZero(inside_p);
    linalg::math_rowvector<arithmetic::Real> t;
    t.resize(d/2);
    for (unsigned i = 0; i < t.size(); ++i)
        t[i].setContext(rc);
    arithmetic::RealContext::UniformRNG urng(rng);
    arithmetic::Real sumsofar(rc), val(rc);
    arithmetic::Real stopbound(rc), lower(rc);
    arithmetic::Real p(rc);
    for (unsigned i = 0; i < total; ++i)
    {
        while (true)
        {
            // Generate a uniformly sampled vector (t[0], ..., t[d/2-1]) which satisfies \sum_{j=0}^i
            // t[j] <= bounding[i] for i=0..d/2-1 and \sum_{j=0}^{d/2-1} t[j] = 1.
            sumsofar = bounding[bounding.size() - 1];
            bool restart = false;
            for (unsigned j = 0; j + 2 < t.size(); ++j)
            {
                // We want that val is distributed in [bounding[bounding.size()-1] - bounding[j],
                // sumsofar] according to the distribution given by (a modified version of) the
                // distribution function integrals[j] with (a modified version of) the density
                // integrals[j-1].
                
                lower = bounding[bounding.size()-1] - bounding[j];
                // Generate uniform random variable on [0, 1]
                urng.randomUniform(t[j]);
                // Evaluate integrals[j] at sumsofar and multiply with uniform random variable on [0,
                // 1]. This is the destination.
                integrals[integrals.size() - 1 - j](sumsofar).evaluate_to(tmp);
                tmp *= t[j];
                // Chose stopping bound
                rc.getEpsilon(stopbound);
                stopbound *= tmp;
                abs(stopbound, stopbound);
                stopbound <<= 4;
                // Solve integrals[j](val) = tmp
                val = sumsofar;
                tmp2 = lower;
                solve_regulafalsi(tmp2, val, integrals[integrals.size() - 1 - j], tmp, 10, rc);
                val += tmp2;
                val >>= 1;
                solve_newton(val, integrals[integrals.size() - 1 - j], integrals[integrals.size() - 2 - j], tmp, stopbound, rc);
                if ((val < lower) || (val > sumsofar))
                {
                    // Apparently, we had numerical(?) problems
                    restart = true;
                    break;
                }
                // Now update sumsofar and t[j].
                t[j] = sumsofar - val;
                sumsofar = val;
            }
            if (restart)
                continue;
            // Now generate t[t.size() - 2] and t[t.size() - 1]
            urng.randomUniform(t[t.size() - 2]);
            tmp = bounding[t.size() - 2] - bounding[bounding.size() - 1];
            if (isNegative(tmp))
                tmp += sumsofar;
            else
                tmp = sumsofar;
            t[t.size() - 2] *= tmp;
            t[t.size() - 1] = sumsofar - t[t.size() - 2];
            // Now check whether a vector belonging to this choice satisfies the pruning conditions
            setZero(sumsofar);
            setOne(p);
            for (unsigned j = 0; j < t.size(); ++j)
            {
//            std::cout << j << ": " << t[j] << " " << p << "\n";
                // We have to determine random x and y such that x^2 + y^2 == t[j], and then check the
                // conditions, i.e. sumsofar + x^2 <= pruning[2 * j] and sumsofar + x^2 + y^2 <=
                // pruning[2 * j + 1]. The latter condition is always satisfied since pruning[2 * j + 1]
                // == bounding[j] and sumsofar + x^2 + y^2 == sumsofar + t[j]. The first condition is
                // only in danger if sumsofar + t[j] is larger than pruning.
                if (sumsofar + t[j] > pruning[2 * j])
                {
                    computeCircularProbability(tmp, t[j], pruning[2 * j] - sumsofar, pi, rc);
                    p *= tmp;
                }
                // Increase sumsofar
                sumsofar += t[j];
                if (sumsofar > pruning[2 * j + 1])
                {
                    // This should not happen, only due to numerical problems...
                    restart = true;
                    break;
                }
            }
            if (restart)
                continue;
            inside_p += p;
            break;
        }
    }
    // Multiply res with inside/total.
    res *= inside_p;
    arithmetic::convert(tmp, total, rc);
    res /= tmp;
    // Return result
    return res;
}

arithmetic::Real computeVolume(const linalg::math_rowvector<arithmetic::Real> & pruning, unsigned d, arithmetic::RandomNumberGenerator & rng, arithmetic::RealContext & rc, bool montecarlo = false)
/*
  Compute that the volume of the set of vectors x that satisfy \sum_{j=0}^i x[j]^2 \le pruning[j] for
  i=0..d-1. We assume that pruning[j] <= 1 for all j.
 */
{
    arithmetic::Real res(rc);
    if (d == 0)
    {
        // In dimension 0, everything is trivial...
        setOne(res);
        return res;
    }
    
    arithmetic::Real pi(rc);
//    arithmetic::convert(pi, 3.14159265358979323846264338327950288l, rc);
    rc.getPi(pi);
    
    if (montecarlo)
    {
        // Do a complete Monte Carlo approach
        arithmetic::RealContext::UniformRNG urng(rng);
        GaussianGenerator<arithmetic::RealContext> gg(rc, urng);
        setZero(res);
        unsigned total = 100000, inside = 0;
        arithmetic::Real tmp(rc), n(rc);
        linalg::math_rowvector<arithmetic::Real> v;
        v.resize(d);
        for (unsigned i = 0; i < v.size(); ++i)
            v[i].setContext(rc);
        for (unsigned i = 0; i < total; ++i)
        {
            // Generate point in unit ball
            while (true)
            {
                // Compute Gaussian distributed point and its norm
                setZero(n);
                for (unsigned j = 0; j < v.size(); ++j)
                {
                    gg.generate(v[j]);
                    square(tmp, v[j]);
                    n += tmp;
                }
                if (!isZero(n))
                    break;
            }
            // Divide by norm
            setOne(tmp);
            sqrt(n, n);
            tmp /= n;
            // Generate radius
            urng.randomUniform(n);
            n = power(n, arithmetic::convert(1, rc) / arithmetic::convert(d, rc));
            tmp *= n;
            // Do all scales
            v *= tmp;
            // Test wether v satisfies conditions
            setZero(n);
            bool is_inside = true;
            for (unsigned j = 0; j < v.size(); ++j)
            {
                square(tmp, v[j]);
                n += tmp;
                if (n > pruning[j])
                {
                    is_inside = false;
                    break;
                }
            }
            if (is_inside)
                ++inside;
        }
        std::cout << d << ": " << inside << "/" << total << " => " << 100.0 * inside / (double)total << "%\n";
        arithmetic::convert(res, inside, rc);
        arithmetic::convert(tmp, total, rc);
        res /= tmp;
        // Multiply res by the volume of an d-dimensional ball of radius 1.
        arithmetic::convert(tmp, (double)d / 2.0, rc);
        power(tmp, pi, tmp);
        res *= tmp;
        arithmetic::convert(tmp, (double)d / 2.0 + 1.0, rc);
        gamma(tmp, tmp);
        res /= tmp;
        return res;
    }
    if (d == 1)
    {
        // In dimension 1, everything is easy :)
        res = sqrt(pruning[0]);
        res <<= 1;
        return res;
    }
    
    if ((d & 1) == 0)
    {
        // The dimension is even.
        
        /*
          We use a hybrid Monte Carlo apporach.
          
          As a bounding volume, we increase pruning[2*i] or pruning[2*i+1] such that both are
          equal. Then, according to Gama, Nguyen, Regev: "Lattice Enumeration Using Extreme
          Pruning", the probability that a vector from the unit ball lies in the volume can be
          computed exactly. For this, note that the distribution of the vector (x[0]^2+x[1]^2,
          x[2]^2+x[3]^2, x[4]^2+x[5]^2, ...) is given by a Dirichlet distribution with parameters
          (1, ..., 1) (a vector of dimension d/2+1), which is a uniform distribution over the set of
          all vectors whose coordinates are non-negative and sum to at most 1.
          
          Thus we compute the probability that an (d/2)-dimensional vector x of non-negative
          coordinates which sum to at most 1 satisfies x[i] <= max(pruning[2*i],
          pruning[2*i+1]). The result is multiplied by the volume of a ball of dimension d: this
          yields the volume of the bounding set.
        */
        
        // Compute bounding in dimension d/2
        bool simple_function = true; // indicates that pruning[2*i] == pruning[2*i+1] for all i
        linalg::math_rowvector<arithmetic::Real> bounding;
        bounding.resize(d / 2);
        for (unsigned i = 0; i < bounding.size(); ++i)
            bounding[i].setContext(rc);
        for (unsigned i = bounding.size(); i > 0; --i)
        {
            if (pruning[2 * i - 2] != pruning[2 * i - 1])
            {
                simple_function = false;
                bounding[i - 1] = std::max(pruning[2 * i - 2], pruning[2 * i - 1]);
            }
            else
                bounding[i - 1] = pruning[2 * i - 2];
            // Make sure that bounding[] is non-decreasing
            if (i < bounding.size())
                if (bounding[i - 1] > bounding[i])
                    bounding[i - 1] = bounding[i];
        }
        
        // Compute probability that if x is sampled from a ball of radius 1, that x[0]^2+x[1]^2 <=
        // bounding[0], x[2]^2+x[3]^2 <= bounding[1], ...
        
        // Compute integrals
        linalg::base_rowvector<realpoly<1, arithmetic::Real> > integrals;
        arithmetic::Real tmp(rc), tmp2(rc);
        integrals.resize(bounding.size() + 1);
        integrals[0] = RealMonomial<arithmetic::Real>(0);
        for (unsigned i = 1; i < integrals.size(); ++i)
        {
            // Integrate
            integral<0>(integrals[i], integrals[i - 1]);
            // Evaluate the integral in bounding[bounding.size()-1]-bounding[bounding.size()-i]
            // using Horner's method. The negative of the result is stored in integrals[i][0].
            tmp = bounding[bounding.size() - 1] - bounding[bounding.size() - i];
            integrals[i](tmp).evaluate_to(tmp2);
            integrals[i][0] = -tmp2;
        }
        // Evaluate last integral at bounding[bounding.size()-1] using Horner's method. Evaluate this into res.
        integrals[integrals.size() - 1](bounding[bounding.size() - 1]).evaluate_to(res);
        // Multiply res by bounding.size()!.
        for (unsigned i = 2; i <= bounding.size(); ++i)
        {
            arithmetic::convert(tmp, i, rc);
            res *= tmp;
        }
        // Multiply res by the volume of an d-dimensional ball of radius 1.
        arithmetic::convert(tmp, (double)d / 2.0, rc);
        power(tmp, pi, tmp);
        res *= tmp;
        arithmetic::convert(tmp, (double)d / 2.0 + 1.0, rc);
        gamma(tmp, tmp);
        res /= tmp;
        // Now res can be multiplied by the probability obtained from the Monte Carlo algorithm.
        if (simple_function)
            // In case of simple functions, there is no need to run the Monte Carlo part...
            return res;
        
        // Apply Monte Carlo algorithm
        unsigned total = 1000;
        arithmetic::Real inside_p(rc);
        setZero(inside_p);
        linalg::math_rowvector<arithmetic::Real> t;
        t.resize(d/2);
        for (unsigned i = 0; i < t.size(); ++i)
            t[i].setContext(rc);
        arithmetic::RealContext::UniformRNG urng(rng);
        arithmetic::Real sumsofar(rc), val(rc);
        arithmetic::Real stopbound(rc), lower(rc);
        arithmetic::Real p(rc);
        for (unsigned i = 0; i < total; ++i)
        {
            while (true)
            {
                // Generate a uniformly sampled vector (t[0], ..., t[d/2-1]) which satisfies
                // \sum_{j=0}^i t[j] <= bounding[i] for i=0..d/2-1.
                sumsofar = bounding[bounding.size() - 1];
                bool restart = false;
                for (unsigned j = 0; j < t.size() - 1; ++j)
                {
                    // We want that val is distributed in [bounding[bounding.size()-1] - bounding[j],
                    // sumsofar] according to the distribution given by (a modified version of) the
                    // distribution function integrals[j] with (a modified version of) the density
                    // integrals[j-1].
                    
                    lower = bounding[bounding.size()-1] - bounding[j];
                    // Generate uniform random variable on [0, 1]
                    urng.randomUniform(t[j]);
                    // Evaluate integrals[j] at sumsofar and multiply with uniform random variable on [0,
                    // 1]. This is the destination.
                    integrals[integrals.size() - 1 - j](sumsofar).evaluate_to(tmp);
                    tmp *= t[j];
                    // Chose stopping bound
                    rc.getEpsilon(stopbound);
                    stopbound *= tmp;
                    abs(stopbound, stopbound);
                    stopbound <<= 4;
                    // Solve integrals[j](val) = tmp
                    val = sumsofar;
                    tmp2 = lower;
                    solve_regulafalsi(tmp2, val, integrals[integrals.size() - 1 - j], tmp, 10, rc);
                    val += tmp2;
                    val >>= 1;
                    solve_newton(val, integrals[integrals.size() - 1 - j], integrals[integrals.size() - 2 - j], tmp, stopbound, rc);
                    if ((val < lower) || (val > sumsofar))
                    {
                        // Apparently, we had numerical(?) problems
                        restart = true;
                        break;
                    }
//            std::cout << "solve " << integrals[integrals.size() - 1 - j] << " for " << tmp << " (error bound " << stopbound << ")\n";
                    // Now update sumsofar and t[j].
                    t[j] = sumsofar - val;
                    sumsofar = val;
                }
                if (restart)
                    continue;
                // Now generate t[t.size() - 1]
                urng.randomUniform(t[t.size() - 1]);
                t[t.size() - 1] *= sumsofar;
                // Now check whether a vector belonging to this choice satisfies the pruning conditions
                setZero(sumsofar);
                setOne(p);
                for (unsigned j = 0; j < t.size(); ++j)
                {
                    // We have to determine random x and y such that x^2 + y^2 == t[j], and then check
                    // the conditions, i.e. sumsofar + x^2 <= pruning[2 * j] and sumsofar + x^2 + y^2 <=
                    // pruning[2 * j + 1]. The latter condition is always satisfied since pruning[2 * j
                    // + 1] == bounding[j] and sumsofar + x^2 + y^2 == sumsofar + t[j]. The first
                    // condition is only in danger if sumsofar + t[j] is larger than pruning.
                    if (sumsofar + t[j] > pruning[2 * j])
                    {
                        computeCircularProbability(tmp, t[j], pruning[2 * j] - sumsofar, pi, rc);
                        p *= tmp;
                    }
                    // Increase sumsofar
                    sumsofar += t[j];
                    if (sumsofar > pruning[2 * j + 1])
                    {
                        // This should not happen, only due to numerical problems...
                        restart = true;
                        break;
                    }
                }
                if (restart)
                    continue;
                inside_p += p;
                break;
            }
        }
//        std::cout << inside_p << "/" << total << " -> " << inside_p / arithmetic::convert(total, rc) << "\n";
        // Multiply res with inside/total.
        res *= inside_p;
        arithmetic::convert(tmp, total, rc);
        res /= tmp;
    }
    else
    {
        // The dimension d is odd.
        
        // Compute bounding in dimension d/2
        bool simple_function = true; // indicates that pruning[2*i] == pruning[2*i+1] for all i
        linalg::math_rowvector<arithmetic::Real> bounding;
        bounding.resize(d / 2);
        for (unsigned i = 0; i < bounding.size(); ++i)
            bounding[i].setContext(rc);
        for (unsigned i = bounding.size(); i > 0; --i)
        {
            if (pruning[2 * i - 2] != pruning[2 * i - 1])
            {
                simple_function = false;
                bounding[i - 1] = std::max(pruning[2 * i - 2], pruning[2 * i - 1]);
            }
            else
                bounding[i - 1] = pruning[2 * i - 2];
            // Make sure that bounding[] is non-decreasing
            if (i < bounding.size())
                if (bounding[i - 1] > bounding[i])
                    bounding[i - 1] = bounding[i];
        }
        
        // Compute probability that if x is sampled from a ball of radius 1, that x[0]^2+x[1]^2 <=
        // bounding[0], x[2]^2+x[3]^2 <= bounding[1], ..., and finally x[d-1] <= pruning[d - 1].
        
        // Compute integrals
        linalg::base_rowvector<realpoly<2, arithmetic::Real> > integrals;
        linalg::base_rowvector<realpoly<1, arithmetic::Real> > integrals2;
        realpoly<1, arithmetic::Real> tmpp, tmpp2, tmpp5;
        realpoly<2, arithmetic::Real> tmpp3, tmpp4;
        // Variable 0 is the variable by which we integrate, and 1 is the last bounding value which
        // is variable in our case. The "usual" polynomial is obtained by substituting both
        // variables by the same one.
        arithmetic::Real tmp(rc), tmp2(rc);
        integrals.resize(bounding.size() + 1);
        integrals2.resize(bounding.size() + 1);
        // On integrals.row(k), we have the integrals where bounding is clipped to bounding[k].
        for (unsigned k = 0; k <= bounding.size(); ++k)
        {
            realpoly<2, arithmetic::Real> current, prev;
            current = RealMonomial<arithmetic::Real>(0, 0); // constant 1
            for (unsigned i = 1; i <= bounding.size(); ++i)
            {
                current.swap(prev);
                // Integrate
                integral<0>(current, prev);
                // The lower integration bound is 0 for bounding.size()-i >= k and
                // bounding[bounding.size()-i] minus the second variable for bounding.size()-i < k.
                if (bounding.size() - i < k)
                {
                    tmpp = RealMonomial<arithmetic::Real>(1);
                    tmpp -= bounding[bounding.size() - i];
                    // Instead of doing
                    //     tmpp2 = RealMonomial<arithmetic::Real>(1);
                    //     current(tmpp)(tmpp2).evaluate_to(tmpp5);
                    //     current[0] -= tmpp5;
                    // we proceed as follows:
                    if (current.degree() >= 0)
                    {
                        tmpp2 = current[current.degree()];
                        for (int i = current.degree() - 1; i >= 0; --i)
                        {
                            tmpp2 *= tmpp;
                            tmpp2 += current[i];
                        }
                        current[0] -= tmpp2;
                    }
                }
            }
            integrals[k].swap(current);
            make_univariate(integrals2[k], integrals[k]);
//            std::cout << k << ": " << integrals[k] << "\n";
//            std::cout << k << ": " << integrals2[k] << "\n";
        }
        // Now let us compute the total volume. We first compute the volume divided by the
        // (d-1)-dimensional volume of a ball. We first integrate by x[d-1], then (at once) by x[0],
        // ..., x[d-2]. The latter part is even-dimensional, whence the volume can be computed
        // directly with the precomputation similarly to above.
        setZero(res);
        // Consider the intervals x[d-1]^2 in [pruning[d-1]-bounding[k],
        // pruning[d-1]-bounding[k-1]] for k = 1 up to bounding.size(), where we assume
        // bounding[bounding.size()] == 0. For this, we need integrals.row(k). Finally, consider the
        // interval x[d-1]^2 in [pruning[d-1]-bounding[0], pruning[d-1]]. For this, we need
        // integrals.row(0).
        linalg::base_rowvector<realpoly<1, arithmetic::Real> > integrals3, integrals3_deriv;
        linalg::base_rowvector<arithmetic::Real> values;
        linalg::base_rowvector<std::pair<arithmetic::Real, arithmetic::Real> > bounds;
        integrals3.resize(bounding.size() + 1);
        integrals3_deriv.resize(bounding.size() + 1);
        values.resize(bounding.size() + 1);
        bounds.resize(bounding.size() + 1);
        for (unsigned j = 0; j < values.size(); ++j)
        {
            values[j].setContext(rc);
            bounds[j].first.setContext(rc);
            bounds[j].second.setContext(rc);
        }
        arithmetic::Real x1s(rc), x2s(rc), x1(rc), x2(rc);
        for (int k = bounding.size(); k >= 0; --k)
        {
            // First, determine squares of integration bounds. The lower bound is
            // sqrt(pruning[d-1]-bounding[k]) for k < bounding.size(), and 0 otherwise.
            // otherwise.
            if ((unsigned)k == bounding.size())
                setZero(x1s);
            else
                x1s = pruning[d - 1] - bounding[k];
            // The upper bound is sqrt(pruning[d-1]-bounding[k-1]) in case k > 0, and
            // sqrt(pruning[d-1]) for k == 0.
            if (k > 0)
                x2s = pruning[d - 1] - bounding[k - 1];
            else
                x2s = pruning[d - 1];
            // Compute x's
            sqrt(x1, x1s);
            sqrt(x2, x2s);
            bounds[k].first = x1;
            bounds[k].second = x2;
            // Check whether the integral part we add will be zero anyway...
            if (x2 != x1)
            {
                // This happens precisely if the integration bounds are not equal
                
                // Evaluate integrals2[k] * bounding.size()! at pruning[d-1] - x[d-1]^2. Then, we
                // integrate this.
                // Substitute and integrate
                tmpp = -RealMonomial<arithmetic::Real>(2);
                tmpp += pruning[d - 1];
                integrals2[k](tmpp).evaluate_to(integrals3_deriv[k]);
                integral<0>(integrals3[k], integrals3_deriv[k]);
                // Evaluate integral function at x1 and x2 and take difference: this is the integral
                integrals3[k](x2).evaluate_to(tmp);
                integrals3[k](x1).evaluate_to(tmp2);
                tmp -= tmp2;
                integrals3[k][0] -= tmp2;
                // Add result to res
                res += tmp;
            }
            // Update values list
            values[k] = res;
        }
        // Multiply res by bounding.size()!.
        for (unsigned i = 2; i <= bounding.size(); ++i)
        {
            arithmetic::convert(tmp, i, rc);
            res *= tmp;
        }
//        std::cout << "values = " << values << "\n";
//        std::cout << "bounds = " << bounds << "\n";
        
        // Now multiply the overall result by the volume of an (d-1)-dimensional ball and then by two.
        arithmetic::convert(tmp, (double)(d - 1) / 2.0, rc);
        power(tmp, pi, tmp);
        res *= tmp;
        tmp2 = tmp;
        arithmetic::convert(tmp, (double)(d - 1) / 2.0 + 1.0, rc);
        gamma(tmp, tmp);
        res /= tmp;
        tmp2 /= tmp;
        res <<= 1;
        tmp2 <<= 1;
        
//        std::cout << "Result = " << res << "\n";
        
        // Now res can be multiplied by the probability obtained from the Monte Carlo algorithm.
        if (simple_function)
            // In case of simple functions, there is no need to run the Monte Carlo part...
            return res;
        
        // Apply Monte Carlo algorithm
        unsigned total = 2000;
        arithmetic::Real inside_p(rc);
        setZero(inside_p);
        linalg::math_rowvector<arithmetic::Real> t;
        t.resize(d/2);
        for (unsigned i = 0; i < t.size(); ++i)
            t[i].setContext(rc);
        arithmetic::RealContext::UniformRNG urng(rng);
        arithmetic::Real sumsofar(rc), val(rc), sumsofar2(rc);
        arithmetic::Real stopbound(rc);
        arithmetic::Real p(rc), lastx(rc), lower(rc);
        for (unsigned i = 0; i < total; ++i)
        {
            // Generate a uniformly sampled vector (t[0], ..., t[d/2-1]) which satisfies
            // \sum_{j=0}^i t[j] <= bounding[i] for i=0..d/2-1.

            while (true)
            {
                // We begin by generating the last coordinate lastx. We generate a non-negative
                // coordinate since the distribution is symmetric with respect to the sign of the last
                // coordinate.
                urng.randomUniform(tmp);
                tmp *= values[0];
//            std::cout << "bait = " << tmp << ", ";
                // Look up in values[] to decide for k.
                unsigned k = bounding.size();
                while ((k > 0) && (values[k] < tmp))
                    --k;
//            std::cout << "k = " << k << ", ";
                // Decrease tmp by values[k + 1]
                if (k + 1 < values.size())
                    tmp -= values[k + 1];
                // Recover lastx using tmp, by solving poly(lastx) = tmp.
                rc.getEpsilon(stopbound);
                stopbound *= tmp;
                stopbound <<= 4;
                tmp2 = bounds[k].first;
                lastx = bounds[k].second;
                solve_regulafalsi(tmp2, lastx, integrals3[k], tmp, 10, rc);
                lastx += tmp2;
                lastx >>= 1;
                solve_newton(lastx, integrals3[k], integrals3_deriv[k], tmp, stopbound, rc);
//            std::cout << "lastx = " << lastx << ", ";
                
                // Now we are left to generate the other coordinates. We proceed as in the even case,
                // but this time we have to compute the polynomials on the fly and cannot precompute
                // them. First generate the polynomial describing the distribution of the inverse of the
                // first coordinate.
                square(tmp, lastx);
                sumsofar = pruning[d - 1] - tmp;
                if (sumsofar > bounding[bounding.size() - 1])
                    sumsofar = bounding[bounding.size() - 1];
                sumsofar2 = sumsofar;
                realpoly<1, arithmetic::Real> poly, deriv;
                // Instead of 
                //     tmpp = RealMonomial<arithmetic::Real>(1);
                //     setZero(tmpp2);
                //     tmpp2[0] = sumsofar;
                //     integrals[k](tmpp)(tmpp2).evaluate_to(poly);
                // we proceed with the following for loop:
                for (int j = integrals[k].degree(); j >= 0; --j)
                    integrals[k][j](sumsofar).evaluate_to(poly[j]);
                derivative<0>(deriv, poly);
//            std::cout << "poly = " << poly << ", ";
                bool restart = false;
                for (unsigned j = 0; j < t.size() - 1; ++j)
                {
                    // We want that val is distributed in [sumsofar2 - bounding[j], sumsofar] for j < k
                    // and [0, sumsofar] for j >= k according to the distribution given by (a modified
                    // version of) the distribution function poly with (a modified version of) the
                    // density deriv.
                    
                    if (j > 0)
                    {
                        // Compute next polynomial on the list
                        poly.swap(deriv);
                        derivative<0>(deriv, poly);
                    }
                    // Compute lower bound
                    if (j < k)
                        lower = sumsofar2 - bounding[j];
                    else
                        setZero(lower);
                    // Generate uniform random variable on [0, 1]
                    urng.randomUniform(t[j]);
                    // Evaluate poly at sumsofar and multiply with uniform random variable on [0,
                    // 1]. This is the destination.
                    poly(sumsofar).evaluate_to(tmp);
                    t[j] *= tmp;
                    // Chose stopping bound
                    rc.getEpsilon(stopbound);
                    stopbound *= t[j];
                    abs(stopbound, stopbound);
                    stopbound <<= 4;
                    // Solve poly(val) = t[j]
                    val = sumsofar;
                    tmp2 = lower;
                    solve_regulafalsi(tmp2, val, poly, t[j], 10, rc);
                    val += tmp2;
                    val >>= 1;
                    solve_newton(val, poly, deriv, t[j], stopbound, rc);
                    if ((val < lower) || (val > sumsofar))
                    {
                        // Apparently, we had numerical(?) problems
                        restart = true;
                        break;
                    }
                    // Now update sumsofar and t[j].
                    t[j] = sumsofar - val;
                    sumsofar = val;
                }
                if (restart)
                    continue;
                // Now generate t[t.size() - 1]
                urng.randomUniform(t[t.size() - 1]);
                t[t.size() - 1] *= sumsofar;
//            std::cout << t << " -> ";
                // Now check whether a vector belonging to this choice satisfies the pruning conditions
                setZero(sumsofar);
                setOne(p);
                for (unsigned j = 0; j < t.size(); ++j)
                {
                    // We have to determine random x and y such that x^2 + y^2 == t[j], and then check
                    // the conditions, i.e. sumsofar + x^2 <= pruning[2 * j] and sumsofar + x^2 + y^2 <=
                    // pruning[2 * j + 1]. The latter condition is always satisfied since pruning[2 * j
                    // + 1] == bounding[j] and sumsofar + x^2 + y^2 == sumsofar + t[j]. The first
                    // condition is only in danger if sumsofar + t[j] is larger than pruning.
                    if (sumsofar + t[j] > pruning[2 * j])
                    {
                        computeCircularProbability(tmp, t[j], pruning[2 * j] - sumsofar, pi, rc);
//                    if (tmp != tmp)
//                    {
//                        std::cout << i << " " << j << " " << tmp << " " << t[j] << " " << pruning[2 * j] - sumsofar << " " << pi << "\n";
//                        std::cout << t << "\n";
//                    }
                        p *= tmp;
                    }
                    // Increase sumsofar
                    sumsofar += t[j];
                    if (sumsofar > pruning[2 * j + 1])
                    {
                        // This should not happen, only due to numerical problems...
                        restart = true;
                        break;
                    }
                }
                if (restart)
                    continue;
//            std::cout << p << "\n";
                inside_p += p;
                break;
            }
        }
//        std::cout << inside_p << "/" << total << " -> " << inside_p / arithmetic::convert(total, rc) << "\n";
        // Multiply res with inside/total.
        res *= inside_p;
        arithmetic::convert(tmp, total, rc);
        res /= tmp;
//        std::cout << "product is " << res << "\n";
    }
    // Return result
    return res;
}

class GSA
{
private:
    double d_q;
    
public:
    GSA(double q)
        : d_q(sqrt(q))
    {
    }
    
    inline void getNorm(arithmetic::Real & dest, unsigned i, arithmetic::RealContext & rc) const
    {
        arithmetic::convert(dest, d_q, rc);
        power(dest, dest, (unsigned long)i);
    }
};

class QVector
{
private:
    const linalg::math_rowvector<long double> & d_qs;
    
public:
    QVector(const linalg::math_rowvector<long double> & qs)
        : d_qs(qs)
    {
    }
    
    inline void getNorm(arithmetic::Real & dest, unsigned i, arithmetic::RealContext & rc) const
    {
        arithmetic::convert(dest, sqrt(d_qs[i]), rc);
    }
};

template<class GSComputer>
class computeRadiusGH
{
public:
    void operator ()(arithmetic::Real & radius, unsigned dimension, arithmetic::RealContext & rc, const GSComputer & gs) const
    {
        arithmetic::Real tmp(rc);
        // Compute gamma function
        arithmetic::convert(radius, dimension, rc);
        radius >>= 1;
        ++radius;
        radius = gamma(radius);
        // Multiply with determinant
        for (unsigned i = 0; i < dimension; ++i)
        {
            gs.getNorm(tmp, i, rc);
            radius *= tmp;
        }
        // Take n-th root
        setOne(tmp);
        arithmetic::Real tmp2(rc);
        arithmetic::convert(tmp2, dimension, rc);
        tmp /= tmp2;
        power(radius, radius, tmp);
        // Divide by sqrt(pi)
        rc.getPi(tmp);
        tmp = sqrt(tmp);
        radius /= tmp;
//        std::cout << "Gaussian heuristics enumeration radius: " << radius << "\n";
    }
};

template<class GSComputer, template<typename> class RadiusComputer>
void computeNodes(linalg::math_rowvector<arithmetic::Real> & destination, const linalg::math_rowvector<arithmetic::Real> & pruning, bool showProgress, arithmetic::RandomNumberGenerator & rng, arithmetic::RealContext & rc, const GSComputer & gs, bool montecarlo = false, const RadiusComputer<GSComputer> & radcomp = RadiusComputer<GSComputer>())
/*
  Computes 1/2 * \sum_{d=0}^{dim-1} R^d * V(d) / prod_{i=dim-1-d}^{dim-1} ||b_i^*||, where V(d) is
  the volume of
      { x \in \R^d \mid \forall j = 0, \dots, d : \sum_{\ell=0}^j x_\ell^2 \le pruning[j] }.
  
  The factor alpha should be such that ||b_i^*|| = \alpha^i ||b_1||. It usually depends more on the
  lattice reduction used than on the lattice.
 */
{
    // Rescale
    destination.resize(pruning.size());
    for (unsigned i = 0; i < destination.size(); ++i)
    {
        destination[i].setContext(rc);
        setZero(destination[i]);
    }
    arithmetic::Real radius(rc), tmp(rc);
    radcomp(radius, pruning.size(), rc, gs);
    for (unsigned i = 1; i <= pruning.size(); ++i)
    {
        if (showProgress)
            std::cout << i << "/" << pruning.size() << "        \r" << std::flush;
        
        // Compute number of nodes
        destination[i - 1] = computeVolume(pruning, i, rng, rc, montecarlo);
//        std::cout << i << "-dimensional volume: " << destination[i - 1] << "\n";
        // Multiply with radius
        arithmetic::convert(tmp, (double)i * 0.5, rc);
        tmp = power(pruning[i - 1], tmp); // take power of square root
        destination[i - 1] *= tmp;
        tmp = power(radius, (unsigned long)i);
        destination[i - 1] *= tmp;
        // Divide by determinant of projected sublattice
        for (unsigned j = pruning.size() - i; j < pruning.size(); ++j)
        {
            gs.getNorm(tmp, j, rc);
            destination[i - 1] /= tmp;
        }
    }
    if (showProgress)
        std::cout << "                                       \r" << std::flush;
}

template<class GSComputer, template<typename> class RadiusComputer>
void computeNodes(TaskManager & tm, linalg::math_rowvector<arithmetic::Real> & destination, const linalg::math_rowvector<arithmetic::Real> & pruning, linalg::base_rowvector<arithmetic::RandomNumberGenerator> & rngs, arithmetic::RealContext & rc, const GSComputer & gs, bool montecarlo = false, const RadiusComputer<GSComputer> & radcomp = RadiusComputer<GSComputer>())
/*
  Computes 1/2 * \sum_{d=0}^{dim-1} R^d * V(d) / prod_{i=dim-1-d}^{dim-1} ||b_i^*||, where V(d) is
  the volume of
      { x \in \R^d \mid \forall j = 0, \dots, d : \sum_{\ell=0}^j x_\ell^2 \le pruning[j] }.
  
  The factor alpha should be such that ||b_i^*|| = \alpha^i ||b_1||. It usually depends more on the
  lattice reduction used than on the lattice.
 */
{
    class Job : public TaskManager::Job
    {
    private:
        arithmetic::Real & d_destination;
        unsigned d_dimension;
        const linalg::math_rowvector<arithmetic::Real> & d_pruning;
        linalg::base_rowvector<arithmetic::RandomNumberGenerator> & d_rngs;
        arithmetic::RealContext & d_rc;
        bool d_montecarlo;
        
    public:
        Job(arithmetic::Real & destination, unsigned dimension, const linalg::math_rowvector<arithmetic::Real> & pruning, linalg::base_rowvector<arithmetic::RandomNumberGenerator> & rngs, arithmetic::RealContext & rc, bool montecarlo)
            : d_destination(destination), d_dimension(dimension), d_pruning(pruning), d_rngs(rngs), d_rc(rc), d_montecarlo(montecarlo)
        {
        }
        
        virtual ~Job()
        {
        }
        
        virtual void run(unsigned thread_no, TaskManager * tm)
        {
            d_destination = computeVolume(d_pruning, d_dimension, d_rngs[thread_no], d_rc, d_montecarlo);
            std::cout << d_dimension << "-dimensional volume: " << d_destination << "\n";
        }
    };
    
    // Rescale
    destination.resize(pruning.size());
    for (unsigned i = 0; i < destination.size(); ++i)
    {
        destination[i].setContext(rc);
        setZero(destination[i]);
    }
    arithmetic::Real radius(rc), tmp(rc);
    radcomp(radius, pruning.size(), rc, gs);
    for (unsigned i = 1; i <= pruning.size(); ++i)
//    for (unsigned i = 35; i <= 45; ++i)
        tm.enqueueJob(new Job(destination[i - 1], i, pruning, rngs, rc, montecarlo));
    tm.waitForDone();
    for (unsigned i = 1; i <= pruning.size(); ++i)
    {
//        std::cout << i << "-dimensional volume: " << destination[i - 1] << "\n";
        // Multiply with radius
        arithmetic::convert(tmp, (double)i * 0.5, rc);
        tmp = power(pruning[i - 1], tmp); // take power of square root
        destination[i - 1] *= tmp;
        tmp = power(radius, (unsigned long)i);
        destination[i - 1] *= tmp;
        // Divide by determinant of projected sublattice
        for (unsigned j = pruning.size() - i; j < pruning.size(); ++j)
        {
            gs.getNorm(tmp, j, rc);
            destination[i - 1] /= tmp;
        }
    }
}

void writePruning(const std::string & fn, const linalg::math_rowvector<arithmetic::Real> & pruning)
{
    std::ofstream f(fn.c_str());
    f << std::setprecision(30) << pruning << "\n";
}

template<class T>
class Buffer
{
private:
    T * d_data;
    unsigned d_size, d_used;
    
    void reserve(unsigned minadd)
    {
        unsigned minsize = d_used + minadd;
        // Make multiple of 256 = 2^8
        minsize += 0xFF;
        minsize &= ~0xFF;
        // Reserve enough memory
        T * d = new T[minsize];
        for (unsigned i = 0; i < d_used; ++i)
            d[i] = d_data[i];
        d_size = minsize;
        delete[] d_data;
        d_data = d;
    }
    
public:
    Buffer()
        : d_data(NULL), d_size(0), d_used(0)
    {
    }
    
    ~Buffer()
    {
        delete[] d_data;
    }
    
    const T * data() const
    {
        return d_data;
    }

    unsigned size() const
    {
        return d_used;
    }
    
    void clear()
    {
        d_used = 0;
    }
    
    void add(const T & c)
    {
        if (d_used + 1 >= d_size)
            reserve(1);
        d_data[d_used++] = c;
    }
    
    void add(const T * c, unsigned n)
    {
        if (d_used + n >= d_size)
            reserve(n);
        while (n)
        {
            d_data[d_used++] = *c++;
            --n;
        }
    }
};

inline bool iswhitespace(char c)
{
    return (c == ' ') || (c == 8) || (c == 10) || (c == 13);
}

inline bool readRealVectorFromStream(std::istream & s, linalg::math_rowvector<arithmetic::Real> & vector, arithmetic::RealContext & rc)
{
    char c;
    c = s.get();
    while (iswhitespace(c))
        c = s.get();
    int dim = -1;
    if (c == 'v')
    {
        c = s.get();
        if (c != 'e')
            return false;
        c = s.get();
        if (c != 'c')
            return false;
        c = s.get();
        if (c != '<')
            return false;
        c = s.get();
        if (!isdigit(c))
            return false;
        dim = 0;
        while (isdigit(c))
        {
            dim *= 10;
            dim += c - '0';
            c = s.get();
        }
        if (c != '>')
            return false;
        c = s.get();
    }
    if (c != '[')
        return false;
    c = s.get();
    std::list<arithmetic::Real> vals;
    Buffer<char> buf;
    while (s)
    {
        while (iswhitespace(c))
            c = s.get();
        if (c == ']')
            break;
        if (c == ',')
            c = s.get();
        while (iswhitespace(c))
            c = s.get();
        // remove number
        buf.clear();
        while (std::isdigit(c) || (c == '-') || (c == '+') || (c == 'e') || (c == 'E') || (c == '.'))
        {
            buf.add(c);
            c = s.get();
        }
        if (buf.size())
        {
            vals.push_back(arithmetic::Real());
            buf.add(0);
            vals.back().setContext(rc);
            arithmetic::convert(vals.back(), buf.data(), rc);
        }
        else
            return false;
    }
    if (dim >= 0)
        if ((unsigned)dim != vals.size())
            return false;
    
    vector.resize(vals.size());
    std::list<arithmetic::Real>::iterator it = vals.begin();
    for (unsigned i = 0; i < vals.size(); ++i, ++it)
    {
        vector[i].setContext(rc);
        vector[i] = *it;
    }
    return true;
}

bool readPruning(const std::string & fn, linalg::math_rowvector<arithmetic::Real> & pruning, arithmetic::RealContext & rc)
{
    std::ifstream f(fn.c_str());
    return readRealVectorFromStream(f, pruning, rc);
}

void parseVectorReal(linalg::math_rowvector<arithmetic::Real> & vec, const std::string & str, arithmetic::RealContext & rc)
{
    linalg::math_rowvector<arithmetic::Integer> output;
    std::istringstream s(str);
    if (!readRealVectorFromStream(s, vec, rc))
        vec.resize(0);
}

void encodeUTF8(std::ostream & s, unsigned c)
{
    if (c < 0x80)
        s << (char)c;
    else
    {
        unsigned bound = 0x80;
        unsigned add = 0, mask = 0x80;
        while (c >= bound)
        {
            bound <<= 6 - 1;
            mask >>= 1;
            mask |= 0x80;
            ++add;
        }
        s << (char)(mask | (c >> (6 * add)));
        while (add)
            s << (char)(0x80 | ((c >> (--add * 6)) & 0x3F));
    }
}

template<class Evaluator, class Scale>
void plot(const Evaluator & eval, const Scale & scale, unsigned columns, unsigned lines)
{
    std::ostringstream out;
    linalg::math_matrix<unsigned> screen;
    screen.resize(lines, columns);
    for (unsigned i = 0; i < screen.rows(); ++i)
        for (unsigned j = 0; j < screen.cols(); ++j)
            screen(i, j) = ' ';
    for (unsigned i = 0; i < screen.cols(); ++i)
    {
        long double y = (1 - scale.toScale(eval(i))) * (lines - 1);
        signed long yy = trunc(y);
        long double frac = y - yy;
        if (yy >= lines)
        {
            yy = 0;
            while (yy < screen.rows())
            {
                screen(yy, i) = 0x2588; // 2588 FULL BLOCK
                ++yy;
            }
        }
        else if (yy >= 0)
        {
            /*
            switch ((unsigned)(frac * 6))
            {
            case 5: screen(yy, i) = '.'; break;
            case 4: screen(yy, i) = ','; break;
            case 3: screen(yy, i) = ':'; break;
            case 2: screen(yy, i) = ':'; break;
            case 1: screen(yy, i) = '^'; break;
            case 0: screen(yy, i) = '\''; break;
            }
            */
            switch ((unsigned)trunc(frac * 8))
            {
            case 7: screen(yy, i) = ' '; break;
            case 6: screen(yy, i) = 0x2581; break;
            case 5: screen(yy, i) = 0x2582; break;
            case 4: screen(yy, i) = 0x2583; break;
            case 3: screen(yy, i) = 0x2584; break;
            case 2: screen(yy, i) = 0x2585; break;
            case 1: screen(yy, i) = 0x2586; break;
            case 0: screen(yy, i) = 0x2587; break;
            }
            while (yy < screen.rows() - 1)
            {
                ++yy;
                screen(yy, i) = 0x2588; // FULL BLOCK
            }
        }
    }
    out << std::setw(10) << "" << " ";
    encodeUTF8(out, 0x2554);
    for (unsigned j = 0; j < screen.cols(); ++j)
        encodeUTF8(out, 0x2550);
    encodeUTF8(out, 0x2557);
    out << "\n";
    for (unsigned i = 0; i < screen.rows(); ++i)
    {
        out << std::right << std::setprecision(3) << std::setw(10) << scale.fromScale(1 - (long double)i / (long double)(lines - 1)) << " ";
        encodeUTF8(out, 0x2551);
        for (unsigned j = 0; j < screen.cols(); ++j)
            encodeUTF8(out, screen(i, j));
        encodeUTF8(out, 0x2551);
        out << "\n";
    }
    out << std::setw(10) << "" << " ";
    encodeUTF8(out, 0x255A);
    unsigned j = 0;
    for (; (j + 1 < columns) && (j < 8); j += 2)
    {
        encodeUTF8(out, 0x2550);
        encodeUTF8(out, 0x256A);
    }
    for (; j + 3 < columns; j += 4)
    {
        encodeUTF8(out, 0x2550);
        encodeUTF8(out, 0x2550);
        encodeUTF8(out, 0x2550);
        encodeUTF8(out, 0x256A);
    }
    for (; j < screen.cols(); ++j)
        encodeUTF8(out, 0x2550);
    encodeUTF8(out, 0x255D);
    out << "\n";
    out << std::setw(10) << "" << "  ";
    unsigned i = 1;
    for (; i < std::min(8u, columns); i += 2)
    {
        out << std::setw(2) << std::right << i + 1;
    }
    i += 2;
    for (; i < columns; i += 4)
    {
        out << std::setw(4) << std::right << i + 1;
    }
    out << "\n";
    std::cout << out.str();
}

class LinearScale
// Returns values in [0, 1)
{
private:
    long double d_min, d_len;
    
public:
    LinearScale(long double min, long double max)
        : d_min(min), d_len(max - min)
    {
    }
    
    
    long double toScale(long double raw) const
    {
        return (raw - d_min) / d_len;
    }
    
    long double fromScale(long double scaled) const
    {
        return scaled * d_len + d_min;
    }
};

class LogarithmicScale
// Returns values in [0, 1)
{
private:
    long double d_min, d_len;
    
public:
    LogarithmicScale(long double min, long double max)
        : d_min(log(min)), d_len(log(max) - log(min))
    {
    }
    
    
    long double toScale(long double raw) const
    {
        return (log(raw) - d_min) / d_len;
    }
    
    long double fromScale(long double scaled) const
    {
        return exp(scaled * d_len + d_min);
    }
};

class VectorEval
{
private:
    const linalg::math_rowvector<arithmetic::Real> & d_v;
    bool d_zeroToNegative;
    
public:
    VectorEval(const linalg::math_rowvector<arithmetic::Real> & v, bool zeroToNegative = false)
        : d_v(v), d_zeroToNegative(zeroToNegative)
    {
    }
    
    long double operator() (unsigned i) const
    {
        long double v = arithmetic::convert<long double>(d_v[i]);
        if (d_zeroToNegative)
            if (v == 0)
                v = -1;
        return v;
    }
};

void plot(const linalg::math_rowvector<arithmetic::Real> & pruning, unsigned lines = 0, bool logarithmic = false)
{
    if (lines == 0)
        lines = 50;
    
    if (logarithmic)
    {
        long double min = 0, max = 0;
        bool first = true;
        for (unsigned i = 0; i < pruning.size(); ++i)
        {
            if (arithmetic::convert<long double>(pruning[i]) == 0)
                continue;
            if (first)
            {
                max = min = arithmetic::convert<long double>(pruning[i]);
                first = false;
            }
            if (arithmetic::convert<long double>(pruning[i]) < min)
                min = arithmetic::convert<long double>(pruning[i]);
            if (arithmetic::convert<long double>(pruning[i]) > max)
                max = arithmetic::convert<long double>(pruning[i]);
        }
        if (first)
        {
            min = 1;
            max = 2;
        }
        if (min <= 0)
            min = 0.00001;
        if (max <= 0)
            min = 0.0001;
        if (max < min)
            max = 1;
        plot(VectorEval(pruning, true), LogarithmicScale(min, max), pruning.size(), lines);
    }
    else
        plot(VectorEval(pruning), LinearScale(0, 1), pruning.size(), lines);
}

template<class GSComputer, template<typename> class RadiusComputer>
arithmetic::Real countNodes(TaskManager & tm, const linalg::math_rowvector<arithmetic::Real> & pruning, bool printAll, bool doPlot, linalg::base_rowvector<arithmetic::RandomNumberGenerator> & rngs, arithmetic::RealContext & rc, const GSComputer & gs, bool montecarlo = false, const RadiusComputer<GSComputer> & radcomp = RadiusComputer<GSComputer>())
{
    linalg::math_rowvector<arithmetic::Real> nodes;
    if (printAll)
        std::cout << "Computing number of nodes...\n";
    computeNodes(tm, nodes, pruning, rngs, rc, gs, montecarlo, radcomp);
    arithmetic::Real res(rc);
    setZero(res);
    for (unsigned i = 0; i < pruning.size(); ++i)
    {
        // Output
        if (printAll && !isZero(nodes[i]))
            std::cout << "  Dimension " << i+1 << ": " << nodes[i] << " nodes\n";
        
        // Add to counter
        res += nodes[i];
    }
    // Divide by 2 and return
    res >>= 1;
    if (doPlot)
        plot(nodes, 0, true);
    return res;
}

template<class GSComputer, template<typename> class RadiusComputer>
arithmetic::Real countNodes(const linalg::math_rowvector<arithmetic::Real> & pruning, bool printAll, bool doPlot, bool showProgress, arithmetic::RandomNumberGenerator & rng, arithmetic::RealContext & rc, const GSComputer & gs, bool montecarlo = false, const RadiusComputer<GSComputer> & radcomp = RadiusComputer<GSComputer>())
{
    linalg::math_rowvector<arithmetic::Real> nodes;
    if (printAll)
        std::cout << "Computing number of nodes...\n";
    computeNodes(nodes, pruning, showProgress, rng, rc, gs, montecarlo, radcomp);
    arithmetic::Real res(rc);
    setZero(res);
    for (unsigned i = 0; i < pruning.size(); ++i)
    {
        // Output
        if (printAll && !isZero(nodes[i]))
            std::cout << "  Dimension " << i+1 << ": " << nodes[i] << " nodes\n";
        
        // Add to counter
        res += nodes[i];
    }
    // Divide by 2 and return
    res >>= 1;
    if (doPlot)
        plot(nodes, 0, true);
    return res;
}

template<class GSComputer, template<typename> class RadiusComputer>
void computeValuation(TaskManager & tm, arithmetic::Real & valuation, const linalg::math_rowvector<arithmetic::Real> & pruning, const arithmetic::Real & reduction_cost,
                      linalg::base_rowvector<arithmetic::RandomNumberGenerator> & rngs, arithmetic::RealContext & rc, const GSComputer & gs, bool montecarlo = false,
                      const RadiusComputer<GSComputer> & radcomp = RadiusComputer<GSComputer>())
{
    class Job : public TaskManager::Job
    {
    private:
        arithmetic::Real & d_destination;
        const linalg::math_rowvector<arithmetic::Real> & d_pruning;
        linalg::base_rowvector<arithmetic::RandomNumberGenerator> & d_rngs;
        arithmetic::RealContext & d_rc;
        bool d_montecarlo;
        
    public:
        Job(arithmetic::Real & destination, const linalg::math_rowvector<arithmetic::Real> & pruning, linalg::base_rowvector<arithmetic::RandomNumberGenerator> & rngs, arithmetic::RealContext & rc, bool montecarlo)
            : d_destination(destination), d_pruning(pruning), d_rngs(rngs), d_rc(rc), d_montecarlo(montecarlo)
        {
        }
        
        virtual ~Job()
        {
        }
        
        virtual void run(unsigned thread_no, TaskManager * tm)
        {
            d_destination = computeProbability(d_pruning, d_rngs[thread_no], d_rc, d_montecarlo);
        }
    };
    
    arithmetic::Real p(rc), one(rc), n(rc), t(rc);
    tm.enqueueJob(new Job(p, pruning, rngs, rc, montecarlo));
    n = countNodes<GSComputer, RadiusComputer>(tm, pruning, false, false, rngs, rc, gs, montecarlo, radcomp); // q^2: 0.920 = LLL, 0.943 = BKZ-10, 0.947 = BKZ-20, 0.951 = BKZ-30, 0.952 = BKZ-40
    setOne(one);
    if (p < one)
        p = one;
    if (n < one)
        n = one;
    n += reduction_cost;
    valuation = n / p;
}

template<class GSComputer, template<typename> class RadiusComputer>
void optimize(TaskManager & tm, linalg::math_rowvector<arithmetic::Real> & pruning, const arithmetic::Real & reduction_cost, bool simple_curve,
              linalg::base_rowvector<arithmetic::RandomNumberGenerator> & rngs, arithmetic::RealContext & rc, const GSComputer & gs, bool montecarlo = false,
              const RadiusComputer<GSComputer> & radcomp = RadiusComputer<GSComputer>())
{
    if (!simple_curve)
    {
        std::cout << "ERROR! Optimizing only implemented for simplified curves at the moment!\n";
        return;
    }
    sanitize(pruning);
    if (pruning.size() < 2)
        return;
    linalg::math_rowvector<arithmetic::Real> best;
    arithmetic::Real best_valuation(rc);
    best.resize(pruning.size());
    for (unsigned i = 0; i < best.size(); ++i)
    {
        best[i].setContext(rc);
        best[i] = pruning[i];
    }
    computeValuation<GSComputer, RadiusComputer>(tm, best_valuation, best, reduction_cost, rngs, rc, gs, montecarlo, radcomp);

    std::cout << "Starting point: valuation = " << best_valuation << ", pruning = " << std::setprecision(10) << pruning << "\n";
    unsigned worse_counter = 0, no_increase_counter = 0;
    while (no_increase_counter < 10000)
    {
        // Try random change
        if (simple_curve)
        {
            unsigned i = rngs[0].random(pruning.size() / 2);
            arithmetic::Real v(rc), m(rc);
            rngs[0].randomUniform(v);
            v <<= 1;
            --v;
            v >>= 6;
//            std::cout << i << " " << v << "\n";
            setOne(m);
            m -= pruning[2 * i];
            if (m < pruning[2 * i])
                m = pruning[2 * i];
            v *= m;
            pruning[2 * i] += v;
            pruning[2 * i + 1] += v;
            for (unsigned j = 2 * i + 2; j < pruning.size(); ++j)
                if (pruning[j] < pruning[j - 1])
                    pruning[j] = pruning[j - 1];
        }
        
        // Compute valuation
        sanitize(pruning);
        if (simple_curve)
            setOne(pruning[pruning.size() - 2]);
        arithmetic::Real val(rc);
        computeValuation<GSComputer, RadiusComputer>(tm, val, pruning, reduction_cost, rngs, rc, gs, montecarlo, radcomp);
//        std::cout << val << " " << pruning << "\n";
        // Compare
        if (val < best_valuation)
        {
            std::cout << "Found new optimum: valuation = " << val << ", pruning = " << std::setprecision(10) << pruning << "\n";
            plot(pruning, 25);
            best_valuation = val;
            best = pruning;
            worse_counter = 0;
            no_increase_counter = 0;
        }
        else
        {
            ++no_increase_counter;
            if ((val > best_valuation) || (val != val))
                if ((++worse_counter == 100) || (val != val))
                {
                    // Go back to last optimium
                    pruning = best;
                    worse_counter = 0;
                }
            if ((no_increase_counter % 100) == 0)
                std::cout << "No increase for " << no_increase_counter << " steps...\n";
        }
    }
    
    pruning = best;
}

void polytest()
{
    poly<2, int> f = Monomial<int>(1, 2) + 3 * Monomial<int>(4, 5);
    std::cout << f << "\n";
    int e = f(4)(2);
    std::cout << e << "\n";
    std::cout << derivative<1>(f) << "\n";
    poly<1, int> g = f(2), h = f(3);
    std::cout << g << " and " << h << "\n";
    poly<1, double> l = f(Monomial<double>(1) - 1)(Monomial<double>(1) - 2);
    std::cout << l << "\n";
    std::cout << derivative<0>(l) << "\n";
}

int main(int argc, char **argv)
{
    arithmetic::initArithmeticThreadAllocators();
    
//    polytest();
//    return 0;
    
    helper::ArgumentParser args(argc, argv);
    const helper::ArgumentParser::Value * v;
    
    if ((args.getValue("help") != NULL) || (args.getValue("h") != NULL))
    {
        help(argv[0]);
        return -1;
    }
    
    unsigned dimension = 120;
    
    if ((v = args.getValue("dimension")) != NULL)
    {
        if (v->isInteger() || (v->getInteger() > 1))
            dimension = v->getInteger();
        else
        {
            std::cerr << "Dimension must be an integer > 1!\n";
            return -1;
        }
    }
    
    arithmetic::RealContext & rc = arithmetic::getThreadRealContext();
    rc.setRealPrecision(256);
    
    linalg::math_rowvector<arithmetic::Real> pruning(dimension);
    arithmetic::Real tmp(rc);
    arithmetic::convert(tmp, dimension, rc);
    for (unsigned i = 0; i < dimension; ++i)
        setOne(pruning[i]);
    
    if ((v = args.getValue("read")) != NULL)
    {
        if (v->isEmpty())
        {
            std::cerr << "Error: argument -read needs input file name (-read=blah)\n";
            return -1;
        }
        std::cout << "Reading '" << v->getText() << "'...\n";
        if (!readPruning(v->getText(), pruning, rc))
        {
            std::cout << "ERROR! Reading failed!\n";
            return -1;
        }
    }
    if ((v = args.getValue("parse")) != NULL)
    {
        linalg::math_rowvector<arithmetic::Real> vec;
        parseVectorReal(vec, v->getText(), rc);
        if (vec.size() == 0)
        {
            std::cerr << "Cannot parse vector!\n";
            return -1;
        }
        dimension = vec.size();
        for (unsigned i = 0; i < dimension; ++i)
            pruning[i] = vec[i];
    }
    
    if ((v = args.getValue("linear")) != NULL)
    {
        std::cout << "Generating linear curve\n";
        for (unsigned i = 0; i < dimension; ++i)
        {
            pruning[i].setContext(rc);
            arithmetic::convert(pruning[i], i + 1, rc);
            pruning[i] /= tmp;
        }
    }

    if ((v = args.getValue("gnrpoly")) != NULL)
    {
        std::cout << "Generating curve from GNR polynomial\n";
        realpoly<1, arithmetic::Real> f;
        /* These values are from the paper. They yield a completely wrong curve which exceeds 1 quickly
           (at around 50) and reaches ~6.4 for 110!
           
           arithmetic::convert(f[0],  9.10e-4, rc);
           arithmetic::convert(f[1],  4.00e-2, rc);
           arithmetic::convert(f[2], -4.00e-3, rc);
           arithmetic::convert(f[3],  2.30e-4, rc);
           arithmetic::convert(f[4], -6.90e-6, rc);
           arithmetic::convert(f[5],  1.21e-7, rc);
           arithmetic::convert(f[6], -1.20e-9, rc);
           arithmetic::convert(f[7],  6.20e-12, rc);
           arithmetic::convert(f[8], -1.29e-14, rc);
        */
        /* These values are taken from the GPU implementation available at
           http://homes.esat.kuleuven.be/~jhermans/gpuenum/index.html */
        arithmetic::convert(f[0],  0.000914465, rc);
        arithmetic::convert(f[1],  0.0400812, rc);
        arithmetic::convert(f[2], -0.00424356, rc);
        arithmetic::convert(f[3],  0.00022931, rc);
        arithmetic::convert(f[4], -6.91288e-06, rc);
        arithmetic::convert(f[5],  1.21218e-07, rc);
        arithmetic::convert(f[6], -1.20165e-09, rc);
        arithmetic::convert(f[7],  6.20066e-12, rc);
        arithmetic::convert(f[8], -1.29185e-14, rc);
        
        arithmetic::Real x(rc), one(rc);
        setOne(one);
        for (unsigned i = 0; i < dimension; ++i)
        {
            // Determine evaluation point
            arithmetic::convert(x, (i + 1) * 110, rc);
            arithmetic::convert(tmp, dimension, rc);
            x /= tmp;
            // Evaluate
            pruning[i] = f(x);
            // Clip
            if (pruning[i] > one)
                pruning[i] = one;
        }
    }
    
    bool work_with_simple_curve = false;
    if ((v = args.getValue("simplify")) != NULL)
    {
        work_with_simple_curve = true;
        if ((v->getText() == "average") || (v->getText() == "avg"))
        {
            std::cout << "Simplifying curve (averaging)...\n";
            for (unsigned i = 0; i + 1 < pruning.size(); i += 2)
            {
                pruning[i] += pruning[i + 1];
                pruning[i] >>= 1;
                pruning[i + 1] = pruning[i];
            }
        }
        else
        {
            std::cout << "Simplifying curve...\n";
            for (unsigned i = 0; i + 1 < pruning.size(); i += 2)
                pruning[i] = pruning[i + 1];
        }
    }
        
    if (!sanitize(pruning))
        std::cout << "WARNING! Pruning curve was not sane!\n";
    
    if ((v = args.getValue("simplify")) != NULL)
        if (((pruning.size() & 1) == 0) && (pruning.size() > 0))
            setOne(pruning[pruning.size() - 2]);
    
    unsigned cores = TaskManager::getNumberOfCores();
    if ((v = args.getValue("maxthreads")) != NULL)
    {
        if (v->isInteger() || (v->getInteger() > 0))
            cores = v->getInteger();
        else
        {
            std::cerr << "Maximal thread number must be an integer > 0!\n";
            return -1;
        }
    }
    if (cores == 0)
        cores = 1;
    TaskManager tm(cores);
    linalg::base_rowvector<arithmetic::RandomNumberGenerator> rngs;
    rngs.resize(cores);
    for (unsigned i = 0; i < rngs.size(); ++i)
        rngs[i].randomizeTime();
    
    if ((v = args.getValue("randomize")) != NULL)
    {
        for (unsigned i = 0; i < rngs.size(); ++i)
            rngs[i].randomizeSeed();
    }
    
    linalg::math_rowvector<long double> qs;
    if ((v = args.getValue("readqs")) != NULL)
    {
        unsigned count = 0;
        qs.resize(dimension);
        for (unsigned i = 0; i < qs.size(); ++i)
            qs[i] = 0;
        const std::list<std::string> & names = args.getNames();
        for (std::list<std::string>::const_iterator namesit = names.begin(); namesit != names.end(); ++namesit)
        {
            std::ifstream f(namesit->c_str());
            linalg::math_matrix<arithmetic::Integer> lattice;
            f >> lattice;
            if (!f)
            {
                std::cout << "Cannot load lattice '" << *namesit << "'!\n";
                return -1;
            }
            if (lattice.rows() == dimension)
            {
                std::cout << "Using lattice " << *namesit << " to determine average q's...\n";
                // Compute normalized GS norms (i.e. divided through length of first vector)
                LatticeReduction lr(lattice);
                lr.setArithmetic(LatticeReduction::A_LongDouble);
                lr.setGramSchmidt(LatticeReduction::G_NumStable); // since the lattices are already reduced, this should be fine
                // lr.setGramSchmidt(LatticeReduction::G_ClassicInteger);
                lr.forceGSRebuild(true);
                qs[0] = lr.getGSSqNormLD(0);
                for (unsigned i = 1; i < qs.size(); ++i)
                    qs[i] += lr.getGSSqNormLD(i) / qs[0];
                qs[0] = ++count;
            }
        }
        if (count == 0)
            qs.resize(0);
        else
        {
            for (unsigned i = 0; i < qs.size(); ++i)
                qs[i] /= (long double)count;
            std::cout << "Average q's: " << qs << "\n";
        }
    }
    
    bool montecarlo = false;
    if ((v = args.getValue("montecarlo")) != NULL)
    {
        montecarlo = true;
    }
    bool onlyprob = false;
    if ((v = args.getValue("onlyprob")) != NULL)
    {
        onlyprob = true;
    }
    bool onlynodes = false;
    if ((v = args.getValue("onlynodes")) != NULL)
    {
        onlynodes = true;
    }
    
    if ((v = args.getValue("optimize")) != NULL)
    {
        std::cout << "Optimizing...\n";
        arithmetic::Real reduction_cost(rc);
        arithmetic::convert(reduction_cost, 10000000, rc);
        if (qs.size())
            optimize<QVector, computeRadiusGH>(tm, pruning, reduction_cost, work_with_simple_curve, rngs, rc, QVector(qs), montecarlo);
        else
            optimize<GSA, computeRadiusGH>(tm, pruning, reduction_cost, work_with_simple_curve, rngs, rc, GSA(0.943), montecarlo);
    }

    arithmetic::Real p(rc), prec(rc), precl2(rc);
    if (!onlynodes)
    {
        p = computeProbability(pruning, rngs[0], rc, montecarlo);
        setOne(prec);
        prec /= p;
        precl2 = log2(prec);
        std::cout << std::setprecision(10) << "Probability: " << p << " (reciprocal: " << prec << ", log_2 = " << precl2 << ")\n";
    }
    arithmetic::Real n(rc), nl2(rc);
    bool showNodes = false;
    if (!onlyprob)
    {
        if ((v = args.getValue("shownodes")) != NULL)
            showNodes = true;
        if (qs.size())
            n = countNodes<QVector, computeRadiusGH>(tm, pruning, showNodes, showNodes, rngs, rc, QVector(qs), montecarlo);
        else
            n = countNodes<GSA, computeRadiusGH>(tm, pruning, showNodes, showNodes, rngs, rc, GSA(0.943), montecarlo); // q^2: 0.920 = LLL, 0.943 = BKZ-10, 0.947 = BKZ-20, 0.951 = BKZ-30, 0.952 = BKZ-40
//        if (qs.size())
//            n = countNodes<QVector, computeRadiusGH>(pruning, showNodes, showNodes, true, rngs[0], rc, QVector(qs), montecarlo);
//        else
//            n = countNodes<GSA, computeRadiusGH>(pruning, showNodes, showNodes, true, rngs[0], rc, GSA(0.943), montecarlo); // q^2: 0.920 = LLL, 0.943 = BKZ-10, 0.947 = BKZ-20, 0.951 = BKZ-30, 0.952 = BKZ-40
        nl2 = log2(n);
        std::cout << "Nodes: ~ " << std::setprecision(10) << n << " (log_2 = " << nl2 << ")\n";
    }
    if (!onlyprob && !onlynodes)
        std::cout << "Total cost: 2^" << precl2 + nl2 << " = " << power(arithmetic::convert(2, rc), precl2 + nl2) << "\n";
    
    if ((v = args.getValue("dump")) != NULL)
    {
        std::cout << "Pruning curve:\n" << std::setprecision(15) << pruning << "\n";
    }
    if ((v = args.getValue("write")) != NULL)
    {
        if (v->isEmpty())
        {
            std::cerr << "Error: argument -write needs output file name (-write=blah)\n";
            return -1;
        }
        std::cout << "Storing as '" << v->getText() << "'...\n";
        writePruning(v->getText(), pruning);
    }
    if ((v = args.getValue("plot")) != NULL)
        plot(pruning, 40);
}

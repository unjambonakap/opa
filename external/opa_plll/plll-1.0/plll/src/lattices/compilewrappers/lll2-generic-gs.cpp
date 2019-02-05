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

#ifdef INTTYPECONTEXT
#ifdef REALTYPECONTEXT

#include "../lll2-internal.hpp"
#include <deque>
#include <cassert>
#ifndef PLLL_CONFIG_NO_GS_CLASSIC
  #include "../gramschmidt-classic.cpp"
#endif
#ifndef PLLL_CONFIG_NO_GS_CLASSICINT
  #include "../gramschmidt-classicint.cpp"
#endif
#ifndef PLLL_CONFIG_NO_GS_GIVENS
  #include "../gramschmidt-givens.cpp"
#endif
#ifndef PLLL_CONFIG_NO_GS_NUMSTABLE
  #include "../gramschmidt-numstable.cpp"
#endif
#include "../lattice.cpp"

namespace plll
{
    class CanDoGivens { };
    class CannotDoGivens { };
    
    template<class RealTypeContext, class IntTypeContext>
    GSInterface<RealTypeContext, IntTypeContext> * createGSInterface(Verbose & verbose, RealTypeContext & rc, IntTypeContext & ic,
                                                                     LatticeReduction::GramSchmidt gs, bool restarts,
                                                                     linalg::math_matrix<typename IntTypeContext::Integer> & A, CanDoGivens)
    {
        switch (gs)
        {
        default:
      #ifndef PLLL_CONFIG_NO_GS_NUMSTABLE
        case LatticeReduction::G_NumStable:
            return new DefaultGSI<RealTypeContext, IntTypeContext, GSNumStableLLL<RealTypeContext, IntTypeContext> >(verbose, rc, ic, A);
      #endif
      #ifndef PLLL_CONFIG_NO_GS_CLASSIC
        case LatticeReduction::G_Classic:
            return new DefaultGSI<RealTypeContext, IntTypeContext, GSClassic<RealTypeContext, IntTypeContext> >(verbose, rc, ic, A, restarts);
      #endif
      #ifndef PLLL_CONFIG_NO_GS_CLASSICINT
        case LatticeReduction::G_ClassicInteger:
            return new DefaultGSI<RealTypeContext, IntTypeContext, GSClassicInteger<RealTypeContext, IntTypeContext> >(verbose, rc, ic, A);
      #endif
      #ifndef PLLL_CONFIG_NO_GS_GIVENS
        case LatticeReduction::G_Givens:
            return new DefaultGSI<RealTypeContext, IntTypeContext, GSGivens<RealTypeContext, IntTypeContext> >(verbose, rc, ic, A, restarts);
      #endif
        }
        assert(!"Invalid GS type!");
        return NULL;
    }
    
    template<class RealTypeContext, class IntTypeContext>
    GSInterface<RealTypeContext, IntTypeContext> * createGSInterface(Verbose & verbose, RealTypeContext & rc, IntTypeContext & ic,
                                                                     LatticeReduction::GramSchmidt gs, bool restarts,
                                                                     linalg::math_matrix<typename IntTypeContext::Integer> & A, CannotDoGivens)
    {
        switch (gs)
        {
        default:
      #ifndef PLLL_CONFIG_NO_GS_NUMSTABLE
        case LatticeReduction::G_NumStable:
            return new DefaultGSI<RealTypeContext, IntTypeContext, GSNumStableLLL<RealTypeContext, IntTypeContext> >(verbose, rc, ic, A);
      #endif
      #ifndef PLLL_CONFIG_NO_GS_CLASSIC
        case LatticeReduction::G_Classic:
            return new DefaultGSI<RealTypeContext, IntTypeContext, GSClassic<RealTypeContext, IntTypeContext> >(verbose, rc, ic, A, restarts);
      #endif
      #ifndef PLLL_CONFIG_NO_GS_CLASSICINT
        case LatticeReduction::G_ClassicInteger:
            return new DefaultGSI<RealTypeContext, IntTypeContext, GSClassicInteger<RealTypeContext, IntTypeContext> >(verbose, rc, ic, A);
      #endif
      #ifndef PLLL_CONFIG_NO_GS_GIVENS
        case LatticeReduction::G_Givens:
            // When assertions are disabled, G_NumStable is used
            assert(!"Givens rotations are not available for this arithmetic!");
      #endif
        }
        assert(!"Invalid GS type!");
        return NULL;
    }
    
    template<class RealTypeContext, class IntTypeContext>
    GSInterface<RealTypeContext, IntTypeContext> * createGSInterface(Verbose & verbose, RealTypeContext & rc, IntTypeContext & ic,
                                                                     LatticeReduction::GramSchmidt gs, bool restarts,
                                                                     linalg::math_matrix<typename IntTypeContext::Integer> & A)
    {
        return createGSInterface<RealTypeContext, IntTypeContext>(verbose, rc, ic, gs, restarts, A,
                                                                  typename helper::SelectFirstType<RealTypeContext::has_squareroot,
                                                                                                   CanDoGivens, CannotDoGivens>::result());
    }
    
    // Instantiate for given arithmetic
    
    template GSInterface<REALTYPECONTEXT, INTTYPECONTEXT> * createGSInterface<REALTYPECONTEXT, INTTYPECONTEXT>(Verbose &, REALTYPECONTEXT &, INTTYPECONTEXT &, LatticeReduction::GramSchmidt,
                                                                                                               bool, linalg::math_matrix<INTTYPECONTEXT::Integer> &);
    
}
#endif
#endif

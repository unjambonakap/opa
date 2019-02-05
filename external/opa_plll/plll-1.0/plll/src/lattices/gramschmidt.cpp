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

#ifndef PLLL_INCLUDE_GUARD__GRAMSCHMIDT_CPP
#define PLLL_INCLUDE_GUARD__GRAMSCHMIDT_CPP

namespace plll
{
/* Interace for GSComputer:
    - Constructor wants arguments:
      1) RealTypeContext & - the real type context
      2) math_matrix<Integer> & A - reference to the matrix (which has to undergo
         the same modifications as the calls to swap(), add(), flip() and trans() indicate)
    - Assume that the rows of A are b_1, ..., b_n. Then the GS orthogonalized basis is b_1^*, ..., b_n^*.
    - Has following getters: (level is set by update(), see below)
      * getCoeff(unsigned i, unsigned j)
        - assumes that level > i > j
        - returns Gram-Schmidt coefficient r_{ij} = <b_i, b_j^*> / <b_j^*, b_j^*> if i > j
        - note that
                                  b_i^* = b_i - \sum_{j=1}^{i-1} r_{ij} b_j^*.
      * getNormSq(unsigned i)
        - assumes that level > i
        - returns ||b_i^*||_2^2 = r_{ii}^2
    - Has following modifiers:
      * Note that the matrix A has to be modified accordingly _before_ calling these modifiers
      * void swap(unsigned i, unsigned j);
        - swap vectors i and j.
        - CURRENTLY only implemented if |i - j| = 1
      * void add(unsigned i, const Integer & m, unsigned j);
        - add m-times the j-th vector to the i-th vector
        - CURRENTLY only implemented if i > j
      * void flip(unsigned i);
        - flip the signs of the i-th vector
      * void trans(unsigned i, unsigned j,
                   const Integer & B00, const Integer & B01, const Integer & B10, const Integer & B11);
        - [ B00 B01 ]   [ ..... < vector i > ..... ]
          [ B10 B11 ] * [ ..... < vector j > ..... ]
        - CURRENTLY not implemented
      * update(unsigned level)
        - extends the current level (if larger than previous level), i.e. applies Gram-Schmidt
          to next basis vectors
*/
}

// We provide four GS computers:

//  * GSGivens<RealTypeContext, IntTypeContext>
//  * GSClassic<RealTypeContext, IntTypeContext>
//  * GSClassicInteger<RealTypeContext, IntTypeContext>
//  * GSNumStableLLL<RealTypeContext, IntTypeContext>

#ifdef PLLL_INTERNAL_CONFIG_HAS_GIVENS
  #include "gramschmidt-givens.cpp"
#endif
#include "gramschmidt-classic.cpp"
#include "gramschmidt-classicint.cpp"
#include "gramschmidt-numstable.cpp"

#endif

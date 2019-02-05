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

#include "../enumimpl.cpp"
#include "../lll2-callback.cpp"

namespace plll
{
    // Explicit template instantiations:
    
    template class GaussianHeuristicEnumBoundSelection<REALTYPECONTEXT, INTTYPECONTEXT>;
    template void StandardEnumBoundSelection<REALTYPECONTEXT, INTTYPECONTEXT>(GaussianFactorComputer &, Lattice<REALTYPECONTEXT, INTTYPECONTEXT> &,
                                                                              REALTYPECONTEXT::Real &);
    template class Enumerator<REALTYPECONTEXT, INTTYPECONTEXT>;
    template Enumerator<REALTYPECONTEXT, INTTYPECONTEXT>::Enumerator(Workspace<REALTYPECONTEXT, INTTYPECONTEXT,
                                                                               Enumerator<REALTYPECONTEXT, INTTYPECONTEXT>,
                                                                               EmptyCallback<REALTYPECONTEXT, INTTYPECONTEXT> > &,
                                                                     unsigned, unsigned,
                                                                     LatticeReduction::EnumCallbackFunction, LatticeReduction::EnumCallbackFunction_LI);
    template Enumerator<REALTYPECONTEXT, INTTYPECONTEXT>::Enumerator(Workspace<REALTYPECONTEXT, INTTYPECONTEXT,
                                                                               Enumerator<REALTYPECONTEXT, INTTYPECONTEXT>,
                                                                               Callback<REALTYPECONTEXT, INTTYPECONTEXT> > &,
                                                                     unsigned, unsigned,
                                                                     LatticeReduction::EnumCallbackFunction, LatticeReduction::EnumCallbackFunction_LI);
}

#endif
#endif

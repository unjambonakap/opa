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

#include <plll/config.hpp>
#include "transform.cpp"
#include "matrixconversion.hpp"
#ifndef PLLL_CONFIG_NO_ARITHMETIC_LONGINT
  #include <plll/arithmetic-nint.hpp>
#endif

namespace plll
{
    class Transform::TransformImpl
    {
    private:
        Transform & d_T;
        
        #define MAKE_CREATE(CONTEXT, NAME)                                                                                                  \
                STD_AUTO_PTR<MatrixConversion<CONTEXT> > d_m_ ## NAME;                                                                     \
                STD_AUTO_PTR<TransformNotifier<CONTEXT> > d_T_ ## NAME;                                                                    \
                                                                                                                                            \
                TransformNotifier<CONTEXT> * createImpl(CONTEXT & ic)                                                                       \
                {                                                                                                                           \
                    d_m_ ## NAME .reset(new MatrixConversion<CONTEXT>(ic, *d_T.d_transform, true, true));                                   \
                    d_T_ ## NAME .reset(new TransformData<CONTEXT>(&d_m_ ## NAME ->matrix(), d_T.d_mode == LatticeReduction::T_Inverse));   \
                    return d_T_ ## NAME .get();                                                                                             \
                }
        #define MAKE_RELEASE(NAME)      \
                d_T_ ## NAME .reset();  \
                d_m_ ## NAME .reset();
        
        #ifndef PLLL_CONFIG_NO_ARITHMETIC_BIGINT
          MAKE_CREATE(arithmetic::IntegerContext, bigint)
        #endif
        #ifndef PLLL_CONFIG_NO_ARITHMETIC_LONGINT
          MAKE_CREATE(arithmetic::NIntContext<long>, longint)
        #endif
            
    public:
        void release()
        {
          #ifndef PLLL_CONFIG_NO_ARITHMETIC_BIGINT
            MAKE_RELEASE(bigint)
          #endif
          #ifndef PLLL_CONFIG_NO_ARITHMETIC_LONGINT
            MAKE_RELEASE(longint)
          #endif
        }
        
        template<class IntTypeContext>
        TransformNotifier<IntTypeContext> * create(IntTypeContext & ic)
        {
            return createImpl(ic);
        }
        
        TransformImpl(Transform & T)
            : d_T(T)
        {
        }
    };
    
    Transform::Transform()
        : d_mode(LatticeReduction::T_Normal), d_transform(), d_impl(new TransformImpl(*this))
    {
    }
    
    Transform::~Transform()
    {
    }
    
    void Transform::reset(const linalg::math_matrix<arithmetic::Integer> & A) // inform that the lattice was changed
    {
        if (d_transform.get())
            enable(d_mode, A);
    }
    
    void Transform::enable(LatticeReduction::Transform mode, const linalg::math_matrix<arithmetic::Integer> & A)
    {
        disable();
        
        d_transform.reset(new linalg::math_matrix<arithmetic::Integer>(A.rows(), A.rows()));
        setUnit<arithmetic::Integer>(*d_transform);
        d_mode = mode;
    }
    
    void Transform::disable()
    {
        d_impl->release();
        d_transform.reset();
    }
    
    bool Transform::isEnabled() const
    {
        return d_transform.get() != NULL;
    }
    
    LatticeReduction::Transform Transform::mode() const
    {
        return d_mode;
    }
    
    const linalg::math_matrix<arithmetic::Integer> * Transform::transform() const
    {
        return d_transform.get();
    }
    
    template<class IntTypeContext>
    TransformNotifier<IntTypeContext> * Transform::createTransformObject(IntTypeContext & ic)
    {
        d_impl->release();
        if (d_transform.get() == NULL)
            return NULL;
        return d_impl->create(ic);
    }
    
    void Transform::releaseTransformObject()
    {
        d_impl->release();
    }
    
    #define SPECIALIZE(CONTEXT) template TransformNotifier<CONTEXT> * Transform::createTransformObject(CONTEXT &);
    
    #ifndef PLLL_CONFIG_NO_ARITHMETIC_BIGINT
      SPECIALIZE(arithmetic::IntegerContext)
    #endif
    #ifndef PLLL_CONFIG_NO_ARITHMETIC_LONGINT
      SPECIALIZE(arithmetic::NIntContext<long>)
    #endif
}

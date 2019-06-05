%include "opa.i"


%{
#include "opa/math/common/base.h"
#include "opa/math/common/bignum.h"
#include "opa/math/common/float.h"
#include "opa/math/common/stats.h"
using namespace opa::math::common;
%}

%{
#include "swig_math_common.h"
%}

%typemap(in) const opa::math::common::bignum& (opa::math::common::bignum tmp){
  $1 = &tmp;
  if (!opa::math::common::PyToBignum($input, $1)) SWIG_fail;
}

%typemap(in) opa::math::common::bignum {
  if (!opa::math::common::PyToBignum($input, &$1)) SWIG_fail;
}
%typemap(typecheck, precedence=SWIG_TYPECHECK_SWIGOBJECT) opa::math::common::bignum, const opa::math::common::bignum& {
  $1 = PyInt_Check($input);
}

%typemap(out) opa::math::common::bignum {
  $result=opa::math::common::PyFromBignum($1);
}

%ignore opa::math::common::operator<<;

%ignore opa::math::common::bitcount_tb;
%include "opa/math/common/base.h"
%include "opa/math/common/float.h"
%include "opa/math/common/bignum.h"
%include "opa/math/common/stats.h"

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

%typemap(out) opa::math::common::bignum {
  $result=opa::math::common::PyFromBignum($1);
}

%ignore opa::math::common::operator<<;

%ignore opa::math::common::bitcount_tb;
%include "opa/math/common/base.h"
%include "opa/math/common/float.h"
%include "opa/math/common/bignum.h"
%include "opa/math/common/stats.h"

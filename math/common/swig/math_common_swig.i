%module opa_math_common_swig

%include "math_common_swig_base.i"


%{
#include "opa/math/common/base.h"
#include "opa/math/common/bignum.h"
#include "opa/math/common/float.h"
#include "opa/math/common/Ring.h"
#include "opa/math/common/Poly.h"
#include "opa/math/common/PolyRing.h"
#include "opa/math/common/Field.h"
#include "opa/math/common/GF_p.h"
#include "opa/math/common/Utils.h"
#include "opa/math/common/Types.h"
#include "opa/math/common/Matrix.h"

using namespace opa::math::common;
%}
namespace opa{
namespace math{
namespace common{
}
}
}

%include "opa/math/common/Ring.h"
%template(Ring_u32) opa::math::common::Ring<u32>;
%include "opa/math/common/Field.h"
%template(Field_u32) opa::math::common::Field<u32>;
%include "opa/math/common/Poly.h"
%template(Poly_u32) opa::math::common::Poly<u32>;
%template(Ring_PU32) opa::math::common::Ring<opa::math::common::Poly<u32> >;
%include "opa/math/common/PolyRing.h"
%template(PolyRing_u32) opa::math::common::PolyRing<u32>;
%include "opa/math/common/GF_p.h"
%include "opa/math/common/Types.h"
%include "opa/math/common/Utils.h"

%nocopyctor;
%include "opa/math/common/Matrix.h"
%template(Matrix_u32) opa::math::common::Matrix<u32>;

%include "opa/math/common/stats.h"

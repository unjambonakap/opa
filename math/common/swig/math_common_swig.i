%module opa_math_common_swig


%include "math_common_swig_base.i"


%{
#include "opa/utils/misc.h"
#include "opa/math/common/base.h"
#include "opa/math/common/float.h"
#include "opa/math/common/Ring.h"
#include "opa/math/common/Poly.h"
#include "opa/math/common/PolyRing.h"
#include "opa/math/common/PolyModRing.h"
#include "opa/math/common/Field.h"
#include "opa/math/common/GF_p.h"
#include "opa/math/common/GF_q.h"
#include "opa/math/common/Utils.h"
#include "opa/math/common/Types.h"
#include "opa/math/common/Matrix.h"
#include "opa/math/common/CyclicCode.h"
#include "opa/math/common/BCHCode.h"
#include "opa/math/common/RSCode.h"
#include "opa/math/common/UtilsGFq.h"

using namespace opa::math::common;
%}
namespace opa{
namespace math{
namespace common{
}
}
}

%include "opa/utils/misc.h"
%include "opa/math/common/Ring.h"
%template(Ring_u32) opa::math::common::Ring<u32>;
%include "opa/math/common/Field.h"
%template(Field_u32) opa::math::common::Field<u32>;
%include "opa/math/common/Poly.h"
%template(Poly_u32) opa::math::common::Poly<u32>;
%template(Poly_Poly_u32) opa::math::common::Poly<opa::math::common::Poly<u32>>;
%template(Ring_PU32) opa::math::common::Ring<opa::math::common::Poly<u32> >;
%template(Ring_P_P_U32) opa::math::common::Ring<opa::math::common::Poly<opa::math::common::Poly<u32> >>;
%template(Field_P_P_U32) opa::math::common::Field<opa::math::common::Poly<opa::math::common::Poly<u32> >>;
%template(Field_PU32) opa::math::common::Field<opa::math::common::Poly<u32> >;
%include "opa/math/common/PolyRing.h"
%template(PolyRingOps_u32) opa::math::common::PolyRingOps<u32>;
%template(PolyRing_u32) opa::math::common::PolyRing<u32>;
%template(PolyRing_u32_F) opa::math::common::PolyRing<u32, opa::math::common::Field<opa::math::common::Poly<u32>>>;

%template(PolyRingOps_Poly_u32) opa::math::common::PolyRingOps<opa::math::common::Poly<u32>>;
%template(PolyRing_Poly_u32) opa::math::common::PolyRing<opa::math::common::Poly<u32>>;

%include "opa/math/common/PolyModRing.h"
%template(PolyModRing_u32) opa::math::common::PolyModRing<u32>;
%template(PolyModRing_u32_F) opa::math::common::PolyModRing<u32, opa::math::common::Field<opa::math::common::Poly<u32>>>;
%template(PolyModField_u32) opa::math::common::PolyModField<u32>;

%include "opa/math/common/GF_p.h"
%include "opa/math/common/GF_q.h"
%include "opa/math/common/Types.h"
%include "opa/math/common/Utils.h"
%template(GF_q_u32) opa::math::common::GF_q<u32>;

%nocopyctor;
%include "opa/math/common/Matrix.h"
%template(Matrix_u32) opa::math::common::Matrix<u32>;

%include "opa/math/common/CyclicCode.h"
%template(CyclicCode_u32) opa::math::common::CyclicCode<u32>;
%template(CyclicCode_P_u32) opa::math::common::CyclicCode<opa::math::common::Poly_u32>;
%include "opa/math/common/BCHCode.h"
%template(BCHCode_u32) opa::math::common::BCHCode<u32>;
%template(BCHCode_P_u32) opa::math::common::BCHCode<opa::math::common::Poly_u32>;

%include "opa/math/common/RSCode.h"
%template(RSCode_u32) opa::math::common::RSCode<u32>;
%template(RSCode_P_u32) opa::math::common::RSCode<opa::math::common::Poly_u32>;
%include "opa/math/common/stats.h"

%include "opa/math/common/UtilsGFq.h"
%template(pack_vector_gfq_u32) opa::math::common::pack_vector_gfq<u32>;
%template(unpack_vector_gfq_u32) opa::math::common::unpack_vector_gfq<u32>;

namespace std {
  %template(v_poly_u32) std::vector<opa::math::common::Poly<u32>>;
}

%module opa_math_common_swig

%{
#include "opa/predef.h"
#include "opa/utils/misc.h"
%}

%include "math_common_swig_base.i"
//%include <std_unique_ptr.i>
%include "opa_std_unique_ptr.i"


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
#include "opa/math/common/RealField.h"
#include "opa/math/common/Utils.h"
#include "opa/math/common/Types.h"
#include "opa/math/common/Matrix.h"
#include "opa/math/common/CyclicCode.h"
#include "opa/math/common/BCHCode.h"
#include "opa/math/common/RSCode.h"
#include "opa/math/common/UtilsGFq.h"
#include "opa/math/common/FFT.h"
#include "opa/math/common/crt_fft.h"

using namespace opa::math::common;
typedef opa::math::common::S64Factors S64Factors;
typedef opa::math::common::BGFactors BGFactors;
typedef opa::math::common::FFT2<std::array<u32, 3>> fft_a3_u32_t ;
%}
namespace opa{
namespace math{
namespace common{
}
}
}
typedef opa::math::common::S64Factors S64Factors;
typedef opa::math::common::BGFactors BGFactors;
%template(a_u32_3) std::array<u32, 3>;
%template(v_a_u32_3) std::vector<std::array<u32, 3>>;
%template(v_float) std::vector<opa::math::common::Float>;

%template(ps64_i) std::pair<s64, int>;
%template(pbg_i) std::pair<opa::math::common::bignum, int>;
%template(p_pu32_pu32) std::pair<opa::math::common::Poly<u32>, opa::math::common::Poly<u32>>;
%template(t_factors_s64) std::vector<std::pair<s64, int>>;
%template(t_factors_bg) std::vector<std::pair<opa::math::common::bignum, int>>;

%include "opa/utils/misc.h"
%include "opa/math/common/Ring.h"
%template(Ring_u32) opa::math::common::Ring<u32>;
%template(Ring_double) opa::math::common::Ring<double>;
%template(Ring_cdouble) opa::math::common::Ring<std::complex<double>>;
%template(Ring_Float) opa::math::common::Ring<opa::math::common::Float>;
%template(Ring_cFloat) opa::math::common::Ring<std::complex<opa::math::common::Float>>;
%include "opa/math/common/Field.h"
%template(Field_u32) opa::math::common::Field<u32>;
%template(Field_double) opa::math::common::Field<double>;
%template(Field_cdouble) opa::math::common::Field<std::complex<double>>;
%template(Field_Float) opa::math::common::Field<opa::math::common::Float>;
%template(Field_cFloat) opa::math::common::Field<std::complex<opa::math::common::Float>>;

%include "opa/math/common/Poly.h"
%template(Poly_u32) opa::math::common::Poly<u32>;
%template(Poly_double) opa::math::common::Poly<double>;
%template(Poly_float) opa::math::common::Poly<opa::math::common::Float>;
%template(Poly_Poly_u32) opa::math::common::Poly<opa::math::common::Poly<u32>>;
%template(Ring_PU32) opa::math::common::Ring<opa::math::common::Poly<u32> >;
%template(Field_PU32) opa::math::common::Field<opa::math::common::Poly<u32> >;
%template(Ring_PD) opa::math::common::Ring<opa::math::common::Poly<double> >;
%template(Field_PD) opa::math::common::Field<opa::math::common::Poly<double> >;

%template(Ring_P_P_U32) opa::math::common::Ring<opa::math::common::Poly<opa::math::common::Poly<u32> >>;
%template(Field_P_P_U32) opa::math::common::Field<opa::math::common::Poly<opa::math::common::Poly<u32> >>;
%template(Field_P) opa::math::common::Field<opa::math::common::Poly<u32> >;
%include "opa/math/common/PolyRing.h"
%template(PolyRingOps_u32) opa::math::common::PolyRingOps<u32>;
%template(PolyRing_u32) opa::math::common::PolyRing<u32>;
%template(PolyRing_u32_F) opa::math::common::PolyRing<u32, opa::math::common::Ring<u32>, opa::math::common::Field<opa::math::common::Poly<u32>>>;
%template(PolyRingOps_d) opa::math::common::PolyRingOps<double>;
%template(PolyRing_d) opa::math::common::PolyRing<double>;
%template(PolyRing_float) opa::math::common::PolyRing<opa::math::common::Float>;

%template(PolyRingOps_Poly_u32) opa::math::common::PolyRingOps<opa::math::common::Poly<u32>>;
%template(PolyRing_Poly_u32) opa::math::common::PolyRing<opa::math::common::Poly<u32>>;

%include "opa/math/common/PolyModRing.h"
%template(PolyModRing_u32) opa::math::common::PolyModRing<u32>;
%template(PolyModRing_u32_F) opa::math::common::PolyModRing<u32, opa::math::common::Ring<u32>, opa::math::common::Field<opa::math::common::Poly<u32>>>;
%template(PolyModField_u32) opa::math::common::PolyModField<u32>;

%template(GF_pT_u32) opa::math::common::GF_pT<u32>;

%include "opa/math/common/RealField.h"
%template(T_RealDouble) opa::math::common::FField<double, double>;
%template(T_RealFloat) opa::math::common::FField<opa::math::common::Float,opa::math::common::Float>;
%template(T_ComplexDouble) opa::math::common::FField<std::complex<double>, double>;
%template(T_ComplexFloat) opa::math::common::FField<std::complex<opa::math::common::Float>,opa::math::common::Float>;

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

%feature("novaluewrapper") std::unique_ptr<opa::math::common::FFTProvider<std::array<u32, 3>>>;

%include "opa/math/common/FFT.h"
%template(fftprovider_u32) opa::math::common::FFTProvider<u32>;
%template(fftprovider_d) opa::math::common::FFTProvider<double>;
%template(fftprovider_float) opa::math::common::FFTProvider<opa::math::common::Float>;
%template(fftprovider_a3_u32) opa::math::common::FFTProvider<std::array<u32, 3>>;
%template(fft_u32) opa::math::common::FFT2<u32>;
%template(fft_dispatcher_u32) opa::math::common::FFT2Dispatcher<u32>;
%template(fft_dispatcher_d) opa::math::common::FFT2Dispatcher<double>;
%template(fft_dispatcher_float) opa::math::common::FFT2Dispatcher<opa::math::common::Float>;
%template(fft_a3_u32) opa::math::common::FFT2<std::array<u32, 3>>;

%template(uptr_fft_a3_u32) std::unique_ptr<opa::math::common::FFTProvider<std::array<u32, 3>>>;
%newobject std::unique_ptr<opa::math::common::FFTProvider<std::array<u32, 3>>>::release;

%typemap(out) std::unique_ptr<opa::math::common::FFTProvider<std::array<u32, 3>>> %{
{
  auto tmpx = new $1_ltype($1.release());
  $result = SWIG_NewPointerObj(tmpx, $&1_descriptor, SWIG_POINTER_OWN);
  }
  %}

wrap_unique_ptr(uptr_fftprovider_u32, opa::math::common::FFTProvider<u32>);
wrap_unique_ptr(uptr_fftprovider_d, opa::math::common::FFTProvider<double>);
wrap_unique_ptr(uptr_fftprovider_float, opa::math::common::FFTProvider<opa::math::common::Float>);
wrap_unique_ptr(uptr_fftdispatcher_d, opa::math::common::FFT2Dispatcher<double>);
wrap_unique_ptr(uptr_fftdispatcher_float, opa::math::common::FFT2Dispatcher<opa::math::common::Float>);

%template(fftreal_d) opa::math::common::FFT_RealF<double>;
%template(fftreal_float) opa::math::common::FFT_RealF<opa::math::common::Float>;


%include "opa/math/common/crt_fft.h"
%template(factprovider_u32) opa::math::common::FactProvider<u32>;
%template(factprovider_d) opa::math::common::FactProvider<double>;
%template(ring_u32_3) opa::math::common::Ring<std::array<u32, 3>>;
%template(field_u32_3) opa::math::common::Field<std::array<u32, 3>>;
%template(crt_u32_3) opa::math::common::CRTField<u32, 3>;
%template(fastpoly_u32) opa::math::common::FastPoly<u32>;
%template(fastpoly_d) opa::math::common::FastPoly<double>;
%template(fastpoly_float) opa::math::common::FastPoly<opa::math::common::Float>;


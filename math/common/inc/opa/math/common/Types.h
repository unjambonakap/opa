#pragma once

#include <opa/math/common/FractionField.h>
#include <opa/math/common/GF_p.h>
#include <opa/math/common/GF_q.h>
#include <opa/math/common/Poly.h>
#include <opa/math/common/PolyRing.h>
#include <opa/math/common/Utils.h>
#include <opa/math/common/Z.h>
#include <opa/math/common/base.h>
#include <opa/math/common/RealField.h>

OPA_NM_MATH_COMMON

typedef opa::math::common::Fraction<bignum> FFloat;
typedef opa::math::common::FractionField<bignum> QField;
typedef FFloat Q;
typedef Float Qf;
typedef Poly<bignum> P_Z;
typedef Poly<P_Z> P_ZZ;
typedef Poly<P_ZZ> P_ZZZ;
typedef OPA_BG Z;

typedef Poly<Q> P_Q;
typedef Poly<P_Q> P_K;
typedef Poly<P_Q> P_QQ;
typedef Poly<P_QQ> P_QQQ;

typedef Poly<Qf> P_Qf;
typedef Poly<P_Qf> P_QQf;
typedef Poly<P_QQf> P_QQQf;

extern GF_p GF2;
extern ZRing Ring_Z;
extern PolyRing<u32> PR_GF2;
extern PolyRing<bignum> PR_Z;
extern QField QF;
extern RealF *QfF;
extern ComplexF *CF;
extern PolyRing<Q> Q_x;
extern PolyRing<P_Q> Q_xy;
extern PolyRing<P_QQ> Q_xyz;
extern PolyRing<Qf> Qf_x;
extern PolyRing<P_Qf> Qf_xy;
extern PolyRing<P_QQf> Qf_xyz;
extern PolyRing<bignum> Z_x;
extern PolyRing<P_Z> Z_xy;
extern PolyRing<P_ZZ> Z_xyz;

void init_math_types();
const GF_p &get_gfp(u32 p);

inline P_Q to_pq(const P_Z &pz) {
  return Q_x.import_vec(import_vec(QF, pz.to_vec()));
}


OPA_NM_MATH_COMMON_END

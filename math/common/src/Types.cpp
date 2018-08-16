#include <opa/math/common/Types.h>
#include <opa/math/common/float.h>

OPA_NAMESPACE(opa, math, common)
GF_p GF2;
ZRing Ring_Z;
PolyRing<u32> PR_GF2;
PolyRing<bignum> PR_Z;
QField QF;
RealF *QfF;
ComplexF *CF;

std::map<u32, GF_p> g_gfps;

PolyRing<FFloat> Q_x;
PolyRing<P_Q> Q_xy;
PolyRing<P_QQ> Q_xyz;
PolyRing<bignum> Z_x;
PolyRing<P_Z> Z_xy;
PolyRing<P_ZZ> Z_xyz;
PolyRing<Qf> Qf_x;
PolyRing<P_Qf> Qf_xy;
PolyRing<P_QQf> Qf_xyz;

void init_math_types() {
  static bool is_init = false;
  if (!is_init) {
    Float_init();
    GF2.init(2);
    PR_GF2.init(&GF2);
    QF.init(&Ring_Z, true /* do_reduce */);
    is_init = true;
    PR_Z.init(&Ring_Z);

    Q_x.init(&QF);
    Q_xy.init(&Q_x);
    Q_xyz.init(&Q_xy);
    Z_x.init(&Ring_Z);
    Z_xy.init(&Z_x);
    Z_xyz.init(&Z_xy);

    QfF = new RealF;
    QfF->init(Float::Float_10pw(FLAGS_opa_float_precision));
    Qf_x.init(QfF);
    Qf_xy.init(&Qf_x);
    Qf_xyz.init(&Qf_xy);

    CF = new ComplexF;
    CF->init(Float::Float_10pw(FLAGS_opa_float_precision));

  }
}
const GF_p &get_gfp(u32 p) {
  if (!g_gfps.count(p)) {
    g_gfps[p].init(p);
  }
  return g_gfps[p];
}

OPA_REGISTER_INIT(init_math_types, init_math_types);

OPA_NAMESPACE_END(opa, math, common)

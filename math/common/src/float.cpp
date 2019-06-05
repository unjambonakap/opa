#include <opa/math/common/float.h>

#include <gmp.h>
#include <mpfr.h>
#include <opa/math/common/bignum.h>
#include <opa/utils/misc.h>
#include <opa_common.h>
DEFINE_int32(opa_float_precision, 100, "");
DEFINE_int32(opa_float_print_precision, 10, "");

using namespace opa::utils;
using namespace std;

OPA_NAMESPACE(opa, math, common)

#define DO_MPFR_OP(op, res)                                                    \
  {                                                                            \
    Float res;                                                                 \
    op;                                                                        \
    return res;                                                                \
  }

struct MpfrWrap : public __mpfr_struct {};
struct MpfrRandState : public __gmp_randstate_struct {};

Float *internal_pi;
Float *internal_e;

MpfrGState Float::g_params;
static bool is_init = false;

void Float_init(int seed) {
  if (!is_init) {
    is_init = true;
    Float::g_params.rand_state = new MpfrRandState();
    Float::g_params.precision = FLAGS_opa_float_precision;
    OPA_TRACES("Initing float", Float::g_params.precision);
    gmp_randinit_default(Float::g_params.rand_state);
    gmp_randseed_ui(Float::g_params.rand_state, seed);

    internal_pi = new Float(Float(-1).acos());
    internal_e = new Float(Float(1).exp());
  }
}
const Float &Float::g_pi() { return *internal_pi; }
const Float &Float::g_e() { return *internal_e; }

std::ostream &operator<<(std::ostream &os, const Float &f) {
  os << f.to_str(FLAGS_opa_float_print_precision);
  return os;
}

void Float::init() {
  OPA_CHECK0(is_init);
  m_state = new MpfrWrap;
  mpfr_init2(m_state, g_params.precision);
}

void Float::release() {
  if (m_state) {
    mpfr_clear(m_state);
    delete m_state;
    m_state = 0;
  }
}

Float::~Float() { release(); }
Float::Float() {
  init();
  reset();
}

Float::Float(s32 v) {
  init();
  from_s32(v);
}

Float::Float(u32 v) {
  init();
  from_u32(v);
}

Float::Float(unsigned long int v) {
  init();
  from_uli(v);
}

Float::Float(double v) {
  init();
  from_double(v);
}

Float::Float(const opa::stolen::StringRef &str) {
  init();
  from_str(str);
}

Float::Float(const Float &peer) {
  init();
  from_Float(peer);
}

Float::Float(const bignum &x) {
  init();
  from_bg(x);
}

Float::Float(const Fraction<bignum> &frac) {
  init();
  from_bg(frac.p);
  Float other(frac.q);
  this->sdiv(other);
}

Float &Float::operator=(const Float &f) {
  from_Float(f);
  return *this;
}

void Float::reset() { from_s32(0); }
void Float::from_str(const opa::stolen::StringRef &str) {
  mpfr_set_str(m_state, str.data(), 10, MPFR_RNDN);
}
void Float::from_double(double v) { mpfr_set_d(m_state, v, MPFR_RNDN); }
void Float::from_s32(s32 v) { mpfr_set_si(m_state, v, MPFR_RNDN); }
void Float::from_u32(u32 v) { mpfr_set_ui(m_state, v, MPFR_RNDN); }
void Float::from_uli(unsigned long int v) {
  mpfr_set_ui(m_state, v, MPFR_RNDN);
}
void Float::from_bg(const bignum &a) {
  mpfr_set_z(m_state, (__mpz_struct *)a.get_internal(), MPFR_RNDN);
}

Float Float::From_bg(const bignum &a) {
  Float res;
  res.from_bg(a);
  return res;
}

void Float::from_Float(const Float &peer) {
  mpfr_set(m_state, peer.m_state, MPFR_RNDN);
}

double Float::to_double() const { return mpfr_get_d(m_state, MPFR_RNDN); }
bignum Float::to_bignum() const { return bignum::fromFloat(*this); }
Fraction<bignum> Float::to_frac(const Float &precision) const {
  Float tmp = *this / precision;
  Fraction<bignum> res;
  res.p = tmp.to_bignum();
  res.q = precision.to_bignum();
  return reduce_fraction(res);
}

std::string Float::to_str(int ndigits) const {
  mpfr_exp_t exp;
  char *res = mpfr_get_str(nullptr, &exp, 10, ndigits, m_state, MPFR_RNDN);
  string x;

  {
    char *c = res;
    if (c[0] == '-') x += "-", ++c;
    x += "0.";
    x += stdsprintf("%se%d", c, exp);
  }

  mpfr_free_str(res);
  return x;
}

Float Float::mul(const Float &peer) const {
  DO_MPFR_OP(mpfr_mul(res.m_state, m_state, peer.m_state, MPFR_RNDN), res);
}
Float Float::div(const Float &peer) const {
  DO_MPFR_OP(mpfr_div(res.m_state, m_state, peer.m_state, MPFR_RNDN), res);
}
Float Float::add(const Float &peer) const {
  DO_MPFR_OP(mpfr_add(res.m_state, m_state, peer.m_state, MPFR_RNDN), res);
}
Float Float::sub(const Float &peer) const {
  DO_MPFR_OP(mpfr_sub(res.m_state, m_state, peer.m_state, MPFR_RNDN), res);
}

bool Float::eq(const Float &peer) const {
  return mpfr_equal_p(m_state, peer.m_state);
}
bool Float::diff(const Float &peer) const {
  return mpfr_lessgreater_p(m_state, peer.m_state);
}
bool Float::lt(const Float &peer) const {
  return mpfr_less_p(m_state, peer.m_state);
}
bool Float::gt(const Float &peer) const {
  return mpfr_greater_p(m_state, peer.m_state);
}
bool Float::le(const Float &peer) const {
  return mpfr_lessequal_p(m_state, peer.m_state);
}
bool Float::ge(const Float &peer) const {
  return mpfr_greaterequal_p(m_state, peer.m_state);
}

Float Float::sqrt() const {
  DO_MPFR_OP(mpfr_sqrt(res.m_state, m_state, MPFR_RNDN), res);
}
Float Float::cos() const {
  DO_MPFR_OP(mpfr_cos(res.m_state, m_state, MPFR_RNDN), res);
}
Float Float::sin() const {
  DO_MPFR_OP(mpfr_sin(res.m_state, m_state, MPFR_RNDN), res);
}
Float Float::tan() const {
  DO_MPFR_OP(mpfr_tan(res.m_state, m_state, MPFR_RNDN), res);
}
Float Float::acos() const {
  DO_MPFR_OP(mpfr_acos(res.m_state, m_state, MPFR_RNDN), res);
}
Float Float::asin() const {
  DO_MPFR_OP(mpfr_acos(res.m_state, m_state, MPFR_RNDN), res);
}
Float Float::atan() const {
  DO_MPFR_OP(mpfr_acos(res.m_state, m_state, MPFR_RNDN), res);
}
Float Float::atan2(const Float &peer) const {
  DO_MPFR_OP(mpfr_atan2(res.m_state, m_state, peer.m_state, MPFR_RNDN), res);
}

Float &Float::smul(const Float &peer) {
  mpfr_mul(m_state, m_state, peer.m_state, MPFR_RNDN);
  return *this;
}
Float &Float::sadd(const Float &peer) {
  mpfr_add(m_state, m_state, peer.m_state, MPFR_RNDN);
  return *this;
}
Float &Float::ssub(const Float &peer) {
  mpfr_sub(m_state, m_state, peer.m_state, MPFR_RNDN);
  return *this;
}
Float &Float::sdiv(const Float &peer) {
  mpfr_div(m_state, m_state, peer.m_state, MPFR_RNDN);
  return *this;
}
Float Float::abs() const {
  DO_MPFR_OP(mpfr_abs(res.m_state, m_state, MPFR_RNDN), res);
}

bool Float::eq(const Float &peer, const Float &eps) const {
  return sub(peer).abs() < eps;
}

Float Float::pow(const Float &pw) const {
  DO_MPFR_OP(mpfr_pow(res.m_state, m_state, pw.m_state, MPFR_RNDN), res);
}

Float Float::neg() const {
  DO_MPFR_OP(mpfr_neg(res.m_state, m_state, MPFR_RNDN), res);
}

Float FloatRng::get_gaussian(const Float &a, const Float &b) const {
  Float res;
  mpfr_nrandom(res.m_state, Float::g_params.rand_state, MPFR_RNDN);
  return res * b + a;
}

Float FloatRng::get_uni(const Float &a, const Float &b) const {
  Float res;
  mpfr_urandom(res.m_state, Float::g_params.rand_state, MPFR_RNDN);
  return res * (b - a) + a;
}

Float Float::exp() const {
  DO_MPFR_OP(mpfr_exp(res.m_state, m_state, MPFR_RNDN), res);
}
Float Float::log_e() const {
  DO_MPFR_OP(mpfr_log(res.m_state, m_state, MPFR_RNDN), res);
}
Float Float::log_2() const {
  DO_MPFR_OP(mpfr_log2(res.m_state, m_state, MPFR_RNDN), res);
}
Float Float::log_10() const {
  DO_MPFR_OP(mpfr_log10(res.m_state, m_state, MPFR_RNDN), res);
}
Float Float::log_base(const Float &base) const {
  return log_e() / base.log_e();
}
Float Float::atanh() const {
  DO_MPFR_OP(mpfr_atanh(res.m_state, m_state, MPFR_RNDN), res);
}

Float Float::tanh() const {
  DO_MPFR_OP(mpfr_tanh(res.m_state, m_state, MPFR_RNDN), res);
}

std::complex<Float> Float::iexp() const {
  return std::complex<Float>(cos(), sin());
}

Float Float::ceil() const { DO_MPFR_OP(mpfr_ceil(res.m_state, m_state), res); }
Float Float::floor() const {
  DO_MPFR_OP(mpfr_floor(res.m_state, m_state), res);
}
Float Float::round() const {
  DO_MPFR_OP(mpfr_round(res.m_state, m_state), res);
}

template <> Float FloatUtil::erfc<Float>(const Float &a) {
  DO_MPFR_OP(mpfr_erfc(res.m_state, a.m_state, MPFR_RNDN), res);
}
template <> Float FloatUtil::erf<Float>(const Float &a) {
  DO_MPFR_OP(mpfr_erf(res.m_state, a.m_state, MPFR_RNDN), res);
}

template <> u32 FloatUtil::cast<Float, u32>(const Float &from) {
  return from.to_bignum().getu32();
}
template <> bignum FloatUtil::cast<Float, bignum>(const Float &from) {
  return from.to_bignum();
}
template <> u32 FloatUtil::cast<double, u32>(const double &from) {
  return Float(from).to_bignum().getu32();
}
template <> bignum FloatUtil::cast<double, bignum>(const double &from) {
  return u32(from);
}

OPA_NAMESPACE_END(opa, math, common)

namespace std {
// clang-format off
opa::math::common::Float sqrt(const opa::math::common::Float &a) { return a.sqrt(); }
opa::math::common::Float cos(const opa::math::common::Float &a) { return a.cos(); }
opa::math::common::Float sin(const opa::math::common::Float &a) { return a.sin(); }
opa::math::common::Float tan(const opa::math::common::Float &a) { return a.tan(); }

opa::math::common::Float atan2(const opa::math::common::Float &a, const opa::math::common::Float &b) { return a.atan2(b); }
opa::math::common::Float acos(const opa::math::common::Float &a) { return a.acos(); }
opa::math::common::Float asin(const opa::math::common::Float &a) { return a.asin(); }
opa::math::common::Float atan(const opa::math::common::Float &a) { return a.atan(); }
opa::math::common::Float abs(const opa::math::common::Float &a) { return a.abs(); }
opa::math::common::Float pow(const opa::math::common::Float &a,
                             const opa::math::common::Float &b){return a.pow(b);}
opa::math::common::Float ceil(const opa::math::common::Float &a){return a.ceil();}
opa::math::common::Float floor(const opa::math::common::Float &a){return a.floor();}
opa::math::common::Float round(const opa::math::common::Float &a){return a.round();}
  }
// clang-format on

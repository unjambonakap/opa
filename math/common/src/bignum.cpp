#include "common/bignum.h"
#include "common/Utils.h"
#include <opa/math/common/FractionField.h>
#include <opa/math/common/Types.h>
#include <opa/math/common/float.h>
#include <opa/utils/string.h>

#include <gmp.h>
#include <mpfr.h>

using namespace std;

OPA_NAMESPACE_DECL3(opa, math, common)

#define DO_MPZ_OP(op, res)                                                                         \
  {                                                                                                \
    bignum res;                                                                                    \
    op;                                                                                            \
    return res;                                                                                    \
  }

struct MpzWrap : public __mpz_struct {};
struct MpzRandState : public __gmp_randstate_struct {};

MpzRandState *g_bignum_rand_state = 0;

void bignum_init(int seed) {
  g_bignum_rand_state = new MpzRandState();
  gmp_randinit_default(g_bignum_rand_state);
  gmp_randseed_ui(g_bignum_rand_state, seed);
}

bignum bignum::froms32(s32 x) {
  bignum a;
  mpz_set_si(a.a, x);
  return a;
}

bignum bignum::froms64(s64 x) { return bignum(x); }

bignum bignum::fromu32(u32 x) {
  bignum a;
  mpz_set_ui(a.a, x);
  return a;
}

bignum bignum::fromFloat(const Float &f) {
  bignum a;
  mpfr_get_z(a.a, (__mpfr_struct *)f.m_state, MPFR_RNDN);
  return a;
}

bignum::bignum(s64 u) {
  init();
  if (u == (s32)u) {
    mpz_set_si(a, u);
    return;
  }

  if (u < 0) {
    mpz_set(a, bignum::fromu64(-u).a);
    *this *= -1;
  } else {
    mpz_set(a, bignum::fromu64(u).a);
  }
}

bignum bignum::fromu64(u64 x) { return bignum::fromu32(x) | bignum::fromu32(x >> 32) << 32; }

void bignum::init() {
  a = new MpzWrap;
  mpz_init(a);
}

bignum::~bignum() {
  mpz_clear(a);
  delete a;
}

u32 bignum::get_size(int base) const { return mpz_sizeinbase(a, base); }

bignum::bignum(const bignum &b) {
  init();
  mpz_set(a, b.a);
}

bignum::bignum() {
  init();
  mpz_set_ui(a, 0);
}

bignum::bignum(const std::string &b, int base) {
  init();
  loadstr(b, base);
}

std::string bignum::str(int base) const {
  char *x = mpz_get_str(0, base, a);
  std::string res(x);
  free(x);
  if (base == 16) {
    if (sgn() >= 0) return "0x" + res;
    return "-0x" + std::string(res.c_str() + 1);
  }
  return res;
}

int bignum::sgn() const { return mpz_sgn(a); }

bignum &bignum::operator=(const bignum &b) {
  mpz_set(a, b.a);
  return *this;
}

bignum bignum::frombytes(opa::stolen::StringRef b) {
  return bignum::fromstr(opa::utils::b2h(b), 16);
}

bignum bignum::fromrbytes(opa::stolen::StringRef b) {
  std::string rb(b.data(), b.size());
  std::reverse(ALL(rb));
  return frombytes(rb);
}

void bignum::loadstr(const char *b, int base) { bad = mpz_set_str(a, b, base) == -1; }

bignum bignum::operator+(const bignum &b) const { return add(b); }

bignum bignum::operator*(const bignum &b) const { return mul(b); }

bignum bignum::operator/(const bignum &b) const { return div(b); }

bignum bignum::operator%(const bignum &b) const { return mod(b); }

bignum bignum::ceil(const bignum &b) const {
  bignum res;
  mpz_cdiv_q(res.a, a, b.a);
  return res;
}

bignum bignum::powm(const bignum &b, const bignum &mod) const {
  bignum res;
  mpz_powm(res.a, a, b.a, mod.a);
  return res;
}
bignum bignum::pow(u32 pw) const { DO_MPZ_OP(mpz_pow_ui(res.a, a, pw), res); }

bignum bignum::neg() const {
  bignum res;
  mpz_neg(res.a, a);
  return res;
}

bignum bignum::rshift(u32 v) const {
  bignum res;
  mpz_div_2exp(res.a, a, v);
  return res;
}

bignum bignum::lshift(u32 v) const {
  bignum res;
  mpz_mul_2exp(res.a, a, v);
  return res;
}

bignum bignum::modplus(const bignum &b) const { return (*this % b + b) % b; }

bignum bignum::mod(const bignum &b) const { DO_MPZ_OP(mpz_fdiv_r(res.a, a, b.a), res); }

void bignum::ediv(const bignum &b, bignum *q, bignum *r) const { mpz_fdiv_qr(q->a, r->a, a, b.a); }

bignum bignum::mul(const bignum &b) const { DO_MPZ_OP(mpz_mul(res.a, a, b.a), res); }

bignum bignum::add(const bignum &b) const { DO_MPZ_OP(mpz_add(res.a, a, b.a), res); }

bignum bignum::sub(const bignum &b) const { DO_MPZ_OP(mpz_sub(res.a, a, b.a), res); }

bignum bignum::div(const bignum &b) const { DO_MPZ_OP(mpz_fdiv_q(res.a, a, b.a), res); }

bignum bignum::sqrt() const { DO_MPZ_OP(mpz_sqrt(res.a, a), res); }

bool bignum::operator<(const bignum &b) const { return mpz_cmp(a, b.a) < 0; }

bool bignum::operator>(const bignum &b) const { return mpz_cmp(a, b.a) > 0; }

bool bignum::operator<=(const bignum &b) const { return mpz_cmp(a, b.a) <= 0; }

bool bignum::operator>=(const bignum &b) const { return mpz_cmp(a, b.a) >= 0; }

bool bignum::operator==(const bignum &b) const { return mpz_cmp(a, b.a) == 0; }

bool bignum::operator!=(const bignum &b) const { return mpz_cmp(a, b.a) != 0; }

void bignum::sneg() { mpz_neg(a, a); }
void bignum::sadd(const bignum &b) { mpz_add(a, a, b.a); }
void bignum::ssub(const bignum &b) { mpz_sub(a, a, b.a); }
void bignum::smul(const bignum &b) { mpz_mul(a, a, b.a); }
void bignum::sdiv(const bignum &b) { mpz_fdiv_q(a, a, b.a); }
void bignum::smod(const bignum &b) { mpz_fdiv_r(a, a, b.a); }
void bignum::srshift(u32 v) { mpz_div_2exp(a, a, v); }
void bignum::slshift(u32 v) { mpz_mul_2exp(a, a, v); }

u32 bignum::get_bit(u32 pos) const { return mpz_tstbit(a, pos); }

void bignum::disp(int base) const { std::cout << str(base) << std::endl; }

bignum bignum::gcd(const bignum &b) const { DO_MPZ_OP(mpz_gcd(res.a, a, b.a), res); }

bignum bignum::egcd(const bignum &b, bignum &u, bignum &v) const {
  DO_MPZ_OP(mpz_gcdext(res.a, u.a, v.a, a, b.a), res);
}

void bignum::ssqrt() { mpz_sqrt(a, a); }

bignum bignum::rand() const {
  OPA_CHECK0(g_bignum_rand_state != 0);
  DO_MPZ_OP(mpz_urandomm(res.a, g_bignum_rand_state, a), res);
}

bignum bignum::rand_signed() const { return this->rand() - (*this) / 2; }

u32 bignum::gets32() const {
  OPA_CHECK0(fitss32());
  return gets32_unsafe();
}

s64 bignum::gets64() const {
  OPA_CHECK0(fitss64());
  return gets64_unsafe();
}

u32 bignum::getu32() const {
  OPA_CHECK0(fitsu32());
  return getu32_unsafe();
}

u64 bignum::getu64() const {
  OPA_CHECK0(fitsu64());
  return getu64_unsafe();
}

bool bignum::fitsu64() const { return (*this >> 64 == 0); }
bool bignum::fitsu32() const { return (*this >> 32 == 0); }
bool bignum::fitss64() const { return (*this >> 32).fitss32(); }
bool bignum::fitss32() const { return mpz_fits_slong_p(a) != 0; }

Float bignum::getfloat() const { return Float::From_bg(*this); }

s32 bignum::gets32_unsafe() const { return mpz_get_si(a); }
s64 bignum::gets64_unsafe() const {
  s64 v = getu64_unsafe();
  return v;
}

u32 bignum::getu32_unsafe() const { return mpz_get_ui(a); }
u64 bignum::getu64_unsafe() const {

  u64 res = 0;
  bignum tmp = bignum::fromu64((1ll << 32) - 1);

  res |= ((u64)((*this >> 32) & tmp).getu32()) << 32;
  res |= (*this & tmp).getu32();
  return res;
}

bignum bignum::inv(const bignum &n) const {
  bignum u, v, d;

  d = egcd(n, u, v) % n;
  if (d != 1) return froms32(-1);
  if (u < 0) return n + u;
  return u;
}

bignum bignum::faste(const bignum &p, const bignum &mod) const { return t_faste(*this, p, mod); }

bignum bignum::operator|(const bignum &b) const { return or_(b); }
bignum bignum::operator^(const bignum &b) const { return xor_(b); }
bignum bignum::operator&(const bignum &b) const { return and_(b); }

bignum bignum::or_(const bignum &b) const { DO_MPZ_OP(mpz_ior(res.a, a, b.a), res); }

bignum bignum::xor_(const bignum &b) const { DO_MPZ_OP(mpz_xor(res.a, a, b.a), res); }

bignum bignum::and_(const bignum &b) const { DO_MPZ_OP(mpz_and(res.a, a, b.a), res); }

bignum bignum::abs() const { DO_MPZ_OP(mpz_abs(res.a, a), res); }

int bignum::legendre(const bignum &p) const { return mpz_legendre(a, p.a); }

std::string bignum::getbytes(int size) const {
  std::string hexstr = str(16);
  if (hexstr.size() & 1) hexstr = "0" + hexstr;

  if (size == -1)
    size = hexstr.size();
  else
    size *= 2;
  size -= hexstr.size();
  OPA_CHECK0(size >= 0);
  hexstr = std::string(size, '0') + hexstr;

  return opa::utils::h2b(hexstr);
}

std::string bignum::getrbytes(int size) const {
  std::string hexstr = str(16);
  if (hexstr.size() & 1) hexstr = "0" + hexstr;

  if (size == -1)
    size = hexstr.size();
  else
    size *= 2;
  size -= hexstr.size();
  OPA_CHECK0(size >= 0);
  hexstr = std::string(size, '0') + hexstr;

  auto tmp = opa::utils::h2b(hexstr);
  std::reverse(ALL(tmp));
  return tmp;
}

bignum bignum::root(u32 b) const { DO_MPZ_OP(mpz_root(res.a, a, b), res); }
bignum bignum::lcm(const bignum &a) const { return *this / gcd(a) * a; }

void bignum::base_repr(const bignum &base, std::vector<bignum> *res) {
  bignum tmp = *this;
  res->clear();
  while (tmp != 0) {
    res->emplace_back(tmp % base);
    tmp /= base;
  }
}

Fraction<bignum> reduce_fraction(const Fraction<bignum> &a) {
  bignum d = a.p.gcd(a.q);
  return QF.import(a.p / d, a.q / d);
}

OPA_NAMESPACE_DECL3_END

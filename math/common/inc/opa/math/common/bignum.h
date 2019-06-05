#pragma once

#include <opa/math/common/base.h>
#include <opa/stolen/StringRef.h>
#include <opa_common.h>

#define OPA_BG opa::math::common::bignum
OPA_NM_MATH_COMMON
class bignum;
OPA_NM_MATH_COMMON_END


OPA_NM_MATH_COMMON
struct MpzWrap;
class Float;
template <class T> class Fraction;

void bignum_init(int seed);
class bignum {

public:
  static bignum identity() { return bignum(1); }
  static bignum froms32(s32 x);
  static bignum fromu64(u64 x);
  static bignum fromu32(u32 x);
  static bignum fromstr(const std::string &b, int base = 16) {
    return bignum(b, base);
  }
  static bignum fromFloat(const Float &f);
  static bignum frombytes(opa::stolen::StringRef b);
  static bignum fromrbytes(opa::stolen::StringRef b);

  ~bignum();
  bignum(const bignum &b);
  bignum();
  bignum(s64 u);
  bignum(const std::string &b, int base = 16);
  std::string str(int base = 16) const;

  void loadstr(const char *b, int base = 16);
  void loadstr(const std::string &b, int base = 16) {
    loadstr(b.c_str(), base);
  }

  u32 get_size(int base = 2) const;
  void base_repr(const bignum &base, std::vector<bignum> *res);

  bignum &operator=(const bignum &b);
  bignum operator+(const bignum &b) const;
  bignum operator*(const bignum &b) const;
  bignum operator/(const bignum &b) const;
  bignum ceil(const bignum &b) const;
  bignum operator%(const bignum &b) const;
  bignum operator-() const { return neg(); }

  bignum operator|(const bignum &b) const;
  bignum operator^(const bignum &b) const;
  bignum operator&(const bignum &b) const;

  void ediv(const bignum &b, bignum *q, bignum *r) const;
  bignum powm(const bignum &b, const bignum &mod) const;
  bignum mod(const bignum &b) const;
  bignum modplus(const bignum &b) const;
  bignum mul(const bignum &b) const;
  bignum add(const bignum &b) const;
  bignum sub(const bignum &b) const;
  bignum div(const bignum &b) const;
  bignum operator>>(u32 v) const { return rshift(v); }
  bignum operator<<(u32 v) const { return lshift(v); }
  bignum root(u32 b) const;

  bignum neg() const;
  bignum or_(const bignum &b) const;
  bignum xor_(const bignum &b) const;
  bignum and_(const bignum &b) const;
  bignum lcm(const bignum &b) const;

  bignum rshift(u32 v) const;
  bignum lshift(u32 v) const;
  bignum pow(u32 pw) const;
  bignum sqrt() const;

  u32 getu32() const;
  u32 gets32() const;
  u64 getu64() const;
  s64 gets64() const;
  bool fitsu64() const;
  bool fitsu32() const;
  bool fitss64() const;
  bool fitss32() const;
  Float getfloat() const;
  u32 getu32_unsafe() const;
  u64 getu64_unsafe() const;
  s32 gets32_unsafe() const;
  s64 gets64_unsafe() const;
  std::string getbytes(int size = -1) const;
  std::string getrbytes(int size = -1) const;

  int legendre(const bignum &p) const;
  template <class T> T get_or(const T &a) const {
    if (*this > bignum(a)) return a;
    return this->auto_conv<T>();
  }

  template <class T> T get() const;

  bool operator<(const bignum &b) const;
  bool operator>(const bignum &b) const;
  bool operator<=(const bignum &b) const;
  bool operator>=(const bignum &b) const;
  bool operator==(const bignum &b) const;
  bool operator!=(const bignum &b) const;

  u32 get_bit(u32 pos) const;

  bignum &operator+=(const bignum &b) {
    sadd(b);
    return *this;
  }
  bignum &operator-=(const bignum &b) {
    ssub(b);
    return *this;
  }
  bignum &operator*=(const bignum &b) {
    smul(b);
    return *this;
  }
  bignum &operator/=(const bignum &b) {
    sdiv(b);
    return *this;
  }
  bignum &operator%=(const bignum &b) {
    smod(b);
    return *this;
  }

  static bignum s_sub(const bignum &a, const bignum &b){
    return a.sub(b);
  }

  void sneg();
  void sadd(const bignum &b);
  void ssub(const bignum &b);
  void smul(const bignum &b);
  void sdiv(const bignum &b);
  void smod(const bignum &b);
  void srshift(u32 v);
  void slshift(u32 v);
  int sgn() const;

  template <class T> T auto_conv() const;

  void ssqrt();
  bignum abs() const;

  bignum gcd(const bignum &b) const;
  bignum egcd(const bignum &b, bignum &u, bignum &v) const;

  bignum inv(const bignum &n) const;
  bignum faste(const bignum &p, const bignum &mod) const;

  void disp(int base = 16) const;
  bignum rand() const;
  bignum rand_signed() const;

  const MpzWrap *get_internal() const {
    return a;
  } // Use as: *(mpz_t *)get_internal()
  OPA_DECL_COUT_OPERATOR(bignum)


    bool bad = false;
private:
  void init();
  MpzWrap *a;
};
template <> inline u32 bignum::auto_conv() const { return getu32(); }
template <> inline s32 bignum::auto_conv() const { return gets32(); }
template <> inline u64 bignum::auto_conv() const { return getu64(); }
template <> inline s64 bignum::auto_conv() const { return gets64(); }
template <> inline bignum bignum::auto_conv() const { return *this; }

template <class T>
void conv_vector_t_to_bignum(const std::vector<T> &input,
                             std::vector<bignum> *output) {
  for (auto &v : input) {
    output->emplace_back(v);
  }
}

template <class T>
void conv_vector_bignum_to_t(const std::vector<bignum> &input,
                             std::vector<T> *output) {
  for (auto &v : input) {
    output->emplace_back(v.auto_conv<T>(v));
  }
}

template <class T>
std::vector<bignum> conv_vector_t_to_bignum(const std::vector<T> &input) {
  std::vector<bignum> res;
  conv_vector_t_to_bignum<T>(input, &res);
  return res;
}

template <class T>
std::vector<T> conv_vector_bignum_to_t(const std::vector<bignum> &input) {
  std::vector<T> res;
  conv_vector_bignum_to_t<T>(input, &res);
  return res;
}

Fraction<bignum> reduce_fraction(const Fraction<bignum> &a);

OPA_NM_MATH_COMMON_END

namespace std {

template <> struct hash<opa::math::common::bignum> {
  typedef opa::math::common::bignum argument_type;
  typedef std::size_t result_type;

  result_type operator()(argument_type const &s) const {
    return s.getu32_unsafe();
  }
};
}

namespace opa {
namespace op {
static OPA_BG operator-(const OPA_BG &a, const OPA_BG &b) { return OPA_BG::s_sub(a, b); }
}
}

using namespace opa::op;

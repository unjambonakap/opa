#pragma once

#include <opa/utils/hash.h>
#include <opa/math/common/FractionField.h>
#include <opa/math/common/base.h>
#include <opa/stolen/StringRef.h>
#include <opa/utils/serialize.h>
#include <opa_common.h>

DECLARE_int32(opa_float_print_precision);
DECLARE_int32(opa_float_precision);

OPA_NM_MATH_COMMON
struct MpfrWrap;
struct MpfrRandState;

struct MpfrGState {
  MpfrRandState *rand_state;
  int precision;
};
class FloatRng;
class bignum;
void Float_init(int seed);

class Float {
  friend class FloatRng;

public:
  static void set_precision(int prec);
  explicit Float();
  Float(s32 v);
  explicit Float(u32 v);
  explicit Float(unsigned long int v);
  explicit Float(const bignum &x);
  explicit Float(const Fraction<bignum> &frac);

  Float(double v);
  explicit Float(const opa::stolen::StringRef &str);
  ~Float();
  Float(const Float &peer);
  Float &operator=(const Float &f);

  Float neg() const;
  Float mul(const Float &peer) const;
  Float div(const Float &peer) const;
  Float add(const Float &peer) const;
  Float sub(const Float &peer) const;
  Float &smul(const Float &peer);
  Float &sadd(const Float &peer);
  Float &ssub(const Float &peer);
  Float &sdiv(const Float &peer);

  bool eq(const Float &peer, const Float &eps) const;
  bool eq(const Float &peer) const;
  bool diff(const Float &peer) const;
  bool lt(const Float &peer) const;
  bool gt(const Float &peer) const;
  bool le(const Float &peer) const;
  bool ge(const Float &peer) const;

  Float operator-() const { return neg(); }
  Float operator*(const Float &peer) const { return mul(peer); }
  Float operator+(const Float &peer) const { return add(peer); }
  Float operator-(const Float &peer) const { return sub(peer); }
  Float operator/(const Float &peer) const { return div(peer); }
  Float &operator*=(const Float &peer) { return smul(peer); }
  Float &operator+=(const Float &peer) { return sadd(peer); }
  Float &operator-=(const Float &peer) { return ssub(peer); }
  Float &operator/=(const Float &peer) { return sdiv(peer); }
  bool operator==(const Float &peer) const { return eq(peer); }
  bool operator!=(const Float &peer) const { return diff(peer); }
  bool operator<=(const Float &peer) const { return le(peer); }
  bool operator<(const Float &peer) const { return lt(peer); }
  bool operator>=(const Float &peer) const { return ge(peer); }
  bool operator>(const Float &peer) const { return gt(peer); }

  std::string to_str(int ndigits = 10) const;
  std::string to_data() const;
  double to_double() const;
  bignum to_bignum() const;
  Fraction<bignum> to_frac(const Float &precision) const;
  void from_Float(const Float &peer);
  void from_str(const opa::stolen::StringRef &str);
  void from_double(double v);
  void from_s32(s32 v);
  void from_u32(u32 v);
  void from_uli(unsigned long int v);
  void from_bg(const bignum &a);

  static Float From_bg(const bignum &a);
  static Float Float_10pw(int pw) { return Float(10).pow(pw); }

  void reset();
  Float sqrt() const;
  Float cos() const;
  Float sin() const;
  Float tan() const;

  Float acos() const;
  Float asin() const;
  Float atan() const;
  Float atan2(const Float &peer) const;
  Float atanh() const;
  Float tanh() const;

  Float exp() const;
  Float log_e() const;
  Float log_2() const;
  Float log_10() const;
  Float log_base(const Float &base) const;
  std::complex<Float> iexp() const;

  Float pow(const Float &pw) const;
  Float abs() const;
  static const Float &g_pi();
  static const Float &g_e();

  Float ceil() const;
  Float round() const;
  Float floor() const;

public:
  static MpfrGState g_params;

  MpfrWrap *m_state;

private:
  void init();
  void release();
};

class FloatRng {
public:
  Float get_gaussian(const Float &a, const Float &b) const;
  Float get_uni(const Float &a, const Float &b) const;
  static FloatRng *get() { return opa::utils::Singleton2<FloatRng>::get(); }
};

class FloatUtil {
public:
  template <class T> static T get_gaussian(const T &a, const T &b);
  template <class T> static const T &get_pi();
  template <class T> static const T &get_e();
  template <class T> static T get_rand_uni(const T &a, const T &b);

  template <class T>
  static bool eq(const std::complex<T> &a, const std::complex<T> &b, const T &eps) {
    return std::abs(a - b) < eps;
  }
  template <class T> static bool eq(const T &a, const T &b, const T &eps) {
    return std::abs(a - b) < eps;
  }

  template <class T> static std::complex<T> iexp(const T &a) {
    return std::complex<T>(std::cos(a), std::sin(a));
  }

  template <class T> static T erfc(const T &a) {
    return opa::utils::Conv2::conv<Float, T>(FloatUtil::erfc<Float>(a));
  }
  template <class T> static T erf(const T &a) {
    return opa::utils::Conv2::conv<Float, T>(FloatUtil::erf<Float>(a));
  }

  template <class T> static std::complex<T> unity_root(int n) {
    T tmp = get_pi<T>() * 2 / n;
    return std::complex<T>(std::cos(tmp), std::sin(tmp));
  }

  template <class T> static bool isz(const T &v, const T &cmp) { return abs(v) <= cmp; }

  template <class From, class To> static To cast(const From &from);
};

template <> inline Float FloatUtil::get_gaussian<Float>(const Float &a, const Float &b) {
  return FloatRng::get()->get_gaussian(a, b);
}
template <> inline Float FloatUtil::get_rand_uni<Float>(const Float &a, const Float &b) {
  return FloatRng::get()->get_gaussian(a, b);
}

template <> inline const Float &FloatUtil::get_pi<Float>() { return Float::g_pi(); }
template <> inline const Float &FloatUtil::get_e<Float>() { return Float::g_e(); }

template <> Float FloatUtil::erfc<Float>(const Float &a);
template <> Float FloatUtil::erf<Float>(const Float &a);

template <>
inline bool FloatUtil::isz<std::complex<double> >(const std::complex<double> &a,
                                                  const std::complex<double> &b) {
  return std::abs(std::real(a)) < std::real(b) && std::abs(std::imag(a)) < std::imag(b);
}

template <>
inline bool FloatUtil::isz<std::complex<Float> >(const std::complex<Float> &a,
                                                 const std::complex<Float> &b) {
  return std::abs(std::real(a)) < std::real(b) && std::abs(std::imag(a)) < std::imag(b);
}

template <> u32 FloatUtil::cast<Float, u32>(const Float &from);
template <> bignum FloatUtil::cast<Float, bignum>(const Float &from);
template <> u32 FloatUtil::cast<double, u32>(const double &from);
template <> bignum FloatUtil::cast<double, bignum>(const double &from);

namespace {
double pi = acos(-1.);
double e = exp(1);
} // namespace

template <> inline double FloatUtil::get_gaussian<double>(const double &a, const double &b) {
  return FloatRng::get()->get_gaussian(a, b).to_double();
}
template <> inline double FloatUtil::get_rand_uni<double>(const double &a, const double &b) {
  return FloatRng::get()->get_uni(a, b).to_double();
}

template <>
inline std::complex<double>
FloatUtil::get_rand_uni<std::complex<double> >(const std::complex<double> &a,
                                               const std::complex<double> &b) {
  return std::complex<double>(get_rand_uni<double>(a.real(), b.real()),
                              get_rand_uni<double>(a.imag(), b.imag()));
}

template <>
inline std::complex<Float>
FloatUtil::get_rand_uni<std::complex<Float> >(const std::complex<Float> &a,
                                              const std::complex<Float> &b) {
  return std::complex<Float>(get_rand_uni<Float>(a.real(), b.real()),
                             get_rand_uni<Float>(a.imag(), b.imag()));
}

template <> inline const double &FloatUtil::get_pi<double>() { return pi; }
template <> inline const double &FloatUtil::get_e<double>() { return e; }

std::ostream &operator<<(std::ostream &os, const Float &f);

template <class T> struct FloatEQCompare {
  T eps;
  FloatEQCompare(const T &eps) : eps(eps) {}
  bool operator()(const std::complex<T> &a, const std::complex<T> &b) const {
    return FloatUtil::eq(a, b, eps);
  }
  bool operator()(const T &a, const T &b) const { return FloatUtil::eq(a, b, eps); }
};

OPA_NM_MATH_COMMON_END

namespace opa {
namespace utils {
template <> inline double Conv2::conv(const opa::math::common::Float &x) { return x.to_double(); }
} // namespace utils
} // namespace opa

namespace std {
template <class T> bool operator<(const std::complex<T> &a, const std::complex<T> &b) {
  if (a.real() != b.real()) return a.real() < b.real();
  return a.imag() < b.imag();
}
} // namespace std

namespace std {

template <> struct hash<opa::math::common::Float> {
  typedef opa::math::common::Float argument_type;
  typedef std::size_t result_type;

  result_type operator()(argument_type const &s) const {
    return std::hash<std::string>{}(s.to_data());
  }
};
template <> struct hash<std::complex<opa::math::common::Float> > {
  typedef std::complex<opa::math::common::Float> argument_type;
  typedef std::size_t result_type;

  result_type operator()(argument_type const &s) const {
    return std::hash<std::pair<opa::math::common::Float, opa::math::common::Float> >()(
      MP(s.real(), s.imag()));
  }
};
} // namespace std

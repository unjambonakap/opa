#pragma once

#include <opa/math/common/Field.h>
#include <opa/math/common/Utils.h>
#include <opa/math/common/base.h>

OPA_NM_MATH_COMMON
template <u32 p> struct MGContextT {
  static constexpr u64 k2_32 = 1ull << 32;
  static constexpr u32 sip_mod = u64inv_egcd(-p & (k2_32 - 1), k2_32);
  static constexpr u32 pmod_2_64 = k2_32 % p * k2_32 % p;
  static constexpr u64 m = (__uint128_t(1) << 64) / p;

  inline u32 mul(u32 x, u32 y) const {
    if (1) {
      s64 a = (s64)x * y;
      uint64_t q = ((__uint128_t(m) * a) >> 64);
      a -= q * p;
      if (a >= p) a -= p;
      return a;
    }
    return raise(reduce(u64(x) * y));
  }

  inline u32 raise(u32 v) const {
    // compute v *2^32 mod p
    return reduce((u64)pmod_2_64 * v);
  }

  inline u32 reduce(u64 v) const {
    uint64_t q = ((__uint128_t(m) * v) >> 64);
    v -= q * p;
    if (v >= p) v -= p;
    return v;

    auto res = (v + u64(m) * p) >> 32;
    if (res >= p) res -= p;
    return res;
  }
};
struct MGContext {
  u32 p;
  u32 sip_mod;
  u32 pmod_2_64;
  MGContext(u32 p) : p(p) {
    constexpr u64 k2_32 = 1ull << 32;
    sip_mod = u64inv_egcd(-p & (k2_32 - 1), k2_32);
    pmod_2_64 = k2_32 % p * k2_32 % p;
  }

  inline u32 mul(u32 a, u32 b) const { return raise(reduce(u64(a) * b)); }

  inline u32 raise(u32 v) const {
    // compute v *2^32 mod p
    return reduce((u64)pmod_2_64 * v);
  }

  inline u32 reduce(u64 v) const {
    // compute v *2^-32 mod p
    u32 m = (u32)v * sip_mod;
    auto res = (v + u64(m) * p) >> 32;
    if (res >= p) res -= p;
    return res;
  }
};

class GF_p : public GF_pT<u32> {
  u32 _p;

public:
  MGContext _mgc;
  GF_p(u32 n, S64Factors factors = {});
  void init(u32 n) {
    OPA_CHECK0(n < (1ull << 31));
    GF_pT<u32>::init(n, n);
    activateU32();
    _p = n;
  }

  virtual u32 import_bg(const bignum &a) const { return importu32((a % _p).getu32()); }
  virtual bignum export_base(const u32 &x) const { return bignum::fromu32(x); }

  inline u32 reduceu32(u32 v) const {
    if (v >= _p) v -= _p;
    return v;
  }

  inline u32 add(const u32 &a, const u32 &b) const override final { return reduceu32(a + b); }
  inline u32 neg(const u32 &a) const override final { return reduceu32(_p - a); }
  inline u32 sub(const u32 &a, const u32 &b) const override final { return reduceu32(a + _p - b); }
  inline u32 mul(const u32 &a, const u32 &b) const override final { return _mgc.mul(a, b); }
  u32 import(const u32 &a) const override final { return a < _p ? a : a % _p; }

  inline u32 inv(const u32 &a) const override final {
    return math::common::u32_faste(a, _p - 2, _p);
  }

  inline bool isZ(const u32 &a) const override final { return a == 0; }

  inline bool isE(const u32 &a) const override final { return a == 1; }

  u32 getZ() const override final { return 0; }

  u32 getE() const override final { return 1; }

  u32 importu32(u32 a) const override final { return a % _p; }

  u32 getRandRaw() const override final { return rng() % _p; }
};

template <u32 N> class GF_pMod : public GF_pT<u32> {

public:
  MGContextT<N> _mgc;
  GF_pMod(S64Factors factors = {}) : GF_pT(N, false) {
    m_factors = to_bgfactors(factors);
    Field<u32>::init(N, N);
    activateU32();
  }

  virtual u32 import_bg(const bignum &a) const { return importu32((a % N).getu32()); }
  virtual bignum export_base(const u32 &x) const { return bignum::fromu32(x); }

  inline u32 reduceu32(u32 v) const {
    if (v >= N) v -= N;
    return v;
  }

  inline u32 add(const u32 &a, const u32 &b) const override final { return reduceu32(a + b); }
  inline u32 neg(const u32 &a) const override final { return reduceu32(N - a); }
  inline u32 sub(const u32 &a, const u32 &b) const override final {
    auto v = s32(a) - b;
    if (v < 0) return v + n;
    return v;
  }
  inline u32 mul(const u32 &a, const u32 &b) const override final { return _mgc.mul(a, b); }
  inline u32 import(const u32 &a) const override final { return a < N ? a : a % N; }

  inline u32 inv(const u32 &a) const override final { return math::common::u32_faste(a, N - 2, N); }
  inline u32 div(const u32 &a, const u32 &b) const override final { return mul(a, inv(b)); }

  inline bool isZ(const u32 &a) const override final { return a == 0; }

  inline bool isE(const u32 &a) const override final { return a == 1; }

  inline u32 getZ() const override final { return 0; }

  inline u32 getE() const override final { return 1; }

  inline u32 importu32(u32 a) const override final { return a % N; }

  u32 getRandRaw() const override final { return rng() % N; }
};

OPA_NM_MATH_COMMON_END

#pragma once

#include <opa/math/common/Field.h>
#include <opa/math/common/Utils.h>
#include <opa/math/common/base.h>

OPA_NM_MATH_COMMON

class GF_p : public GF_pT<u32> {
  u32 _p;

public:
  GF_p(u32 n, S64Factors factors = {});
  void init(u32 n) {
    GF_pT<u32>::init(n, n);
    activateU32();
    _p = n;
  }

  virtual u32 import_bg(const bignum &a) const { return importu32((a % _p).getu32()); }
  virtual bignum export_base(const u32 &x) const { return bignum::fromu32(x); }

  inline u32 mul(const u32 &a, const u32 &b) const override final { return (u64)a * b % _p; }

  inline u32 add(const u32 &a, const u32 &b) const override final { return ((u64)a + b) % _p; }

  inline u32 inv(const u32 &a) const override final {
    return math::common::u32_faste(a, _p - 2, _p);
  }

  inline u32 neg(const u32 &a) const override final { return ((u64)_p - a) % _p; }

  inline bool isZ(const u32 &a) const override final { return a == 0; }

  inline bool isE(const u32 &a) const override final { return a == 1; }

  u32 import(const u32 &a) const override final { return ((u64)a % _p + _p) % _p; }
  inline u32 sub(const u32 &a, const u32 &b) const override final { return (_p + (u64)a - b) % _p; }

  u32 getZ() const override final { return 0; }

  u32 getE() const override final { return 1; }

  u32 importu32(u32 a) const override final { return a % _p; }

  u32 getRandRaw() const override final { return rng() % _p; }
};

template <u32 N> class GF_pMod : public GF_pT<u32> {

public:
  GF_pMod(S64Factors factors = {}) : GF_pT(N, false) {
    m_factors = to_bgfactors(factors);
    Field<u32>::init(N, N);
    activateU32();
  }


  virtual u32 import_bg(const bignum &a) const { return importu32((a % N).getu32()); }
  virtual bignum export_base(const u32 &x) const { return bignum::fromu32(x); }

  inline u32 mul(const u32 &a, const u32 &b) const override final { return (u64)a * b % N; }

  inline u32 add(const u32 &a, const u32 &b) const override final { return ((u64)a + b) % N; }
  inline u32 sub(const u32 &a, const u32 &b) const override final { return (N + (u64)a - b) % N; }

  inline u32 inv(const u32 &a) const override final { return math::common::u32_faste(a, N - 2, N); }

  inline u32 neg(const u32 &a) const override final { return ((u64)N - a) % N; }

  inline bool isZ(const u32 &a) const override final { return a == 0; }

  inline bool isE(const u32 &a) const override final { return a == 1; }

  u32 import(const u32 &a) const override final { return ((u64)a % N + N) % N; }

  u32 getZ() const override final { return 0; }

  u32 getE() const override final { return 1; }

  u32 importu32(u32 a) const override final { return a % N; }

  u32 getRandRaw() const override final { return rng() % N; }
};

OPA_NM_MATH_COMMON_END

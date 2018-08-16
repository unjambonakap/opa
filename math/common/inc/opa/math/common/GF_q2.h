#pragma once

#include <opa/math/common/Field.h>
#include <opa/math/common/Poly.h>
#include <opa/math/common/PolyRing.h>
#include <opa/math/common/UtilsPoly.h>

OPA_NAMESPACE(opa, math, common)
template <int N>
class PolyF2 {
 public:
  static constexpr int SZ = (N + 64) / 64;
  PolyF2(){
    memset(tb, 0, sizeof(tb));

  }

  void sadd(const PolyF2<N> &other) {
    REP(i, SZ)
    tb[i] ^= other.tb[i];
  }

  bool isZ() const {
    u64 has = 0;
    REP(i, SZ) has |= tb[i];
    return has == 0;
  }

  bool isE() const {
    u64 has = 0;
    REP(i, SZ - 1) has |= tb[i + 1];
    return has == 0 && tb[0] == 1;
  }

  u64 tb[SZ];
};

template <int N>
class GF_q2 : public Field<PolyF2<N>> {
 public:
  void init(const PolyF2<N> &mod_poly) {
    m_mod_poly = mod_poly;
    qm1 = bignum(1).lshift(N) - 1;
  }

  virtual PolyF2<N> inv(const PolyF2<N> &a) const override {
    auto res = this->faste(a, qm1 - 1);
    return res;
  }

  virtual PolyF2<N> mul(const PolyF2<N> &a, const PolyF2<N> &b) const override {
    auto res=get_poly();
  }

  virtual PolyF2<N> add(const PolyF2<N> &a, const PolyF2<N> &b) const override {
    auto res = a;
    res.sadd(b);
    return res;
  }

  virtual PolyF2<N> neg(const PolyF2<N> &a) const override { return a; }

  virtual bool isZ(const PolyF2<N> &a) const override { return a.isZ(); }
  virtual bool isE(const PolyF2<N> &a) const override { return a.isE(); }
  virtual PolyF2<N> import(const PolyF2<N> &a) const override { return a; }

  virtual PolyF2<N> getZ() const override { return get_poly(); }
  virtual PolyF2<N> getE() const override {
    auto res = get_poly();
    res.tb[0] = 1;
    return res;
  }

  PolyF2<N> get_poly() const {
    PolyF2<N> res;
    return res;
  }

 private:
  bignum qm1;
  PolyF2<N> m_mod_poly;
  PolyF2<N> m_cache[2*N];
};

OPA_NAMESPACE_END(opa, math, common)

#pragma once

#include <opa/math/common/Ring.h>
#include <opa/math/common/Utils.h>

OPA_NAMESPACE_DECL3(opa, math, common)
extern bignum zrand_bound;

class ZRing : public Ring<OPA_BG> {
public:
  ZRing() : Ring<OPA_BG>(-1, -1) {}
  virtual bool compareRank(const OPA_BG &a, const OPA_BG &b) const {
    return a.abs() < b.abs();
  }

  virtual bool lt(const OPA_BG &a, const OPA_BG &b) const {
    return a < b;
  }

  virtual OPA_BG mul(const OPA_BG &a, const OPA_BG &b) const;
  virtual OPA_BG add(const OPA_BG &a, const OPA_BG &b) const;
  virtual OPA_BG neg(const OPA_BG &a) const;

  virtual bool isInv(const OPA_BG &a) const;
  virtual OPA_BG inv(const OPA_BG &a) const;

  virtual bool isZ(const OPA_BG &a) const;
  virtual bool isE(const OPA_BG &a) const;
  virtual OPA_BG import(const OPA_BG &a) const;

  virtual OPA_BG getZ() const;
  virtual OPA_BG getE() const;
  virtual OPA_BG importu32(u32 a) const;
  virtual bool ediv(const OPA_BG &a, const OPA_BG &b, OPA_BG *q,
                    OPA_BG *r) const;

  virtual OPA_BG getRand() const;

  virtual OPA_BG abs(const OPA_BG &a) const {
    return a.abs();
  }

  virtual void norm_fraction(OPA_BG &a, OPA_BG &b) const {
    if (b < 0) a *= -1, b *= -1;
  }
  virtual ~ZRing() {}
};

OPA_NAMESPACE_DECL3_END

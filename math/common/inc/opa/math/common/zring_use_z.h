#pragma once

#include <opa/math/common/Utils.h>
#include <opa/math/common/Poly.h>
#include <opa/math/common/Ring.h>
#include <opa/math/common/Field.h>
#include <opa/math/common/bignum.h>

OPA_NAMESPACE_DECL3(opa, math, common)

class ZRing : public Ring<OPA_BG> {
private:
public:
  ZRing() : Ring<OPA_BG>(-1, -1) {}
  virtual ~ZRing() {}

  virtual bool isInv(const OPA_BG &a) const { return a == 1 || a == -1; }
  virtual OPA_BG inv(const OPA_BG &a) const {
    assert(isInv(a));
    return a;
  }
  virtual bool compareRank(const OPA_BG &a, const OPA_BG &b) const {
    return a.abs() < b.abs();
  }
  virtual bool ediv(const OPA_BG &a, const OPA_BG &b, OPA_BG *q,
                    OPA_BG *r) const {

    if (b == 0) return false;

    if (q) *q = a / b;
    if (r) *r = a % b;
    return true;
  }

  virtual OPA_BG mul(const OPA_BG &a, const OPA_BG &b) const { return a * b; }

  virtual OPA_BG add(const OPA_BG &a, const OPA_BG &b) const { return a + b; }

  virtual OPA_BG neg(const OPA_BG &a) const { return -a; }

  virtual bool isZ(const OPA_BG &a) const { return a == 0; }

  virtual bool isE(const OPA_BG &a) const { return a == 1; }

  virtual OPA_BG getZ() const { return OPA_BG(); }

  virtual OPA_BG getE() const { return OPA_BG(1); }
  virtual OPA_BG getRand() const { return OPA_BG(2).lshift(100).rand(); }

  OPA_BG import(const OPA_BG &p) const { return p; }
};

OPA_NAMESPACE_DECL3_END

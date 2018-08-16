#pragma once

#include <opa/math/common/Ring.h>
#include <opa/math/common/Utils.h>

OPA_NAMESPACE_DECL3(opa, math, common)

class Zn_BG : public Ring<OPA_BG> {
public:
  Zn_BG(const OPA_BG &n) { init(n); }
  Zn_BG() {}
  void init(const OPA_BG &n) { Ring<OPA_BG>::init(n, n); }
  virtual bool compareRank(const OPA_BG &a, const OPA_BG &b) const {
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

  ~Zn_BG() {}
};

OPA_NAMESPACE_DECL3_END

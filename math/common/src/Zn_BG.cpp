#include "common/Zn_BG.h"

OPA_NAMESPACE_DECL3(opa, math, common)

OPA_BG Zn_BG::mul(const OPA_BG &a, const OPA_BG &b) const {
  return a * b % getSize();
}
OPA_BG Zn_BG::add(const OPA_BG &a, const OPA_BG &b) const {
  return (a + b) % getSize();
}
OPA_BG Zn_BG::neg(const OPA_BG &a) const { return (getSize() - a) % getSize(); }

bool Zn_BG::isInv(const OPA_BG &a) const { return getSize().gcd(a) == 1; }
OPA_BG Zn_BG::inv(const OPA_BG &a) const {
  OPA_BG u, v;
  getSize().egcd(a, u, v);
  return v;
}

bool Zn_BG::isZ(const OPA_BG &a) const { return a == 0; }
bool Zn_BG::isE(const OPA_BG &a) const { return a == 1; }
OPA_BG Zn_BG::import(const OPA_BG &a) const {
  return (getSize() + a % getSize()) % getSize();
}

OPA_BG Zn_BG::getZ() const { return 0; }
OPA_BG Zn_BG::getE() const { return 1; }

bool Zn_BG::ediv(const OPA_BG &a, const OPA_BG &b, OPA_BG *q, OPA_BG *r) const {
  if (isZ(b))
    return false;
  OPA_BG d_a = a.gcd(getSize());
  OPA_BG d_b = b.gcd(getSize());
  if (d_a % d_b != 0)
    return false;
  if (q)
    *q = a / d_b * this->inv(b / d_b) % getSize();
  if (r)
    *r = getZ();
  return true;
}
OPA_BG Zn_BG::importu32(u32 a) const { return OPA_BG::fromu64(a) % getSize(); }

OPA_BG Zn_BG::getRand() const {
  OPA_BG x = getSize().rand();
  return x;
}

OPA_NAMESPACE_DECL3_END

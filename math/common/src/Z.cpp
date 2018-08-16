#include "common/Z.h"

OPA_NAMESPACE_DECL3(opa, math, common)

bignum zrand_bound = bignum(2).lshift(100);

OPA_BG ZRing::mul(const OPA_BG &a, const OPA_BG &b) const { return a * b; }
OPA_BG ZRing::add(const OPA_BG &a, const OPA_BG &b) const { return a + b; }
OPA_BG ZRing::neg(const OPA_BG &a) const { return -a; }

bool ZRing::isInv(const OPA_BG &a) const { return a == 1 || a == -1; }
OPA_BG ZRing::inv(const OPA_BG &a) const {
  OPA_CHECK(isInv(a), a);
  return a;
}

bool ZRing::isZ(const OPA_BG &a) const { return a == 0; }
bool ZRing::isE(const OPA_BG &a) const { return a == 1; }
OPA_BG ZRing::import(const OPA_BG &a) const { return a; }

OPA_BG ZRing::getZ() const { return 0; }
OPA_BG ZRing::getE() const { return 1; }

bool ZRing::ediv(const OPA_BG &a, const OPA_BG &b, OPA_BG *q, OPA_BG *r) const {
  if (b == 0) return false;
  if (q)
    *q = a / b;
  if (r)
    *r = a % b;
  return true;
}
OPA_BG ZRing::importu32(u32 a) const { return OPA_BG::fromu64(a); }
OPA_BG ZRing::getRand() const { return zrand_bound.rand_signed(); }

OPA_NAMESPACE_DECL3_END

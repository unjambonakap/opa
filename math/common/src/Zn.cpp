#include "common/Zn.h"

OPA_NAMESPACE_DECL3(opa, math, common)

u32 Zn::mul(const u32 &a, const u32 &b) const {
  return (u64)a * b % getSizeU32();
}
u32 Zn::add(const u32 &a, const u32 &b) const {
  return ((u64)a + b) % getSizeU32();
}
u32 Zn::neg(const u32 &a) const {
  return ((u64)getSizeU32() - a) % getSizeU32();
}

bool Zn::isInv(const u32 &a) const { return gcd((u64)getSizeU32(), a) == 1; }
u32 Zn::inv(const u32 &a) const {
  u32 u, v;
  u32 res = this->_egcd(getSizeU32(), a, getE(), getZ(), getZ(), getE(), u, v);
  return v;
}

bool Zn::isZ(const u32 &a) const { return a == 0; }
bool Zn::isE(const u32 &a) const { return a == 1; }
u32 Zn::import(const u32 &a) const {
  return ((u64)getSizeU32() + a % getSizeU32()) % getSizeU32();
}

u32 Zn::getZ() const { return 0; }
u32 Zn::getE() const { return 1; }

bool Zn::ediv(const u32 &a, const u32 &b, u32 *q, u32 *r) const {
  if (isZ(b))
    return false;
  if (q)
    *q = a / b;
  if (r)
    *r = a % b;
  return true;
}

u32 Zn::div(const u32 &a, const u32 &b) const {
  if (isZ(a)) return getZ();

  u32 ga = _gcd(getSizeU32(), a);
  u32 gb = _gcd(getSizeU32(), b);
  OPA_CHECK(ga%gb==0, a, b, ga, gb);

  Zn tmp(getSizeU32() / gb);
  u32 res = mul(a / ga, mul(tmp.inv(b / gb), (ga / gb)));
  OPA_CHECK(mul(b,res)==a, b, res, a, mul(b,res));
  return res;
}

u32 Zn::importu32(u32 a) const { return a % getSizeU32(); }

u32 Zn::getRand() const {
  u32 x = rng() % getSizeU32();
  return x;
}

bool Zn::compareRank(const u32 &a, const u32 &b) const {
  return _gcd(getSizeU32(), a) < _gcd(getSizeU32(), b);
}

OPA_NAMESPACE_DECL3_END

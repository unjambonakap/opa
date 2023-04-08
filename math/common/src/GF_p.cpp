#include "common/Utils.h"
#include <common/GF_p.h>

OPA_NAMESPACE_DECL3(opa, math, common)

void GF_p::init(u32 n) {
  Field<u32>::init(n, n);
  activateU32();
}

GF_p::GF_p(u32 n, S64Factors factors) : GF_pT(n, false) { 
  m_factors = to_bgfactors(factors);
  init(n); }

u32 GF_p::mul(const u32 &a, const u32 &b) const {
  return (u64)a * b % getSizeU32();
}

u32 GF_p::add(const u32 &a, const u32 &b) const {
  return ((u64)a + b) % getSizeU32();
}

u32 GF_p::inv(const u32 &a) const {
  return math::common::u32_faste(a, getSizeU32() - 2, getSizeU32());
}

u32 GF_p::neg(const u32 &a) const {
  return ((u64)getSizeU32() - a) % getSizeU32();
}

bool GF_p::isZ(const u32 &a) const { return a == 0; }

bool GF_p::isE(const u32 &a) const { return a == 1; }

u32 GF_p::import(const u32 &a) const {
  return ((u64)a % getSizeU32() + getSizeU32()) % getSizeU32();
}

u32 GF_p::getZ() const { return 0; }

u32 GF_p::getE() const { return 1; }

u32 GF_p::importu32(u32 a) const { return a % getSizeU32(); }

u32 GF_p::getRandRaw() const { return rng() % getSizeU32(); }

OPA_NAMESPACE_DECL3_END

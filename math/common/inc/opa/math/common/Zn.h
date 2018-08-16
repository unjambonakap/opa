#pragma once

#include <opa/math/common/Ring.h>
#include <opa/math/common/Utils.h>

OPA_NAMESPACE_DECL3(opa, math, common)

class Zn : public Ring<u32> {
public:
  Zn(u32 n) : Ring<u32>(n, n) { activateU32(); }
  virtual bool compareRank(const u32 &a, const u32 &b) const;

  virtual u32 mul(const u32 &a, const u32 &b) const;
  virtual u32 add(const u32 &a, const u32 &b) const;
  virtual u32 neg(const u32 &a) const;

  virtual bool isInv(const u32 &a) const;
  virtual u32 inv(const u32 &a) const;

  virtual bool isZ(const u32 &a) const;
  virtual bool isE(const u32 &a) const;
  virtual u32 import(const u32 &a) const;

  virtual u32 getZ() const;
  virtual u32 getE() const;
  virtual u32 importu32(u32 a) const;
  virtual bool ediv(const u32 &a, const u32 &b, u32 *q, u32 *r) const;
  virtual u32 div(const u32 &a, const u32 &b) const;

  virtual u32 getRand() const;

  ~Zn() {}
};

OPA_NAMESPACE_DECL3_END

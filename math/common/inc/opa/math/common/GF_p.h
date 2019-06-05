#pragma once

#include <opa/math/common/Field.h>
#include <opa/math/common/Utils.h>
#include <opa/math/common/base.h>

OPA_NM_MATH_COMMON

class GF_p : public Field<u32> {
public:
  GF_p(u32 n);
  GF_p();
  void init(u32 n);

  virtual u32 mul(const u32 &a, const u32 &b) const;
  virtual u32 add(const u32 &a, const u32 &b) const;
  virtual u32 inv(const u32 &a) const;
  virtual u32 neg(const u32 &a) const;
  virtual bool isZ(const u32 &a) const;
  virtual bool isE(const u32 &a) const;
  virtual u32 import(const u32 &a) const;
  virtual u32 import_bg(const bignum &a) const {
    return importu32((a % getSizeU32()).getu32());
  }
  virtual u32 getRandRaw() const;
  virtual bignum export_base(const u32 &x) const {
    return bignum::fromu32(x);
  }

  virtual u32 getZ() const;
  virtual u32 getE() const;
  virtual u32 importu32(u32 a) const;
};

OPA_NM_MATH_COMMON_END

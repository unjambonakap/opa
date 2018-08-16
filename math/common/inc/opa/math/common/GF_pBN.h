#pragma once

#include <opa/math/common/Utils.h>
#include <opa/math/common/Field.h>
#include <opa/math/common/bignum.h>

OPA_NAMESPACE_DECL3(opa, math, common)

class GF_pBN : public Field<bignum> {
public:
  GF_pBN();
  GF_pBN(const bignum &n) { this->init(n); }
  void init(const bignum &n);

  virtual bignum mul(const bignum &a, const bignum &b) const;
  virtual bignum add(const bignum &a, const bignum &b) const;
  virtual bignum inv(const bignum &a) const;
  virtual bignum neg(const bignum &a) const;
  virtual bool isZ(const bignum &a) const;
  virtual bool isE(const bignum &a) const;
  virtual bignum import(const bignum &a) const;
  virtual bignum getRandRaw() const;
  virtual bignum import_bg(const bignum &a) const { return import(a); }

  virtual bignum getZ() const;
  virtual bignum getE() const;
  virtual bignum importu32(u32 a) const;
};

OPA_NAMESPACE_DECL3_END

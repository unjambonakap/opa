#include <common/GF_pBN.h>

OPA_NAMESPACE_DECL3(opa, math, common)

void GF_pBN::init(const bignum &n) {
  GF_pT<bignum>::init(n, n);
}

bignum GF_pBN::mul(const bignum &a, const bignum &b) const {
  return a * b % getSize();
}

bignum GF_pBN::add(const bignum &a, const bignum &b) const {
  return (a + b) % getSize();
}

bignum GF_pBN::inv(const bignum &a) const { return a.inv(getSize()); }

bignum GF_pBN::neg(const bignum &a) const {
  return (getSize() - a) % getSize();
}

bool GF_pBN::isZ(const bignum &a) const { return a == 0; }

bool GF_pBN::isE(const bignum &a) const { return a == 1; }

bignum GF_pBN::import(const bignum &a) const {
  return (a + getSize()) % getSize();
}

bignum GF_pBN::getZ() const { return 0; }

bignum GF_pBN::getE() const { return 1; }

bignum GF_pBN::importu32(u32 a) const { return bignum(a) % getSize(); }

bignum GF_pBN::getRandRaw() const { return getSize().rand(); }

OPA_NAMESPACE_DECL3_END

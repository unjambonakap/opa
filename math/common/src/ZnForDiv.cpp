#include "common/ZnForDiv.h"

OPA_NAMESPACE_DECL3(opa, math, common)

NumForDiv ZnForDiv::mul(const NumForDiv &a, const NumForDiv &b) const {
  NumForDiv res = getZ();;
  res.rem = zn.mul(a.rem, b.rem);
  REP (i, a.factors.size())
    res.factors[i] =a.factors[i] + b.factors[i];
  return res;
}

NumForDiv ZnForDiv::add(const NumForDiv &a, const NumForDiv &b) const {
  OPA_CHECK0(false);
  NumForDiv res;
  return res;
  // REP (i, a.factors.size()) {
  //  int mx = std::max(a.factors[i], b.factors[i]);
  //  res.factors.push_back(std::min(a.factors[i], b.factors[i]));
  //}
  // return ((u64)a + b) % getSizeU32();
}
NumForDiv ZnForDiv::neg(const NumForDiv &a) const {
  OPA_CHECK0(false);
  return a;
}

bool ZnForDiv::isInv(const NumForDiv &a) const {
  if (a.rem == 0) return false;
  for (auto &f : a.factors)
    if (f != 0) return false;
  return true;
}

NumForDiv ZnForDiv::inv(const NumForDiv &a) const { 
  NumForDiv res = getZ();
  res.rem = zn.inv(a.rem);
  return res;
}

bool ZnForDiv::isZ(const NumForDiv &a) const {
  if (a.rem != 0) return false;
  return isInv(a);
}

bool ZnForDiv::isE(const NumForDiv &a) const {
  if (a.rem != 1) return false;
  for (auto &f : a.factors)
    if (f != 0) return false;
  return true;
}

NumForDiv ZnForDiv::import(const NumForDiv &a) const {
  return a;
}

NumForDiv ZnForDiv::getZ() const {
  NumForDiv res;
  res.rem = 0;
  res.factors.resize(m_factors.size(), 0);
  return res;
}
NumForDiv ZnForDiv::getE() const {
  NumForDiv res;
  res.rem = 1;
  res.factors.resize(m_factors.size(), 0);
  return res;
}

bool ZnForDiv::ediv(const NumForDiv &a, const NumForDiv &b, NumForDiv *q,
                    NumForDiv *r) const {
  OPA_CHECK0(0);
  return false;
}

NumForDiv ZnForDiv::div(const NumForDiv &a, const NumForDiv &b) const {
  if (isZ(a)) return getZ();
  NumForDiv res = getZ();;
  REP (i, a.factors.size()) {
    res.factors[i] =a.factors[i] - b.factors[i];
    OPA_CHECK0(res.factors.back() >= 0);
  }
  res.rem = zn.div(a.rem, b.rem);
  return res;
}

NumForDiv ZnForDiv::importu32(u32 a) const {
  a %= getSizeU32();
  NumForDiv res = getZ();
  REP (i, m_factors.size()) {
    while (a % m_factors[i].first == 0)
      res.factors[i]++, a /= m_factors[i].first;
  }
  res.rem = a;
  return res;
}

NumForDiv ZnForDiv::getRand() const { return importu32(rng()); }

bool ZnForDiv::compareRank(const NumForDiv &a, const NumForDiv &b) const {
  OPA_CHECK0(false);
  return false;
}

u32 ZnForDiv::to_u32(const NumForDiv &a) const {
  u32 res = a.rem;
  REP (i, a.factors.size()) {
    OPA_CHECK(a.factors[i] >= 0, a);
    REP (j, a.factors[i]) { res = zn.mul(res, m_factors[i].first); }
  }
  return res;
}

OPA_NAMESPACE_DECL3_END

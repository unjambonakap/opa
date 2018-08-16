#pragma once

#include <opa/math/common/Ring.h>
#include <opa/math/common/Utils.h>
#include <opa/math/common/Zn.h>

OPA_NAMESPACE_DECL3(opa, math, common)

struct NumForDiv {
  std::vector<int> factors;
  u32 rem;
  OPA_DECL_COUT_OPERATOR2(NumForDiv, a.rem, a.factors);
  OPA_DECL_EQ_OPERATOR(NumForDiv, factors, rem);
};

class ZnForDiv : public Ring<NumForDiv> {
public:
  ZnForDiv(u32 n) : Ring<NumForDiv>(n, n), zn(n) {
    activateU32();
    m_factors = factor_small(n);
    OPA_DISP0(m_factors);
  }
  virtual bool compareRank(const NumForDiv &a, const NumForDiv &b) const;

  virtual NumForDiv mul(const NumForDiv &a, const NumForDiv &b) const;
  virtual NumForDiv add(const NumForDiv &a, const NumForDiv &b) const;
  virtual NumForDiv neg(const NumForDiv &a) const;

  virtual bool isInv(const NumForDiv &a) const;
  virtual NumForDiv inv(const NumForDiv &a) const;

  virtual bool isZ(const NumForDiv &a) const;
  virtual bool isE(const NumForDiv &a) const;
  virtual NumForDiv import(const NumForDiv &a) const;

  virtual NumForDiv getZ() const;
  virtual NumForDiv getE() const;
  virtual NumForDiv importu32(u32 a) const;
  virtual bool ediv(const NumForDiv &a, const NumForDiv &b, NumForDiv *q,
                    NumForDiv *r) const;
  virtual NumForDiv div(const NumForDiv &a, const NumForDiv &b) const;

  virtual NumForDiv getRand() const;
  u32 to_u32(const NumForDiv &a) const;

  ~ZnForDiv() {}

  S64Factors m_factors;
  Zn zn;
};

OPA_NAMESPACE_DECL3_END

#pragma once

#include <opa/math/common/Utils.h>
#include <opa/math/common/base.h>
#include <opa/math/common/bignum.h>
#include <opa/utils/serialize.h>

OPA_NM_MATH_COMMON

class SmallestDivHelper {
public:
  bignum compute(const BGFactors &factors, const bignum &lb);

  struct State {
    int pos;
    std::vector<bignum> factors;
    bignum v;
    bool operator>(const State &other) const { return v > other.v; }
  };

  State res;
};

static std::vector<bignum> factor_known(const bignum &x, const BGFactors &factors){
  std::vector<bignum> res;
  bignum tmp = x;
  for (const auto &factor :factors){
    while (tmp%factor.first == 0) {
      tmp /= factor.first;
      res.push_back(factor.first);
    }
  }

  return res;
}

class SmoothPrimeOrderHelper {
public:
  typedef std::pair<u64, S64Factors> Entry;

  SmoothPrimeOrderHelper(int lb, int ub, const std::vector<int> &prime_set);
  void build() const;
  const std::vector<Entry> &data() const {
    build();
    return m_data;
  }

  int m_lb;
  int m_ub;
  std::vector<int> m_prime_set;
  mutable std::vector<Entry> m_data;
  mutable bool m_init = false;
};

extern SmoothPrimeOrderHelper spoh_for_gfp_split;

OPA_NM_MATH_COMMON_END

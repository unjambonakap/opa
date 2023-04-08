#pragma once

#include <opa/math/common/PolyModRing.h>
#include <opa/math/common/PolyRing.h>
#include <opa/math/common/Ring.h>
#include <opa/math/common/UtilsPoly.h>
#include <opa/math/common/base.h>

OPA_NM_MATH_COMMON

template <typename T> class Ring_col : public PolyModRing<T> {
public:
  typedef Poly<T> U;

private:
  std::vector<std::pair<Poly<T>, int> > m_factors;

public:
  Ring_col(const Field<T> *base_field, const U &p) { init(base_field, p); }
  Ring_col() {}
  ~Ring_col() {}

  void init(const Field<T> *base_field, const U &p) {
    PolyRing<T> tmp_pr(base_field);
    m_factors = tmp_pr.factor2(p);
    order = 1;
    for (auto &e : m_factors) {
      bignum p = t_faste_u32(base_field->getSize(), e.first.deg());
      order *= (p - 1) * p.pow(e.second - 1);
    }
    PolyModRing<T>::init(base_field, p);
  }
  bignum order;
};

OPA_NM_MATH_COMMON_END

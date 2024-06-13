#pragma once

#include <opa_common.h>
#include <opa/math/common/UtilsPoly.h>
#include <opa/math/common/Types.h>
#include <opa/math/common/FractionField.h>
#include <opa/math/common/Matrix.h>
#include <opa/math/common/GF_pBN.h>
#include <opa/math/common/GF_q.h>
#include <opa/math/common/Zn_BG.h>

OPA_NAMESPACE(opa, math, adv)
typedef OPA_MATH::bignum bignum;
template<class T>
class Group {
  public:
  virtual T e() const = 0;
  virtual T mul(const T &a, const T &b) const = 0;
  virtual T inv(const T &a) const = 0;


};

template<class LA_t, LG_t> 
class Lie {
  public:
  Group<LG_t> lie_group;


  virtual LG_t exp(const LA_t) const = 0;


};

OPA_NAMESPACE_END(opa, math, adv)

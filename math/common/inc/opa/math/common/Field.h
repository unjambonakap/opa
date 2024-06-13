#pragma once

#include <opa/math/common/Ring.h>
#include <opa/math/common/Utils.h>
#include <opa/math/common/base.h>

OPA_NM_MATH_COMMON

template <class T> class Field : public Ring<T> {
public:
  void init(const bignum &size, const bignum &car) {
    Ring<T>::init(size, car);
    init();
  }

  void init() {
    if (m_factors.empty() && this->getSize() > 0 && this->getSize() < (1ll << 31) - 1)
      this->compute_order_factors();
    m_hasPrimElem = false;
  }

  Field(const bignum &size, const bignum &car) { init(size, car); }
  Field() {}

  virtual ~Field() {}

  // interface

  virtual T inv(const T &a) const = 0;
  virtual T mul(const T &a, const T &b) const = 0;
  virtual T add(const T &a, const T &b) const = 0;
  virtual T neg(const T &a) const = 0;

  virtual bool isZ(const T &a) const = 0;
  virtual bool isE(const T &a) const = 0;
  virtual T import(const T &a) const = 0;

  virtual T getZ() const = 0;
  virtual T getE() const = 0;
  virtual T getRandRaw() const {
    OPA_ABORT0(true);
    return T();
  }

  // EOI

  virtual bool isInv(const T &a) const { return !isZ(a); }
  virtual bool ediv(const T &a, const T &b, T *q, T *r) const {
    if (isZ(b)) return false;
    if (q) *q = mul(a, inv(b));
    if (r) *r = getZ();
    return true;
  }

  virtual bool compareRank(const T &a, const T &b) const { return !isZ(a) < !isZ(b); }

  virtual T getPrimitiveElem() const {
    if (m_hasPrimElem) return m_cachedPrimElem;

    while (true) {
      T x = getRandRaw();

      if (is_prim_elem(x)) {
        m_cachedPrimElem = x;
        m_hasPrimElem = true;
        return x;
      }
    }
  }

  virtual T getNthRoot(u32 n) const {
    OPA_CHECK0((this->getSize() - 1) % n == 0);
    return this->faste(getPrimitiveElem(), (this->getSize() - 1) / n);
  }

  virtual T getRand() const {
    bignum x = this->getSize().rand();
    if (x == 0) return getZ();
    return this->faste(getPrimitiveElem(), x - 1);
  }

  virtual bignum getOrderOf(const T &a) const {
    OPA_CHECK0(!isZ(a));

    bignum cur = this->getSize() - 1;
    for (auto f : m_factors) {
      REP (j, f.ND) {
        if (!isE(this->faste(a, cur / f.ST))) break;
        cur /= f.ST;
      }
    }
    return cur;
  }

  virtual bool isField() const { return true; }

  void set_order_factors(const BGFactors &factors) { m_factors = factors; }

  const BGFactors &get_factors() const {
    if (m_factors.size() == 0) {
      compute_order_factors();
    }
    return m_factors;
  }

  void compute_order_factors() const {
    if (m_factors.size()) return;
    m_factors = factor_large(this->getSize() - 1);
  }

  bool is_prim_elem(const T &x) const {
    if (isZ(x)) return false;

    for (auto u : m_factors) {
      T nx = this->faste(x, this->getSize() / u.ST);
      if (isE(nx)) return false;
    }
    return true;
  }

  void set_prim_elem(const T &x_) {
    T x = this->import(x_);
    OPA_CHECK0(is_prim_elem(x));
    m_cachedPrimElem = x;
    m_hasPrimElem = true;
  }

public:
  mutable BGFactors m_factors;
  mutable T m_cachedPrimElem;
  mutable bool m_hasPrimElem;
};

template <class T> class GF_pT : public Field<T> {
public:
  T n;

  GF_pT() {}
  GF_pT(T n, bool init = true) : n(n) {
    if (init) this->init(bignum(n), bignum(n));
  }
  template <class U> T importT(const U &a) const { return T(a % U(n)); }
};

OPA_NM_MATH_COMMON_END

#pragma once

#include <opa/math/common/GF_p.h>
#include <opa/math/common/PolyRing.h>
#include <opa/math/common/Ring.h>
#include <opa/math/common/Utils.h>
#include <opa/math/common/base.h>
#include <opa/math/common/bignum.h>

OPA_NM_MATH_COMMON

template <class T> class Poly;

template <class T, class BaseType = Ring<Poly<T> > >
class PolyModRing : public PolyRing<T, BaseType> {
public:
  typedef PolyRing<T, BaseType> PolyRingBase;
  Poly<T> m_mod;
  PolyRingBase m_pur_pr;

public:
  typedef T XType;
  typedef Poly<T> PT;

  void init(const Ring<T> *ring, const Poly<T> &mod, bignum order = -1) {
    bignum size = -1;
    if (ring->getCar() != -1) size = ring->getCar().pow(mod.deg());
    m_pur_pr.init(ring);
    m_mod = to_pur(mod);
    PolyRing<T, BaseType>::init(ring, size);
  }

  PolyModRing() {}
  PolyModRing(const Ring<T> *ring, const Poly<T> &mod) { init(ring, mod); }

  const PolyRingBase *pur_pr() const { return m_pur_pr; }

  virtual bool isInv(const Poly<T> &a) const {
    OPA_DISP0(m_pur_pr.gcd(to_pur(a), m_mod), a, m_mod);
    Poly<T> tmp =  m_pur_pr.gcd(to_pur(a), m_mod);
    return (tmp.deg() == 0 && this->m_ring->isInv(tmp[0]));
  }

  virtual Poly<T> inv(const Poly<T> &a) const {
    Poly<T> u, v;
    m_pur_pr.egcd(to_pur(a), m_mod, u, v);
    return sto_ipur(u);
  }

  virtual Poly<T> mul(const Poly<T> &a, const Poly<T> &b) const {
    Poly<T> res;
    res = m_pur_pr.mulmod(to_pur(a), to_pur(b), m_mod);
    return sto_ipur(res);
  }

  bool ediv(const Poly<T> &a, const Poly<T> &b, Poly<T> *q, Poly<T> *r) const {
    if (b.deg() < 0) return false;
    if (!isInv(b)) return false;

    if (q) {
      *q = this->mul(b, this->inv(b));
    }
    if (r) *r = this->getZ();
    return true;
  }

  Poly<T> import_base(OPA_BG num) const {
    return import(m_pur_pr.import_base(num));
  }

  Poly<T> sto_pur(Poly<T> &res) const {
    res.unsafe_change_ring(&m_pur_pr);
    return res;
  }

  Poly<T> sto_ipur(Poly<T> &res) const {
    res.unsafe_change_ring(this);
    return res;
  }

  Poly<T> to_pur(const Poly<T> &t) const {
    Poly<T> res = t;
    res.unsafe_change_ring(&m_pur_pr);
    return res;
  }

  Poly<T> to_ipur(const Poly<T> &t) const {
    Poly<T> res = t;
    res.unsafe_change_ring(this);
    return res;
  }

  Poly<T> &reduce_pur(Poly<T> &p) const {
    if (p.deg() >= m_mod.deg()) {
      p = m_pur_pr.mod(sto_pur(p), m_mod);
      sto_ipur(p);
    }
    return p;
  }

  Poly<T> &reduce_ipur(Poly<T> &p) const {
    if (p.deg() >= m_mod.deg()) {
      sto_pur(p);
      reduce_pur(p);
    }
    return p;
  }

  Poly<T> import(const Poly<T> &p) const {
    Poly<T> tmp = m_pur_pr.import(p);
    return reduce_pur(tmp);
  }

  Poly<T> import(const std::vector<T> &px, bool rev = false) const {
    Poly<T> tmp = m_pur_pr.import(px, rev);
    return reduce_pur(tmp);
  }

  Poly<T> &set1(Poly<T> &p, int deg, const T &get) const {
    PolyRingBase::set1(p, deg, get);
    return reduce_ipur(p);
  }

  virtual Poly<T> getRand() const { return this->rand(m_mod.deg() - 1); }

  Poly<T> xpwv(int d, const T &v) const {
    PT tmp = PolyRingBase::xpwv(d, v);
    return reduce_ipur(tmp);
  }

  Poly<T> xpw(int d) const {
    PT tmp = PolyRingBase::xpw(d);
    return reduce_ipur(tmp);
  }
  Poly<T> getRandRaw() const { return this->randDim(this->m_mod.deg()); }
};

template <class T>
class PolyModField : public PolyModRing<T, Field<Poly<T> > > {
public:
  void init(const Ring<T> *ring, const Poly<T> &mod) {
    PolyModRing<T, Field<Poly<T> > >::init(ring, mod);
  }

  virtual Poly<T> inv(const Poly<T> &a) const {
    Poly<T> tmpa = this->to_pur(a);
    Poly<T> u, v, d;
    d = this->m_pur_pr.egcd(a, this->m_mod, u, v);
    OPA_CHECK(this->m_pur_pr.isE(d), a, this->m_mod);
    return this->sto_ipur(u);
  }

  virtual bool ediv(const Poly<T> &a, const Poly<T> &b, Poly<T> *q,
                    Poly<T> *r) const {
    return Field<Poly<T> >::ediv(a, b, q, r);
  }
};

OPA_NM_MATH_COMMON_END

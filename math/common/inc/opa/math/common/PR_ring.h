#pragma once

#include <opa/math/common/Poly.h>
#include <opa/math/common/PolyRing.h>
#include <opa/math/common/Ring.h>
#include <opa/math/common/UtilsPoly.h>
#include <opa/math/common/base.h>

OPA_NM_MATH_COMMON

template <typename T> class PR_ring : public Ring<Poly<T> > {
public:
  typedef Poly<T> U;

protected:
  const Field<T> *m_base_field;
  PolyRing<T> m_pr;
  U m_modPoly;
  bignum m_order;

public:
  void init(const Field<T> *base_field, const U &p, const bignum &order);

  virtual const void *get_underlying_ring() const override {
    return m_base_field;
  }
  virtual const void *get_poly_ring() const override {
    return nullptr; /*todo: reenable */
  }

  PR_ring(const Field<T> *base_field, const U &p, const bignum &order);
  PR_ring() {}
  ~PR_ring() {}
  U theta() const { return m_pr->x(); }

  U from_base_field(const T &a) const {
    return m_pr.constant(m_base_field->import(a));
  }

  OPA_BG to_base(const U &p) const { return m_pr.to_base(p); }
  U import_base(const OPA_BG &num) const { return m_pr.import_base(num); }
  virtual U mul(const U &a, const U &b) const;
  virtual U add(const U &a, const U &b) const;
  virtual U inv(const U &a) const;
  virtual bool isInv(const U &a) const { OPA_CHECK0(false); }

  virtual U neg(const U &a) const;

  virtual bool isZ(const U &a) const;
  virtual bool isE(const U &a) const;
  virtual U import(const U &a) const;

  virtual U import_vec(const std::vector<T> &a) const {
    return this->import(this->m_pr.import(a));
  }

  virtual U importu32(u32 a) const;
  virtual int get_poly_pos() const override { return m_pr.get_poly_pos(); }

  virtual U getZ() const;
  virtual U getE() const;
  virtual U getRandRaw() const;
  virtual bool ediv(const U &a, const U &b, U *q, U *r) const {
    OPA_CHECK0(false);
  }

  virtual bool compareRank(const U &a, const U &b) const { OPA_CHECK0(false); }

  const Field<T> *getBaseField() const { return m_base_field; }
  const U &getModPoly() const { return m_modPoly; }
  OPA_ACCESSOR_R(PolyRing<T>, m_pr, pr);
  OPA_ACCESSOR_R(bignum, m_order, order);
};

template <class T> Poly<T> PR_ring<T>::getRandRaw() const {
  return import(m_pr.randDim(m_modPoly.deg()));
}

template <class T>
PR_ring<T>::PR_ring(const Field<T> *base_field, const Poly<T> &p,
                    const bignum &order) {
  this->init(base_field, p, order);
}

template <class T>
void PR_ring<T>::init(const Field<T> *base_field, const Poly<T> &p,
                      const bignum &order) {
  m_modPoly = p;
  m_base_field = base_field;

  m_order = order;
  m_pr.init(base_field);
  m_pr.set_cur_ring(this);
  Ring<Poly<T> >::init(t_faste_u32(base_field->getSize(), p.deg()),
                       base_field->getCar().getu32());
}

template <class T>
Poly<T> PR_ring<T>::mul(const Poly<T> &a, const Poly<T> &b) const {
  return m_pr.mulmod(a, b, m_modPoly);
}

template <class T>
Poly<T> PR_ring<T>::add(const Poly<T> &a, const Poly<T> &b) const {
  return m_pr.addmod(a, b, m_modPoly);
}

template <class T> Poly<T> PR_ring<T>::inv(const Poly<T> &a) const {
  OPA_CHECK0(order() != 0);
  Poly<T> res = m_pr.faste(a, m_order - 1, m_modPoly);
  // OPA_DISP0(m_order, res, a, m_modPoly, m_pr.mulmod(res, a, m_modPoly));
  OPA_CHECK0(m_pr.mulmod(res, a, m_modPoly) == getE());
  return res;
}

template <class T> Poly<T> PR_ring<T>::neg(const Poly<T> &a) const {
  return m_pr.neg(a);
}

template <class T> bool PR_ring<T>::isZ(const Poly<T> &a) const {
  return m_pr.isZ(a);
}

template <class T> bool PR_ring<T>::isE(const Poly<T> &a) const {
  return m_pr.isE(a);
}

template <class T> Poly<T> PR_ring<T>::getZ() const { return m_pr.getZ(); }

template <class T> Poly<T> PR_ring<T>::getE() const { return m_pr.getE(); }

template <class T> Poly<T> PR_ring<T>::import(const Poly<T> &a) const {
  return m_pr.mod(m_pr.import(a), m_modPoly);
}

template <class T> Poly<T> PR_ring<T>::importu32(u32 a) const {
  return m_pr.constant(m_base_field->importu32(a));
}

OPA_NM_MATH_COMMON_END

#pragma once

#include <opa/math/common/Field.h>
#include <opa/math/common/Poly.h>
#include <opa/math/common/PolyModRing.h>
#include <opa/math/common/PolyRing.h>
#include <opa/math/common/UtilsPoly.h>

OPA_NAMESPACE_DECL3(opa, math, common)

template <typename T> class GF_q : public PolyModField<T> {
private:
  const Field<T> *m_base_field;
  Poly<T> m_modPoly;
  bignum m_qm1;
  std::vector<Poly<T> > pw_cache;

public:
  void construct(const Poly<T> &p) {
    m_modPoly = p;
    m_qm1 = this->getSize() - 1;
  }

  GF_q(const Field<T> *base_field, int m) { this->init(base_field, m); }

  GF_q(const Field<T> *base_field, const Poly<T> &p) {
    this->init(base_field, p);
  }

  void init(const Field<T> *base_field, int m) {
    puts("try ot find ireed");
    Poly<T> irred = find_primitive_poly(*base_field, m);
    this->init(base_field, irred);
  }

  void init(const Field<T> *base_field, const Poly<T> &p) {
    PolyModField<T>::init(base_field, p);
    construct(p);
    REP (i, 2 * this->m_mod.size() - 1) { pw_cache.pb(this->import(this->xpw(i))); }
  }

  T get_norm(const Poly<T> &a) const {
    return this->faste(a, m_qm1 / (this->getCar() - 1)).get_safe(0);
  }

  GF_q() {}
  ~GF_q() {}

  Poly<T> theta() const { return this->x(); }

  Poly<T> from_base_field(const T &a) const {
    return this->constant(m_base_field->import(a));
  }

  virtual Poly<T> mul(const Poly<T> &a, const Poly<T> &b) const {
    return PolyModField<T>::mul(a, b);
    Poly<T> res = this->getZ();

    /*
    REP (i, a.size() + b.size() - 1) {
      T v = m_base_field->getZ();
      FOR (j, std::max(0, i - b.deg()), std::min(a.size(), i)) {
        v = m_base_field->add(v, m_base_field->mul(b[i - j], a[j]));
      }

      Poly<T> tmp = pw_cache[i];
      m_pr.sadd(res, m_pr.smulc(tmp, v));
    }
    return res;
    */
  }

  virtual Poly<T> inv(const Poly<T> &a) const {
    Poly<T> res = this->faste(a, m_qm1 - 1, m_modPoly);
    OPA_CHECK0(this->mul(res, a) == this->getE());
    return res;
  }

  const Field<T> *getBaseField() const {
    return (const Field<T> *)this->get_underlying_ring();
  }
  const Poly<T> &getModPoly() const { return this->m_mod; }
};

OPA_NAMESPACE_DECL3_END

#pragma once

#include <opa/math/common/Field.h>

OPA_NAMESPACE_DECL3(opa, math, common)

template <class T> class Fraction {

public:
  typedef Fraction<T> SelfType;
  typedef Field<SelfType> SelfField;

  Fraction<T>() {}
  Fraction<T>(const T &p, const T &q, const SelfField *field) {
    this->field = field;
    this->set(p, q);
  }

  Fraction<T>(const T &p, const SelfField *field) {
    this->field = field;
    this->set(p);
  }

  void set(const T &p) {
    this->p = p;
    this->q = field->getE().q;
  }

  void set(const T &p, const T &q) {
    this->p = p;
    this->q = q;
  }
  Fraction<T> clone() const { return Fraction<T>(p, q, field); }
  Fraction<T> create(const T &p, const T &q) const {
    return Fraction<T>(p, q, field);
  }
  Fraction<T> create(const T &p) const { return Fraction<T>(p, field); }

  bool operator==(const SelfType &x) const { return p == x.p && q == x.q; }
  bool operator<(const SelfType &x) const { return this->field->lt(*this, x); }
  SelfType operator-() const { return field->neg(*this); }

#define OPA_DEFINE_OP(typ, op, func)                                           \
  typ operator op(const typ &a) const { return func(*this, a); }
  OPA_DEFINE_OP(SelfType, +, this->field->add);
  OPA_DEFINE_OP(SelfType, -, this->field->sub);
  OPA_DEFINE_OP(SelfType, *, this->field->mul);
  OPA_DEFINE_OP(SelfType, /, this->field->div);

#define OPA_DEFINE_CMP_OPS_FROM_LT_AND_EQ(typ)                                 \
  bool operator!=(const typ &a) const { return !(*this == a); }                \
  bool operator>(const typ &a) const { return !(*this <= a); }                 \
  bool operator<=(const typ &a) const { return (*this == a || *this < a); }    \
  bool operator>=(const typ &a) const { return !(*this < a); }

  OPA_DEFINE_CMP_OPS_FROM_LT_AND_EQ(SelfType);

  friend std::ostream &operator<<(std::ostream &os, const SelfType &a) {
    if (a.field == nullptr){
      os << "fraction_not_init";
      return os;
    }
    if (a.q == a.field->getE().q)
      os << a.p;
    else {
      os << "(";
      os << a.p;
      os << " / ";
      os << a.q;
      os << ")";
    }
    return os;
  }
  T p;
  T q;
  const SelfField *field = nullptr;
};

template <class T> class FractionField : public Field<Fraction<T> > {

public:
  typedef Fraction<T> QT;
  typedef T BaseType;
  FractionField() {}

  FractionField(const Ring<T> *base_ring, bool do_reduce = true) {
    this->init(base_ring, do_reduce);
  }

  void init(const Ring<T> *base_ring, bool do_reduce = true) {
    Field<QT>::init(-1, base_ring->getCar());
    m_base_ring = base_ring;
    m_do_reduce = do_reduce;
  }

  const Ring<T> *get_base_ring() const { return m_base_ring; }

  virtual QT getPrimitiveElem() const {
    assert(0);
    return QT();
  }

  virtual QT inv(const QT &a) const { return reduce(create(a.q, a.p)); }

  virtual QT mul(const QT &a, const QT &b) const {
    QT res =
      this->create(m_base_ring->mul(a.p, b.p), m_base_ring->mul(a.q, b.q));
    return reduce(res);
  }

  virtual QT add(const QT &a, const QT &b) const {
    T p, q;
    q = m_base_ring->mul(a.q, b.q);
    p =
      m_base_ring->add(m_base_ring->mul(a.p, b.q), m_base_ring->mul(b.p, a.q));
    QT res = this->create(p, q);
    return reduce(res);
  }

  virtual QT neg(const QT &a) const {
    return this->import(m_base_ring->neg(a.p), a.q);
  }

  virtual QT importu32(u32 a) const {
    T p = m_base_ring->importu32(a);
    return this->import(p, m_base_ring->getE());
  }

  virtual QT import(const QT &a) const { return this->import(a.p, a.q); }

  virtual QT create(const T &p, const T &q) const {
    QT res = QT(p, q, this);
    return res;
  }

  virtual QT importq(const T &q) const {
    return this->import(m_base_ring->getE(), q);
  }

  virtual QT import(const T &p) const {
    return this->import(p, m_base_ring->getE());
  }

  virtual QT import(const T &p, const T &q) const {
    QT res = this->create(m_base_ring->import(p), m_base_ring->import(q));
    return reduce(res);
  }

  virtual bool isZ(const QT &a) const { return m_base_ring->isZ(a.p); }
  virtual bool isE(const QT &a) const {
    return m_base_ring->isE(a.p) && m_base_ring->isE(a.q);
  }

  virtual QT getZ() const {
    return create(m_base_ring->getZ(), m_base_ring->getE());
  }
  virtual QT getE() const {
    return create(m_base_ring->getE(), m_base_ring->getE());
  }
  virtual int get_poly_pos() const { return m_base_ring->get_poly_pos(); }

  QT import_integer(const T &a) const { return create(a, m_base_ring->getE()); }

  bool is_integer(const QT &a) const { return m_base_ring->isE(a.q); }
  T to_base_or_fail(const QT &a) const {
    OPA_CHECK(is_integer(a), a);
    return a.p;
  }
  T project(const QT &a) const { return to_base_or_fail(a); }

  T integer_part(const QT &a) const { return m_base_ring->div(a.p, a.q); }
  T floor(const QT &a) const { return m_base_ring->div(a.p, a.q); }
  T ceil(const QT &a) const {
    return m_base_ring->div(m_base_ring->add(a.p, m_base_ring->add(a.q, m_base_ring->getE())),
                            a.q);
  }

  QT rational_part(const QT &a) const {
    return create(m_base_ring->mod(a.p, a.q), a.q);
  }

  std::vector<QT> import_vec(const std::vector<T> &a) const {
    std::vector<QT> res;
    for (auto &x : a) res.push_back(this->import_integer(x));
    return res;
  }

  std::vector<T> to_base_or_fail(const std::vector<QT> &a) const {
    std::vector<T> res;
    for (auto &x : a) res.push_back(this->to_base_or_fail(x));
    return res;
  }

  virtual QT getRand() const  {
    return reduce(create(m_base_ring->getRand(), m_base_ring->getRandNZ()));
  }

  virtual bool lt(const QT &a, const QT &b) const  {
    T l = m_base_ring->mul(a.p, b.q);
    T r = m_base_ring->mul(a.q, b.p);
    return m_base_ring->lt(l, r);
  }

  QT reduce(const QT &f) const {
    if (!m_do_reduce) return f;
    if (m_base_ring->isZ(f.p)) return getZ();

    T d = m_base_ring->gcd(f.p, f.q);
    T p = m_base_ring->div(f.p, d);
    T q = m_base_ring->div(f.q, d);

    if (m_base_ring->isInv(q)) {
      T iq = m_base_ring->inv(q);
      q = m_base_ring->getE();
      p = m_base_ring->mul(iq, p);
    }
    m_base_ring->norm_fraction(p, q);
    return create(p, q);
  }

  QT operator()(const T &a) const { return import(a); }
  QT operator()(const T &a, const T &b) const { return import(a, b); }

private:
  bool m_do_reduce = true;
  const Ring<T> *m_base_ring = nullptr;
};

OPA_NAMESPACE_DECL3_END

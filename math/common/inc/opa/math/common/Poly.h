#pragma once

#include <opa/math/common/PolyRing.h>
#include <opa/math/common/Ring.h>
#include <opa/math/common/Utils.h>
#include <opa/math/common/UtilsRing.h>
#include <opa/math/common/base.h>

OPA_NM_MATH_COMMON

template <class T> class PolyRingOps;

std::string get_poly_letter(int id);

template <class T> class Poly {
public:
  // It cannot be PolyRing<T, *> because of second template argument
  typedef Ring<Poly<T>> RingType;

private:
  friend class PolyRing<T>;
  friend class PolyRing<T, Field<Poly<T> > >;
  std::vector<T> poly;
  const RingType *m_ring = nullptr;

public:
  void init(const RingType *ring);

  const Ring<T> *get_underlying_ring() const {
    return static_cast<const Ring<T> *>(m_ring->get_underlying_ring());
  }

  const PolyRingOps<T> *get_poly_ring() const {
    return static_cast<const PolyRingOps<T> *>(m_ring->get_poly_ring());
  }
  void unsafe_change_ring(const RingType *nring) { m_ring = nring; }

  Poly() { m_ring = 0; }
  Poly(const RingType *ring) { init(ring); }

  //Poly(u32 a) {
  //  OPA_ASSERT0(a == 0 || a == 1);
  //  if (a == 1) poly.pb(1);
  //}
  Poly(const RingType *ring, u32 a) {
    init(ring);
    OPA_ASSERT0(a == 0 || a == 1);
    if (a == 1) poly.pb(get_underlying_ring()->getE());
  }

  Poly<T> &operator=(const Poly<T> &a) {
    copy_poly(a);
    return *this;
  }

  T operator()(const T &a) const;

  void copy_poly(const Poly<T> &a);

  Poly(const Poly<T> &a) { copy_poly(a); }
  Poly(const RingType *ring, const std::vector<T> &px) {
    init(ring);
    poly = px;
  }
  ~Poly() {}

  const RingType *get_ring() const { return m_ring; }
  typename std::vector<T>::iterator begin() { return poly.begin(); }
  typename std::vector<T>::iterator end() { return poly.end(); }
  typename std::vector<T>::const_iterator begin() const { return poly.begin(); }
  typename std::vector<T>::const_iterator end() const { return poly.end(); }

  int deg() const { return poly.size() - 1; }
  int size() const { return poly.size(); }

  T lc() const {
    if (poly.size() == 0) return get_underlying_ring()->getZ();
    return poly.back();
  }

  T get_safe(int a) const {
    OPA_ASSERT0(a >= 0);
    if (a <= deg()) return poly[a];
    return get_underlying_ring()->getZ();
  }
  const T &get(int a) const {
    OPA_ASSERT0(a >= 0 && a < poly.size());
    return poly[a];
  }
  const T &operator[](int a) const { return get(a); }
  T &operator[](int a) { return get(a); }

  T &get(int a) {
    OPA_ASSERT0(a >= 0 && a < poly.size());
    return poly[a];
  }

  void normalize() {
    get_poly_ring()->normalize(*this);
  }

  // forcing the user to normalize
  T &get_force(int a) {
    if (a >= poly.size()) {
      poly.resize(a + 1, get_underlying_ring()->getZ());
    }
    return get(a);
  }

  const T &operator&(int a) const { return get(a); }
  T &operator&(int a) { return get(a); }

  const std::vector<T> &toVector() const { return poly; }
  const std::vector<T> &to_vec() const { return poly; }
  std::vector<T> to_vec(int sz) const {
    return get_poly_ring()->toVector(*this, sz);
  }
  std::vector<T> &vec_unsafe() { return poly; }

  std::vector<T> toVector(int sz) const {
    return get_poly_ring()->toVector(*this, sz);
  }

  void disp() const { std::cout << *this << std::endl; }
  void mulx(const T &v) {
    T pw = get_underlying_ring()->getE();
    for (auto &x : poly){
      x = get_underlying_ring()->mul(x, pw);
      pw = get_underlying_ring()->mul(pw, v);
    }
  }

  bool operator==(const Poly<T> &x) const { return poly == x.poly; }
  bool operator!=(const Poly<T> &x) const { return poly != x.poly; }
  bool operator<(const Poly<T> &x) const {
    if (deg() != x.deg()) return deg() < x.deg();
    return poly < x.poly;
  }

  friend std::ostream &operator<<(std::ostream &os, const Poly<T> &a) {
    if (!a.m_ring) {
      os << "0";
      return os;
    }
    const Ring<T> *ur = a.get_underlying_ring();

    OPA_CHECK0(ur != nullptr);

    os << "(";
    assert(a.m_ring);
    std::string poly_letter = get_poly_letter(ur->get_poly_pos());
    for (int i = 0; i < a.size(); ++i) {
      if (ur->isZ(a.poly[i])) continue;
      os << a.poly[i] << '*' << poly_letter << "^" << i;
      if (i != a.size() - 1) os << ' ';
    }
    os << ")";
    return os;
  }

  Poly<T> operator+(const Poly<T> &b) const;
  Poly<T> operator-(const Poly<T> &b) const;
  Poly<T> operator*(const Poly<T> &b) const;
  Poly<T> operator/(const Poly<T> &b) const;
  Poly<T> operator%(const Poly<T> &b) const;
  Poly<T> operator-() const;
  Poly<T> operator+(const T &b) const;
  Poly<T> operator-(const T &b) const;
  Poly<T> operator*(const T &b) const;
  Poly<T> operator/(const T &b) const;
  Poly<T> monic() const;
  Poly<T> powm(const bignum &v) const;
  Poly<T> extract_pw(const std::vector<int> &lst) const;
  Poly<T> derivative() const;
  T linear_root() const;

  T cont() const;
  Poly<T> pp() const;
  void pp_and_cont(Poly<T> *pp, T *cont) const;
};

template <class T> Poly<T> get_unity(const Poly<T> &a) {
  return a.get_ring()->getE();
}

OPA_NM_MATH_COMMON_END

#include <opa/math/common/Poly.h>

OPA_NM_MATH_COMMON

template <class T> void Poly<T>::init(const RingType *ring) { m_ring = ring; }

template <class T> Poly<T> Poly<T>::operator+(const Poly<T> &b) const {
  return m_ring->add(*this, b);
}

template <class T> Poly<T> Poly<T>::operator-(const Poly<T> &b) const {
  return m_ring->sub(*this, b);
}

template <class T> Poly<T> Poly<T>::operator*(const Poly<T> &b) const {
  return m_ring->mul(*this, b);
}

template <class T> Poly<T> Poly<T>::operator/(const Poly<T> &b) const {
  return m_ring->div(*this, b);
}

template <class T> Poly<T> Poly<T>::operator%(const Poly<T> &b) const {
  return m_ring->mod(*this, b);
}

template <class T> Poly<T> Poly<T>::operator-() const {
  return m_ring->neg(*this);
}

template <class T> Poly<T> Poly<T>::operator+(const T &b) const {
  return get_poly_ring()->addc(*this, b);
}

template <class T> Poly<T> Poly<T>::operator-(const T &b) const {
  return get_poly_ring()->subc(*this, b);
}
template <class T> Poly<T> Poly<T>::operator*(const T &b) const {
  return get_poly_ring()->mulc(*this, b);
}

template <class T> Poly<T> Poly<T>::operator/(const T &b) const {
  return get_poly_ring()->divc(*this, b);
}

template <class T> Poly<T> Poly<T>::powm(const bignum &v) const {
  return m_ring->faste(*this, v);
}
template <class T> T Poly<T>::operator()(const T &a) const {
  return get_poly_ring()->eval(*this, a);
}

template <class T> Poly<T> Poly<T>::monic() const {
  return get_poly_ring()->monic(*this);
}

template <class T> Poly<T> Poly<T>::derivative() const {
  return get_poly_ring()->derivative(*this);
}

template <class T> T Poly<T>::linear_root() const {
  return get_poly_ring()->linear_root(*this);
}

template <class T> void Poly<T>::copy_poly(const Poly<T> &a) {
  this->poly = a.poly;
  this->m_ring = a.m_ring;
}

template <class T>
Poly<T> Poly<T>::extract_pw(const std::vector<int> &lst) const {
  Poly<T> res = m_ring->getZ();
  for (auto &x : lst) {
    res = res + get_poly_ring()->xpwv(x, get_safe(x));
  }
  return res;
}

template <class T> T Poly<T>::cont() const {
  T res;
  pp_and_cont(nullptr, &res);
  return res;
}

template <class T> Poly<T> Poly<T>::pp() const {
  Poly<T> res;
  pp_and_cont(&res, nullptr);
  return res;
}

template <class T>
void Poly<T>::pp_and_cont(Poly<T> *res_pp, T *res_cont) const {
  auto base_ring = this->get_underlying_ring();
  T cont = gcd_list(*base_ring, *this);
  if (res_cont != nullptr) *res_cont = cont;
  if (res_pp != nullptr) {
    *res_pp = *this / cont;
  }
}

OPA_NM_MATH_COMMON_END

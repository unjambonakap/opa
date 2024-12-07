#pragma once

#include <opa/math/common/Utils.h>
#include <opa/math/common/bignum.h>
#include <opa/utils/misc.h>
#include <opa/utils/range.h>

OPA_NM_MATH_COMMON
template <class T> T get_unity(const T &a) { return 1; }

template <typename T> struct GroupDesc {
  opa::utils::BinaryOp<T> mul;
  opa::utils::BinaryOp<T> div;
  T e;

  T multb(const std::vector<T> &a) const { return QQ::fold_left(a, e, mul); }
};

template <class T> class Ring : public virtual utils::Initable {
  bignum m_size;
  bignum m_car;
  u32 m_u32size;
  bool is_u32size;
  T m_me;

protected:
  void activateU32() {
    is_u32size = true;
    m_u32size = m_size.getu32();
  }

public:
  typedef T Type;
  u32 getSizeU32() const {
    OPA_CHECK0(is_u32size);
    return m_u32size;
  }

  virtual ~Ring() {}
  void init(const bignum &size, const bignum &car) {
    if (this->is_init()) return;
    Initable::init();
    m_size = size;
    m_car = car;
    is_u32size = false;
  }

  virtual const void *get_underlying_ring() const {
    OPA_CHECK0(0);
    return nullptr;
  }
  virtual const void *get_poly_ring() const {
    OPA_CHECK0(0);
    return nullptr;
  }
  GroupDesc<T> get_add_group() const {
    return GroupDesc<T>{ .mul = STD_FUNC2(this->add),
                         .div = STD_FUNC2(this->sub),
                         .e = this->getZ() };
  }

  GroupDesc<T> get_mul_group() const {
    return GroupDesc<T>{ .mul = STD_FUNC2(this->mul),
                         .div = STD_FUNC2(this->div),
                         .e = this->getE() };
  }

  Ring() {}
  Ring(const bignum &size, const bignum &car) { init(size, car); }
  virtual bignum export_base(const T &x) const {
    OPA_CHECK0(0);
    return bignum();
  }

  const bignum &getSize() const { return m_size; }
  bignum getCar() const { return m_car; }
  virtual bignum compute_order(const T &a, const BGFactors &factors);
  virtual T faste(T a, bignum p) const {
    T x = getE();

    for (; p > 0; p.srshift(1)) {
      if (p.get_bit(0)) x = mul(x, a);
      a = mul(a, a);
    }
    return x;
  }

  virtual T faste(T a, bignum p, const T &pmod) const {
    T x = getE();

    for (; p > 0; p.srshift(1)) {
      if (p.get_bit(0)) x = mod(mul(x, a), pmod);
      a = mod(mul(a, a), pmod);
    }
    return x;
  }

  virtual T div(const T &a, const T &b) const {
    T res, r;
    bool ok = ediv(a, b, &res, &r);
    OPA_CHECK(ok, a, b);
    return res;
  }

  virtual T mod(const T &a, const T &b) const {
    T res;
    bool ok = ediv(a, b, 0, &res);
    OPA_CHECK(ok, a, b);
    return res;
  }

  virtual T gcd(const T &a, const T &b) const;
  virtual T lcm(const T &a, const T &b) const {
    T d = this->gcd(a, b);
    return this->mul(this->div(a, d), b);
  }

  virtual T egcd(const T &a, const T &b, T &u, T &v) const;
  virtual T egcd2(const T &a, const T &b, T &u, T &v, T &u2, T &v2, bool &gcd_on_a) const {
    if (compareRank(a, b)) return egcd2(b, a, v, u, v2, u2, gcd_on_a ^= 1);
    T res = _egcd2(a, b, getE(), getZ(), getZ(), getE(), u, v, u2, v2, gcd_on_a);
    return res;
  }

  T sgn_v(const T &val, int pw_minus1) const { return pw_minus1 % 2 == 0 ? val : this->neg(val); }
  virtual bool isField() const { return false; }
  virtual T getRandRaw() const {
    OPA_CHECK0(0);
    return T();
  }

  virtual T getRandNZ() const {
    T tmp;
    while (true) {
      tmp = getRand();
      if (!isZ(tmp)) return tmp;
    }
  }

  T dot(const std::vector<T> &a, const std::vector<T> &b) const {
    T res = getZ();
    OPA_CHECK_EQ0(a.size(), b.size());
    REP (i, a.size()) res = this->add(res, this->mul(a[i], b[i]));
    return res;
  }

  // interface
  virtual bool ediv(const T &a, const T &b, T *q, T *r) const = 0;
  virtual bool compareRank(const T &a,
                           const T &b) const = 0; // return rank(a) < rank(b)

  virtual T mulv(const std::vector<T> &a) const { return this->mul(a); }

  virtual T mul(const std::vector<T> &a) const { return this->get_mul_group().multb(a); }
  virtual T add(const std::vector<T> &a) const { return this->get_add_group().multb(a); }

  virtual T mul(const T &a, const T &b) const = 0;
  virtual T add(const T &a, const T &b) const = 0;
  virtual T neg(const T &a) const = 0;
  virtual T sub(const T &a, const T &b) const { return add(a, neg(b)); }
  virtual bool isInv(const T &a) const = 0;
  virtual T inv(const T &a) const = 0;

  virtual bool isZ(const T &a) const = 0;
  virtual bool isE(const T &a) const = 0;
  virtual bool eq(const T &a, const T &b) const { return a == b; }
  virtual bool lt(const T &a, const T &b) const {
    OPA_CHECK0(0);

    return false;
  }

  virtual T import(const T &a) const = 0;

  virtual T getZ() const = 0;
  virtual T getE() const = 0;
  virtual T getME() const { return this->neg(this->getE()); }
  virtual T importu32(u32 a) const {
    OPA_CHECK0(false);
    return T();
  }
  virtual T import_bg(const bignum &a) const { return importu32(a.getu32()); }
  virtual T getRand() const {
    OPA_CHECK0(false);
    return T();
  }
  virtual int get_poly_pos() const {
    return 0;
  } // return id of poly variable
    // eoi

  virtual int sel_for_stab(const std::vector<T> &tb) const {
    REP (i, tb.size())
      if (this->isInv(tb[i])) return i;
    return -1;
  }

  virtual void norm_fraction(T &a, T &b) const {}
  virtual T abs(const T &a) const {
    OPA_CHECK0(false);
    return a;
  }
  virtual T sqrt(const T &a) const {
    OPA_CHECK0(false);
    return a;
  }

protected:
  T _gcd(const T &a, const T &b) const;
  T _egcd(const T &a, const T &b, const T &ua, const T &va, const T &ub, const T &vb, T &u,
          T &v) const;
  T _egcd2(const T &a, const T &b, const T &ua, const T &va, const T &ub, const T &vb, T &u, T &v,
           T &u2, T &v2, bool &gcd_on_a) const;
};

template <class T> T Ring<T>::gcd(const T &a, const T &b) const {
  if (compareRank(a, b)) return _gcd(b, a);
  return _gcd(a, b);
}

template <class T> T Ring<T>::_gcd(const T &a, const T &b) const {
  T xa = a;
  T xb = b;
  while (!isZ(xb)) {
    T r;
    bool res = ediv(xa, xb, 0, &r);
    // std::cout << b << std::endl;
    OPA_CHECK(res, xa, xb, res, a, b);
    xa = xb;
    xb = r;
  }
  return xa;
}

template <class T> T Ring<T>::egcd(const T &a, const T &b, T &u, T &v) const {
  if (compareRank(a, b)) return egcd(b, a, v, u);
  T res = _egcd(a, b, getE(), getZ(), getZ(), getE(), u, v);

  if (this->isInv(res)) {
    T tmp = this->inv(res);
    res = this->mul(tmp, res);
    u = this->mul(tmp, u);
    v = this->mul(tmp, v);
  }
  return res;
}

template <class T>
T Ring<T>::_egcd(const T &a, const T &b, const T &ua, const T &va, const T &ub, const T &vb, T &u,
                 T &v) const {
  if (isZ(b)) {
    u = ua;
    v = va;
    return a;
  }

  T r, q;
  bool res = ediv(a, b, &q, &r);
  assert(res);
  q = neg(q);

  return _egcd(b, r, ub, vb, add(ua, mul(ub, q)), add(va, mul(vb, q)), u, v);
}

template <class T>
T Ring<T>::_egcd2(const T &a, const T &b, const T &ua, const T &va, const T &ub, const T &vb, T &u,
                  T &v, T &u2, T &v2, bool &gcd_on_a) const {
  if (isZ(b)) {
    u = ua;
    v = va;
    return a;
  }

  T r, q;
  bool res = ediv(a, b, &q, &r);
  assert(res);
  q = neg(q);

  T d = _egcd2(b, r, ub, vb, add(ua, mul(ub, q)), add(va, mul(vb, q)), u, v, u2, v2, gcd_on_a ^= 1);
  if (d == b) {
    u2 = ua;
    v2 = va;
  }
  return d;
}

template <class T> bignum Ring<T>::compute_order(const T &a, const BGFactors &factors) {
  bignum cur = getSize() - 1;
  assert(isE(this->faste(a, cur)));

  for (const auto &factor : factors) {
    REP (i, factor.ND) {
      bignum tmp = cur / factor.ST;
      if (isE(this->faste(a, tmp)))
        cur = tmp;
      else
        break;
    }
  }
  return cur;
}

OPA_NM_MATH_COMMON_END

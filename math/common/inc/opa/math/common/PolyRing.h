#pragma once

#include <opa/math/common/GF_p.h>
#include <opa/math/common/Ring.h>
#include <opa/math/common/Utils.h>
#include <opa/math/common/base.h>
#include <opa/math/common/bignum.h>

OPA_NM_MATH_COMMON
constexpr bool kPolyRev = true;

template <class T> class Poly;

template <class T> class PolyRingOps {
public:
  virtual ~PolyRingOps<T>() {}

  virtual Poly<T> addc(const Poly<T> &p, const T &a) const = 0;
  virtual Poly<T> &normalize(Poly<T> &a) const = 0;

  virtual Poly<T> subc(const Poly<T> &p, const T &a) const = 0;
  virtual Poly<T> mulc(const Poly<T> &p, const T &a) const = 0;
  virtual Poly<T> divc(const Poly<T> &p, const T &a) const = 0;
  virtual std::vector<T> toVector(const Poly<T> &p, int sz) const = 0;
  virtual Poly<T> derivative(const Poly<T> &p) const = 0;
  virtual Poly<T> monic(const Poly<T> &p) const = 0;
  virtual T eval(const Poly<T> &p, const T &a) const = 0;
  virtual Poly<T> import_base(OPA_BG num) const = 0;
  virtual Poly<T> import(const Poly<T> &p) const = 0;
  virtual Poly<T> import(const std::vector<T> &px,
                         bool rev = !kPolyRev) const = 0;

  virtual Poly<T> xpwv(int d, const T &v) const = 0;
  virtual Poly<T> xpw(int d) const = 0;
  virtual Poly<T> get_poly() const = 0;
  virtual T linear_root(const Poly<T> &p) const = 0;
};

template <class T, class BaseType = Ring<Poly<T> > >
class PolyRing : public BaseType, public PolyRingOps<T> {
protected:
  friend class Poly<T>;
  const Ring<T> *m_ring;

public:
  int poly_rand_deg = -1;
  typedef T XType;
  typedef Poly<T> PT;
  typedef PolyRing<T, BaseType> CurType;
  virtual const Ring<T> *underlying_ring() const { return m_ring; }

  virtual const void *get_underlying_ring() const  { return m_ring; }
  virtual const void *get_poly_ring() const  {
    return (void *)get_poly_ring_impl();
  }
  const PolyRingOps<T> *get_poly_ring_impl() const { return this; }

  void init(const Ring<T> *ring, bignum size = -1) {
    m_ring = ring;
    BaseType::init(size, m_ring->getCar());
    OPA_CHECK0(m_ring != nullptr);
  }

  PolyRing() {}
  PolyRing(const Ring<T> *ring) { init(ring); }

  virtual ~PolyRing<T, BaseType>() {}
  virtual Poly<T> get_poly() const  { return Poly<T>(this); }

  virtual bool isInv(const Poly<T> &a) const  {
    return a.deg() == 0 && underlying_ring()->isInv(a[0]);
  }
  virtual Poly<T> inv(const Poly<T> &a) const  {
    return constant(underlying_ring()->inv(a[0]));
  }

  virtual bool compareRank(const Poly<T> &a, const Poly<T> &b) const  {
    return a.size() < b.size();
  }

  virtual Poly<T> mulmoddeg(const Poly<T> &a, const Poly<T> &b, int deg) const {
    Poly<T> res = get_poly();
    if (isZ(a) || isZ(b)) return res;

    res.poly.resize(std::max(0, a.size() + b.size() - 1), m_ring->getZ());
    for (int i = 0; i < a.size(); ++i)
      for (int j = 0; j < b.size(); ++j) {
        if (i + j >= deg) continue;
        res.poly[i + j] =
          m_ring->add(res.poly[i + j], m_ring->mul(a.poly[i], b.poly[j]));
      }
    return normalize(res);
  }

  virtual Poly<T> mul(const Poly<T> &a, const Poly<T> &b) const  {
    Poly<T> res = get_poly();
    if (isZ(a) || isZ(b)) return res;

    res.poly.resize(std::max(0, a.size() + b.size() - 1), m_ring->getZ());
    for (int i = 0; i < a.size(); ++i)
      for (int j = 0; j < b.size(); ++j)
        res.poly[i + j] =
          m_ring->add(res.poly[i + j], m_ring->mul(a.poly[i], b.poly[j]));
    return normalize(res);
  }

  Poly<T> mul_fft(const PT &a, const PT &b) const {
    int n = a.deg() + b.deg() + 1;
    // OPA_CHECK0(m_ring->
    OPA_CHECK0(false);
    return Poly<T>();
  }

  virtual Poly<T> add(const Poly<T> &a, const Poly<T> &b) const {
    if (a.size() < b.size()) return add(b, a);
    Poly<T> res(a);
    for (int i = 0; i < b.size(); ++i)
      res.get(i) = m_ring->add(res.get(i), b.get(i));
    return normalize(res);
  }

  virtual Poly<T> sub(const Poly<T> &a, const Poly<T> &b) const {
    return add(a, this->neg(b));
  }

  virtual Poly<T> &sadd(Poly<T> &a, const Poly<T> &b) const {
    a.poly.resize(std::max(a.size(), b.size()), m_ring->getZ());
    REP (i, b.size())
      a.get(i) = m_ring->add(a.get(i), b.get(i));
    normalize(a);
    return a;
  }

  virtual Poly<T> add(const Poly<T> &a, const Poly<T> &b, int shift) const  {
    Poly<T> res(a);
    if (b.deg() == -1) return res;
    if (b.deg() + shift > a.deg())
      res.poly.resize(b.deg() + shift + 1, m_ring->getZ());
    for (int i = 0; i < b.size(); ++i)
      res.get(i + shift) = m_ring->add(res.get(i + shift), b.get(i));
    return normalize(res);
  }

  virtual Poly<T> neg(const Poly<T> &a) const  {
    Poly<T> res(a);
    for (int i = 0; i < res.size(); ++i) res.get(i) = m_ring->neg(res.get(i));
    return res;
  }

  virtual bool isZ(const Poly<T> &a) const  { return a.deg() == -1; }

  virtual bool isE(const Poly<T> &a) const  {
    return a.deg() == 0 && m_ring->isE(a.get(0));
  }

  virtual Poly<T> getZ() const  { return get_poly(); }

  virtual Poly<T> getE() const  {
    Poly<T> res = get_poly();
    res.poly.push_back(m_ring->getE());
    return res;
  }

  virtual Poly<T> getRand() const  {
    OPA_CHECK0(poly_rand_deg != -1);
    return randDim(poly_rand_deg + 1);
  }

  virtual int get_poly_pos() const { return 1 + m_ring->get_poly_pos(); }


  virtual Poly<T> importu32(u32 v) const  {
    return constant(m_ring->importu32(v));
  }

  virtual bignum export_base(const Poly<T> &x) const {
    OPA_BG res = 0;
    REPV(i, x.size()){
      res = res * m_ring->getSize() + m_ring->export_base(x.get(i));
    }
    return res;

  }
  virtual Poly<T> import_base(OPA_BG num) const  {
    Poly<T> res = get_poly();
    OPA_BG per_elem = m_ring->getSize();
    while (num > 0) {
      OPA_BG cur = num % per_elem;
      res.poly.push_back(m_ring->importu32(cur.getu32()));
      num /= per_elem;
    }
    return normalize(res);
  }

  virtual Poly<T> import(const Poly<T> &p) const  {
    Poly<T> res = get_poly();
    for (int i = 0; i < p.size(); ++i)
      res.poly.push_back(m_ring->import(p.poly[i]));
    return normalize(res);
  }

  virtual std::vector<Poly<T>> import_vecs(const std::vector<Poly<T>> &plist) const  {
    std::vector<Poly<T>> res;
    for (auto &p : plist) res.push_back(this->import(p));
    return res;
  }


  Poly<T> import_vec(const std::vector<T> &px, bool rev = !kPolyRev) const {
    return import(px, rev);
  }

  virtual Poly<T> import(const std::vector<T> &px,
                         bool rev = !kPolyRev) const  {
    Poly<T> res = get_poly();
    res.poly = px;
    if (rev) std::reverse(ALL(res.poly));
    return normalize(res);
  }

  std::vector<Poly<T> > import_consts(const std::vector<T> &consts) const {
    std::vector<Poly<T> > res;
    for (auto &e : consts) {
      res.emplace_back(this->constant(e));
    }
    return res;
  }

  Poly<T> &normalize(Poly<T> &p) const {
    OPA_CHECK0(m_ring != nullptr);
    while (p.size() && m_ring->isZ(p.poly.back())) p.poly.pop_back();
    return p;
  }

  Poly<T> &selfrev(Poly<T> &p) const {
    for (int i = 0; i < p.size() / 2; ++i)
      std::swap(p.poly[i], p.poly[p.size() - 1 - i]);
    return normalize(p);
  }

  Poly<T> rev(const Poly<T> &p) const {
    Poly<T> res = p;
    this->selfrev(res);
    return res;
  }

  Poly<T> &set1(Poly<T> &p, int deg, const T &get) const {
    if (p.size() <= deg) p.poly.resize(deg + 1, m_ring->getZ());
    p.poly[deg] = get;
    normalize(p);
    return p;
  }

  Poly<T> &sresize(Poly<T> &p, int size) const {
    p.poly.resize(size, m_ring->getZ());
    normalize(p);
    return p;
  }

  Poly<T> resize(const Poly<T> &p, int size) const {
    Poly<T> res = p;
    this->sresize(res, size);
    return res;
  }

  Poly<T> &mul1(Poly<T> &p, int deg, const T &u) const {
    if (p.size() > deg) p.poly[deg] = m_ring->mul(p.poly[deg], u);
    normalize(p);
    return p;
  }

  Poly<T> &add1(Poly<T> &p, int deg, const T &u) const {
    if (p.size() <= deg) p.poly.resize(deg + 1, m_ring->getZ());
    p.poly[deg] = m_ring->add(p.poly[deg], u);
    normalize(p);
    return p;
  }

  virtual std::vector<T> toVector(const Poly<T> &p, int sz) const  {
    std::vector<T> res = p.toVector();
    res.resize(sz, m_ring->getZ());
    return res;
  }

  Poly<T> mulmod(const Poly<T> &a, const Poly<T> &b,
                 const Poly<T> &pmod) const {
    return this->mod(mul(a, b), pmod);
  }

  Poly<T> addmod(const Poly<T> &a, const Poly<T> &b,
                 const Poly<T> &pmod) const {
    return this->mod(add(a, b), pmod);
  }

  std::vector<std::pair<Poly<T>, int> > factor2(Poly<T> a,
                                                int deg_lim = -1) const {
    std::vector<Poly<T> > factors = factor(a, deg_lim);
    std::map<Poly<T>, int> tmp;
    for (auto &e : factors) {
      ++tmp[e];
    }
    std::vector<std::pair<Poly<T>, int> > res;
    for (auto &e : tmp) {
      res.push_back(e);
    }
    return res;
  }

  bool isIrred(const Poly<T> &a) const {
    std::vector<Poly<T> > res = factor(a);
    return res.size() == 1;
  }

  Poly<T> rand(int deg) const {
    assert(deg >= 0);

    std::vector<T> tb(deg + 1);
    tb[deg] = m_ring->getRandNZ();
    for (int i = 0; i < deg; ++i) tb[i] = m_ring->getRand();
    return import(tb);
  }

  Poly<T> randDim(int dim) const {
    assert(dim > 0);

    std::vector<T> tb(dim);
    for (int i = 0; i < dim; ++i) tb[i] = m_ring->getRand();
    return import(tb);
  }

  Poly<T> x() const {
    Poly<T> res = get_poly();
    res.poly.push_back(m_ring->getZ());
    res.poly.push_back(m_ring->getE());
    return res;
  }

  virtual Poly<T> xpwv(int d, const T &v) const  {
    Poly<T> res = xpw(d);
    res.poly[d] = v;
    normalize(res);
    return res;
  }

  virtual Poly<T> xpw(int d) const  {
    Poly<T> res = get_poly();
    if (d < 0) return res;
    res.poly.resize(d + 1, m_ring->getZ());
    res.poly[d] = m_ring->getE();
    return res;
  }

  Poly<T> E() const {
    Poly<T> res = get_poly();
    res.poly.push_back(m_ring->getE());
    return res;
  }

  Poly<T> constant(const T &a) const {
    Poly<T> res = get_poly();
    res.poly.push_back(a);
    normalize(res);
    return res;
  }

  virtual T eval(const Poly<T> &p, const T &a) const  {
    T res = m_ring->getZ();
    T x = m_ring->getE();
    for (int i = 0; i < p.size(); ++i) {
      res = m_ring->add(res, m_ring->mul(x, p.get(i)));
      x = m_ring->mul(x, a);
    }
    return res;
  }

  bool is_monic(const Poly<T> &p) const {
    if (p.deg() == -1) return false;
    return m_ring->isE(p[p.deg()]);
  }

  virtual Poly<T> addc(const Poly<T> &p, const T &a) const  {
    Poly<T> res = p;
    return saddc(res, a);
  }

  Poly<T> &saddc(Poly<T> &p, const T &a) const {
    if (isZ(p)) {
      p = this->constant(a);
    } else {
      p[0] = m_ring->add(p[0], a);
      normalize(p);
    }
    return p;
  }

  Poly<T> &ssubc(Poly<T> &p, const T &a) const {
    return saddc(p, m_ring->neg(a));
  }

  virtual Poly<T> subc(const Poly<T> &p, const T &a) const  {
    Poly<T> res = p;
    return ssubc(res, a);
  }

  virtual Poly<T> mulc(const Poly<T> &p, const T &a) const  {
    Poly<T> res = p;
    return smulc(res, a);
  }

  Poly<T> &smulc(Poly<T> &p, const T &a) const {
    REP (i, p.deg() + 1)
      p.get(i) = m_ring->mul(a, p.get(i));
    return normalize(p);
  }

  Poly<T> &sdivc(Poly<T> &p, const T &a) const {
    REP (i, p.deg() + 1)
      p.get(i) = m_ring->div(p.get(i), a);
    return p;
  }

  virtual Poly<T> divc(const Poly<T> &p, const T &a) const  {
    Poly<T> res = p;
    return sdivc(res, a);
  }

  Poly<T> monic(const Poly<T> &p) const {
    Poly<T> res = p;
    return smonic(res);
  }

  Poly<T> &smonic(Poly<T> &p) const {
    OPA_CHECK0(p.deg() >= 0);
    OPA_CHECK(m_ring->isInv(p.lc()), p);
    return smulc(p, m_ring->inv(p.lc()));
  }

  virtual T linear_root(const Poly<T> &p) const {
    OPA_CHECK0(p.deg() == 1);
    Poly<T> t2 = this->monic(p);
    return m_ring->neg(t2[0]);
  }

  virtual Poly<T> derivative(const Poly<T> &p) const  {
    Poly<T> res = xpw(p.deg() - 1);
    REP (i, p.deg()) {
      res[i] = m_ring->mul(m_ring->importu32(i + 1), p[i + 1]);
    }
    return normalize(res);
  }

  bool ediv(const Poly<T> &a, const Poly<T> &b, Poly<T> *q, Poly<T> *r) const {
    if (b.deg() < 0) return false;
    Poly<T> qq = get_poly();
    Poly<T> rr = a;

    if (0 && get_poly_pos() == 2) {
      puts("\nDO EDIV >>");
      std::cout << a << std::endl;
      std::cout << b << std::endl;
    }

    if (a.deg() >= b.deg()) {
      qq.poly.resize(a.deg() - b.deg() + 1, m_ring->getZ());

      // todo: shit here
      // T ibc = m_ring->getZ();
      // if (m_ring->isInv(b.lc())) {
      //  ibc = m_ring->inv(b.lc());
      //}
      // T ibc = m_ring->inv(b.get(b.deg()));

      while (rr.deg() >= b.deg()) {
        T e;
        if (!m_ring->ediv(rr.lc(), b.lc(), &e, nullptr)) return false;
        int shift = rr.deg() - b.deg();

        int old = rr.deg();
        rr = add(rr, mulc(b, m_ring->neg(e)), shift);
        OPA_CHECK(old > rr.deg(), old, rr, a, b);
        qq.poly[shift] = e;
      }
    }

    if (q) *q = qq;
    if (r) *r = rr;
    return true;
  }

  void equalDegFactor(Poly<T> a, int deg, std::vector<Poly<T> > &res) const {
    bignum q = m_ring->getSize();

    if (a.deg() == deg) {
      res.push_back(a);
      return;
    }

    bignum pw;
    if (q.get_bit(0)) {
      pw = 1;
      for (int i = 0; i < deg; ++i) pw = pw * q;
      pw = (pw - 1) / 2;
    } else {
      pw = 0;
      for (bignum i = 1; i <= q; i.slshift(1), pw += 1)
        ;
      pw *= deg;
      pw -= 1;
    }

    // :(
    bool special_case = deg == 1 && a.deg() == 2;

    while (1) {

      Poly<T> spe_root;
      Poly<T> a_special = a;

      if (special_case) {
        // create artificial third factor
        spe_root = import({ m_ring->getRandRaw(), m_ring->getE() });
        a_special = a * spe_root;
      }

      Poly<T> b = this->rand(a_special.deg() - 1);

      if (m_ring->getSize().get_bit(0)) {
        b = this->faste(b, pw, a);
        b = add(b, getE());
      } else {
        Poly<T> tmp = b;
        REP (i, pw.getu32() - 1) {
          tmp = this->mod(mul(tmp, tmp), a_special);
          b = add(b, tmp);
        }
      }

      b = this->gcd(a, b);
      if (b.deg() > 0 && b.deg() < a.deg()) {
        Poly<T> q = this->div(a, b);
        equalDegFactor(q, deg, res);
        equalDegFactor(b, deg, res);
        return;
      }
    }
  }

  std::vector<Poly<T> > factor(Poly<T> a, int deg_lim = -1) const {
    std::vector<Poly<T> > res;
    Poly<T> rem = x();
    if (deg_lim == -1) deg_lim = a.deg();
    deg_lim = std::min(deg_lim, a.deg() / 2);

    for (int i = 1; i <= deg_lim; ++i) {
      rem = this->faste(rem, m_ring->getSize(), a);

      Poly<T> tmp = add(rem, neg(x()));

      while (1) {
        tmp = this->gcd(a, tmp);

        if (tmp.deg() <= 0) break;

        equalDegFactor(tmp, i, res);
        a = this->div(a, tmp);
      }
    }

    if (a.deg() > 0) res.push_back(a);
    return res;
  }

  virtual PT gcd(const PT &a, const PT &b) const {
    if (compareRank(a, b)) return gcd(b, a);
    if (isZ(b) || m_ring->isInv(b.lc())) {
      if (b.deg() == 0) return getE();
      return BaseType::gcd(a, b);
    }
    PT a2, b2;
    T ca, cb;
    a.pp_and_cont(&a2, &ca);
    b.pp_and_cont(&b2, &cb);
    T cd = m_ring->gcd(ca, cb);
    PT d = _gcd(a2, b2);
    // OPA_DISP0(cd, d);
    return d * cd;
  }

  PT _gcd(const PT &a, const PT &b) const {
    if (isZ(b)) return a;
    if (b.deg() == 0) return getE();
    T mc = m_ring->faste(b.lc(), a.deg() - b.deg() + 1);
    PT na = a * mc;

    PT r, q;
    bool res = ediv(na, b, &q, &r);
    r = r.pp();
    OPA_CHECK(res, a, b, na);
    q = neg(q);

    PT d = _gcd(b, r);
    // OPA_DISP0(a, b, r, d);
    return d;
  }

  virtual PT egcd(const PT &a, const PT &b, PT &u, PT &v) const {
    if (isZ(b) || m_ring->isInv(b.lc())) return BaseType::egcd(a, b, u, v);
    OPA_CHECK0(false); // ufd shit
    PT u2, v2;
    bool ignore = false;
    return egcd2(a, b, u, v, u2, v2, ignore);
  }

  struct GCDContext {
    PT u, v, u2, v2;
    T coeff_mul;
    T gcd_mul;
    T final_mul;
    PT res;
    bool gcd_on_a;
    const CurType *ct;

    void go(const PT &a, const PT &b, const PT &ua, const PT &va, const PT &ub,
            const PT &vb) {
      // OPA_DISP0(a, b, ct->isZ(b));

      if (b.deg() == 0) {
        res = ct->getE();
        final_mul = b[0];
        u = ub;
        v = vb;
        u2 = ua;
        v2 = va;
        return;
      }

      if (ct->isZ(b)) {
        final_mul = ct->m_ring->getE();
        u = ua;
        v = va;
        res = a;
        return;
      }
      T mc = ct->m_ring->faste(b.lc(), a.deg() - b.deg() + 1);
      PT na = a * mc;

      PT r, q;
      coeff_mul = ct->m_ring->mul(coeff_mul, mc);

      bool res = ct->ediv(na, b, &q, &r);
      OPA_CHECK(res, a, b, na);
      OPA_CHECK_EQ(b * q + r, na, b, q, r, na);

      gcd_on_a ^= 1;
      go(b, r, ub, vb, ua * mc - ub * q, va * mc - vb * q);
      if (this->res == b) {
        u2 = ua;
        v2 = va;
      }
    }
  };

  virtual GCDContext egcd2_ctx(const PT &a, const PT &b) const {
    // TODO: a version with coefficient reduction directly integrated into
    // determinant computation, or use a fraction field
    GCDContext ctx;
    ctx.coeff_mul = m_ring->getE();
    ctx.ct = this;

    ctx.go(a, b, getE(), getZ(), getZ(), getE());
    return ctx;
  }

  virtual PT egcd2(const PT &a, const PT &b, PT &u, PT &v, PT &u2, PT &v2,
                   bool &gcd_on_a) const {
    if (isZ(b) || m_ring->isInv(b.lc()))
      return BaseType::egcd2(a, b, u, v, u2, v2, gcd_on_a);
    if (compareRank(a, b)) return egcd2(b, a, v, u, v2, u2, gcd_on_a ^= 1);
    GCDContext ctx = egcd2_ctx(a, b);
    u = ctx.u;
    v = ctx.v;
    u2 = ctx.u2;
    v2 = ctx.v2;
    gcd_on_a ^= ctx.gcd_on_a;
    return ctx.res;
  }
};

template <class T> using PolyRingBaseField = PolyRing<T, Field<Poly<T> > >;

OPA_NM_MATH_COMMON_END

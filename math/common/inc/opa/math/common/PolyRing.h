#pragma once

#include <opa/math/common/GF_p.h>
#include <opa/math/common/Ring.h>
#include <opa/math/common/Utils.h>
#include <opa/math/common/base.h>
#include <opa/math/common/bignum.h>
#include <unordered_map>

OPA_NM_MATH_COMMON
constexpr bool kPolyRev = true;

template <class T, class BaseRing = Ring<T> > class Poly;

template <class T, class BaseRing = Ring<T> > class PolyRingOps {
public:
  typedef Poly<T, BaseRing> PT;
  virtual ~PolyRingOps() {}

  virtual PT addc(const PT &p, const T &a) const = 0;
  virtual PT &normalize(PT &a) const = 0;

  virtual PT &set1(PT &p, int deg, const T &get) const = 0;
  virtual PT subc(const PT &p, const T &a) const = 0;
  virtual PT pointwise_mul(const PT &p, const PT &a) const = 0;
  virtual PT mulc(const PT &p, const T &a) const = 0;
  virtual PT divc(const PT &p, const T &a) const = 0;
  virtual std::vector<T> toVector(const PT &p, int sz) const = 0;
  virtual PT derivate(const PT &p) const = 0;
  virtual PT integrate(const PT &p) const = 0;
  virtual PT monic(const PT &p) const = 0;
  virtual T eval(const PT &p, const T &a) const = 0;
  virtual PT eval2(const PT &p, const PT &x) const = 0;
  virtual PT import_base(OPA_BG num) const = 0;
  virtual PT import(const PT &p) const = 0;
  virtual PT import(const std::vector<T> &px, bool rev = !kPolyRev) const = 0;
  virtual PT xma(const T &v) const = 0;
  virtual PT xpa(const T &v) const = 0;
  virtual PT divxpw(const PT &p, u32 xpw) const = 0;
  virtual PT mulxpw(const PT &p, u32 xpw) const = 0;
  virtual PT mulx(const PT &p, const T &v) const = 0;
  virtual int get_xpw(const PT &p) const = 0;

  virtual PT resize(const PT &p, int size) const = 0;
  virtual PT rev(const PT &p, int target_size = 0) const = 0;
  virtual PT xpwv(int d, const T &v) const = 0;
  virtual PT xpw(int d) const = 0;
  virtual PT get_poly() const = 0;
  virtual PT interpolate(const std::vector<T> &xl, const std::vector<T> &yl) const = 0;
  virtual T linear_root(const PT &p) const = 0;
};

template <class T, class BaseRing = Ring<T>, class BaseType = Ring<Poly<T, BaseRing> > >
class PolyRing : public BaseType, public PolyRingOps<T, BaseRing> {
protected:
  friend class Poly<T, BaseRing>;
  const BaseRing *m_ring;

public:
  int poly_rand_deg = -1;
  typedef T XType;
  typedef Poly<T, BaseRing> PT;
  typedef PolyRing<T, BaseRing, BaseType> CurType;
  virtual const BaseRing *underlying_ring() const { return m_ring; }

  virtual const void *get_underlying_ring() const { return m_ring; }
  virtual const void *get_poly_ring() const { return (void *)get_poly_ring_impl(); }
  const PolyRingOps<T, BaseRing> *get_poly_ring_impl() const { return this; }

  void init(const BaseRing *ring, bignum size = -1) {
    m_ring = ring;
    BaseType::init(size, m_ring->getCar());
    OPA_CHECK0(m_ring != nullptr);
  }

  PolyRing() {}
  PolyRing(const BaseRing *ring) { init(ring); }

  virtual ~PolyRing() {}
  virtual PT get_poly() const { return PT(this); }

  virtual int get_xpw(const PT &p) const {
    REP (i, p.size())
      if (!m_ring->isZ(p[i])) return i;
    return -1;
  }

  virtual bool isInv(const PT &a) const { return a.deg() == 0 && underlying_ring()->isInv(a[0]); }
  virtual PT inv(const PT &a) const { return constant(underlying_ring()->inv(a[0])); }

  virtual bool compareRank(const PT &a, const PT &b) const { return a.size() < b.size(); }

  virtual PT mulmoddeg(const PT &a, const PT &b, int deg) const {
    PT res = get_poly();
    if (isZ(a) || isZ(b)) return res;

    res.poly.resize(std::max(0, a.size() + b.size() - 1), m_ring->getZ());
    for (int i = 0; i < a.size(); ++i)
      for (int j = 0; j < b.size(); ++j) {
        if (i + j >= deg) continue;
        res.poly[i + j] = m_ring->add(res.poly[i + j], m_ring->mul(a.poly[i], b.poly[j]));
      }
    return normalize(res);
  }
  virtual PT pointwise_mul(const PT &p, const PT &a) const {
    PT res = p;
    REP (i, std::min(p.size(), a.size())) res[i] = m_ring->mul(res[i], a[i]);
    return res;
  }

  virtual PT mul(const PT &a, const PT &b) const {
    PT res = get_poly();
    if (isZ(a) || isZ(b)) return res;

    res.poly.resize(std::max(0, a.size() + b.size() - 1), m_ring->getZ());
    for (int i = 0; i < a.size(); ++i)
      for (int j = 0; j < b.size(); ++j)
        res.poly[i + j] = m_ring->add(res.poly[i + j], m_ring->mul(a.poly[i], b.poly[j]));
    return normalize(res);
  }

  PT mul_fft(const PT &a, const PT &b) const {
    int n = a.deg() + b.deg() + 1;
    // OPA_CHECK0(m_ring->
    OPA_CHECK0(false);
    return PT();
  }

  virtual PT add(const PT &a, const PT &b) const {
    if (a.size() < b.size()) return add(b, a);
    PT res(a);
    for (int i = 0; i < b.size(); ++i) res.get(i) = m_ring->add(res.get(i), b.get(i));
    return normalize(res);
  }

  virtual PT sub(const PT &a, const PT &b) const { return add(a, this->neg(b)); }

  virtual PT &sadd(PT &a, const PT &b) const {
    a.poly.resize(std::max(a.size(), b.size()), m_ring->getZ());
    REP (i, b.size()) a.get(i) = m_ring->add(a.get(i), b.get(i));
    normalize(a);
    return a;
  }

  virtual PT add(const PT &a, const PT &b, int shift) const {
    PT res(a);
    if (b.deg() == -1) return res;
    if (b.deg() + shift > a.deg()) res.poly.resize(b.deg() + shift + 1, m_ring->getZ());
    for (int i = 0; i < b.size(); ++i)
      res.get(i + shift) = m_ring->add(res.get(i + shift), b.get(i));
    return normalize(res);
  }

  virtual PT neg(const PT &a) const {
    PT res(a);
    for (int i = 0; i < res.size(); ++i) res.get(i) = m_ring->neg(res.get(i));
    return res;
  }

  virtual bool isZ(const PT &a) const { return a.deg() == -1; }

  virtual bool isE(const PT &a) const { return a.deg() == 0 && m_ring->isE(a.get(0)); }

  virtual PT getZ() const { return get_poly(); }

  virtual PT getE() const {
    PT res = get_poly();
    res.poly.push_back(m_ring->getE());
    return res;
  }

  virtual PT getRand() const {
    OPA_CHECK0(poly_rand_deg != -1);
    return randDim(poly_rand_deg + 1);
  }

  virtual int get_poly_pos() const { return 1 + m_ring->get_poly_pos(); }

  virtual PT importu32(u32 v) const { return constant(m_ring->importu32(v)); }

  virtual bignum export_base(const PT &x) const {
    OPA_BG res = 0;
    REPV (i, x.size()) {
      res = res * m_ring->getSize() + m_ring->export_base(x.get(i));
    }
    return res;
  }
  virtual PT import_base(OPA_BG num) const {
    PT res = get_poly();
    OPA_BG per_elem = m_ring->getSize();
    while (num > 0) {
      OPA_BG cur = num % per_elem;
      res.poly.push_back(m_ring->importu32(cur.getu32()));
      num /= per_elem;
    }
    return normalize(res);
  }

  virtual PT import(const PT &p) const {
    PT res = get_poly();
    for (int i = 0; i < p.size(); ++i) res.poly.push_back(m_ring->import(p.poly[i]));
    return normalize(res);
  }

  virtual std::vector<PT> import_vecs(const std::vector<PT> &plist) const {
    std::vector<PT> res;
    for (auto &p : plist) res.push_back(this->import(p));
    return res;
  }

  PT import_vec(const std::vector<T> &px, bool rev = !kPolyRev) const { return import(px, rev); }

  virtual PT import(const std::vector<T> &px, bool rev = !kPolyRev) const {
    PT res = get_poly();
    res.poly = px | STD_TSFX(m_ring->import(x)) | STD_VEC;
    if (rev) std::reverse(ALL(res.poly));
    return normalize(res);
  }

  std::vector<PT> import_consts(const std::vector<T> &consts) const {
    std::vector<PT> res;
    for (auto &e : consts) {
      res.emplace_back(this->constant(e));
    }
    return res;
  }

  PT &normalize(PT &p) const {
    OPA_CHECK0(m_ring != nullptr);
    while (p.size() && m_ring->isZ(p.poly.back())) p.poly.pop_back();
    return p;
  }

  PT &selfrev(PT &p, int target_size = -1) const {
    if (target_size != -1) p.poly.resize(target_size, m_ring->getZ());
    for (int i = 0; i < p.size() / 2; ++i) std::swap(p.poly[i], p.poly[p.size() - 1 - i]);
    return normalize(p);
  }

  PT rev(const PT &p, int target_size = -1) const override {
    PT res = p;
    this->selfrev(res, target_size);
    return res;
  }

  PT &set1(PT &p, int deg, const T &get) const override {
    if (p.size() <= deg) p.poly.resize(deg + 1, m_ring->getZ());
    p.poly[deg] = get;
    normalize(p);
    return p;
  }

  PT &sresize(PT &p, int size) const {
    p.poly.resize(size, m_ring->getZ());
    normalize(p);
    return p;
  }

  PT resize(const PT &p, int size) const override {
    PT res = p;
    this->sresize(res, size);
    return res;
  }

  PT xpa(const T &v) const override { return this->xma(m_ring->neg(v)); }

  PT xma(const T &v) const override {
    auto res = this->x();
    res[0] = this->m_ring->neg(v);
    return res;
  }

  PT interpolate(const std::vector<T> &xl, const std::vector<T> &yl) const override {
    std::unordered_set<T> sx(ALL(xl));
    OPA_CHECK0(sx.size() == xl.size()); // no dups
    OPA_CHECK0(xl.size() == yl.size());
    PT mulp = this->getE();
    REP (i, xl.size()) mulp = this->mul(mulp, this->xma(xl[i]));

    PT res = this->getZ();
    REP (i, xl.size()) {
      auto h = mulp / this->xma(xl[i]);
      res = res + h * (m_ring->div(yl[i], h(xl[i])));
    }
    return res;
  }

  PT &mul1(PT &p, int deg, const T &u) const {
    if (p.size() > deg) p.poly[deg] = m_ring->mul(p.poly[deg], u);
    normalize(p);
    return p;
  }

  PT &add1(PT &p, int deg, const T &u) const {
    if (p.size() <= deg) p.poly.resize(deg + 1, m_ring->getZ());
    p.poly[deg] = m_ring->add(p.poly[deg], u);
    normalize(p);
    return p;
  }

  virtual std::vector<T> toVector(const PT &p, int sz) const {
    std::vector<T> res = p.toVector();
    res.resize(sz, m_ring->getZ());
    return res;
  }

  PT mulmod(const PT &a, const PT &b, const PT &pmod) const { return this->mod(mul(a, b), pmod); }

  PT addmod(const PT &a, const PT &b, const PT &pmod) const { return this->mod(add(a, b), pmod); }

  std::vector<std::pair<PT, int> > factor2(PT a, int deg_lim = -1) const {
    std::vector<PT> factors = factor(a, deg_lim);
    std::map<PT, int> tmp;
    for (auto &e : factors) {
      ++tmp[e];
    }
    std::vector<std::pair<PT, int> > res;
    for (auto &e : tmp) {
      res.push_back(e);
    }
    return res;
  }

  bool isIrred(const PT &a) const {
    std::vector<PT> res = factor(a);
    return res.size() == 1;
  }

  PT rand(int deg) const {
    assert(deg >= 0);

    std::vector<T> tb(deg + 1);
    tb[deg] = m_ring->getRandNZ();
    for (int i = 0; i < deg; ++i) tb[i] = m_ring->getRand();
    return import(tb);
  }

  PT randDim(int dim) const {
    assert(dim > 0);

    std::vector<T> tb(dim);
    for (int i = 0; i < dim; ++i) tb[i] = m_ring->getRand();
    return import(tb);
  }

  PT x() const {
    PT res = get_poly();
    res.poly.push_back(m_ring->getZ());
    res.poly.push_back(m_ring->getE());
    return res;
  }

  virtual PT xpwv(int d, const T &v) const {
    PT res = xpw(d);
    res.poly[d] = v;
    normalize(res);
    return res;
  }

  virtual PT xpw(int d) const {
    PT res = get_poly();
    if (d < 0) return res;
    res.poly.resize(d + 1, m_ring->getZ());
    res.poly[d] = m_ring->getE();
    return res;
  }

  PT E() const {
    PT res = get_poly();
    res.poly.push_back(m_ring->getE());
    return res;
  }

  PT constant(const T &a) const {
    PT res = get_poly();
    res.poly.push_back(a);
    normalize(res);
    return res;
  }

  std::vector<T> eval_batch(const PT &p, const std::vector<T> &tb) const {
    std::vector<T> res(tb.size());
    REP (i, tb.size()) res[i] = this->eval(p, tb[i]);
    return res;
  }

  virtual PT eval2(const PT &p, const PT &x) const override {
    PT res = this->getZ();
    PT cur = this->getE();

    for (int i = 0; i < p.size(); ++i) {
      res = res + this->mulc(cur, p[i]);
      cur = cur * x;
    }
    return res;
  }

  virtual T eval(const PT &p, const T &a) const override {
    T res = m_ring->getZ();
    T x = m_ring->getE();
    for (int i = 0; i < p.size(); ++i) {
      res = m_ring->add(res, m_ring->mul(x, p.get(i)));
      x = m_ring->mul(x, a);
    }
    return res;
  }

  bool is_monic(const PT &p) const {
    if (p.deg() == -1) return false;
    return m_ring->isE(p[p.deg()]);
  }

  virtual PT addc(const PT &p, const T &a) const {
    PT res = p;
    return saddc(res, a);
  }

  PT &saddc(PT &p, const T &a) const {
    if (isZ(p)) {
      p = this->constant(a);
    } else {
      p[0] = m_ring->add(p[0], a);
      normalize(p);
    }
    return p;
  }

  PT &ssubc(PT &p, const T &a) const { return saddc(p, m_ring->neg(a)); }

  virtual PT subc(const PT &p, const T &a) const {
    PT res = p;
    return ssubc(res, a);
  }

  virtual PT mulc(const PT &p, const T &a) const {
    PT res = p;
    return smulc(res, a);
  }

  PT &smulc(PT &p, const T &a) const {
    REP (i, p.deg() + 1) p.get(i) = m_ring->mul(a, p.get(i));
    return normalize(p);
  }

  PT &sdivc(PT &p, const T &a) const {
    REP (i, p.deg() + 1) p.get(i) = m_ring->div(p.get(i), a);
    return p;
  }

  virtual PT divc(const PT &p, const T &a) const {
    PT res = p;
    return sdivc(res, a);
  }

  virtual PT mulxpw(const PT &p, u32 xpw) const override {
    std::vector<T> res(p.poly.size() + xpw, this->underlying_ring()->getZ());
    REP (i, p.poly.size()) res[i + xpw] = p.poly[i];
    return this->import_vec(res);
  }

  virtual PT mulx(const PT &p, const T &v) const {
    // x => x*v
    std::vector<T> res(p.poly.size());
    T x = this->underlying_ring()->getE();
    REP (i, p.poly.size()) {
      res[i] = this->underlying_ring()->mul(p.poly[i], x);
      x = this->underlying_ring()->mul(x, v);
    }
    return this->import_vec(res);
  }

  virtual PT divxpw(const PT &p, u32 xpw) const override {
    std::vector<T> res;
    FOR (i, xpw, p.poly.size()) res.pb(p[i]);
    return this->import_vec(res);
  }

  PT combine_lh(const PT xl, const PT &xh, int l) const {
    std::vector<T> res = xl.poly;
    res.resize(l + xh.poly.size(), m_ring->getZ());
    REP (i, xh.poly.size()) res[i + l] = xh[i];
    return import_vec(res);
  }

  PT monic(const PT &p) const {
    PT res = p;
    return smonic(res);
  }

  PT &smonic(PT &p) const {
    OPA_CHECK0(p.deg() >= 0);
    OPA_CHECK(m_ring->isInv(p.lc()), p);
    return smulc(p, m_ring->inv(p.lc()));
  }

  virtual T linear_root(const PT &p) const {
    OPA_CHECK0(p.deg() == 1);
    PT t2 = this->monic(p);
    return m_ring->neg(t2[0]);
  }

  virtual PT integrate(const PT &p) const {
    if (p.deg() < 0) return this->getZ();

    PT res = xpw(p.deg() + 1);
    REP (i, p.deg()) {
      res[i + 1] = m_ring->div(p[i], m_ring->importu32(i + 1));
    }
    return normalize(res);
  }

  virtual PT derivate(const PT &p) const {
    PT res = xpw(p.deg() - 1);
    REP (i, p.deg()) {
      res[i] = m_ring->mul(m_ring->importu32(i + 1), p[i + 1]);
    }
    return normalize(res);
  }

  bool ediv(const PT &a, const PT &b, PT *q, PT *r) const {
    if (b.deg() < 0) return false;
    PT qq = get_poly();
    PT rr = a;

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

  void equalDegFactor(PT a, int deg, std::vector<PT> &res) const {
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

      PT spe_root;
      PT a_special = a;

      if (special_case) {
        // create artificial third factor
        spe_root = import({ m_ring->getRandRaw(), m_ring->getE() });
        a_special = a * spe_root;
      }

      PT b = this->rand(a_special.deg() - 1);

      if (m_ring->getSize().get_bit(0)) {
        b = this->faste(b, pw, a);
        b = add(b, getE());
      } else {
        PT tmp = b;
        REP (i, pw.getu32() - 1) {
          tmp = this->mod(mul(tmp, tmp), a_special);
          b = add(b, tmp);
        }
      }

      b = this->gcd(a, b);
      if (b.deg() > 0 && b.deg() < a.deg()) {
        PT q = this->div(a, b);
        equalDegFactor(q, deg, res);
        equalDegFactor(b, deg, res);
        return;
      }
    }
  }

  std::vector<PT> factor(PT a, int deg_lim = -1) const {
    std::vector<PT> res;
    PT rem = x();
    if (deg_lim == -1) deg_lim = a.deg();
    deg_lim = std::min(deg_lim, a.deg() / 2);

    for (int i = 1; i <= deg_lim; ++i) {
      rem = this->faste(rem, m_ring->getSize(), a);

      PT tmp = add(rem, neg(x()));

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

    void go(const PT &a, const PT &b, const PT &ua, const PT &va, const PT &ub, const PT &vb) {
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

  virtual PT egcd2(const PT &a, const PT &b, PT &u, PT &v, PT &u2, PT &v2, bool &gcd_on_a) const {
    if (isZ(b) || m_ring->isInv(b.lc())) return BaseType::egcd2(a, b, u, v, u2, v2, gcd_on_a);
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

template <class T> using PolyRingBaseField = PolyRing<T, Ring<T>, Field<Poly<T> > >;

OPA_NM_MATH_COMMON_END

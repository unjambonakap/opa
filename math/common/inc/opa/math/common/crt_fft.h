#pragma once

#include <memory>
#include <opa/math/common/FFT.h>
#include <opa/math/common/GF_p.h>
#include <opa/math/common/PolyModRing.h>
#include <opa/math/common/PolyRing.h>
#include <opa/math/common/Ring.h>
#include <opa/math/common/UtilsPoly.h>
#include <opa/math/common/base.h>
#include <opa_common_base.h>
#include <span>

OPA_NM_MATH_COMMON

template <class F, class T, typename LiftT> class CRT_T {
public:
  const F &m_f;
  std::vector<LiftT> coeffs;
  LiftT n;

  CRT_T(const F &f, const std::vector<const GF_pT<T> *> &pFields) : m_f(f) {
    n = f.getE();
    for (auto &x : pFields) n = n * x->n; // do operations on N

    REP (i, pFields.size()) {
      LiftT nx = n / pFields[i]->n;
      LiftT tv = pFields[i]->inv(pFields[i]->importT(nx));
      coeffs.push_back(nx * tv);
    }
  }

  T solve(const std::vector<T> &vec) const {
    OPA_CHECK0(vec.size() == coeffs.size());

    LiftT res = m_f.getZ();
    REP (i, vec.size()) {
      res += vec[i] * coeffs[i];
    }
    res %= n;
    return m_f.importT(res);
  }
};

template <class TE, int N> class CRTFieldFFTProvider;

template <class TE, int N> class CRTField : public Field<std::array<TE, N> > {

public:
  typedef typename std::conditional<std::is_same<TE, u32>::value, s128, bignum>::type TLift;
  typedef std::array<TE, N> T;
  typedef CRT_T<GF_pT<TE>, TE, TLift> CRT_Type;

private:
  CRT_Type m_crt;
  std::vector<const GF_pT<TE> *> fields;

public:
  const GF_pT<TE> &target;

  CRTField(std::vector<const GF_pT<TE> *> fields, const GF_pT<TE> &target)
      : fields(fields), target(target), m_crt(target, fields) {
    assert(fields.size() == N);
  }

  T import(const TE &v) const {
    T res;
    REP (i, N) res[i] = fields[i]->importT(v);
    return res;
  }

  std::vector<T> import_vec(const std::vector<TE> &v) const {
    std::vector<T> res(v.size());
    REP (i, v.size()) res[i] = this->import(v[i]);
    return res;
  }

  std::vector<TE> lower_vec(const std::vector<T> &v) const {
    std::vector<TE> res(v.size());
    REP (i, v.size()) res[i] = this->lower(v[i]);
    return res;
  }

  TE lower(const T &v) const {
    std::vector<TE> tb(ALL(v));
    return m_crt.solve(tb);
  }

  virtual T getNthRoot(u32 n) const override {
    T v;
    REP (i, N) v[i] = fields[i]->getNthRoot(n);
    return v;
  }

  std::unique_ptr<FFTProvider<T> > make_fft(int npw) const {
    return std::unique_ptr<FFTProvider<T> >(new FFT2Dispatcher<T>(this, npw));
  }

  std::unique_ptr<FFTProvider<TE> > make_ffte_dispatcher(int npw) const {
    auto tmp = [&](int npw_req) { return this->make_fft_provider(npw_req).release(); };
    return std::unique_ptr<FFTProvider<TE> >(new FFT2Dispatcher<TE>(npw, tmp));
  }

  std::unique_ptr<FFTProvider<TE> > make_fft_provider(int npw) const {
    auto tmp = new CRTFieldFFTProvider<TE, N>(*this, npw);
    return std::unique_ptr<FFTProvider<TE> >(tmp);
  }

#define MAKE_ZARY_OP(op_name)                                                                      \
  virtual T op_name() const override {                                                             \
    T res;                                                                                         \
    REP (i, N) res[i] = fields[i]->op_name();                                                      \
    return res;                                                                                    \
  }
#define MAKE_UNARY_OP(op_name)                                                                     \
  virtual T op_name(const T &a) const override {                                                   \
    T res;                                                                                         \
    REP (i, N) res[i] = fields[i]->op_name(a[i]);                                                  \
    return res;                                                                                    \
  }
#define MAKE_BINARY_OP(op_name)                                                                    \
  virtual T op_name(const T &a, const T &b) const override {                                       \
    T res;                                                                                         \
    REP (i, N) res[i] = fields[i]->op_name(a[i], b[i]);                                            \
    return res;                                                                                    \
  }
  MAKE_UNARY_OP(inv);
  MAKE_UNARY_OP(import);
  MAKE_UNARY_OP(neg);
  MAKE_ZARY_OP(getZ);
  MAKE_ZARY_OP(getE);
  MAKE_BINARY_OP(mul);
  MAKE_BINARY_OP(add);

  virtual T importu32(u32 v) const override {
    T res;
    REP (i, N) res[i] = fields[i]->importu32(v);
    return res;
  }

  virtual bool isZ(const T &a) const override {
    REP (i, N)
      if (!fields[i]->isZ(a[i])) return false;
    return true;
  }
  virtual bool isE(const T &a) const override {
    REP (i, N)
      if (!fields[i]->isE(a[i])) return false;
    return true;
  }
};

template <class T, int N> class CRTFieldFFTProvider : public FFTProvider<T> {
  const CRTField<T, N> &crtf;
  typedef typename CRTField<T, N>::T TLift;
  std::unique_ptr<FFTProvider<TLift> > _fft;

public:
  CRTFieldFFTProvider(const CRTField<T, N> &crtf, int npw) : crtf(crtf), _fft(crtf.make_fft(npw)) {
    this->configure(npw, &crtf.target);
  }

  std::vector<T> fft(const std::vector<T> &tb) const override {
    auto a = crtf.import_vec(tb);
    return crtf.lower_vec(_fft->fft(a));
  }

  std::vector<T> ifft(const std::vector<T> &tb, bool normalize = true) const override {
    auto a = crtf.import_vec(tb);
    return crtf.lower_vec(_fft->ifft(a, normalize));
  }

  std::vector<T> mul(const std::vector<T> &a, const std::vector<T> &b) const override {

    return crtf.lower_vec(_fft->mul(crtf.import_vec(a), crtf.import_vec(b)));
  }
};

GF_p f1{ 998244353, { { 2, 23 }, { 7, 1 }, { 17, 1 } } };
GF_p f2{ 880803841, { { 2, 23 }, { 3, 1 }, { 5, 1 }, { 7, 1 } } };
GF_p f3{ 754974721, { { 2, 23 }, { 3, 2 }, { 5, 1 } } };

CRTField<u32, 3> crtf1(const GF_p &target) {
  auto res = CRTField<u32, 3>({ &f1, &f2, &f3 }, target);
  return res;
}

template <class T> class FactProvider {

private:
  mutable std::vector<T> _fact;

public:
  typedef Ring<T> R;
  const R &ring;

  FactProvider(const R &ring) : ring(ring) { _fact.pb(ring.getE()); }

  T ifact(u32 i) const { return ring.inv(fact(i)); }
  T fact(u32 i) const {
    while (_fact.size() <= i) {
      _fact.pb(ring.mul(_fact.back(), ring.importu32(_fact.size() + 1)));
    }
    return _fact[i];
  }
};

template <class T> class FastPoly {

  const PolyRing<T> &pr;
  const FFTProvider<T> &provider;
  const FactProvider<T> &fact_provider;

public:
  typedef Poly<T> PT;
  typedef FastPoly<T> This;

  const Ring<T> &ring() const { return *pr.underlying_ring(); }

  FastPoly(const PolyRing<T> &pr, const FFTProvider<T> &provider)
      : pr(pr), provider(provider), fact_provider(*pr.underlying_ring()) {}

  PT mul(const PT &a, const PT &b) const {
    return pr.import_vec(provider.mul(a.to_vec(), b.to_vec()));
  }

  PT mulxmod(const PT &a, const PT &b, int xpw) const { return pr.resize(this->mul(a, b), xpw); }

  PT invxmod(const PT &x, int n) const { return _invxmod(x, log2_high_bit(n - 1) + 1).resize(n); }

  PT _invxmod(const PT &x, int degpw) const {
    OPA_CHECK0(ring().isE(x[0]));
    if (degpw == 0) return pr.getE();
    int ndegpw = degpw - 1;
    int l = 1 << ndegpw;
    PT xl = pr.resize(x, l);
    PT xh = pr.divxpw(x, l);

    PT yl = _invxmod(xl, ndegpw);
    PT c = pr.divxpw(mul(xl, yl), l);
    // xl * yl = 1 + c * x^l
    // (xl + x^l xh ) * (yl + y^l yh) = 1 (x^2l)
    // => c + xh*yl + xl*yh = 0 (x^)
    // => yh = -yl * (c + xh * yl)

    PT yh = this->mulxmod(-yl, c + this->mulxmod(xh, yl, l), l);
    return pr.combine_lh(yl, yh, l);
  }

  PT mod(const PT &a, const PT &b) const { return a - this->mul(this->div(a, b), b); }

  PT div(const PT &a, const PT &b) const {
    OPA_CHECK0(!pr.isZ(b));
    if (pr.isZ(a)) return a;
    if (a.deg() < b.deg()) return pr.getZ();

    PT sb = b.monic();
    PT rb = sb.rev();
    PT ra = a.rev();
    int dq = a.deg() - b.deg();
    pr.sresize(rb, dq + 1);
    pr.sresize(ra, dq + 1);
    int pw = log2_high_bit(dq + 1) + 1;

    PT rq = mul(ra, this->_invxmod(rb, pw));
    pr.sresize(rq, dq + 1);

    PT q = rq.rev(dq + 1);
    return pr.mulc(q, ring().inv(b.lc()));
  }

  struct MultiEvaluator {
    std::vector<std::vector<PT> > polys;
    std::vector<T> res;
    const std::vector<T> &pts;
    const This &fp;
    int nd;

    MultiEvaluator(const This &fp, const std::vector<T> &pts) : fp(fp), pts(pts) {
      nd = log2_high_bit(pts.size() - 1) + 1;
      polys.resize(nd + 1);
      REP (i, 1 << nd)
        polys[nd].pb(i >= pts.size()
                       ? fp.pr.getE()
                       : fp.pr.import_vec({ fp.ring().neg(pts[i]), fp.ring().getE() }));

      REPV (i, nd + 1) {
        if (i == 1) break;
        auto &nxt = polys[i - 1];
        auto &cur = polys[i];
        int n = cur.size();
        REP (j, n / 2) nxt.emplace_back(fp.mul(cur[2 * j], cur[2 * j + 1]));
      }

      res.resize(pts.size());
    }

    void solve(const PT &poly, int depth = 0, int pos = 0) {
      if (poly.deg() <= 4) {
        pos <<= (nd - depth);
        REP (i, 1 << (nd - depth)) {
          int idx = pos + i;
          if (idx < pts.size()) res[idx] = poly(pts[idx]);
        }
        return;
      }
      int p0 = 2 * pos;
      int p1 = 2 * pos + 1;
      auto a0 = fp.mod(poly, polys[depth + 1][p0]);
      auto a1 = fp.mod(poly, polys[depth + 1][p1]);
      solve(a0, depth + 1, p0);
      solve(a1, depth + 1, p1);
    }
  };

  struct InterpolatorRes {
    Poly<T> res;
    Poly<T> xprod;
  };

  InterpolatorRes _interpolate(std::span<const T> x, std::span<const T> y) const {
    if (x.size() == 1) return InterpolatorRes{ pr.constant(y[0]), pr.xma(x[0]) };
    int n = x.size();
    int mid = n / 2;
    auto xl = x.first(mid);
    auto yl = y.first(mid);
    auto xh = x.subspan(mid);
    auto yh = y.subspan(mid);

    auto rl = this->_interpolate(xl, yl);

    auto xhv = std::vector(ALL(xh));
    auto ev1 = this->eval(rl.res, xhv);
    auto ev2 = this->eval(rl.xprod, xhv);

    std::vector<T> nyh(ALL(yh));
    auto ring = pr.underlying_ring();
    REP (i, nyh.size()) nyh[i] = ring->div(ring->sub(yh[i], ev1[i]), ev2[i]);

    auto rh = this->_interpolate(xh, nyh);

    auto res =
      InterpolatorRes{ rl.res + this->mul(rh.res, rl.xprod), this->mul(rl.xprod, rh.xprod) };
    return res;
  }

  PT interpolate(const std::vector<T> &xl, const std::vector<T> &yl) const {
    std::unordered_set<T> sx(ALL(xl));
    OPA_CHECK0(sx.size() == xl.size()); // no dups
    OPA_CHECK0(xl.size() == yl.size());

    if (xl.size() == 0) return {};

    return _interpolate(std::span(ALL(xl)), std::span(ALL(yl))).res;
  }

  std::vector<T> eval(const PT &a, const std::vector<T> &pts) const {
    if (pts.size() == 0) return {};
    if (pts.size() == 1) return { a(pts[0]) };

    MultiEvaluator mev(*this, pts);
    mev.solve(a);
    return mev.res;
  }

  PT exp_dumb(const PT &p, int n) const {
    auto &r = this->ring();
    OPA_CHECK0(r.isZ(p.get_safe(0)));

    PT res = this->pr.getE();
    PT x = p.resize(n);
    T v = r.getE();
    FOR (i, 1, n) {
      v = r.mul(v, r.importu32(i));
      if (this->ring().isZ(v)) continue;
      res = res + this->pr.mulc(x, r.inv(v));
      x = this->mul(x, p).resize(n);
    }
    return res;
  }

  T exp_eval(const PT &p, int v) const {
    auto x = p.get_safe(v);

    if (ring().isZ(x)) return x;
    return ring().mul(x, fact_provider.fact(v));
  }

  PT exp(const PT &p, int n) const {

    auto g = pr.getE();
    auto c2 = pr.importu32(2);
    auto f = p.resize(2) + ring().getE();
    auto dh = pr.derivate(p);
    int curn = 2;
    while (curn < n) {
      int nn = 2 * curn - 1;
      auto df = pr.derivate(f);
      g = this->mul(g, c2 - this->mul(g, f).resize(n));
      auto tmp = this->mul(f, this->mul(g, dh.resize(nn).resize(nn))).resize(nn) -
                 this->mul(df, g).resize(nn);
      f = f + this->mul(f, pr.integrate(tmp)).resize(nn);
      curn = nn;
    }
    return f.resize(n);
  }

  PT sqrt(const PT &p, int n) const {
    auto &r = this->ring();
    auto pp = p.resize(n);
    OPA_CHECK0(r.isE(p.get_safe(0)));
    auto f = pr.getE();
    auto g = pr.getE();
    int tn = 1;
    auto c2 = pr.importu32(2);
    auto i2 = r.inv(r.importu32(2));

    while (tn < n) {

      int nn = tn * 2;
      g = this->mul(g, c2 - this->mul(g, f).resize(n));
      auto ff = this->mul(f, f);
      auto tp = pp.resize(nn);
      auto diff = pr.mulc(tp - ff, i2);
      f = f + this->mul(diff, g).resize(nn);
      tn = nn;
    }
    return f.resize(n);
  }

  PT log(const PT &p, int n) const {
    auto &r = this->ring();
    auto pp = p.resize(n);
    OPA_CHECK0(!r.isZ(p.get_safe(0)));
    auto ip = this->invxmod(pp, n);
    auto dp = pr.derivate(pp);
    return pr.integrate(this->mul(dp, ip)).resize(n);
  }

  PT log_dumb(const PT &p, int n) const {
    auto &r = this->ring();
    OPA_CHECK0(!r.isZ(p.get_safe(0)));
    auto c = p[0];
    auto q = p;
    q[0] = r.getZ();
    auto ic = r.inv(r.neg(c));
    auto cx = r.neg(r.getE());

    PT res = this->pr.getZ();
    PT x = q.resize(n);
    T v = r.getE();
    FOR (i, 1, n) {
      cx = r.mul(cx, ic);
      auto v = r.div(cx, r.importu32(i));
      auto tmp = this->pr.mulc(x, v);
      res = res + tmp;
      x = this->mul(x, q).resize(n);
    }
    return res;
  }
};

OPA_NM_MATH_COMMON_END

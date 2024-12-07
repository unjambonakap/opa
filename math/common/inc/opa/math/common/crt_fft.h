#pragma once

#include <memory>
#include <opa/math/common/FFT.h>
#include <opa/math/common/GF_p.h>
#include <opa/math/common/PolyModRing.h>
#include <opa/math/common/PolyRing.h>
#include <opa/math/common/Ring.h>
#include <opa/math/common/UtilsPoly.h>
#include <opa/math/common/base.h>
#include <opa/utils/contest.h>
#include <opa/utils/range.h>
#include <opa_common_base.h>
#include <span>

OPA_NM_MATH_COMMON

template <class T> std::vector<T> &operator+=(std::vector<T> &a, const std::vector<T> &b) {
  a.reserve(a.size() + b.size());
  for (auto &x : b) a.pb(x);
  return a;
}

template <class T> std::vector<T> operator+(const std::vector<T> &a, const std::vector<T> &b) {
  std::vector<T> res = a;
  res.reserve(a.size() + b.size());
  for (auto &x : b) res.pb(x);
  return res;
}

static constexpr int kMulGroup = true;
template <class R, bool MulGroup, class T = typename R::Type> class RingOpsProvider {
public:
  const R &r;
  T e;
  RingOpsProvider(const R &r) : r(r) { e = MulGroup ? r.getE() : r.getZ(); }
  T mul(const T &a, const T &b) const { return MulGroup ? r.mul(a, b) : r.add(a, b); }
  T div(const T &a, const T &b) const { return MulGroup ? r.div(a, b) : r.sub(a, b); }
  T inv(const T &a) const { return MulGroup ? r.inv(a) : r.neg(a); }
};

template <class R, class T = typename R::Type> class RingOps {
public:
  const R &r;
  T e;
  RingOps(const R &r) : r(r) { e = r.getE(); }
  T mul(const T &a, const T &b) const { return r.mul(a, b); }
  T add(const T &a, const T &b) const { return r.add(a, b); }
  T div(const T &a, const T &b) const { return r.div(a, b); }

  T agg_mul(const std::vector<T> &v) const { return QQ::fold_left(v, r.getE(), STD_FUNC2(r.mul)); }
  T agg_add(const std::vector<T> &v) const { return QQ::fold_left(v, r.getZ(), STD_FUNC2(r.add)); }
  T vec_mul(const std::vector<T> &l, const std::vector<T> &r) const {
    return LQ::zip(l, r) | STD_TSFTUPLE(r.mul(p0, p1)) | STD_VEC;
  }
  T vec_add(const std::vector<T> &l, const std::vector<T> &r) const {
    return LQ::zip(l, r) | STD_TSFTUPLE(r.add(p0, p1)) | STD_VEC;
  }
};

template <class T> std::vector<T> one_hot(T v, int pos, int n, T tz = {}) {
  std::vector<T> res(n, tz);
  res[pos] = v;
  return res;
}

template <class T, class BaseRing> struct MatrixGroup {
  typedef Matrix<T, BaseRing> MT;
  typedef MatrixGroup<T, BaseRing> This;
  MT e;
  MT mul(const MT &a, const MT &b) const { return b * a; }
  MT inv(const MT &a) const { return a.inverse(); }
  MT div(const MT &a, const MT &b) const { return a * b.inverse(); }
};

VV_t(int) mat_ids(int n, int m) {
  return opa::utils::product<int>({ STD_RANGE(0, n) | STD_VEC, STD_RANGE(0, m) | STD_VEC });
}

template <class T, class BaseRing> VV_t(int) mat_ids(const Matrix<T, BaseRing> &m) {
  return mat_ids(m.n, m.m);
}

template <class T, typename F, typename R> void init_mat(Matrix<T, R> &m, F &&fill) {
  STD_FOREACHX(mat_ids(m), m(x[0], x[1]) = fill(x[0], x[1]));
}
template <typename F, typename R, class T = typename R::Type>
Matrix<T, R> create_mat(const R *r, int n, int m, F &&fill) {
  Matrix<T, R> mat(r, n, m);
  init_mat(mat, std::forward<F &&>(fill));
  return mat;
}

template <class T, typename F, typename R> Matrix<T, R> omap(const Matrix<T, R> &m, F &&tsf) {
  auto res = m.clone();
  STD_FOREACHX(mat_ids(m), res(x[0], x[1]) = tsf(m(x[0], x[1])));
  return res;
}

template <typename F, class T = typename std::invoke_result<F, int, int>::type>
VV_t(T) ids_vec(int n, int m, F &&tsf) {
  auto res = VV_t(T)(n, V_t(T)(m));
  STD_FOREACHX(mat_ids(n, m), res[x[0]][x[1]] = tsf(x[0], x[1]));
  return res;
}

template <class T, typename F, typename R, typename U = typename std::invoke_result<F, T>::type>
std::vector<std::vector<U> > omap_vec(const Matrix<T, R> &m, F &&tsf) {
  auto res = VV_t(U)(m.n, V_t(U)(m.m));
  STD_FOREACHX(mat_ids(m), res[x[0]][x[1]] = tsf(m(x[0], x[1])));
  return res;
}
template <class T, typename F, typename R, typename RU,
          typename U = typename std::invoke_result<F, T>::type>
Matrix<U, RU> omap(const Matrix<T, R> &m, RU *u, F &&tsf) {
  auto res = Matrix<U, RU>(u, m.n, m.m);
  STD_FOREACHX(mat_ids(m), res(x[0], x[1]) = tsf(m(x[0], x[1])));
  return res;
}

template <typename T, typename OpsProvider> struct CumulativeQuery {
  const OpsProvider &g;
  std::vector<T> q;
  std::vector<T> iq;

  T query(int a, int b) const {
    OPA_CHECK(a >= 0 && b >= 0, a, b);
    OPA_CHECK(a < q.size() && b + 1 < q.size(), a, b, q.size());
    return g.mul(q[b + 1], iq[a]);
  }
};
template <typename T, typename OpsProvider>
static CumulativeQuery<T, OpsProvider> MakeCumulativeQuery(const OpsProvider &g,
                                                           const std::vector<T> &data) {
  auto q = data | LQ::partial_sum(STD_FUNC2(g.mul)) | STD_VEC;
  q.insert(q.begin(), g.e);
  return CumulativeQuery<T, OpsProvider>{
    .g = g,
    .q = q,
    .iq = q | STD_TSFX(g.inv(x)) | STD_VEC,
  };
}

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

  std::unique_ptr<IFFTProvider<T> > make_fft(int npw) const {
    return std::unique_ptr<IFFTProvider<T> >(new FFT2Dispatcher<T>(this, npw));
  }

  std::unique_ptr<IFFTProvider<TE> > make_ffte_dispatcher(int npw) const {
    auto tmp = [&](int npw_req) { return this->make_fft_provider(npw_req).release(); };
    return std::unique_ptr<IFFTProvider<TE> >(new FFT2Dispatcher<TE>(npw, tmp));
  }

  std::unique_ptr<IFFTProvider<TE> > make_fft_provider(int npw) const {
    auto tmp = new CRTFieldFFTProvider<TE, N>(*this, npw);
    return std::unique_ptr<IFFTProvider<TE> >(tmp);
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
  std::unique_ptr<IFFTProvider<TLift> > _fft;

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

CRTField<u32, 3> crtf1(const GF_pT<u32> &target) {
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
      _fact.pb(ring.mul(_fact.back(), ring.importu32(_fact.size())));
    }
    return _fact[i];
  }
};

template <class T, class Op> T faste_t(T a, u64 p, Op mul, T id) {
  T x = id;
  for (; p != 0; p >>= 1, a = mul(a, a))
    if (p & 1) x = mul(x, a);
  return x;
}

template <class T, class BaseRing = Ring<T> > class FastPoly {

  typedef PolyRing<T, BaseRing> PR;
  const IFFTProvider<T> &provider;
  std::unique_ptr<FactProvider<T> > owned_fact_provider;

public:
  const FactProvider<T> *fact_provider;
  typedef Poly<T, BaseRing> PT;
  typedef FastPoly<T, BaseRing> This;
  const PR &pr;

  const BaseRing &ring() const { return *pr.underlying_ring(); }

  FastPoly(const PR &pr, const IFFTProvider<T> &provider) : pr(pr), provider(provider) {
    owned_fact_provider = std::make_unique<FactProvider<T> >(*pr.underlying_ring());
    fact_provider = owned_fact_provider.get();
  }

  PT mul(const PT &a, const PT &b) const {
    if (a.size() == 0 || b.size() == 0) return pr.getZ();
    return pr.import_vec(provider.mul(a.to_vec(), b.to_vec()));
  }

  PT tower_prod(std::vector<PT> lst) const {
    if (lst.size() == 1) return lst[0];
    auto mid = lst.begin() + lst.size() / 2;
    return this->mul(tower_prod(V_t(PT)(lst.begin(), mid)), tower_prod(V_t(PT)(mid, lst.end())));
  }

  PT pow(const PT &a, int pw) const {
    auto pts = provider.fft(a.to_vec(a.deg() * pw + 1));
    auto fa = pts | STD_TSFX(fact_provider->ring.faste(x, pw)) | STD_VEC;
    return pr.import_vec(provider.ifft(fa));
  }

  PT mulxmod(const PT &a, const PT &b, int xpw) const { return pr.resize(this->mul(a, b), xpw); }

  PT invxmod(const PT &x, int n) const { return _invxmod(x, log2_high_bit(n - 1) + 1).resize(n); }

  PT _invxmod(const PT &x, int degpw) const {
    OPA_CHECK0(ring().isInv(x[0]));
    if (degpw == 0) return pr.constant(ring().inv(x[0]));
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
      s64 cost_stupid = poly.size() * (1ull << (nd - depth));
      double l2 = std::log2(poly.size());

      if (poly.deg() < 4 || cost_stupid < l2 * l2 * poly.deg()) {
        pos <<= (nd - depth);
        REP (i, 1 << (nd - depth)) {
          int idx = pos + i;
          if (idx < pts.size()) res[idx] = fp.pr.eval(poly, pts[idx]);
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
    PT res;
    PT xprod;
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
    if (pts.size() == 1) return { pr.eval(a, pts[0]) };
    if (a.deg() == 0) return std::vector<T>(pts.size(), a.get_safe(0));
    double cost = std::log2(pts.size());
    cost = cost * cost;
    if (a.deg() < cost) return pts | STD_TSFX(pr.eval(a, x)) | STD_VEC;

    MultiEvaluator mev(*this, pts);
    mev.solve(a);
    return mev.res;
  }

  PT shift_poly(const PT &a, const T &x) const {
    // a(t) => a(x+t)
    auto pts = STD_RANGE(0, a.deg() + 1) | STD_VECT(u32);
    auto yv = this->eval(a, pts);
    return this->interpolate(pts | STD_TSFY(this->ring().sub(y, x)) | STD_VEC, yv);
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

  PT exp2normal(int n) const {
    std::vector<T> res(n);
    REP (i, n) res[i] = fact_provider->fact(i);
    return pr.import_vec(res);
  }

  PT normal2exp(int n) const {
    std::vector<T> res(n);
    REP (i, n) res[i] = fact_provider->ifact(i);
    return pr.import_vec(res);
  }

  PT normal2exp(const PT &p) const { return pr.pointwise_mul(p, normal2exp(p.size())); }
  PT exp2normal(const PT &p) const { return pr.pointwise_mul(p, exp2normal(p.size())); }

  T exp_eval(const PT &p, int v) const {
    auto x = p.get_safe(v);

    if (ring().isZ(x)) return x;
    return ring().mul(x, fact_provider->fact(v));
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

  std::vector<T> linear_rec_kth(const std::vector<T> &f, const PT &g, u64 k, int m = -1) const {
    // g annihilator
    // https://codeforces.com/blog/entry/115696
    auto gr = g.rev();
    int n = g.deg();
    if (m == -1) m = f.size() - n;
    int needf = m + n;
    if (f.size() < needf) {
      auto h = this->mulxmod(pr.import_vec(f), g, n);
      auto nf = this->mulxmod(h, this->invxmod(g, needf), needf);
      return linear_rec_kth(nf.to_vec(needf), g, k, m);
    }
    if (f.size() > needf) {
      auto nf = f;
      nf.resize(needf);
      return linear_rec_kth(nf, g, k, m);
    }
    auto xs = faste_t(pr.x(), k, STD_FUNCXY(this->mod(this->mul(x, y), gr)), pr.getE());
    auto res = this->mul(xs, pr.import_vec(f, kPolyRev)).to_vec(m + 2 * n);
    res.resize(m + n);
    std::reverse(ALL(res));
    res.resize(m);
    return res;
  }

  std::vector<T> linear_rec_kth_poly(const std::vector<T> &f, u64 k, int m = -1) const {
    u32 n = f.size() - 1;
    if (m == -1) m = n;
    if (n == 0) return {};
    if (QQ::all_of(f, STD_FUNCX(x == f[0]))) return STD_REPEAT(f[0], m) | STD_VEC;
    if (k <= n) {
      auto fv = std::vector<T>(f.begin() + k, f.end());
      return fv + linear_rec_kth_poly(f, n + 1, m - k);
    }
    auto ro = RingOpsProvider<BaseRing, kMulGroup>(ring());
    auto cq1 = MakeCumulativeQuery(ro, STD_RANGE(1u, n + 1) | STD_VECT(T));
    auto cqx = MakeCumulativeQuery(ro, STD_RANGE(k - n, k + m + n + 1) | STD_VECT(T));

    auto &r = ring();
    auto fp = STD_RANGE(0, n + 1) |
              STD_TSFX(r.mulv({ f[x], cq1.iq[x], cq1.iq[n - x], r.sgn_v(r.getE(), n - x) })) |
              STD_VEC;
    auto gp = STD_RANGE(k - n, k + m + 1) | STD_TSFX(r.inv(x)) | STD_VEC;
    OPA_DISP0(fp.size(), gp.size());
    auto px = this->mul(pr.import_vec(fp), pr.import_vec(gp));
    auto res = STD_RANGE(0, m) | STD_TSFX(r.mul(px.get_safe(n + x), cqx.query(x, x + n))) | STD_VEC;
    return res;
  }
};

template <class T, class BaseRing>
std::vector<T> partial_fraction_expansion(const opa::math::common::FastPoly<T, BaseRing> &fp,
                                          const Poly<T, BaseRing> &target,
                                          const std::vector<T> &roots) {
  // ret ki, ki/(x-ri) = target

  auto gfp = fp.ring();
  T prod = QQ::fold_left(roots, (T)1, STD_FUNCXY(gfp.mul(x, gfp.neg(y))));

  auto qinv =
    fp.interpolate(LQ::concat(roots, STD_SINGLEV(T(0))) | STD_VEC,
                   LQ::concat(STD_REPEAT((T)0, roots.size()), STD_SINGLEV(prod)) | STD_VEC);
  auto qd = fp.pr.derivate(qinv);
  auto qde = fp.eval(qd, roots);
  auto te = fp.eval(target, roots);
  auto tmp = LQ::zip(te, qde) | STD_TSFTUPLE(gfp.div(p0, p1)) | STD_VEC;

  return tmp;
}

template <class T, class BaseRing> struct PRecursirveSeq {
  Poly<T, BaseRing> p0;
  Poly<T, BaseRing> pconst;
  std::vector<Poly<T, BaseRing> > coeffs;
  std::vector<T> c0;
  int n() const { return c0.size(); }
  PRecursirveSeq<T, BaseRing> shift(T t0) {
    auto pr = p0.get_poly_ring();
    auto nx = pr->xpa(t0);
    return PRecursirveSeq<T, BaseRing>{ .p0 = pr->eval2(p0, nx),
                                        .pconst = pr->eval2(pconst, nx),
                                        .coeffs = coeffs | STD_TSFX(pr->eval2(x, nx)) | STD_VEC,
                                        .c0 = c0 };
  }
};

template <class TIn, class TState, class TFunc>
std::vector<TState> tsf_accumulate(const std::vector<TIn> &in, const TState &s0,
                                   const TFunc &tfunc) {
  struct Data {
    TState state;
  };

  Data data{ .state = s0 };
  return in | STD_TSFX(data.state = tfunc(data.state, x)) | STD_VEC;
}

template <class T, class BaseRing> struct PRecursirveSeqComp {
  typedef Poly<T, BaseRing> PT;
  typedef PolyRing<T, BaseRing> PR;
  typedef Matrix<T, BaseRing> MT;
  typedef Matrix<PT, PR> MPT;

  static constexpr int t0 = 1;
  struct BlockEntry {
    int log2_B;
    int step;
    std::vector<MT> mats;
    MPT mpt;
  };

  struct PrecompData {
    std::vector<MT> mats;
    int step;
    int log2_B;
    std::vector<u32> eval_pos;
  };

  BlockEntry mat_be, norm_be;
  bool has_const;

  std::vector<T> mat2state(const std::vector<T> &ms) const {
    return ms | LQ::reverse | LQ::drop(has_const ? 1 : 0) | STD_VEC;
  }

  std::vector<T> state2mat(const std::vector<T> &state) const {
    auto cur = state | LQ::reverse | STD_VEC;
    if (has_const) cur.pb(f->getE());
    return cur;
  }

  BlockEntry compute_be(const MPT &mat, int max_bound) const {

    auto &pr = this->fp->pr;

    auto mids = mat_ids(mat);
    int deg = QQ::max(mids | STD_TSFX(mat(x[0], x[1]).deg()));
    OPA_DISP0(mat);

    auto shift = [&](const std::vector<MT> &mats, int s) -> std::vector<MT> {
      VVV_t(T) vals = ids_vec(
        mat.n, mat.m,
        STD_FUNCAB(fp->linear_rec_kth_poly(mats | STD_TSFY(y(a, b)) | STD_VEC, s, mats.size())));
      std::vector<MT> res(mats.size(), MT(f, mat.n, mat.m));
      STD_FOREACHX(STD_RANGE(0, mats.size()), init_mat(res[x], STD_FUNCAB(vals[a][b][x])));
      return res;
    };

    std::vector<MT> cur =
      STD_RANGE(0, deg + 1) | STD_TSFX(omap(mat, f, STD_FUNCY(this->fp->pr.eval(y, x)))) | STD_VEC;

    int log2_B = 0;
    while ((1ull << log2_B) * cur.size() < max_bound) {
      ++log2_B;
      auto s1 = shift(cur, cur.size());
      auto s2 = shift(cur, 2 * cur.size());
      auto s3 = shift(cur, 3 * cur.size());
      cur += s1 + s2 + s3;
      cur = STD_RANGE(0, cur.size() / 2) | STD_TSFX(cur[2 * x + 1] * cur[2 * x]) | STD_VEC;
      cur.pop_back();
      OPA_DISP0(cur.size(), cur.size() * (1ull << (log2_B)));
    }
    OPA_DISP0(cur.size(), 1 << log2_B, cur.size() * (1ull << log2_B));
    return BlockEntry{ .log2_B = log2_B, .step = 1 << log2_B, .mats = cur, .mpt = mat };
  }

  std::vector<T> eval_be(const BlockEntry &be, std::vector<T> cur, int t) const {
    int tcur = 0;
    STD_FOREACHX(STD_RANGE(0, be.mats.size()) |
                   LQ::filter(STD_FUNCY(1ull * (y + 1) * be.step <= t)),
                 ({ tcur += be.step, cur = be.mats[x].eval(cur); }));
    OPA_DISP0(tcur, t, cur, be.mats.size(), t / be.step);
    for (; tcur < t; ++tcur) cur = omap(be.mpt, f, STD_FUNCY(fp->pr.eval(y, tcur))).eval(cur);
    return cur;
  }

  opa::math::common::FastPoly<T, BaseRing> *fp;
  BaseRing *f;
  PRecursirveSeq<T, BaseRing> seq;

  T next(const std::vector<T> &vals, u32 t) const {
    auto cur = f->getZ();
    T ft = f->importu32(t);

    auto pvals = seq.coeffs | STD_TSFX(fp->pr.eval(x, ft)) | STD_VEC;
    auto p0v = fp->pr.eval(seq.p0, ft);
    auto r = f->add(f->dot(vals, pvals), fp->pr.eval(seq.pconst, ft));
    return f->div(r, p0v);
  }

  std::vector<T> run_update(std::vector<T> vals, u32 st, u32 nt) const {
    FOR (t, st, nt) {
      vals.pb(
        this->next(std::span{ vals }.subspan(vals.size() - seq.n()) | LQ::reverse | STD_VEC, t));
    }
    return vals;
  }

  MPT seq2mat(const PRecursirveSeq<T, BaseRing> &seq) const {
    auto &pr = this->fp->pr;
    int sz = seq.n() + (has_const ? 1 : 0);
    auto mat = MPT(&pr, sz, sz);
    init_mat(mat, [&](int i, int j) -> PT {
      if (i == seq.n()) return j == seq.n() ? pr.getE() : pr.getZ();
      if (i == 0) {
        if (j == seq.n()) return pr.mul(seq.pconst, seq.p0);
        return seq.coeffs[j];
      }
      return i == j + 1 ? seq.p0 : pr.getZ();
    });
    return mat;
  }

  void setup(int maxbound) {

    auto &pr = this->fp->pr;

    has_const = seq.pconst != pr.getZ();
    mat_be = compute_be(seq2mat(seq), maxbound);
    auto mat_norm = MPT(&pr, 1, 1);
    mat_norm(0, 0) = seq.p0;
    norm_be = compute_be(mat_norm, maxbound);
  }

  T compute(u32 pos) {
    if (pos < seq.n()) return seq.c0[pos];
    pos -= seq.n() - 1;
    auto a1 = eval_be(mat_be, state2mat(seq.c0), pos);
    auto a2 = eval_be(norm_be, { f->getE() }, pos);
    return f->div(a1[0], a2[0]);
  }

  /*
  if (0) {
    T c2 = f->importu32(2);
    int B = 10;

    OPA_DISP0("A2.5");
    REP (i, precomp.log2_B) {
      auto m2 = omap(mat, [&](const auto &a) { return fp->shift_poly(a, 1 << i); });
      mat = m2 * mat;
    }

    OPA_DISP0("A3.0");
    if (0) {
      int poly_size = B * (seq.p0.deg() + QQ::max(seq.coeffs | STD_TSFX(x.deg())) + 1);
      auto mats = STD_RANGE(1, poly_size + B + 1) |
                  STD_TSFX(omap(mat, f, STD_FUNCY(fp->pr.eval(y, x)))) | STD_VEC;
      auto mg = MatrixGroup<T, BaseRing>{ .e = mats[0].identity() };
      auto mcq = MakeCumulativeQuery(mg, mats);
      auto xvs = STD_RANGE(1, poly_size + 1) | STD_VECT(T);
      auto yvs = xvs | STD_TSFX(mcq.query(x - 1, x - 1 + B - 1)) | STD_VEC;
      init_mat(mat, STD_FUNCAB(fp->interpolate(xvs, yvs | STD_TSFX(x(a, b)) | STD_VEC)));
    }

    OPA_DISP0("A3");
    int tx = t0;
    while (tx < maxbound) {
      precomp.eval_pos.pb(tx);
      tx += precomp.step;
    }
    OPA_DISP0("A4");
    auto evals = omap_vec(mat, [&](const PT &px) { return fp->eval(px, precomp.eval_pos); });

    auto ids = STD_RANGE(0, precomp.eval_pos.size()) | STD_VEC;
    precomp.mats =
      ids | STD_TSFX(create_mat(f, mat.n, mat.m, STD_FUNCAB(evals[a][b][x]))) | STD_VEC;
  }
  */
  // precompute done
};

OPA_NM_MATH_COMMON_END

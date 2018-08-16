#ifndef _H_OPA_MATH_COMMON_UTILS_POLY
#define _H_OPA_MATH_COMMON_UTILS_POLY

#include <opa/math/common/FFT.h>
#include <opa/math/common/Field.h>
#include <opa/math/common/Matrix.h>
#include <opa/math/common/Poly.h>
#include <opa/math/common/PolyRing.h>
#include <opa/math/common/bignum.h>
#include <opa/math/common/utils_num.h>
#include <opa/utils/DataStruct.h>
#include <opa/utils/misc.h>

OPA_NAMESPACE_DECL3(opa, math, common)

template <class T> Poly<T> findIrred(const Field<T> &f, int deg);

template <class T>
bool checkIrredDummy(const Field<T> &field, const Poly<T> &a);

template <class T> Poly<T> findIrred(const Field<T> &f, int deg) {
  PolyRing<T> pr(&f);
  while (1) {
    Poly<T> p = pr.rand(deg);
    if (pr.isIrred(p)) {
      if (f.getSize() <= 4) assert(checkIrredDummy(f, p));
      return p;
    }
  }
}

template <class T> Poly<T> find_primitive_poly(const Field<T> &f, int deg) {
  PolyRing<T> pr(&f);
  while (1) {
    Poly<T> tmp = pr.import(findIrred(f, deg));
    bignum pw = f.getSize().pow(deg) - 1;
    BGFactors factors = factor_large(pw);

    bool ok = 1;
    for (auto &factor : factors) {
      if (pr.isE(pr.faste(pr.x(), pw / factor.ST, tmp))) {
        ok = 0;
        break;
      }
    }
    if (ok) return tmp;
  }
}

template <class T>
bool go(const PolyRing<T> &pr, const Poly<T> &a, std::vector<T> &lst) {
  Poly<T> p = pr.import(lst);
  // printf("CHECKING >> "); p.disp();
  // if (p.deg()>=1){printf(">>>> RES "); pr.mod(a,p).disp();}
  if (p.deg() >= 1)
    if (pr.isZ(pr.mod(a, p))) {
      puts("BIG FAIL");
      p.disp();
      return true;
    }
  if (lst.size() >= a.size() / 2) return false;
  lst.push_back(0);
  for (int i = 0; i < pr.underlying_ring()->getSizeU32(); ++i) {
    lst.back() = i;
    if (go(pr, a, lst)) return true;
  }
  lst.pop_back();
  return false;
}

template <class T>
bool checkIrredDummy(const Field<T> &field, const Poly<T> &a) {
  int sz = a.size() / 2;
  PolyRing<T> pr(&field);
  std::vector<T> tb;
  return !go(pr, a, tb);
}

template <class T> std::vector<T> genRandVector(int n, const Field<T> &field) {
  std::vector<T> res(n, field.getZ());
  for (int i = 0; i < n; ++i) res[i] = field.getRand();
  return res;
}

template <class T>
std::vector<T> genRandWeight(int n, int weight, const Field<T> &field) {
  std::vector<T> res(n, field.getZ());
  for (int i = 0; i < n && weight; ++i) {
    int m = (n - i) * weight;
    int k = rng() % m;
    if (k < m * weight / (n - i)) {
      res[i] = field.getRandNZ();
      --weight;
    }
  }
  return res;
}

template <class T, typename RealType> class PolyFFTHelper_WithReal {
public:
  typedef Poly<T> PT;
  const Field<T> *m_f;
  FFT_Real<RealType> m_fft;

  void init(const Field<T> *f, int maxd = -1, RealType thresh = 1e-9) {
    bignum order = f->getSize() - 1;
    if (maxd == -1) maxd = order.getu32();
    m_f = f;

    m_fft.init(maxd, thresh);
  }

  std::vector<RealType> to_realtype(const std::vector<T> &a) const {
    std::vector<RealType> res;
    for (auto &x : a) res.emplace_back(x);
    return res;
  }

  std::vector<T> from_realtype(const std::vector<RealType> &a) const {
    std::vector<T> res;
    for (auto &x : a)
      res.push_back(m_f->import(FloatUtil::cast<RealType, T>(round(x))));
    return res;
  }

  std::vector<T> multiply(const std::vector<T> &a,
                          const std::vector<T> &b) const {
    std::vector<RealType> aa, bb;
    aa = to_realtype(a);
    bb = to_realtype(b);

    auto ia = m_fft.fft(aa);
    auto ib = m_fft.fft(bb);

    REP (i, ia.size()) { ia[i] = m_fft.cf.mul(ia[i], ib[i]); }
    auto r = m_fft.ifft(ia);
    return from_realtype(r);
  }
};

template <class T> class PolyFFTPlanner {
public:
  typedef Poly<T> PT;
  class FFTExecutor;

  struct ExecutionData {
    std::vector<T> input;
    std::vector<T> output;
    bool inverse = false;
  };

  class FFTExecutor : public utils::Initable {
  public:
    FFTExecutor(const Field<T> *f, const T &root, int order) {
      m_f = f;
      m_root = root;
      m_iroot = f->inv(root);
      m_order = order;
    }
    virtual ~FFTExecutor() {}
    virtual void execute(ExecutionData *ex_data) const {
      this->do_execute(ex_data);
    }

    virtual std::vector<T> multiply(const std::vector<T> &a,
                                    const std::vector<T> &b) const {
      ExecutionData ed1, ed2, ied;
      ed1.input = a;
      ed2.input = b;

      ed1.input.resize(this->m_order, this->m_f->getZ());
      ed2.input.resize(this->m_order, this->m_f->getZ());

      this->execute(&ed1);
      this->execute(&ed2);

      REP (i, this->m_order)
        ied.input.push_back(this->m_f->mul(ed1.output[i], ed2.output[i]));
      ied.inverse = true;
      this->execute(&ied);

      T iv = this->m_f->inv(this->m_f->importu32(this->m_order));
      REP (i, ied.output.size()) {
        ied.output[i] = this->m_f->mul(ied.output[i], iv);
      }
      return ied.output;
    }

    virtual void do_execute(ExecutionData *ex_data) const = 0;

    std::vector<T> compute_roots(const Field<T> *f, const T &root, int num) {
      std::vector<T> roots;
      T cur = f->getE();
      REP (i, num) {
        roots.push_back(cur);
        cur = f->mul(cur, root);
      }
      return roots;
    }

    T root(bool inv) const { return inv ? m_iroot : m_root; }

    const Field<T> *m_f;
    T m_root;
    T m_iroot;
    int m_order;
  };

  class SplitFFTExecutor : public FFTExecutor {
  public:
    struct SubData {
      u32 l;
      const FFTExecutor *subex;
    };

    SplitFFTExecutor(const Field<T> *f, const T &root, const SubData &d1,
                     const SubData &d2)
        : FFTExecutor(f, root, d1.l * d2.l) {
      m_d1 = d1;
      m_d2 = d2;
      N1 = m_d1.l;
      N2 = m_d2.l;
      m_roots = this->compute_roots(f, root, N1 * N2);
      m_iroots = this->compute_roots(f, this->m_iroot, N1 * N2);
    }

    virtual void do_execute(ExecutionData *ex_data) const override {
      // ex_data->output = m_fft.fft(ex_data.input);
      ex_data->output.resize(N1 * N2);
      OPA_CHECK0(m_d1.subex->m_root == this->m_f->faste(this->m_root, N2));
      OPA_CHECK0(m_d2.subex->m_root == this->m_f->faste(this->m_root, N1));

      std::vector<std::vector<T> > A;
      A.resize(N2, std::vector<T>(N1));
      REP (i, N1) {
        ExecutionData ed;
        ed.inverse = ex_data->inverse;
        REP (j, N2)
          ed.input.push_back(ex_data->input[i + j * N1]);
        m_d2.subex->execute(&ed);
        REP (j, N2)
          A[j][i] = this->m_f->mul(ed.output[j], root(i * j, ex_data->inverse));
      }

      /*
         p0 = a0 + a3 x3 + a6 x6
         p1 = a1 + a4 x3 + a7 x6
         p2 = a2 + a5 x3 + a8 x6

         bi = p0(wi) + wi * p1(wi) + w2i * p2(wi)
         b(n3 +i) = p0(wi) + wi * wn3 * p1(wi) + w2i * w(2n3) * p2(wi)
DFT1:
  p0, p1
DFT2:
         */

      REP (i, N2) {
        ExecutionData ed;
        ed.input = A[i];
        ed.inverse = ex_data->inverse;
        m_d1.subex->execute(&ed);
        REP (j, N1)
          ex_data->output[j * N2 + i] = ed.output[j];
      }

      if (0) {
        SimpleFFTExecutor tmp(this->m_f, this->m_root, this->m_order);
        ExecutionData e2 = *ex_data;
        tmp.execute(&e2);
        OPA_CHECK(e2.output == ex_data->output, e2.output, ex_data->output,
                  e2.input, ex_data->input);
      }
    }

    T root(int i, bool inverse) const {
      return inverse ? m_iroots[i] : m_roots[i];
    }

    int N1, N2;
    std::vector<T> m_roots;
    std::vector<T> m_iroots;
    SubData m_d1;
    SubData m_d2;
  };

  class CRTFFTExecutor : public FFTExecutor {
  public:
    struct SubData {
      std::vector<const FFTExecutor *> subex;
    };

    static constexpr int lim_simple_fft = 17;

    CRTFFTExecutor(const Field<T> *f, const T &root, int order,
                   const SubData &subdata)
        : FFTExecutor(f, root, order) {
      m_subdata = subdata;
      std::vector<bignum> primes;
      for (auto &x : subdata.subex) {
        primes.emplace_back(x->m_f->getSize());
      }
      crt.init(primes);
    }

    virtual std::vector<T> multiply(const std::vector<T> &a,
                                    const std::vector<T> &b) const {

      // ex_data->output = m_fft.fft(ex_data.input);
      int N = a.size() + b.size() - 1;
      std::vector<std::vector<bignum> > tb(
        N, std::vector<bignum>(m_subdata.subex.size(), 0));

      REP (i, m_subdata.subex.size()) {
        auto cur = m_subdata.subex[i];

        std::vector<T> na, nb;
        REP (j, a.size())
          na.push_back(cur->m_f->import(a[j]));
        REP (j, b.size())
          nb.push_back(cur->m_f->import(b[j]));
        auto cres = cur->multiply(na, nb);

        if (0) {
          SimpleFFTExecutor tmp(cur->m_f, cur->m_root, cur->m_order);
          auto ctmp = tmp.multiply(na, nb);
          OPA_CHECK(ctmp == cres, ctmp, cres, cur->m_f->getSizeU32(),
                    cur->m_root, cur->m_order);
        }

        REP (j, N)
          tb[j][i] = cres[j];
      }

      std::vector<T> res(N);
      REP (i, N)
        res[i] = this->m_f->import_bg(crt.solve(tb[i]));
      return res;
    }

    virtual void do_execute(ExecutionData *ex_data) const override {
      OPA_CHECK0(false);
    }

    CRT_Batched crt;
    SubData m_subdata;
  };

  class SimpleFFTExecutor : public FFTExecutor {
  public:
    SimpleFFTExecutor(const Field<T> *f, const T &root, int order)
        : FFTExecutor(f, root, order) {
      m_fft.init(f, this->m_order, root);
    }

    virtual void do_execute(ExecutionData *ex_data) const override {
      if (ex_data->inverse) {
        ex_data->output = m_fft.ifft(ex_data->input, false);
      } else {
        ex_data->output = m_fft.fft(ex_data->input);
      }
    }

    FFT<T> m_fft;
  };

  class FFTExecutorPW2 : public FFTExecutor {
  public:
    FFTExecutorPW2(const Field<T> *f, const T &root, int pw2)
        : FFTExecutor(f, root, 1 << pw2) {
      m_fft.init(f, pw2, root);
    }

    virtual void do_execute(ExecutionData *ex_data) const override {
      if (ex_data->inverse) {
        ex_data->output = m_fft.ifft(ex_data->input, false);
      } else {
        ex_data->output = m_fft.fft(ex_data->input);
      }
    }

    FFT2<T> m_fft;
  };

  static constexpr bool kWantDft = true;
  static constexpr bool kForceSimple = true;
  static constexpr bool kCheckSplit = true;

  void init(const Field<T> *f, int order = -1,
            bool force_simple = !kForceSimple, bool want_dft = kWantDft) {
    m_f = f;
    m_force_simple = force_simple;

    if (want_dft) {
      if (order == -1) order = (m_f->getSize() - 1).getu32();
      m_factors = f->get_factors();
      SmallestDivHelper helper;
      helper.compute(m_factors, order);
      order = helper.res.v.getu32();
    }

    T root = f->getPrimitiveElem();
    m_order = order;
    OPA_CHECK0(order != 0);
    m_root = m_f->faste(root, (m_f->getSize() - 1) / order);

    main_executor = build(m_f, m_root, order, m_factors, kCheckSplit, want_dft);
  }

  FFTExecutor *build(const Field<T> *f, const T &root, int order,
                     const BGFactors &factors, bool check_split,
                     bool want_dft) {
    OPA_CHECK(!want_dft || f->getOrderOf(root) == order, root, order);

    check_split &=
      want_dft && !m_force_simple && (order > CRTFFTExecutor::lim_simple_fft);

    if (check_split) {

      std::vector<bignum> order_factors = factor_known(order, factors);
      bool is_smooth = true;
      for (auto factor : order_factors) {
        if (factor > CRTFFTExecutor::lim_simple_fft) {
          is_smooth = false;
          break;
        }
      }

      if (is_smooth && order_factors.size() != 1) {
        bignum border(order);
        std::sort(ALL(order_factors));

        int pw2_end = 0;
        bignum curdiv = 1;
        for (; pw2_end < order_factors.size() && order_factors[pw2_end] == 2;
             ++pw2_end)
          curdiv *= 2;

        FFTExecutor *last;
        if (1 || pw2_end == 0) {
          curdiv = order_factors[0];
          last = build(f, f->faste(root, border / curdiv), curdiv.getu32(), {},
                       !kCheckSplit, kWantDft);
          pw2_end = 1;
        } else {
          last = this->build_pw2(f, f->faste(root, border / curdiv), pw2_end);
        }

        for (int i = pw2_end; i < order_factors.size(); ++i) {
          FFTExecutor *next =
            build(f, f->faste(root, border / order_factors[i]),
                  order_factors[i].getu32(), {}, !kCheckSplit, kWantDft);

          last = executors.add(new SplitFFTExecutor(
            f, f->faste(root, border / (curdiv * order_factors[i])),
            { order_factors[i].getu32(), next }, { curdiv.getu32(), last }));
          curdiv *= order_factors[i];
        }
        return last;
      }
    }

    if (want_dft) {
      return executors.add(new SimpleFFTExecutor(f, root, order));
    }

    bignum bound = f->getSize() * f->getSize() * order;
    std::vector<SmoothPrimeOrderHelper::Entry> crt_fact =
      this->get_crt_split(bound, order);

    typename CRTFFTExecutor::SubData subdata;
    for (auto &crt_entry : crt_fact) {
      Field<T> *nf = fields.add(new GF_p(crt_entry.first));
      SmallestDivHelper helper;

      auto nfactors = to_bgfactors(crt_entry.second);
      helper.compute(nfactors, order);
      int norder = helper.res.v.getu32();
      T nroot = nf->getNthRoot(norder);

      subdata.subex.push_back(
        this->build(nf, nroot, norder, nfactors, kCheckSplit, kWantDft));
    }
    return executors.add(new CRTFFTExecutor(f, root, order, subdata));
  }

  FFTExecutor *build_pw2(const Field<T> *f, const T &root, int pw2) {
    return executors.add(new FFTExecutorPW2(f, root, pw2));
  }

  std::vector<SmoothPrimeOrderHelper::Entry> get_crt_split(const bignum &bound_,
                                                           u32 order) {
    std::vector<SmoothPrimeOrderHelper::Entry> res;
    bignum bound = bound_;
    const auto &smooth_data = spoh_for_gfp_split.data();
    int pos = std::lower_bound(
                ALL(smooth_data),
                SmoothPrimeOrderHelper::Entry{ bound.get_or<s32>(1e9), {} }) -
              smooth_data.begin();

    const int maxcnd = 100;
    pos = std::min<int>(pos, smooth_data.size() - maxcnd);
    std::vector<pii> cnds;

    for (int i = pos; i < smooth_data.size() && i < pos + maxcnd; ++i) {
      int order_cost = 0;
      for (auto &e : smooth_data[i].second) order_cost += e.second;
      cnds.emplace_back(-order_cost, i);
    }

    std::sort(ALL(cnds));
    for (auto &e : cnds) {
      bound /= smooth_data[e.second].first;
      res.push_back(smooth_data[e.second]);
      if (bound == 0) break;
    }
    OPA_CHECK0(bound == 0);
    return res;
  }

  void execute(ExecutionData *ex_data) const {
    main_executor->execute(ex_data);
  }

  std::vector<T> fft(const std::vector<T> &a) const {
    ExecutionData ed;
    ed.input = a;

    ed.input.resize(m_order, m_f->getZ());
    main_executor->execute(&ed);
    return ed.output;
  }

  std::vector<T> ifft(const std::vector<T> &a) const {
    ExecutionData ed;
    ed.input = a;
    ed.inverse = true;

    ed.input.resize(m_order, m_f->getZ());
    main_executor->execute(&ed);

    T iv = m_f->inv(m_f->importu32(m_order));
    REP (i, ed.output.size()) { ed.output[i] = m_f->mul(ed.output[i], iv); }
    return ed.output;
  }

  bool m_force_simple;
  const Field<T> *m_f;
  T m_root;
  u32 m_order;
  BGFactors m_factors;
  FFTExecutor *main_executor;
  utils::ObjContainer<FFTExecutor> executors;
  utils::ObjContainer<Field<T> > fields;
};

template <class T> T resultant(const Poly<T> &a, const Poly<T> &b) {
  int m = a.deg();
  int n = b.deg();
  Matrix<T> res(a.get_underlying_ring(), n + m, n + m);
  REP (i, n)
    REP (j, m + 1)
      res(i, i + j) = a[m - j];
  REP (i, m)
    REP (j, n + 1)
      res(n + i, i + j) = b[n - j];
  OPA_DISP0(res);
  return res.get_det_row_echelon();
}

template <class T>
Poly<T> make_squarefree(const Poly<T> &poly, Poly<T> *rem = nullptr) {
  auto pr = poly.get_poly_ring();
  Poly<T> f_dot = pr->derivative(poly);
  Poly<T> common = poly.get_ring()->gcd(poly, f_dot);
  Poly<T> res = poly / common;
  OPA_DISP0(poly, f_dot, common, common * f_dot, res);
  if (rem != nullptr) *rem = common;
  return res;
}

template <class T> bool is_squarefree(const Poly<T> &poly) {
  auto pr = poly.get_poly_ring();
  OPA_TRACE("Is squarefree ", poly);
  Poly<T> f_dot = pr->derivative(poly);
  OPA_TRACE("Is squarefree ", poly, f_dot);
  Poly<T> common = poly.get_ring()->gcd(poly, f_dot);
  OPA_TRACE("Done squarefree", poly, f_dot, common);
  return common.deg() == 0;
}

template <typename PR, typename U>
Poly<Poly<U> > import_change_var(const PR &pr, const Poly<U> &p) {
  Poly<Poly<U> > res = pr.get_poly();
  for (auto &e : p) {
    res.vec_unsafe().push_back(
      ((const PolyRingOps<U> *)(pr.underlying_ring()->get_poly_ring()))
        ->import(std::vector<U>{ e }));
  }
  return res;
}

template <typename R, typename T>
std::vector<typename R::Type> import_vec(const R &r, const std::vector<T> &a) {
  std::vector<typename R::Type> res(a.size());
  REP (i, a.size())
    res[i] = r.import(a[i]);
  return res;
}

template <typename U, typename T>
Poly<U> import_poly(const PolyRing<U> &pr, const Poly<T> &a) {
  Poly<U> conv = pr.get_poly();
  for (auto &x : a) conv.vec_unsafe().emplace_back(x);
  return conv;
}

template <typename U, typename T>
Poly<U> import_poly(const PolyRingBaseField<U> &pr, const Poly<T> &a) {
  Poly<U> conv = pr.get_poly();
  for (auto &x : a) conv.vec_unsafe().emplace_back(x);
  return conv;
}

template <typename U>
Poly<Poly<U> > poly_switch_vars(const Poly<Poly<U> > &src) {
  const PolyRingOps<Poly<U> > &pr_xy = *src.get_poly_ring();
  const Ring<Poly<U> > &pr_x_ring = *src.get_underlying_ring();
  const PolyRingOps<U> &pr_x =
    *(const PolyRingOps<U> *)pr_x_ring.get_poly_ring();

  Poly<Poly<U> > res = pr_xy.get_poly();
  REP (ydeg, src.size()) {
    const auto &cur_x = src[ydeg];
    REP (xdeg, cur_x.size()) {
      auto &res_x = res.get_force(xdeg);
      res_x = res_x + pr_x.xpwv(ydeg, cur_x[xdeg]);
    }
  }

  res.normalize();
  return res;
}

template <typename BaseU, typename FracU>
BaseU force_poly_to_frac_base(const Poly<FracU> &frac_poly,
                              Poly<BaseU> *base_poly) {
  std::vector<BaseU> denoms;
  for (auto &x : frac_poly) denoms.push_back(x.q);
  const auto &base_ring = *base_poly->get_underlying_ring();
  BaseU d = lcm_list(base_ring, denoms);
  for (auto &x : frac_poly)
    base_poly->vec_unsafe().push_back(
      base_ring.mul(x.p, base_ring.div(d, x.q)));
  return d;
}

template <typename T> T resultant2(const Poly<T> &ia, const Poly<T> &ib) {
  const auto &ring = *ia.get_underlying_ring();
  if (ia.deg() < ib.deg()) {
    T res = resultant2(ib, ia);
    if ((ib.deg() & 1) && (ia.deg() & 1)) res = ring.neg(res);
    return res;
  }

  Poly<T> a, b;
  T ca, cb;
  ia.pp_and_cont(&a, &ca);
  ib.pp_and_cont(&b, &cb);
  T g = ring.getE();
  T h = ring.getE();
  T s = ring.getE();
  T t = ring.mul(ring.faste(ca, b.deg()), ring.faste(cb, a.deg()));

  while (1) {
    int alpha = a.deg() - b.deg();
    if ((a.deg() & 1) && (b.deg() & 1)) s = ring.neg(s);
    Poly<T> a2 = a * ring.faste(b.lc(), alpha + 1);
    Poly<T> r = a2 % b;
    a = b;
    b = r / ring.mul(g, ring.faste(h, alpha));

    g = a.lc();
    h = ring.div(ring.faste(g, alpha), ring.faste(h, alpha - 1));
    if (b.deg() <= 0) break;
  }
  h = ring.div(ring.faste(b.lc(), a.deg()), ring.faste(h, a.deg() - 1));
  return ring.mul(s, ring.mul(t, h));
}

template <class T> T resultant_slow(const Poly<T> &a, const Poly<T> &b) {
  int m = a.deg();
  int n = b.deg();
  Matrix<T> res(a.get_underlying_ring(), n + m, n + m);
  REP (i, n)
    REP (j, m + 1)
      res(i, i + j) = a[m - j];
  REP (i, m)
    REP (j, n + 1)
      res(n + i, i + j) = b[n - j];
  return res.get_det_slow();
}

OPA_NAMESPACE_DECL3_END

#endif

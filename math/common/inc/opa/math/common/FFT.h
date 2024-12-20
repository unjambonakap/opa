#pragma once

#include <opa/math/common/Field.h>
#include <opa/math/common/RealField.h>

OPA_NM_MATH_COMMON

template <class T> class IFFTProvider {

protected:
  void configure0(int npw, const Field<T> *f) {
    m_npw = npw;
    m_n = 1 << npw;
    m_f0 = f;
  }

public:
  int m_npw;
  int m_n;
  const Field<T> *m_f0;

  virtual std::vector<T> fft(const std::vector<T> &tb) const = 0;
  virtual std::vector<T> ifft(const std::vector<T> &tb, bool normalize = true) const = 0;

  virtual std::vector<T> mul(const std::vector<T> &a, const std::vector<T> &b) const {

    OPA_CHECK0(a.size() + b.size() - 1 <= m_n);
    auto ia = fft(a);
    auto ib = fft(b);
    REP (i, m_n) ia[i] = m_f0->mul(ia[i], ib[i]);
    auto res = ifft(ia);
    return res;
  }
};

template <class T, class F> class FFTProviderF : public IFFTProvider<T> {

protected:
  void configure(int npw, const F *f) {
    m_f = f;
    this->configure0(npw, f);
  }

public:
  const F *m_f;
};

template <class T> class FFTProvider : public FFTProviderF<T, Field<T> > {};

template <class T, class F> class FFT2F : public FFTProviderF<T, F> {
  typedef std::vector<T> vec;
  std::vector<T> w, dw, y, dy;

public:
  FFT2F() {}
  FFT2F(const F *f, int npw) { this->init(f, npw); }

  void init(const F *f, int npw) { this->init(f, npw, f->getNthRoot(1 << npw)); }

  void init(const F *f, int npw, const T &nth_root) {
    this->configure(npw, f);
    m_w = nth_root;
    m_iw = this->m_f->inv(m_w);
    m_normalizer = f->inv(f->importu32(this->m_n));

    m_iwl.resize(this->m_npw + 1);
    m_wl.resize(this->m_npw + 1);
    m_wl[0] = m_w;
    m_iwl[0] = m_iw;
    REP (i, this->m_npw) {
      m_wl[i + 1] = f->mul(m_wl[i], m_wl[i]);
      m_iwl[i + 1] = f->mul(m_iwl[i], m_iwl[i]);
    }

    m_order.resize(this->m_n);
    REP (i, this->m_n) {
      int x = 0;
      REP (j, this->m_npw) x |= (i >> j & 1) << (this->m_npw - 1 - j);
      m_order[i] = x;
    }
  }

  std::vector<T> fft(const std::vector<T> &tb) const override {
    OPA_CHECK0(tb.size() <= this->m_n);
    std::vector<T> res;
    _do_fft(m_wl, tb, res);
    return res;
  }

  std::vector<T> ifft(const std::vector<T> &tb, bool normalize = true) const override {
    OPA_CHECK0(tb.size() <= this->m_n);
    std::vector<T> res;
    _do_fft(m_iwl, tb, res);
    if (normalize) {
      REP (i, this->m_n) res[i] = this->m_f->mul(res[i], m_normalizer);
    }
    return res;
  }

private:
  void _do_fft(const std::vector<T> &wl, const std::vector<T> &tb, std::vector<T> &res) const {
    res.resize(this->m_n);
    REP (i, this->m_n) res[m_order[i]] = i >= tb.size() ? this->m_f->getZ() : tb[i];
    REP (k, this->m_npw) {
      for (int j = 0; j < this->m_n; j += 1 << (k + 1)) {
        int step = 1 << k;
        const T &curpw = wl[this->m_npw - 1 - k];
        T tmp = this->m_f->getE();

        REP (i, step) {
          T b = this->m_f->mul(tmp, res[j + i + step]);
          tmp = this->m_f->mul(tmp, curpw);
          res[j + i + step] = this->m_f->sub(res[j + i], b);
          res[j + i] = this->m_f->add(res[j + i], b);
        }
      }
    }
  }

  T m_w, m_iw;
  T m_normalizer;

  std::vector<int> m_order;
  std::vector<T> m_wl;
  std::vector<T> m_iwl;
};
template <class T> class FFT2 : public FFT2F<T, Field<T> > {

public:
  FFT2() {}
  FFT2(const Field<T> *f, int npw) : FFT2F<T, Field<T> >(f, npw) {}
};

template <class T> class FFT2Dispatcher : public IFFTProvider<T> {
  std::vector<std::unique_ptr<IFFTProvider<T> > > _ffts;

public:
  typedef std::function<IFFTProvider<T> *(int)> MakeProvider;
  FFT2Dispatcher(const Field<T> *f, int npw) {
    FOR (i, 4, npw + 1) {
      _ffts.emplace_back(new FFT2<T>(f, i));
    }
  }

  FFT2Dispatcher(int npw, const MakeProvider &provider) {
    FOR (i, 4, npw + 1) {
      _ffts.emplace_back(provider(i));
    }
  }
  template <class Provider> FFT2Dispatcher(int npw, const Provider &provider) {
    FOR (i, 4, npw + 1) {
      _ffts.emplace_back(provider(i));
    }
  }

  const IFFTProvider<T> &sel(int a, int b = 1) const {

    int s = std::max(1, a) + std::max(1, b) - 2;
    int sz = s < 0 ? 0 : log2_high_bit(s) + 1;
    for (auto &e : _ffts)
      if (e->m_npw >= sz) return *e;
    OPA_CHECK0(false);
  }

  std::vector<T> fft(const std::vector<T> &tb) const override { return sel(tb.size()).fft(tb); }

  std::vector<T> ifft(const std::vector<T> &tb, bool normalize = true) const override {
    return sel(tb.size()).ifft(tb, normalize);
  }

  std::vector<T> mul(const std::vector<T> &a, const std::vector<T> &b) const override {
    return sel(a.size(), b.size()).mul(a, b);
  }
};

template <class T> class FFT_RealF : public IFFTProvider<T> {
public:
  typedef std::complex<T> CT;
  FFT2<CT> _fft;
  FField<CT, T> ff;
  FFT_RealF(const Field<T> *f, int npw) { _fft.init(&ff, npw); }

  virtual std::vector<T> fft(const std::vector<T> &tb) const override { OPA_CHECK0(false); }
  virtual std::vector<T> ifft(const std::vector<T> &tb, bool normalize = true) const override {
    OPA_CHECK0(false);
  }

  std::vector<CT> lift(const std::vector<T> &a) const {
    std::vector<CT> res;
    for (auto x : a) res.pb(CT(x));
    return res;
  }
  std::vector<T> lower(const std::vector<CT> &a) const {
    std::vector<T> res;
    for (auto x : a) res.pb(std::real(x));
    return res;
  }

  virtual std::vector<T> mul(const std::vector<T> &a, const std::vector<T> &b) const {
    return this->lower(_fft.mul(this->lift(a), this->lift(b)));
  }

  static std::unique_ptr<FFT2Dispatcher<T> > make_dispatcher(const Field<T> *f, int npw) {
    return std::make_unique<FFT2Dispatcher<T> >(
      npw, [&](int req) { return (FFTProvider<T> *)new FFT_RealF(f, req); });
  }
};

template <class T> class FFT_Real {
public:
  typedef std::complex<T> CT;
  void init(int n, const T &thresh = 1e-9) {
    this->n = n;
    this->pw = log2_high_bit(n - 1) + 1;
    int nn = 1 << pw;
    this->thresh = thresh;

    fft_v.init(&cf, pw, FloatUtil::unity_root<T>(nn));
  }

  std::vector<CT> fft(const std::vector<T> &a) const {
    std::vector<CT> aa;
    REP (i, a.size()) aa.emplace_back(a[i]);
    return fft_v.fft(aa);
  }

  std::vector<T> ifft(const std::vector<CT> &a) const {
    auto tmp = fft_v.ifft(a);
    std::vector<T> res;
    REP (i, tmp.size()) res.push_back(std::real(tmp[i]));
    while (res.size() > 0 && std::abs(res.back()) < thresh) res.pop_back();
    return res;
  }

  ComplexField<T> cf;
  FFT2<CT> fft_v;
  T thresh;

  int n;
  int pw;
};

template <class T> class FFT {
  typedef std::vector<T> vec;
  int m_n;
  const Field<T> *m_field;
  T m_w, m_iw;

  void _fft(const std::vector<T> &in, std::vector<T> &out, const T &w) const {
    OPA_CHECK0(in.size() <= m_n);

    T v = m_field->getE();
    int mx = std::min<int>(in.size(), m_n);
    REP (i, m_n) {
      T res = m_field->getZ();
      REPV (j, mx) res = m_field->add(in[j], m_field->mul(res, v));

      v = m_field->mul(v, w);
      out.push_back(res);
    }
  }

public:
  void init(const Field<T> *f, int n, const T &nthRoot) {
    m_field = f;
    m_n = n;
    m_w = nthRoot;
    m_iw = m_field->inv(m_w);
  }

  std::vector<T> fft(const std::vector<T> &tb) const {
    std::vector<T> ans;
    _fft(tb, ans, m_w);
    return ans;
  }

  std::vector<T> ifft(const std::vector<T> &tb, bool norm = true) const {
    std::vector<T> ans;
    T in = m_field->inv(m_field->import(m_n));

    _fft(tb, ans, m_iw);
    if (norm) {
      for (int i = 0; i < m_n; ++i) ans[i] = m_field->mul(ans[i], in);
    }
    return ans;
  }
};

OPA_NM_MATH_COMMON_END

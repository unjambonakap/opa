#pragma once

#include <opa/math/common/Field.h>
#include <opa/math/common/RealField.h>

OPA_NAMESPACE_DECL3(opa, math, common)

template <class T> class FFT2 {
  typedef std::vector<T> vec;

public:
  void init(const Field<T> *f, int npw, const T &nth_root) {
    m_f = f;
    m_npw = npw;
    m_n = 1 << npw;
    m_w = nth_root;
    m_iw = m_f->inv(m_w);
    m_normalizer = f->inv(f->importu32(m_n));

    m_iwl.resize(m_npw + 1);
    m_wl.resize(m_npw + 1);
    m_wl[0] = m_w;
    m_iwl[0] = m_iw;
    REP (i, m_npw) {
      m_wl[i + 1] = f->mul(m_wl[i], m_wl[i]);
      m_iwl[i + 1] = f->mul(m_iwl[i], m_iwl[i]);
    }

    m_order.resize(m_n);
    REP (i, m_n) {
      int x = 0;
      REP (j, m_npw)
        x |= (i >> j & 1) << (m_npw - 1 - j);
      m_order[i] = x;
    }
  }

  std::vector<T> fft(const std::vector<T> &tb) const {
    OPA_CHECK0(tb.size() <= m_n);
    std::vector<T> res;
    _do_fft(m_wl, tb, res);
    return res;
  }

  std::vector<T> ifft(const std::vector<T> &tb, bool normalize = true) const {
    OPA_CHECK0(tb.size() <= m_n);
    std::vector<T> res;
    _do_fft(m_iwl, tb, res);
    if (normalize) {
      REP (i, m_n)
        res[i] = m_f->mul(res[i], m_normalizer);
    }
    return res;
  }

private:
  void _do_fft(const std::vector<T> &wl, const std::vector<T> &tb,
               std::vector<T> &res) const {
    res.resize(m_n);
    REP (i, m_n)
      res[m_order[i]] = i >= tb.size() ? m_f->getZ() : tb[i];
    REP (k, m_npw) {
      for (int j = 0; j < m_n; j += 1 << (k + 1)) {
        int step = 1 << k;
        const T &curpw = wl[m_npw - 1 - k];
        T tmp = m_f->getE();

        REP (i, step) {
          T b = m_f->mul(tmp, res[j + i + step]);
          tmp = m_f->mul(tmp, curpw);
          res[j + i + step] = m_f->sub(res[j + i], b);
          res[j + i] = m_f->add(res[j + i], b);
        }
      }
    }
  }

  T m_w, m_iw;
  T m_normalizer;
  const Field<T> *m_f;
  int m_npw;
  int m_n;

  std::vector<int> m_order;
  std::vector<T> m_wl;
  std::vector<T> m_iwl;
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
    REP (i, a.size())
      aa.emplace_back(a[i]);
    return fft_v.fft(aa);
  }

  std::vector<T> ifft(const std::vector<CT> &a) const {
    auto tmp = fft_v.ifft(a);
    std::vector<T> res;
    REP (i, tmp.size())
      res.push_back(std::real(tmp[i]));
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
      REPV (j, mx)
        res = m_field->add(in[j], m_field->mul(res, v));

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

OPA_NAMESPACE_DECL3_END

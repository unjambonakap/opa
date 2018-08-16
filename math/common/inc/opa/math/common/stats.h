#pragma once

#include <opa/math/common/bignum.h>
#include <opa/math/common/float.h>

OPA_NAMESPACE_DECL3(opa, math, common)

bignum nchoosek(int n, int k);
Float erf_binomial(const Float &prob, const Float &nsamples,
                   const Float &obs_below);
Float erf_binomial_small(const Float &prob, int nsamples, int obs_below);
Float erf_at_conf(const Float &prob, const Float &conf);

template <class T, int N> class DataWhitener {
public:
  typedef DataWhitener<T, N> SelfType;

  SelfType &add(const T &v) {
    m_data.push_back(v);
    return *this;
  }
  SelfType &add(const std::vector<T> &v) {
    for (auto &x : v) this->add(x);
    return *this;
  }

  void compute(bool fixed_dir = true) {
    m_mean = T(0);
    for (auto &x : m_data) m_mean += x;
    m_mean /= m_data.size();
    // two phases for numerical precision
    T temp_sum = T(0);
    for (auto &x : m_data) temp_sum += (x - m_mean) * (x - m_mean);

    if (fixed_dir) {
      double length = 0;
      REP (i, N)
        length += temp_sum[i];
      m_std = T(sqrt(length / (N * m_data.size())));

    } else {
      REP (i, N)
        m_std[i] = std::sqrt(temp_sum[i] / m_data.size());
    }
    m_data.clear();
  }

  T remap(const T &v) const { return (v - m_mean) / m_std; }
  T unmap(const T &v) const { return v * m_std + m_mean; }

  std::vector<T> remap(const std::vector<T> &v) const {
    std::vector<T> res;
    res.reserve(v.size());
    for (auto &x : v) res.push_back(remap(x));
    return res;
  }

  std::vector<T> unmap(const std::vector<T> &v) const {
    std::vector<T> res;
    res.reserve(v.size());
    for (auto &x : v) res.push_back(unmap(x));
    return res;
  }

private:
  std::vector<T> m_data;
  T m_mean = T(0.);
  T m_std = T(0.);
};

OPA_NAMESPACE_DECL3_END

#pragma once

#include <opa/utils/base.h>

namespace opa {
namespace utils {

template <class T>
std::vector<std::vector<T> > product(const std::vector<std::vector<T> > &tb, int k = -1) {
  if (k == -1) k = tb.size();
  if (k == 0) return { {} };
  std::vector<std::vector<T> > res;
  for (auto e : product(tb, k - 1)) {
    for (auto x : tb[k - 1]) {
      res.pb(e);
      res.back().pb(x);
    }
  }
  return res;
}

template <class T> std::vector<std::vector<T> > product2(const std::vector<T> &tb, int k) {
  if (k == 0) return { {} };
  std::vector<std::vector<T> > res;
  for (auto e : product2(tb, k - 1)) {
    for (auto x : tb) {
      res.pb(e);
      res.back().pb(x);
    }
  }

  return res;
}

template <class T> std::vector<std::vector<T> > combinations2(const std::vector<T> &tb, int K) {
  int n = tb.size();
  if (K > n) return {};
  if (K == 0) {

    return { {} };
  } else {

    if (K == tb.size()) {
      return { tb };
    }

    std::vector<std::vector<T> > res;
    auto copy = tb;
    REP (i, n - K + 1) {
      auto last = copy.back();
      copy.pop_back();

      for (auto x : combinations2<T>(copy, K - 1)) {
        res.pb(x);
        res.back().pb(last);
      }
    }
    return res;
  }
}
template <class T, int K> std::vector<std::array<T, K> > combinations(const std::vector<T> &tb) {
  int n = tb.size();
  if (K > n) return {};
  if constexpr (K == 0) {

    return { {} };
  } else {

    if (K == tb.size()) {
      std::array<T, K> u;
      REP (i, K) u[i] = tb[i];
      return { u };
    }

    std::vector<std::array<T, K> > res;
    auto copy = tb;
    REP (i, n - K + 1) {
      auto last = copy.back();
      copy.pop_back();

      for (auto x : combinations<T, K - 1>(copy)) {
        res.emplace_back();
        REP (i, K - 1) res.back()[i] = x[i];
        res.back()[K - 1] = last;
      }
    }
    return res;
  }
}

} // namespace utils
} // namespace opa

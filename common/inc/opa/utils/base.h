#pragma once

#include <glib/gtl/map_util.h>
#include <glib/hash/hash.h>
#include <opa/stolen/StringRef.h>
#include <opa/utils/misc.h>
#include <opa/utils/serialize.h>
#include <opa/utils/string.h>

OPA_NAMESPACE(opa, utils)
constexpr double misc_base_eps = 1e-9;

template <typename T> class Range {

public:
  static Range<T> Build_range(const T &from, const T &to, int n) {
    Range<T> res;
    res.build_range(from, to, n);
    return res;
  }
  static Range<T> StepRange(const T &from, const T &to, int step) {
    Range<T> res;
    for (T x = from; x < to; x += step) res.m_tb.push_back(x);
    return res;
  }

  void build_range(const T &from, const T &to, int n) {
    T step = (to - from) / (n - 1);
    T cur = from;

    REP (i, n) {
      m_tb.pb(cur);
      cur += step;
    }
  }
  OPA_ACCESSOR_R(std::vector<T>, m_tb, tb);

private:
  std::vector<T> m_tb;
};

static Range<int> range(int a, int b) { return Range<int>::StepRange(a, b, 1); }

struct Range1D {
  double low, high;

  double get_length() const { return high - low; }
  double get_middle() const { return (high + low) / 2; }
  double is_in(double p, double eps = misc_base_eps) const {
    return p > low - eps && p < high + eps;
  }
  double is_in_fast(double p) const { return p >= low && p <= high; }

  double dist(double p) const {
    if (is_in(p)) return 0.;
    return p < low ? low - p : p - high;
  }
  bool is_almost_empty(const double eps = misc_base_eps) const {
    return high - low < eps;
  }

  Range<double> linspace(int n) const {
    return Range<double>::Build_range(low, high, n);
  }

  static Range1D MakeRange(double a, double b) {
    if (a > b) std::swap(a, b);
    return Range1D{ a, b };
  }
};
template <class T> void sort_pair(std::pair<T, T> &a) {
  if (a.second < a.first) std::swap(a.first, a.second);
}

template <
  template <typename T, typename... Ts> class Container, typename Functor,
  typename T, // <-- This is the one we'll override in the return container
  typename U = typename std::result_of<Functor(T)>::type, typename... Ts>
Container<U> transform_container(const Container<T, Ts...> &c, Functor &&f) {
  Container<U> ret;
  std::transform(std::begin(c), std::end(c), std::inserter(ret, std::end(ret)),
                 f);
  return ret;
}

template <
  template <typename T, typename... Ts> class Container, typename Functor,
  typename T, // <-- This is the one we'll override in the return container
  typename U = typename std::result_of<Functor(T, T)>::type, typename... Ts>
Container<U, Ts...> transform_container_binary(const Container<T, Ts...> &a,
                                               const Container<T, Ts...> &b,
                                               Functor &&f) {
  Container<U, Ts...> ret;
  REP (i, a.size()) { ret.push_back(f(a[i], b[i])); }
  return ret;
}

template <template <typename T, typename... Ts> class Container, typename T,
          typename... Ts>
T get_or(const Container<T, Ts...> &container, int pos, const T &other) {
  if (container.size() <= pos) return other;
  return container[pos];
}

template <class Container, class T = typename Container::mapped_type>
std::vector<T> get_values(const Container &container) {
  std::vector<T> res;
  for (auto &x : container) res.push_back(x.second);
  return res;
}

template <class Container, class T = typename Container::value_type>
bool contained_sorted(const Container &source, const Container &ref) {
  for (auto it = source.begin(), it2 = ref.begin(); it != source.end(); ++it) {
    while (it2 != ref.end() && *it2 < *it) ++it2;
    if (it2 == ref.end() || *it2 != *it) return false;
    ++it2;
  }
  return true;
}

template <class Container, class T = typename Container::value_type>
std::vector<T> remove_sorted(const Container &source, const Container &bad) {
  std::vector<T> res;
  for (auto it = source.begin(), it2 = bad.begin(); it != source.end(); ++it) {
    while (it2 != bad.end() && *it2 < *it) ++it2;
    if (it2 != bad.end() && *it2 == *it) {
      ++it2;
      continue;
    }
    res.push_back(*it);
  }
  return res;
}

template <class T> class ScopedPush {
public:
  ScopedPush(T &target, const T &nv) : target(target) {
    old = target;
    target = nv;
  }
  ~ScopedPush() { target = old; }

  T old;
  T &target;
};

OPA_NAMESPACE_END(opa, utils)

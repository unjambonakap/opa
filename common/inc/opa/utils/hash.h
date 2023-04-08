
#pragma once

#include <absl/hash/hash.h>
#include <opa/stolen/StringRef.h>
#include <opa/utils/serialize.h>
#include <opa/utils/string.h>

namespace std {
template <class A, class B> struct hash<std::pair<A, B> > {
  size_t operator()(const std::pair<A, B> &v) const {
    return absl::HashOf(v.first, v.second);
  }
};
template <class T> struct hash<std::vector<T> > {
  size_t operator()(const std::vector<T> &v) const {
    size_t curhash = 0;
    for (auto &x : v) {
      curhash = absl::HashOf(curhash, std::hash<T>()(x));
    }
    return curhash;
  }
};
template <class T, int N> struct hash<std::array<T, N> > {
  size_t operator()(const std::array<T, N> &v) const {
    size_t curhash = 0;
    for (auto &x : v) {
      curhash = absl::HashOf(curhash, std::hash<T>()(x));
    }
    return curhash;
  }
};
}

OPA_NAMESPACE(opa, utils)
namespace hash_utils {
// Recursive case (hash up to Index)
template <class Tuple, size_t Index = std::tuple_size<Tuple>::value - 1>
struct _hash_helper {
  static size_t hash(Tuple const &tuple) {
    size_t prev = _hash_helper<Tuple, Index - 1>::hash(tuple);
    using TypeForIndex = typename std::tuple_element<Index, Tuple>::type;
    size_t thisHash = std::hash<TypeForIndex>()(std::get<Index>(tuple));
    return absl::HashOf(prev, thisHash);
  }
};

// Base case (hash 0th element)
template <class Tuple> struct _hash_helper<Tuple, 0> {
  static size_t hash(Tuple const &tuple) {
    using TypeForIndex = typename std::tuple_element<0, Tuple>::type;
    return std::hash<TypeForIndex>()(std::get<0>(tuple));
  }
};

// Recursive case (elements equal up to Index)
template <class Tuple, size_t Index = std::tuple_size<Tuple>::value - 1>
struct _eq_helper {
  static bool equal(Tuple const &a, Tuple const &b) {
    auto aValue = std::get<Index>(a);
    auto bValue = std::get<Index>(b);
    using TypeForIndex = typename std::tuple_element<Index, Tuple>::type;
    if (!std::equal_to<TypeForIndex>(aValue, bValue)) return false;
    bool prev = _eq_helper<Tuple, Index - 1>::equal(a, b);
    return prev;
  }
};

// Base case (0th elements equal)
template <class Tuple> struct _eq_helper<Tuple, 0> {
  static bool equal(Tuple const &a, Tuple const &b) {
    using TypeForIndex = typename std::tuple_element<0, Tuple>::type;
    auto &aValue = std::get<0>(a);
    auto &bValue = std::get<0>(b);
    return std::equal_to<TypeForIndex>()(aValue, bValue);
  }
};

template <typename... TT> struct hash {
  size_t operator()(std::tuple<TT...> const &tt) const {
    return _hash_helper<std::tuple<TT...> >::hash(tt);
  }
};

template <typename... TT> struct equal_to {
  bool operator()(std::tuple<TT...> const &a,
                  std::tuple<TT...> const &b) const {
    return _eq_helper<std::tuple<TT...> >::equal(a, b);
  }
};
}

template <typename Subclass, typename ReturnValue, typename Context,
          typename... Args>
class MemoizeHelper {
public:
  typedef std::tuple<Args...> TupleInput;

  void init(const Context *context) { this->context = context; }

  template <typename FuncType>
  const ReturnValue &get(const FuncType &func, const Args &... args) const {
    TupleInput input(args...);
    OPA_CHECK0(context != nullptr);
    auto res = store.emplace(input, nullptr);
    if (res.second) {
      retv_pool.emplace_back();
      ReturnValue *cur = &retv_pool.back();
      res.first->second = cur;
      *cur =
        ((const Subclass *)this)->template real_get<FuncType>(func, args...);
      return *cur;
    }
    return *res.first->second;
  }

  template <typename FuncType>
  ReturnValue real_get(const FuncType &func, const Args &... args) const;

  const Context *context = nullptr;

  mutable std::deque<ReturnValue> retv_pool;
  mutable std::unordered_map<TupleInput, ReturnValue *,
                             ::opa::utils::hash_utils::hash<Args...> >
    store;
};

template <typename ReturnValue, typename Context, typename... Args>
class MemoizeHelperFunc
  : public MemoizeHelper<MemoizeHelperFunc<ReturnValue, Context, Args...>,
                         ReturnValue, Context, Args...> {
public:
  template <typename FuncType>
  ReturnValue real_get(const FuncType &func, const Args &... args) const {
    return func(*this->context, args...);
  }
};

template <typename ReturnValue, typename Context, typename... Args>
class MemoizeHelperClass
  : public MemoizeHelper<MemoizeHelperClass<ReturnValue, Context, Args...>,
                         ReturnValue, Context, Args...> {
public:
  template <typename Method>
  ReturnValue real_get(const Method &method, const Args &... args) const {
    return (this->context->*method)(args...);
  }
};

#define CALL_mem(self_obj_type, st, func, ...)                                 \
  st.get(&self_obj_type::func, __VA_ARGS__)
#define CALL(self_obj_type, st, func, ...)                                     \
  CALL_mem(self_obj_type, st, func, __VA_ARGS__)
#define CALL1(self_obj_type, name, ...)                                        \
  CALL(self_obj_type, name##_mem, name##_, __VA_ARGS__)
#define CALL1_OUTSIDE(obj, name, ...)                                          \
  (obj).CALL(&decltype(obj), objname##_mem, name##_, __VA_ARGS__)

OPA_NAMESPACE_END(opa, utils)

namespace std {
template <typename... TT> struct hash<std::tuple<TT...> > {
  size_t operator()(std::tuple<TT...> const &tt) const {
    return opa::utils::hash_utils::_hash_helper<std::tuple<TT...> >::hash(tt);
  }
};
}

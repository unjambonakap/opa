#pragma once

#include <opa/math/common/base.h>

OPA_NM_MATH_COMMON

template <class T> class Operations {
public:
  typedef std::function<T(const T &)> import_func_t;
  struct OperationsParams {
    import_func_t import_func;
  } OperationsParams m_params;

  Operations(const OperationsParams &params) : m_params(params) {}
#define BINARY_OP                                                              \
  (a, b, op) {                                                                 \
    T _a = import(a);                                                          \
    T _b = import(b);                                                          \
    return import(_a op _b);                                                   \
  }

  T add(const T &a, const T &b) const { BINARY_OP(a, b, +); }
  T sub(const T &a, const T &b) const { BINARY_OP(a, b, -); }
  T mul(const T &a, const T &b) const { BINARY_OP(a, b, *); }
  T neg(const T &a) const { return import(-a); }
#undef BINARY_OP

private:
  T import(const T &a) const { return m_params.import_func(a); }
};

OPA_NM_MATH_COMMON_END

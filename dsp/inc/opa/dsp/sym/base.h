#pragma once

#include <opa_common.h>
#include <opa/dsp/inc.h>

OPA_NAMESPACE(opa, dsp)

template <class R, class... Args>
class Function : public opa::utils::Initable {
 public:
  virtual R eval(const Args &... args) const = 0;
  R operator()(const Args &... args) const { return eval(args...); }

  Function<R, Args...> *base() { return this; }
};

template <class R, class... Args>
using FunctionSptr = std::shared_ptr<Function<R, Args...>>;

template <class R, class... Args>
class FuncFromFunc : public Function<R, Args...> {
 public:
  typedef std::function<R(const Args &... args)> FuncType;
  void init(FuncType func) { m_func = func; }
  virtual R eval(const Args &... args) const override {
    return m_func(args...);
  }

  static FunctionSptr<R, Args...> Make(FuncType func) {
    auto x = new FuncFromFunc<R, Args...>();
    x->init(func);
    return FunctionSptr<R, Args...>(x->base());
  }

  FuncType m_func;
};

template <class R, class... Args>
class ResultMixer : public Function<R, Args...> {
 public:
  struct Params {
    FunctionSptr<R, Args...> f1, f2;
    FunctionSptr<R, R, R> mixer;
  };

  void init(const Params &params) { m_params = params; }
  virtual R eval(const Args &... args) const override {
    return m_params.mixer->eval(m_params.f1->eval(args...),
                                m_params.f2->eval(args...));
  }

  static FunctionSptr<R, Args...> Make(const Params &params) {
    auto x = new ResultMixer<R, Args...>();
    x->init(params);
    return FunctionSptr<R, Args...>(x->base());
  }

  Params m_params;
};

template <class R, class... Args>
FunctionSptr<R, Args...> operator+(FunctionSptr<R, Args...> f1,
                                   FunctionSptr<R, Args...> f2) {
  return ResultMixer<R, Args...>::Make(
      {f1, f2, FuncFromFunc<R, R, R>::Make([](const R &a,
                                              const R &b) { return a + b; })});
}
template <class R, class... Args>
FunctionSptr<R, Args...> operator*(FunctionSptr<R, Args...> f1,
                                   FunctionSptr<R, Args...> f2) {
  return ResultMixer<R, Args...>::Make(
      {f1, f2, FuncFromFunc<R, R, R>::Make([](const R &a,
                                              const R &b) { return a * b; })});
}

template <class R, class A, class A2>
FunctionSptr<R, A> compose(FunctionSptr<R, A2> f1, FunctionSptr<A2, A> f2) {
  return FuncFromFunc<R, A>::Make([f1, f2](const A &a) {
    return f1->eval(f2->eval(a));
  });
}

template <typename R, typename A, typename B>
FunctionSptr<R, A> fix2(FunctionSptr<R, A, B> f1, const B &b) {
  return FuncFromFunc<R, A>::Make([f1, b](const A &a) {
    return f1->eval(a, b);
  });
}

template <typename R, typename A, typename B>
FunctionSptr<R, B> fix1(FunctionSptr<R, A, B> f1, const A &a) {
  return FuncFromFunc<R, B>::Make([f1, a](const B &b) {
    return f1->eval(a, b);
  });
}

template <typename R, typename A, typename B>
FunctionSptr<R, A, B> ignore2(FunctionSptr<R, A> f1) {
  return FuncFromFunc<R, A, B>::Make([f1](const A &a, const B &b) {
    return f1->eval(a);
    ;
  });
}
template <typename R, typename A, typename B>
FunctionSptr<R, A, B> ignore1(FunctionSptr<R, B> f1) {
  return FuncFromFunc<R, A, B>::Make([f1](const A &a, const B &b) {
    return f1->eval(b);
    ;
  });
}

template <typename R, typename... Args>
FunctionSptr<R, Args...> conj(FunctionSptr<R, Args...> f1) {
  return FuncFromFunc<R, Args...>::Make([f1](const Args &... args) {
    return std::conj(f1->eval(args...));
  });
}

OPA_NAMESPACE_END(opa, dsp)

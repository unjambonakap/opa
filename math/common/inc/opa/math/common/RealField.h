#pragma once

#include <opa/math/common/float.h>
#include <opa/math/common/Field.h>
#include <opa/math/common/FractionField.h>
#include <opa/math/common/bignum.h>

OPA_NAMESPACE_DECL3(opa, math, common)

template <class T, class U> class FField : public Field<T> {
public:
  FField(U eps = 0) : Field<T>(-1, -1) { this->eps = eps; }
  // interface
  void init(U eps) { this->eps = eps; }

  virtual T inv(const T &a) const override { return T(1) / a; }
  virtual T mul(const T &a, const T &b) const override { return a * b; }
  virtual T add(const T &a, const T &b) const override { return a + b; }
  virtual T neg(const T &a) const override { return -a; }

  virtual bool isZ(const T &a) const override { return eq(a, getZ()); }
  virtual bool isE(const T &a) const override { return eq(a, getE()); }
  virtual T import(const T &a) const override { return a; }
  virtual T importu32(u32 a) const override { return T((s32)a); }

  virtual T getZ() const override { return T(0); }
  virtual T getE() const override { return T(1); }
  virtual T getRandRaw() const override { assert(0); }
  virtual T getRand() const override { return FloatUtil::get_rand_uni<T>(T(0), T(1)); }
  virtual bool eq(const T &a, const T &b) const override {
    return std::abs(a - b) <= eps;
  }
  U eps;
};

template <class T> using RealField = FField<T, T>;
template <class T> using ComplexField = FField<std::complex<T>, T>;

typedef RealField<double> RealDouble;
typedef RealField<Float> RealF;
typedef FractionField<bignum> RealBG;
typedef ComplexField<Float> ComplexF;
typedef ComplexField<double> ComplexDouble;

OPA_NAMESPACE_DECL3_END

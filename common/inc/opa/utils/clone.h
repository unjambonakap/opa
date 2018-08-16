#pragma once

#include <opa_common.h>

OPA_NAMESPACE_DECL2(opa, utils)

template <typename T> struct Clonable {

  virtual T *clone() const = 0;

protected:
  Clonable() {}
  Clonable(const Clonable &) {}
  virtual Clonable &operator=(const Clonable &) { return *this; }
  virtual ~Clonable() {}
};

#define OPA_DECL_CLONE(basetyp, typ)                                           \
  virtual basetyp *clone() const override { return new typ(*this); }

template <typename T> struct Duplicatable {

  virtual T *duplicate() const = 0;

protected:
  Duplicatable() {}
  Duplicatable(const Duplicatable &) {}
  virtual Duplicatable &operator=(const Duplicatable &) { return *this; }
  virtual ~Duplicatable() {}
};

#define OPA_DECL_DUPLICATE(basetyp, typ)                                       \
  virtual basetyp *duplicate() const override { return new typ(); }

OPA_NAMESPACE_DECL2_END

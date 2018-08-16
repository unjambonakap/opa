#pragma once

#include <opa/math/common/rng.h>
#include <opa/math/common/Types.h>
#include <opa/utils/base.h>

#define OPA_MATH_CO                                                            \
  namespace opa {                                                              \
  namespace math {                                                             \
  namespace co {
#define OPA_MATH_CO_END                                                        \
  }                                                                            \
  }                                                                            \
  }

OPA_MATH_CO

using namespace opa::math::common;

OPA_MATH_CO_END

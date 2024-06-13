#include "common/Utils.h"
#include <common/GF_p.h>

OPA_NAMESPACE_DECL3(opa, math, common)

GF_p::GF_p(u32 n, S64Factors factors) : GF_pT(n, false) {
  m_factors = to_bgfactors(factors);
  init(n);
}

OPA_NAMESPACE_DECL3_END

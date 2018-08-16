#include "opa/utils/debug.h"

OPA_NAMESPACE_DECL2(opa, utils)

std::string DebugCtx::flush_and_str() {
  FEV (it, cur_ctx)
    (*it)->flush1();
  return OPA_STREAM_STR(node);
}

OPA_NAMESPACE_DECL2_END

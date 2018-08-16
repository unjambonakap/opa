#include <opa_common.h>
#include "jit.h"

using namespace asmjit;

OPA_NAMESPACE(opa, crypto, la)

uintptr_t FunctionCaller::call() {
  switch (args.size()) {
    OPA_CASE(0, ret = ((f0)(func_ptr))(););
    OPA_CASE(1, ret = ((f1)(func_ptr))(args[0]););
    OPA_CASE(2, ret = ((f2)(func_ptr))(args[0], args[1]););
    OPA_CASE(3, ret = ((f3)(func_ptr))(args[0], args[1], args[2]););
  default:
    OPA_CHECK0(0);
  }
  return ret;
}

asmjit::VarType JitUtil::sizeToType(int size) const {
  if (size <= 8)
    return kVarTypeUInt8;
  if (size <= 16)
    return kVarTypeUInt16;
  if (size <= 32)
    return kVarTypeUInt32;
  if (size <= 64)
    return kVarTypeUInt64;
  OPA_CHECK0(0);
}
OPA_NAMESPACE_END(opa, crypto, la)

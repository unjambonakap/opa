#include <opa_common.h>
#include <opa/crypto/la/context.h>
#include <opa/crypto/la/graph.h>
#include "jit.h"

OPA_NAMESPACE(opa, crypto, la)
CipherContext::~CipherContext() {}

CipherContext::CipherContext() { init(); }
void CipherContext::init() {
  m_jitutil.reset(new JitUtil());
  /*
  m_runtime.reset(new asmjit::JitRuntime());
  m_compiler.reset(new asmjit::X86Compiler(m_runtime.get()));
  */
}
void CipherContext::start_evaluate(CipherNode *node) const {

  if (!is_debug())
    return;
  append(opa::utils::SPrintf("Evaluate iv=%s: %s", node->iv.str().c_str(),
                             node->block->desc().c_str()));
  ++depth;
}

void CipherContext::end_evaluate(CipherNode *node) const {
  if (!is_debug())
    return;
  --depth;
  append(opa::utils::SPrintf("Result ov=%s", node->ov.str().c_str()));
}
opa::utils::Singleton<CipherContext> CipherContext::singleton;

OPA_NAMESPACE_END(opa, crypto, la)

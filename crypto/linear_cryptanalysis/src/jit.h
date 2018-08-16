#pragma once
#include <opa_common.h>
#define ASMJIT_BUILD_X64
#include <asmjit/asmjit.h>

#include <opa/crypto/la/base.h>

OPA_NAMESPACE(opa, crypto, la)

class FunctionCaller {
public:
  typedef uintptr_t (*f0)(void);
  typedef uintptr_t (*f1)(uintptr_t);
  typedef uintptr_t (*f2)(uintptr_t, uintptr_t);
  typedef uintptr_t (*f3)(uintptr_t, uintptr_t, uintptr_t);

  uintptr_t call();

  uintptr_t func_ptr;
  uintptr_t ret;
  std::vector<uintptr_t> args;
};

struct JitContext {
  asmjit::X86Compiler *c;
  std::vector<std::unique_ptr<asmjit::X86GpVar> > vars;
};

struct JitBuilder {
public:
  std::vector<asmjit::X86GpVar *> input_vars;
  asmjit::X86GpVar *key_var = nullptr;
  asmjit::X86GpVar *output_var = nullptr;
  asmjit::X86Compiler *c = nullptr;

  JitBuilder() {}
  JitBuilder(asmjit::X86Compiler *c) { this->c = c; }
  JitBuilder(asmjit::X86Compiler *c,
             const std::vector<asmjit::X86GpVar *> &input_vars,
             asmjit::X86GpVar *key_var, asmjit::X86GpVar *output_var) {
    this->c = c;
    this->input_vars = input_vars;
    this->key_var = key_var;
    this->output_var = output_var;
  }
};
class JitUtil {
public:
  asmjit::VarType sizeToType(int size) const;
};

OPA_NAMESPACE_END(opa, crypto, la)

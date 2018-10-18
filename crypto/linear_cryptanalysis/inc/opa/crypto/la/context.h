#pragma once

#include <opa_common.h>
#include <opa/utils/misc.h>
#include <opa/threading/dispatcher.h>
#include <opa/math/common/fast_gf2.h>

namespace asmjit {
class X86Compiler;
class JitRuntime;
}

OPA_NAMESPACE_DECL3(opa, crypto, la)
class CipherNode;
class CipherBlock;
class JitUtil;

class CipherContext : public opa::utils::Initable {
public:
  OPA_ACCESSOR(bool, m_is_debug, is_debug);
  virtual void init() override;
  CipherContext();
  virtual ~CipherContext();

  void reset() {
    dbg().clear();
    depth = 0;
  }

  void start_evaluate(CipherNode *node) const;

  void append(const std::string &msg) const { m_dbg += get_pad() + msg + "\n"; }

  void end_evaluate(CipherNode *node) const;
  OPA_ACCESSOR(std::string, m_dbg, dbg);
  OPA_ACCESSOR_R(int, m_rel_lim, rel_lim);
  asmjit::JitRuntime &runtime() const { return *m_runtime; }
  asmjit::X86Compiler *compiler() const { return m_compiler.get(); }
  OPA_ACCESSOR_R(JitUtil, *m_jitutil, jit);
  std::vector<u64> real_key;

  template <class T>
  std::unique_ptr<T> instanciate(bool from_classstore = false) const {
    T *res;
    auto cs = opa::utils::ClassStore::instance.get();
    if (!from_classstore && cs->has<T>())
      res = cs->get2<T>();
    else
      res = new T;
    res->set_context(this);
    return std::unique_ptr<T>(res);
  }
  opa::threading::Dispatcher *dispatcher = nullptr;

  opa::math::common::BitVec target_key;
  static opa::utils::Singleton<CipherContext> singleton;

private:
  std::string get_pad() const { return std::string(2 * depth, ' '); }

  UPTR(JitUtil) m_jitutil;
  mutable UPTR(asmjit::JitRuntime) m_runtime;
  mutable UPTR(asmjit::X86Compiler) m_compiler;

  int m_rel_lim = 15;
  mutable int depth = 0;
  mutable std::string m_dbg;
  bool m_is_debug = false;
};

OPA_NAMESPACE_DECL3_END

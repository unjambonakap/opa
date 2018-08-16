#pragma once
#include <opa_common.h>

#include <opa/crypto/la/graph.h>
#include <opa/crypto/la/feistel.h>
#include <opa/crypto/la/blocks.h>

OPA_NAMESPACE(opa, crypto, la)

class DESBlock : public CipherGraph {
public:
  static constexpr int blk_size = 32;
  static constexpr int keysize = 48;

  DESBlock() {}
  virtual void init() override;
  virtual void init_graph() override;
  virtual u64 fast_eval2(u64 a, u64 b) const override;
  OPA_TGEN_CIPHERBLOCK0()

public:
  std::array<SboxDesc, 8> m_sbox;
  ExpansionDesc m_ed;
  PermutationDesc m_perm;
  UPTR(PermutationBlock) m_permb;
  UPTR(SboxBlock) m_sboxb[8];
  UPTR(ExpansionBlock) m_eb;
  UPTR(XorBlock) m_xor;
};

class DES {};

OPA_NAMESPACE_END(opa, crypto, la)

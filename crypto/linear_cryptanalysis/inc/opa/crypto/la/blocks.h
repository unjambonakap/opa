
#pragma once
#include <opa/crypto/la/base.h>
#include <opa_common.h>

OPA_NAMESPACE_DECL3(opa, crypto, la)

class SboxDesc : public opa::utils::Initable {
public:
  // output_size is in bits
  virtual void init(u32 output_size, const std::vector<u8> &mp) {
    opa::utils::Initable::init();
    m_output_size = output_size;
    m_input_size = log2_high_bit(mp.size());
    // OPA_DISP("sbox desc >> ", m_output_size, m_input_size);
    this->mp = mp;
  }

  OPA_ACCESSOR_R(u32, m_input_size, input_size);
  OPA_ACCESSOR_R(u32, m_output_size, output_size);

  u8 get(u8 e) const {
    check_init();
    return mp[e];
  }
  const u8 *ptr() const { return &mp[0]; }

  static SboxDesc Rand(int in_bitsize, int out_bitsize);
  static SboxDesc RandWeak(int in_bitsize, int out_bitsize, int nweak, double mx_v);

  u32 m_output_size;
  u32 m_input_size;
  std::vector<u8> mp;
};

// input i goes to perm[i]
class PermutationDesc : public opa::utils::Initable {
public:
  virtual void init(const std::vector<u8> &perm, bool is_inv = false) {
    opa::utils::Initable::init();
    m_size = perm.size();
    m_perm = perm;
    m_iperm.resize(size());
    REP (i, m_size)
      m_iperm[perm[i]] = i;
    if (is_inv)
      std::swap(m_perm, m_iperm);
  }
  u8 get(u8 pos) const { return m_perm[pos]; }
  u8 iget(u8 pos) const { return m_iperm[pos]; }
  u64 get2(u64 v) const;

  u16 *ptr() const {
    if (m_tb.size() == 0) {
      m_tb.resize(1 << size());
      REP (i, m_tb.size()) {
        u16 a = 0;
        REP (j, size())
          if (i >> j & 1)
            a |= 1 << get(j);
        m_tb[i] = a;
      }
    }
    return &m_tb[0];
  }

  OPA_ACCESSOR_R(u32, m_size, size);
  std::vector<u8> m_perm;
  std::vector<u8> m_iperm;
  mutable std::vector<u16> m_tb;
  u32 m_size;
};

class ExpansionDesc : public opa::utils::Initable {
public:
  virtual void init(const std::vector<u8> &mp_) {
    opa::utils::Initable::init();
    this->mp = mp_;
  }

  u8 get_bit(u32 pos, const BitVec &v) const { return v.get(mp[pos]); }

  u64 get(u64 v) const {
    u64 res = 0;
    REP (i, mp.size()) { res |= (v >> mp[i] & 1) << i; }
    return res;
  }

  u32 size() const { return mp.size(); }

  u32 from(u32 a) const { return mp[a]; }

  std::vector<u8> mp;
};

class ExpansionBlock : public CipherBlock {
public:
  virtual void init(u32 input_size, const ExpansionDesc &ed) {
    ed.check_init();
    m_ed = ed;

    CipherBlock::init(Params()
                        .call(CallInfo(input_size, 0, false))
                        .input_size(input_size)
                        .output_size(ed.size())
                        .fast_eval(true));
    desc() = opa::utils::SPrintf("ExpansionBlock(sz=%d)", input_size);
  }

  virtual void do_evaluate(const BitVec &iv, const BitVec &kv,
                           BitVec &ov) const override {
    REP (i, ov.size())
      ov.set(i, m_ed.get_bit(i, iv));
  }

  virtual void do_get_relations(Relations &rels) const override;
  virtual void do_get_relations_diff(Relations &rels) const override;

  virtual u64 fast_eval1(u64 a) const override { return m_ed.get(a); }

  virtual void setup_jit(const JitBuilder &builder) const override {
    OPA_CHECK0(false);
  }

private:
  ExpansionDesc m_ed;
};

class AddBlock : public CipherBlock {
public:
  virtual void init(u32 sz) {
    bool can_fast = sz <= 64;
    CipherBlock::init(Params()
                        .call(CallInfo(sz, 0, true))
                        .input_size(sz * 2)
                        .output_size(sz)
                        .fast_eval(can_fast));
    desc() = opa::utils::SPrintf("AddBlock(sz=%d)", sz);
    mask = (1 << sz) - 1;
    this->sz = sz;
  }

  virtual void do_evaluate(const BitVec &iv, const BitVec &kv,
                           BitVec &ov) const override;

  virtual void do_get_relations(Relations &rels) const override;
  virtual void do_get_relations_diff(Relations &rels) const override;

  virtual void do_get_pre_fail(const Basis &out, Basis &res) const override;

  virtual void setup_jit(const JitBuilder &builder) const override;

  virtual u64 fast_eval2(u64 a, u64 b) const override { return a + b & mask; }

  u64 mask;
  u32 sz;
};

class AddConstBlock : public CipherBlock {
public:
  virtual void init(u64 v, u32 sz) {
    bool can_fast = sz <= 63;
    mask = 1ull << sz;

    CipherBlock::init(
      Params().call(CallInfo(sz, 0, true)).input_size(sz).output_size(sz));
    desc() = opa::utils::SPrintf("AddConstBlock(sz=%d, v=%d)", sz, v);
    m_v = v;
  }
  virtual void do_evaluate(const BitVec &iv, const BitVec &kv,
                           BitVec &ov) const override;

  virtual void do_get_relations(Relations &rels) const override;
  virtual void do_get_relations_diff(Relations &rels) const override;

  virtual void do_get_pre_fail(const Basis &out, Basis &res) const override;

  virtual void setup_jit(const JitBuilder &builder) const override;

  virtual u64 fast_eval1(u64 a) const override { return a + m_v & mask; }

  u64 mask;
  u64 m_v;
};

class SboxBlock : public CipherBlock {
public:
  virtual void init(const SboxDesc &sbox) {
    sbox.check_init();
    CipherBlock::init(Params()
                        .call(CallInfo(sbox.input_size(), 0, true))
                        .input_size(sbox.input_size())
                        .output_size(sbox.output_size())
                        .fast_eval(true));
    desc() = opa::utils::SPrintf("SboxBlock(input=%d,output=%d)", input_size(),
                                 output_size());
    this->sbox = sbox;
  }

  virtual u64 fast_eval1(u64 a) const override { return sbox.get(a); }

  virtual void do_evaluate(const BitVec &iv, const BitVec &kv,
                           BitVec &ov) const override {
    u64 tmp = iv.to<u64>();
    ov.from(sbox.get(tmp));
  }

  virtual void do_get_relations(Relations &rels) const override;
  virtual void do_get_relations_diff(Relations &rels) const override;

  virtual void setup_jit(const JitBuilder &builder) const override;

  SboxDesc sbox;
};

class XorBlock : public CipherBlock {
public:
  virtual void init(u32 sz) {
    m_blk_size = sz;
    CipherBlock::init(Params()
                        .call(CallInfo(sz, 0, true))
                        .input_size(2 * sz)
                        .output_size(sz)
                        .fast_eval(true));
    desc() = opa::utils::SPrintf("XorBlock(%d)", sz);
  }
  virtual void do_evaluate(const BitVec &iv, const BitVec &kv,
                           BitVec &ov) const override {
    u64 a = iv.to<u64>(0, m_blk_size);
    u64 b = iv.to<u64>(m_blk_size, m_blk_size);
    ov.from(a ^ b);
  }

  virtual void do_get_relations(Relations &rels) const override;
  virtual void do_get_relations_diff(Relations &rels) const override;

  virtual void setup_jit(const JitBuilder &builder) const override;

private:
  u32 m_blk_size;
};

class PermutationBlock : public CipherBlock {
public:
  virtual void init(const PermutationDesc &perm, bool fast_eval = true) {
    perm.check_init();
    CipherBlock::init(Params()
                        .call(CallInfo(perm.size(), 0, true))
                        .input_size(perm.size())
                        .output_size(perm.size())
                        .fast_eval(fast_eval));
    desc() = opa::utils::SPrintf("PermutationBlock(%d)", perm.size());
    m_perm = perm;
  }
  virtual void do_evaluate(const BitVec &iv, const BitVec &kv,
                           BitVec &ov) const override {
    REP (i, m_perm.size())
      ov.set(m_perm.get(i), iv.get(i));
  }

  virtual u64 fast_eval1(u64 a) const override { return m_perm.ptr()[a]; }

  virtual void do_get_relations(Relations &rels) const override;
  virtual void do_get_relations_diff(Relations &rels) const override;

  virtual void setup_jit(const JitBuilder &builder) const override;

private:
  PermutationDesc m_perm;
};

class RotationBlock : public PermutationBlock {
public:
  virtual void init(int blk_size, int shift) {
    PermutationDesc pdesc;
    shift %= blk_size;
    if (shift < 0)
      shift += blk_size;

    std::vector<u8> tb(blk_size);
    REP (i, blk_size)
      tb[i] = (i + shift) % blk_size;
    pdesc.init(tb);
    PermutationBlock::init(pdesc);
    desc() = opa::utils::SPrintf("RotationBlock(%d,%d)", blk_size, shift);
  }

private:
  virtual void init(const PermutationDesc &perm, bool fast_eval) override {}
};
OPA_NAMESPACE_DECL3_END

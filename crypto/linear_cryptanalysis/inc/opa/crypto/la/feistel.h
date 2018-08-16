#pragma once
#include <opa_common.h>

#include <opa/crypto/la/graph.h>
#include <opa/or/best_first_search.h>
#include <opa/or/beam_search.h>
#include <opa/threading/runner.h>

OPA_NAMESPACE(opa, crypto, la)

class XorBlock;

// Structure representing differential or linear path.
struct RelState : public opa::utils::ProtobufParams {
  RangeCoverage l;
  RangeCoverage r;
  RangeCoverage lin;
  RangeCoverage rin;

  int round;
  double cost;

  bool operator<(const RelState &x) const {
    if (round != x.round)
      return round < x.round;
    if (l != x.l)
      return l < x.l;
    if (r != x.r)
      return r < x.r;
    return cost < x.cost;
  }
  OPA_TGEN_IMPL(round, l, r, lin, rin, cost);

private:
  friend std::ostream &operator<<(std::ostream &os, const RelState &r) {
    os << RAW_OPA_DISP_VARS(r.round, r.cost, r.l, r.r, r.lin, r.rin);
    return os;
  }
};

class FeistelBlock : public CipherGraph,
                     public std::enable_shared_from_this<FeistelBlock> {
public:
  struct Params : public CipherBlock::Params {
    SPTR(CipherBlock) round;
    int nround;
    bool mid_xor;
    Params() {}
    Params(SPTR(CipherBlock) round, int nround, bool mid_xor = true) {
      this->round = round;
      this->nround = nround;
      this->mid_xor = mid_xor;
    }
  };

  virtual void init(const Params &params);

  virtual void init_graph() override;

  virtual void do_get_relations(Relations &rels) const override;
  virtual u64 fast_eval_tb(u64 a, const u64 *b) const override;

  BitVec get_round_key(const BitVec &kv, int r) const;
  void undo_last_round(const CipherData &data, const BitVec &kv,
                       CipherData *output_data) const;
  void undo_last_round_one(const CipherData::IOPair &pair, const BitVec &kv,
                           CipherData::IOPair *output_pair) const;

  RangeCoverage range_lr_to_one(const RangeCoverage &l, const RangeCoverage &r) const;
  BitVec range_lr_to_bitvec(const RangeCoverage &l, const RangeCoverage &r) const;

  void u64_to_lr(u64 lr, u32 &l, u32 &r) const;
  u64 lr_to_u64(u32 l, u32 r) const;

  OPA_ACCESSOR(SPTR(CipherBlock), m_params.round, round);
  OPA_ACCESSOR_R(int, m_params.nround, nround);
  OPA_ACCESSOR_R(u32, m_blk_size, blk_size);
  OPA_ACCESSOR_R(bool, m_params.mid_xor, mid_xor);
  OPA_ACCESSOR_R(Params, m_params, params);

private:
  Params m_params;
  UPTR(XorBlock) xorblock;
  UPTR(IdBlock) idb;
  u32 m_blk_size;
};

OPA_NAMESPACE_END(opa, crypto, la)

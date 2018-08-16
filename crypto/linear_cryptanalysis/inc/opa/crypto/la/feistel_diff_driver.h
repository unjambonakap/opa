#pragma once
#include <opa_common.h>

#include <opa/crypto/la/feistel_solver.h>
#include <opa/utils/time.h>

OPA_NAMESPACE(opa, crypto, la)

struct RequiredDiff {
  BitVec diff;
  int count;
  RequiredDiff() {}
  RequiredDiff(const BitVec &diff, int count) : diff(diff), count(count) {}
  void gen_cipher_data(const BitVec &key, const CipherBlock *target,
                       CipherData *out_cipher_data) const {
    target->gen_rand_diff(diff, count, key, out_cipher_data);
  }
};

struct RequiredDiffs {
  std::vector<RequiredDiff> entries;
  void gen_cipher_data(const BitVec &key, const CipherBlock *target,
                       CipherData *out_cipher_data) const {
    for (const auto &entry : entries)
      entry.gen_cipher_data(key, target, out_cipher_data);
  }
};

class FeistelDiffSolver : public FeistelSolver {
public:
  RequiredDiffs required_data_to_diffs(const RequiredData &required_data) const;
  FeistelDiffSolver() : FeistelSolver(true) {}

  OPA_DECL_DUPLICATE(FeistelSolver, FeistelDiffSolver);

protected:
  Result rec(const RelState &s, int parid);
  virtual void setup_solver_rels() override;
  virtual void data_to_outrel(const RelState &s, const CipherData &data,
                              OutRel &out) const override;
  virtual double get_score_count(const FindKeysContext *ctx, int pos, u64 key,
                                 int round_id) const override;

private:
  double eval_one_outrel(const OutRel &out_rel, u64 key, u64 check_mask) const;
  Result add_rel(const RelState &rel_);
  void fill_rels();
  std::map<RelState, int> m_rec_cache;
  vi m_rec_par;

  std::vector<RangeCoverage> m_buckets;
  std::vector<RelState> m_final_rels;

  opa::utils::TimeAlloc m_alloc;
  opa::utils::KBestContainer<double, RelState> m_best_rels;
};

OPA_NAMESPACE_END(opa, crypto, la)

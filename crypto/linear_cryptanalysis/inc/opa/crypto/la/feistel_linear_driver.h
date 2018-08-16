#pragma once
#include <opa_common.h>
#include <opa/crypto/la/feistel_solver.h>

#include <opa/crypto/la/graph.h>
#include <opa/or/best_first_search.h>
#include <opa/or/beam_search.h>
#include <opa/threading/runner.h>

OPA_NAMESPACE(opa, crypto, la)

class FeistelBlock;
class XorBlock;

class FeistelLinearSolver : public FeistelSolver {
public:
  FeistelSolverResult solve(const CipherData &data);
  FeistelLinearSolver() : FeistelSolver(false) {}

  OPA_DECL_DUPLICATE(FeistelSolver, FeistelLinearSolver);

protected:
  Result rec(const RelState &s, int parid);
  virtual void setup_solver_rels() override;
  virtual void data_to_outrel(const RelState &s, const CipherData &data,
                              OutRel &out) const override;
  virtual double get_score_count(const FindKeysContext *ctx, int pos, u64 key,
                                 int round_id) const override;

private:
  void conv_one_data(const RelState &s, const CipherData::IOPair &entry,
                     OutRel::InputLinearRel &out_linear_rel) const;

  Result add_rel(const RelState &rel_);
  int evaluate_cost(std::vector<int> *order = nullptr);
  void fill_rels();

  std::set<RangeCoverage> seen;
  std::map<RelState, int> m_rec_cache;

  vi m_rec_par;
  std::vector<Entry> m_pending_entries;
};

OPA_NAMESPACE_END(opa, crypto, la)

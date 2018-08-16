#pragma once
#include <opa_common.h>

#include <opa/crypto/la/graph.h>
#include <opa/or/best_first_search.h>
#include <opa/or/beam_search.h>
#include <opa/threading/runner.h>
#include <opa/utils/clone.h>
#include <opa/crypto/la/feistel.h>

OPA_NAMESPACE(opa, crypto, la)

struct FeistelSolverResult {
  std::vector<BitVec> keys;
};

struct FeistelSolverSharedContext {
  FeistelSolverResult *result;
  bool should_stop = false;
};

struct FeistelSolverSecondPhaseData {
  FeistelSolverSharedContext *shared_context;
  BitVec key;
};

class FeistelSolver : public opa::utils::Initable,
                      public opa::utils::Duplicatable<FeistelSolver> {
public:
  struct FindKeysContext;
  struct RecState {
    // Not owned;
    FindKeysContext *context = nullptr;
    u64 vbase = 0;
    u64 vkey = 0;
    int round = 0;
    double cost = 1.;
    RecState() {}
    RecState(FindKeysContext *context, int round, u64 a, u64 b) {
      this->context = context;
      vbase = a;
      vkey = b;
      this->round = round;
    }
    void update_score(double prob, double expected) {
      double mul = prob / expected;
      const double kFrac = 0.0;
      // cost = pow(cost, kFrac) * pow(mul, 1. - kFrac);
      cost *= mul;
    }

    bool operator<(const RecState &a) const { OPA_LT_OP(a, cost); }
    OPA_DECL_COUT_OPERATOR2(RecState, a.round, a.vbase, a.vkey, a.cost);
  };

  struct FindKeysContext {
    // Not owned;
    const CipherData *cipher_data;
    const std::vector<OutRel> *out_rels;
    std::map<std::pair<int, u64>, double> cached_scores;
    opa::OR::BeamSearch<RecState> search;
    const FeistelSolverSecondPhaseData *second_phase_data;
  };

  // Represent one relation.
  struct Entry {
    // Vectors we must now to be able to evaluate the relation.
    Basis basis;
    RelState rel;

    std::vector<u64> key_to_rel;
    std::vector<u64> bit_to_cachekey; // finding orthogonal vectors might help
                                      // remove this
    std::set<int> used_in_basis;
    std::vector<u64> bit_to_key; // only for additional keybits
    int bit_pos;
    double cum_score;

    Entry() {}

    // These are the two data that must be set when calling finish_rels()
    Entry(const RelState &r, const Basis &b) {
      basis = b;
      rel = r;
    }
  };

  struct Step {
    std::vector<int> checks;
    int main;
  };

  struct RequiredDataPerRound {
    std::vector<RelState> rels;
  };

  struct RequiredData {
    std::vector<RequiredDataPerRound> data_per_round;
  };

  struct Result {
    bool more;
    Result(bool more = true) { this->more = more; }
  };
  typedef std::function<Result(const RelState &rel)> cb_found_rel_t;
  FeistelSolver(bool diff_mode) : m_diff_mode(diff_mode) {}

  virtual void init(FeistelBlock *feistel);
  virtual void init(int nround, FeistelBlock *feistel,
                    const CipherBlock *target);
  const RequiredData &first_phase_find_rels();
  void second_phase_entry(const CipherData &cipher_data,
                          FeistelSolverResult *output_result);
  void
  second_phase_find_keys(const CipherData &cipher_data,
                         const FeistelSolverSecondPhaseData &second_phase_data);

  double test_r2(int nr, const RelState &rel, bool diff = false) const;
  double test_r2_internal(int nr, const RelState &rel, bool diff = false) const;

protected:
  virtual void setup_solver_rels() = 0;
  virtual void data_to_outrel(const RelState &s, const CipherData &data,
                              OutRel &out) const = 0;

  virtual double get_score_count(const FindKeysContext *ctx, int pos, u64 key,
                                 int round_id) const = 0;

  void evaluate_entry_real_cost(Entry &entry);

  bool add_rel(const RelState &rel);
  void find_keys_start(FindKeysContext *context);
  void find_keys_rec(const RecState &s);
  void find_keys_end(const RecState &s);
  u64 key_to_base(u64 k);
  u64 base_to_key(u64 b);
  bool filter_res(const FindKeysContext *context,
                  std::vector<std::pair<double, RecState> > &res,
                  const Entry &entry) const;
  double get_score(FindKeysContext *context, int data_pos, const RecState &s,
                   u64 cache_key, int round_id) const;
  void finish_rels(const std::vector<Entry> &entries, int main_size,
                   const std::vector<RelState> &rels);

  int N; // key size

  int nround;
  const CipherBlock *target;
  FeistelBlock *feistel;

private:
  bool m_diff_mode;
  opa::math::common::Matrix<u32> inv;
  std::vector<u64> real_key;
  std::vector<Entry> m_entries;
  std::vector<Step> m_steps;
  std::vector<Basis> m_incremental_basis;
  std::vector<RelState> m_rels;

  std::shared_ptr<FeistelSolver> m_next_solver;
  std::shared_ptr<FeistelBlock> m_next_blk;
  std::shared_ptr<RequiredData> m_required_data;
};

OPA_NAMESPACE_END(opa, crypto, la)

#pragma once
#include <opa/crypto/la/algo.h>
#include <opa/crypto/stream/common.h>
#include <opa/crypto/stream/utils.h>
#include <opa/math/common/Matrix.h>
#include <opa/math/common/algebra/basis.h>
#include <opa/math/common/stats.h>
#include <opa/utils/csv.h>
#include <opa/utils/misc.h>
#include <opa_common.h>

OPA_NAMESPACE(opa, crypto, stream)

class SolverPlan : public opa::utils::Initable {
public:
  void init(int input_len, const std::set<int> &known_data);

  void set_final_rels(int pos, const std::vector<SolverRel> &rels);

  void find_one_rels(StepDescription &step, const SolverState &s, u64 lim,
                     int zero_ub, const KeyRelsList &key_rels);

  void find_rels(StepDescription &step, const std::vector<SolverState> &to_use,
                 u64 last_count, int rem);

  void setup_step(StepDescription &step, int rem, int nfix);

  u64 get_num_of_rels(const std::vector<pii> &used) const;

  void build_plan();

  void eval_rels(const StepDescription &step,
                 const opa::math::common::BitVec &key,
                 std::vector<int> &out_data, u64 &out_tot) const;

  OPA_ACCESSOR_R(std::vector<la::Relation>, m_relations, filtered_relations);
  int nsteps() const { return m_steps.size(); }
  OPA_ACCESSOR_R(std::vector<StepDescription>, m_steps, steps);

  int m_input_len;
  std::vector<la::Relation> m_relations;
  std::vector<StepDescription> m_steps;
  RelsStore m_rels_store;
  std::set<int> m_known_data;
};

OPA_NAMESPACE_END(opa, crypto, stream)

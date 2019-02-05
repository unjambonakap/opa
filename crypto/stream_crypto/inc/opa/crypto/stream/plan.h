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


class SolverPlan : public opa::utils::BaseStorable {
public:
  OPACS_GETTER_BY_CL(SolverPlan);
  void init(int input_len);

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

  int nsteps() const { return m_steps.size(); }
  OPA_ACCESSOR_R(std::vector<StepDescription>, m_steps, steps);

  int m_input_len;
  std::vector<StepDescription> m_steps;
  RelsStore m_rels_store;
  int lim_nfix = 31;
  int n_bruteforce = 4;


  OPA_TGEN_IMPL(m_input_len, m_steps, m_rels_store);
};

OPA_NAMESPACE_END(opa, crypto, stream)

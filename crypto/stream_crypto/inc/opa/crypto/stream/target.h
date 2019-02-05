#pragma once
#include <opa/crypto/stream/common.h>
#include <opa/utils/DataStruct.h>
#include <opa/crypto/stream/plan.h>
#include <opa/crypto/stream/solver.h>
#include <opa/math/common/Matrix.h>

OPA_NAMESPACE(opa, crypto, stream)

class SboxLfsr : public opa::utils::Initable {
public:
  struct Params {
    std::vector<LFSR_GF2_small> lfsrs;
    la::SboxDesc sbox;
  };

  virtual void init(const Params &params);
  void reset();
  u32 get();
  OPA_ACCESSOR(std::vector<LFSR_GF2_small>, m_params.lfsrs, lfsrs);
  OPA_ACCESSOR(la::SboxDesc, m_params.sbox, sbox);
  void set_rand();
  opa::math::common::BitVec get_state() const;

private:
  Params m_params;
  opa::utils::BitDataStream m_bs;
  int m_n;
};
typedef std::unordered_map<int, int> ObservedData;

class SboxLfsrSolver {
public:
  struct Params {
    SboxLfsr lfsr;
    opa::math::common::BitVec ans_input;
    bool check = false;
  };

  void init(const Params &params);
  void test_rel(const la::Relation &rel);

  // dont like the no const here, but I would need mutable stuff for rels_store
  void load_plan_from_data(const std::set<int> &known_data);

  std::vector<la::Relation> get_used_relations() const;

  opa::math::common::BitVec solve(const ObservedData &obs_data);

  Params m_params;
  SPTR(SolverPlan) m_plan;
  SPTR(SolverPlan) m_own_plan;

  void configure_output_bits(const std::set<int> &data);

  std::vector<int> m_output_bit_list;
  std::map<int, int> m_output_bit_rmp;
  std::vector<int> m_input_pos;
  std::vector<opa::math::common::MulMatrixF2<false> > m_tb;
  void check_key(const opa::math::common::BitVec &key,
                 const ObservedData &obs_data) const;

  int m_ni;
  int m_nlfsr;
  int m_no;
  int m_input_len;
};

void get_sbox_lfsr_plan(const SboxLfsr &lfsr, const std::set<int> &known_data,
                        SolverPlan *plan);
void solve_sbox_lfsr(const SboxLfsr &lfsr, const ObservedData &obs_data,
                     const SolverPlan &plan);

OPA_NAMESPACE_END(opa, crypto, stream)

#pragma once
#include <opa/crypto/stream/common.h>
#include <opa/crypto/la/algo.h>
#include <opa/math/common/fast_gf2.h>
#include <opa/math/common/algebra/basis.h>
#include <opa/math/common/Matrix.h>

OPA_NAMESPACE(opa, crypto, stream)

class SolverPlan;
class StepDescription;

class Solver {

public:
  struct Params {
    opa::math::common::BitVec ans_input;
    bool check = false;
    const SolverPlan *plan;
    int input_len;
  };

  void init(const Params &params);

  opa::math::common::BitVec solve();

  void check_key(const opa::math::common::BitVec &key);

  void test_rel(const la::Relation &rel);

  bool solve_step(const StepDescription &step, opa::math::common::BitVec &key);

  Params m_params;
};

OPA_NAMESPACE_END(opa, crypto, stream)

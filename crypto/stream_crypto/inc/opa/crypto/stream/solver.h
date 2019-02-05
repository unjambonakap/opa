#pragma once
#include <opa/crypto/la/algo.h>
#include <opa/crypto/stream/common.h>
#include <opa/crypto/stream/plan.h>
#include <opa/math/common/Matrix.h>
#include <opa/math/common/algebra/basis.h>
#include <opa/math/common/fast_gf2.h>
#include <opa/threading/auto_job.h>
#include <opa/threading/runner.h>

OPA_NAMESPACE(opa, crypto, stream)

struct SolverWork : public opa::utils::ProtobufParams {
  opa::math::common::BitVecRepr input;
  OPA_TGEN_IMPL(input);
};

struct SolverRes : public opa::utils::ProtobufParams {
  SolverRes(){}
  SolverRes(opa::math::common::BitVecRepr res, double score):res(res), score(score){}

  opa::math::common::BitVecRepr res;
  double score;
  OPA_TGEN_IMPL(res, score);
  OPA_DECL_LT_OPERATOR(SolverRes, score);
};

typedef std::vector<SolverRes> SolverResList;

class Solver : public opa::threading::AutoJob<SolverWork, SolverResList> {
public:
  struct Params : public opa::utils::ProtobufParams {
    opa::math::common::BitVec ans_input;
    mutable opa::math::common::BitVecRepr ans_input_repr;
    bool check = false;
    std::shared_ptr<SolverPlan> plan;
    int input_len = -1;

    virtual void before_store() const { ans_input_repr = ans_input.to_repr(); }
    virtual void after_load() { ans_input = ans_input_repr.to_bitvec(); }
    OPA_TGEN_IMPL(ans_input_repr, check, plan, input_len);
  };

  void init(const Params &params);

  SolverResList solve(opa::math::common::BitVec &key);
  SolverResList solve();
  void solve_rec(int pos, const opa::math::common::BitVec &key, SolverResList &res);
  opa::math::common::BitVec get_best() const;
  opa::math::common::BitVec solve_best();

  void test_rel(const la::Relation &rel);

  virtual void auto_worker_do_work(const SolverWork &data,
                                   SolverResList &out_res) override;

  virtual void auto_server_get_work() override;

  virtual void
  auto_server_set_work_result(const SolverResList &res,
                              opa::threading::DataId data_id) override {
    reslist.insert(reslist.end(), ALL(res));
  }

  virtual bool server_want_more_results() const override {
    return true;
  }

  SolverResList reslist;
  Params m_params;

  OPA_CLOUDY_JOB_DECL2(Solver);
  OPA_TGEN_IMPL(m_params);
};

OPA_NAMESPACE_END(opa, crypto, stream)

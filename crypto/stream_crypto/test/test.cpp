#include <gtest/gtest.h>
#include <opa_common.h>
#include <opa/crypto/stream/plan.h>
#include <opa/crypto/stream/utils.h>
#include <opa/threading/dispatcher.h>
#include <opa/threading/runner.h>

using namespace std;
using namespace opa::crypto::stream;
using namespace opa::math::common;
using namespace opa::threading;

Dispatcher *g_dispatcher;

TEST(Plan, Quadrisection) {
  RelsStore store;
  int input_len = 42;
  int ni = 1 << 15;
  int nx = 4;
  store.init(input_len);
  store.add_new_typ(0);
  REP (i, ni) {
    SolverRel rel;
    rel.v = BitVec::rand(input_len);
    store.add_new_rel(rel);
  }
  KeyRelsList key_rels = store.to_key_rels(input_len);
  SolverState state;
  state.tb = { MP(0, nx) };

  StepDescription desc;
  StepHelper::Params params(&store, &key_rels, &state, input_len, &desc);
  params.dispatcher = g_dispatcher;
  StepHelper helper(params);
  helper.solve(-1, 0);

  int nrels = desc.rels.size();
  u64 tot_rels = get_count_for_sizes(state.tb, store.get_sizes());
  u64 expected = tot_rels >> input_len;
  OPA_DISP0(tot_rels, input_len, nrels, expected);
  puts("LA");
}

GTEST_API_ int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  opa::init::opa_init(argc, argv);
  Runner::Build();
  Runner runner;
  runner.run_both();
  Dispatcher *dispatcher = runner.dispatcher();
  g_dispatcher = dispatcher;
  return RUN_ALL_TESTS();
}

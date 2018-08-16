#include <opa/crypto/stream/solver.h>
#include <opa/crypto/stream/target.h>
#include <opa/crypto/stream/plan.h>

using namespace std;
using namespace opa::math::common;
using namespace opa::crypto::la;

OPA_NAMESPACE(opa, crypto, stream)

void Solver::init(const Params &params) { m_params = params; }

BitVec Solver::solve() {
  BitVec key(m_params.input_len);

  for (int i = 0; i < m_params.plan->nsteps(); ++i) {
    OPA_CHECK0(solve_step(m_params.plan->steps()[i], key));
  }

  puts("SOLUTION >>> >+++++====J");
  key.disp();
  m_params.ans_input.disp();
  return key;
}

bool Solver::solve_step(const StepDescription &step, BitVec &key) {
  int t, need;
  int nfix = step.fixed_input.size();
  puts("");
  OPA_DISP("SOLVE ROUND >> ", nfix, step.fixed_input);
  // in a relation:
  // [0,rem[ -> bits we try to nullify
  // [rem, nbits=rem+nfix[ = we'll bruteforce those, find most likely
  // [nbits,  input_len[ -> key bits that are either already known or we
  // know
  // the relation will have zeroes there
  // [input_len, input_len+output_len[ -> output bits
  OPA_DISP0(step.rels.size());

  u64 nfixpow = 1ull << nfix;

  vector<int> tb(nfixpow, 0);

  u64 tot;
  m_params.plan->eval_rels(step, key, tb, tot);
  OPA_DISP0(tot);
  do_walsh(tb.data(), nfix);

  s64 best = -1;
  REP64(i, nfixpow) {
    if (best == -1 || abs(tb[best]) < abs(tb[i]))
      best = i;
  }

  OPA_DISP0("bEST >> ", best, abs(tb[best]), key);
  if (m_params.check) {
    u32 correct_key = step.get_fixed_key_from_bitvec(m_params.ans_input);
    OPA_DISP0(correct_key, abs(tb[correct_key]), m_params.ans_input.str());
    OPA_CHECK0(correct_key == best);
  }

  step.update_bitvec_from_fixed_key(best, key);
  return true;
}

OPA_NAMESPACE_END(opa, crypto, stream)

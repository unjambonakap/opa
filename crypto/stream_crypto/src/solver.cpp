#include <opa/crypto/stream/plan.h>
#include <opa/crypto/stream/solver.h>
#include <opa/crypto/stream/target.h>
#include <opa/threading/runner.h>

using namespace std;
using namespace opa::math::common;
using namespace opa::crypto::la;

OPA_NAMESPACE(opa, crypto, stream)

void Solver::init(const Params &params) {
  puts("");
  m_params = params;
  if (m_params.check) {
    int cnt = 0;
    auto &relstore = m_params.plan->m_rels_store;
    double expected = 0;
    for (auto &step : m_params.plan->steps()) {
      OPA_DISP("on step ", step.rels.size());
      SolverRel cur;
      cur.init(m_params.input_len);

      for (auto &relids : step.rels) {
        double bias = 1;
        for (auto &relid : relids) {

          auto &rel = relstore.rels[relid];
          bias *=
            relstore.get_type_bias(relstore.get_rel_type(relid));
          cur.sadd(rel);
        }
        cnt += cur.check(m_params.ans_input);
        expected += bias_to_prob(bias);
      }
      OPA_DISP0(cnt, int(expected));
    }
  }
}

SolverResList Solver::solve() {
  BitVec key(m_params.input_len);
  auto &step0 = m_params.plan->steps()[0];
  OPA_DISP0(step0.n_bruteforce);
  if (step0.n_bruteforce != -1) {
    SolverResList res;

    REP (i, 1 << step0.n_bruteforce) {
      key.from(i, m_params.input_len - step0.n_bruteforce, step0.n_bruteforce);
      auto tmp = solve(key);
      res.insert(res.end(), ALL(tmp));
    }
    puts("REturn");
    return res;
  } else
    return solve(key);
}

BitVec Solver::solve_best() {
  reslist = solve();
  puts("SOLUTION >>> >+++++====J");
  BitVec key = get_best();
  key.disp();
  if (m_params.check) {
    m_params.ans_input.disp();
  }
  return key;
}

BitVec Solver::get_best() const {
  OPA_CHECK0(reslist.size() > 0);
  return std::max_element(ALL(reslist))->res.to_bitvec();
}

SolverResList Solver::solve(opa::math::common::BitVec &key) {

  //OPA_DISP("start solve ", key);
  SolverResList res;
  for (auto &step : m_params.plan->steps()) {
    if (step.n_bruteforce != -1) continue;
    if (m_params.check) {
      int last_found = step.fixed_input.back() + 1;
      if (m_params.ans_input.extract(last_found) != key.extract(last_found)) {
        //OPA_DISP("Differing ", last_found,
        //         m_params.ans_input.extract(last_found),
        //         key.extract(last_found));
        return res;
      }
    }
    this->solve_rec(0, key, res);
  }
  return res;
}

void Solver::solve_rec(int pos, const opa::math::common::BitVec &key,
                       SolverResList &res) {
  if (pos == m_params.plan->steps().size()) {
    res.emplace_back(
      SolverRes{ key.to_repr(), m_params.plan->m_rels_store.get_score(key) });
    return;
  }
  auto &step = m_params.plan->steps()[pos];
  if (step.n_bruteforce != -1) return solve_rec(pos + 1, key, res);

  int t, need;
  int nfix = step.fixed_input.size();
  puts("");
  OPA_DISP("SOLVE ROUND >> ", nfix, step.fixed_input, key);
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

  opa::utils::KBestContainer<int, int> container;
  container.init(3);
  s64 best = -1;
  REP64(i, nfixpow) { container.add(abs(tb[i]), i); }
  std::vector<int> bests;
  container.get_items(bests);
  bool found = false;
  for (auto &cnd : bests) {
    if (m_params.check) {
      u32 correct_key = step.get_fixed_key_from_bitvec(m_params.ans_input);
      OPA_DISP0(correct_key, abs(tb[correct_key]), m_params.ans_input.str(),
                cnd, abs(tb[cnd]));
      if (correct_key != cnd) continue;
      found = 1;
    }
    BitVec nkey = key;
    step.update_bitvec_from_fixed_key(cnd, nkey);
    solve_rec(pos + 1, nkey, res);
  }

  if (m_params.check) {
    OPA_CHECK0(found);
  }
}

void Solver::auto_worker_do_work(const SolverWork &data,
                                 SolverResList &out_res) {
  BitVec key = data.input.to_bitvec();
  out_res = solve(key);
}

void Solver::auto_server_get_work() {
  auto &step0 = m_params.plan->steps()[0];

  BitVec key(m_params.input_len);
  OPA_DISP0(step0.n_bruteforce, key.size(), m_params.input_len);

  REP (i, 1 << step0.n_bruteforce) {
    key.from<u32>(i, m_params.input_len - step0.n_bruteforce);
    SolverWork work;
    OPA_DISP("sending ", i, key, key.size());
    work.input = key.to_repr();
    bool more;
    cb()(work, more);
    if (!more) break;
  }
}

OPA_CLOUDY_REGISTER_BASE(Solver);
OPA_CLOUDY_JOB_IMPL(Solver);

OPA_NAMESPACE_END(opa, crypto, stream)

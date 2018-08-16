#include <opa/crypto/stream/plan.h>
#include <gsl/gsl_poly.h>
#include <opa/threading/runner.h>

using namespace std;
using namespace opa::math::common;

OPA_NAMESPACE(opa, crypto, stream)

double bias_to_prob(double bias) { return (1 + bias) / 2; }
double prob_to_bias(double prob) { return 2 * prob - 1; }

// N: number of rels
// M: number total of keys.
double get_success_prob(u64 N, double p0, u64 M, int rem) {
  N /= (1ull << rem);

  double p1 = 1. / 2;
  double u0 = N * p0;
  double u1 = N * p1;
  double v0 = N * p0 * (1 - p0);
  double v1 = N * p1 * (1 - p1);

  double lim;
  {
    // neyman pearson test
    double a = 1. / 2 / v1 - 1. / 2 / v0;
    double b = u0 / v0 - u1 / v1;
    double c =
      u1 * u1 / 2 / v1 - u0 * u0 / 2 / v0 + 1. / 2 * log(v1 / v0) - log(M - 1);
    double s1 = 0, s2 = 0;
    int res = gsl_poly_solve_quadratic(a, b, c, &s1, &s2);
    if (res != 2)
      return 0;
    lim = s1;
    if (s1 < u1 || s1 > u0)
      lim = s2;
    // OPA_DISP0(N, p0, s1, s2, lim, u0, u1);
    // OPA_CHECK0(lim >= u1 && lim <= u0);
    lim = min(s1, s2);
  }
  double success_prob =
    1. / 2 * FloatUtil::erfc<double>((lim - u0) / sqrt(2 * v0));
  // OPA_DISP0(N, p0, M, lim, success_prob);
  return success_prob;
}

u64 get_num_event_for_success_prob(u64 ub, double prob, int nfix, int rem,
                                   const double thresh) {
  u64 T = 1;
  u64 H = ub;
  u64 best = ub;
  while (T <= H) {
    u64 M = (T + H) / 2;
    double tmp = get_success_prob(M, prob, 1ull << nfix, rem);
    if (tmp > thresh) {
      best = M, H = M - 1;
    } else
      T = M + 1;
  }
  return best;
}

void SolverPlan::init(int input_len, const std::set<int> &known_data) {
  opa::utils::Initable::init();
  m_known_data = known_data;
  m_input_len = input_len;
  m_rels_store.init(input_len);
}

void SolverPlan::find_one_rels(StepDescription &step, const SolverState &s,
                               u64 lim, int zero_ub,
                               const KeyRelsList &key_rels) {
  StepHelper::Params params(&m_rels_store, &key_rels, &s, zero_ub, &step);
  params.dispatcher = opa::threading::Runner::EasyDispatcher();

  StepHelper helper(params);
  helper.solve(lim);
}

void SolverPlan::find_rels(StepDescription &step,
                           const vector<SolverState> &to_use, u64 last_count,
                           int rem) {

  KeyRelsList key_rels = m_rels_store.to_key_rels(rem);

  REP (i, to_use.size()) {
    u64 lim = 0;
    if (i == to_use.size() - 1)
      lim = last_count;
    if (1) {
      if (i != to_use.size() - 1)
        continue;
      lim = -1;
    }
    find_one_rels(step, to_use[i], lim, rem, key_rels);
  }
}

void SolverPlan::setup_step(StepDescription &step, int rem, int nfix) {

  priority_queue<SolverState> q;
  q.push(SolverState({ pii(0, 0) }, 1));

  u64 count = 0;
  double tot_prob = 1;
  const double k_prob_threshold = 0.95;
  vector<SolverState> to_use;

  bool first = 1;
  bool ok = 0;
  while (!q.empty()) {
    SolverState s = q.top();
    q.pop();

    if (!first) {
      u64 nv = get_num_of_rels(s.tb);
      u64 ncount = count + nv;
      double s_prob = bias_to_prob(s.cost);

      double nprob =
        tot_prob * (1. * count / ncount) + (s_prob * 1. * nv / ncount);

      // not the first dummy push
      to_use.pb(s);

      double success_prob = get_success_prob(ncount, nprob, 1ull << nfix, rem);
      OPA_DISP("ON state >> ", s.tb, s.cost, ncount, nprob, nfix, rem,
               success_prob);
      if (success_prob > k_prob_threshold) {
        u64 want = get_num_event_for_success_prob(ncount, nprob, nfix, rem,
                                                  k_prob_threshold);
        OPA_CHECK0(want >= count);
        want = want - count;
        OPA_DISP("want >> ", want, ncount, nv / want, want >> rem);
        want *= 2;
        find_rels(step, to_use, want, rem);
        ok = 1;
        break;
      }

      count = ncount;
      tot_prob = nprob;
    }
    first = 0;

    FOR (i, s.tb.back().ST, m_rels_store.get_num_types()) {
      SolverState ns = s;
      if (i == ns.tb.back().ST) {
        ++ns.tb.back().ND;
      } else {
        ns.tb.pb(MP(i, 1));
      }
      ns.cost *= m_rels_store.get_type_bias(i);
      q.push(ns);
    }
  }
  OPA_CHECK0(ok);
}

u64 SolverPlan::get_num_of_rels(const vector<pii> &used) const {
  bignum res = 1;
  for (auto &e : used) {
    res *= nchoosek(m_rels_store.get_num_rels_by_type(e.ST), e.ND);
  }
  return res.getu64();
}

void SolverPlan::build_plan() {
  m_steps.clear();
  /*
     At each step, we will fix bits [rem, rem+nfix[.
     The bits from [rem+nfix, inf[ are already found.
     We need to find rels where [0, rem[ are zeroed.
     */

  int rem = m_input_len;
  while (rem > 0) {
    int lim_nfix = 31;
    int nfix = min(lim_nfix, rem);
    OPA_DISP("ON STEP >> ", rem);
    rem -= nfix;
    m_steps.emplace_back();
    REP (i, nfix)
      m_steps.back().fixed_input.pb(rem + i);
    setup_step(m_steps.back(), rem, nfix);
  }
}

void StepDescription::update_bitvec_from_fixed_key(
  u32 fixed_key, opa::math::common::BitVec &out_key) const {
  REP (i, this->fixed_input.size())
    out_key.set(fixed_input[i], fixed_key >> i & 1);
}

u32 StepDescription::get_fixed_key_from_bitvec(
  const opa::math::common::BitVec &key) const {
  u32 fixed_key = 0;
  REP (i, fixed_input.size())
    fixed_key |= s64(key.get(fixed_input[i])) << i;
  return fixed_key;
}

void SolverPlan::eval_rels(const StepDescription &step,
                           const opa::math::common::BitVec &key,
                           std::vector<int> &out_data, u64 &out_tot) const {
  OPA_CHECK0(step.fixed_input.size() <= 31);
  // first 31 bits are the keys. 32th is the constant
  vector<u32> evaled_rels(m_rels_store.get_num_rels());
  out_tot = 0;
  OPA_DISP("EVALING", key.str());

  REP (i, evaled_rels.size()) {
    const auto &rel = m_rels_store.rels[i];
    u32 v = rel.c;
    REP (j, m_input_len) {
      if (rel.v.get(j))
        v ^= key.get(j);
    }
    u32 w = step.get_fixed_key_from_bitvec(rel.v);
    evaled_rels[i] = w | (v << 31);
  }

  REP (i, step.rels.size()) {
    // IMPORTANT: bits of key that are not yet set should be at 0
    // Otherwise this won't work
    auto &rel_ids = step.rels[i];
    u32 state = 0;
    for (auto &j : rel_ids) {
      OPA_CHECK0(j >= 0 && j < evaled_rels.size());
      state ^= evaled_rels[j];
    }

    if (0) {
      SolverRel tmp(m_input_len);
      for (auto &j : rel_ids) {
        tmp.sadd(m_rels_store.rels[j]);
      }
      OPA_DISP(i, tmp.v.str(), rel_ids);
      REP (j, step.fixed_input[0]) { OPA_CHECK0(tmp.v.get(j) == 0); }
    }

    u32 w = state & ~(1ull << 31);

    int v = state >> 31;
    ++out_tot;
    v = 2 * v - 1;
    out_data[w] += v;
  }
}

OPA_NAMESPACE_END(opa, crypto, stream)

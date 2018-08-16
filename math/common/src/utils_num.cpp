#include <opa/math/common/utils_num.h>
#include <opa/or/best_first_search.h>

OPA_NAMESPACE_DECL3(opa, math, common)

bignum SmallestDivHelper::compute(const BGFactors &factors,
                                  const bignum &lb) {
  typedef OR::BestFirstSearch<State, std::greater<State> > Search_t;

  Search_t search;
  Search_t::Params params;
  bool found = true;
  params.search_params.func = [&](OR::Search<State> *self, const State &state) {
    if (state.v >= lb) {
      this->res = state;
      self->stop();
      return;
    }
    if (state.pos == factors.size()) return;

    State ns = state;
    ++ns.pos;
    REP (i, factors[state.pos].second + 1) {
      search.add(ns);
      ns.v *= factors[state.pos].first;
      ns.factors.push_back(factors[state.pos].first);
    }
  };

  OPA_CHECK0(found == true);
  State s0;
  s0.pos = 0;
  s0.v = 1;
  search.init(params);
  search.add(s0);
  search.start();
  return res.v;
}

SmoothPrimeOrderHelper spoh_for_gfp_split(1000, 1e9, { 2, 3, 5, 7 });

SmoothPrimeOrderHelper::SmoothPrimeOrderHelper(
  int lb, int ub, const std::vector<int> &prime_set) {
  m_lb = lb;
  m_ub = ub;
  m_prime_set = prime_set;
}

void SmoothPrimeOrderHelper::build() const {
  if (m_init) return;
  m_init = true;
  struct State {
    int pos;
    u64 v;
    S64Factors factors;
  };

  typedef OR::Search_dfs<State> Search_t;

  Search_t search;
  Search_t::SearchParams params;
  params.func = [&](OR::Search<State> *self, const State &state) {
    if (state.pos == this->m_prime_set.size()) {
      if (state.v >= this->m_lb) {
        u64 v2 = state.v + 1;
        if (!testPrime(v2)) return;
        this->m_data.emplace_back(v2, state.factors);
      }
      return;
    }

    State ns = state;
    ++ns.pos;
    search.add(ns);

    ns.factors.emplace_back(this->m_prime_set[state.pos], 0);
    while (1) {
      ns.v *= this->m_prime_set[state.pos];
      ++ns.factors.back().second;
      if (ns.v >= this->m_ub) break;
      search.add(ns);
    }
  };

  State s0;
  s0.pos = 0;
  s0.v = 1;
  search.init(params);
  search.add(s0);
}

OPA_NAMESPACE_DECL3_END

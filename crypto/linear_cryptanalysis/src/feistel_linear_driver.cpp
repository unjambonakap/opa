#include <opa/crypto/la/feistel_linear_driver.h>
#include <opa/crypto/la/feistel.h>
#include <opa/crypto/la/blocks.h>
#include <opa/math/common/Utils.h>

using namespace std;

OPA_NAMESPACE(opa, crypto, la)

FeistelSolverResult FeistelLinearSolver::solve(const CipherData &data) {
  const RequiredData &required_data = this->first_phase_find_rels();
  FeistelSolverResult result;
  for (auto &r : required_data.data_per_round.back().rels) {
    puts("Testing REL >> ");
    OPA_DISP0(r);
    this->test_r2(10, r);
  }
  OPA_DISP0(opa::utils::Join("\n", required_data.data_per_round.back().rels));
  this->second_phase_entry(data, &result);
  return result;
}

void FeistelLinearSolver::conv_one_data(
  const RelState &s, const CipherData::IOPair &entry,
  OutRel::InputLinearRel &out_linear_rel) const {

  out_linear_rel.first = entry.second.extract(feistel->blk_size()).to<u64>();
  out_linear_rel.second =
    entry.first.dot(s.lin + s.rin.shift(feistel->blk_size())) ^
    entry.second.extract(0, feistel->blk_size()).dot(s.l) ^
    entry.second.extract(feistel->blk_size()).dot(s.r);
}

void FeistelLinearSolver::data_to_outrel(const RelState &s,
                                         const CipherData &data,
                                         OutRel &out_rel) const {
  out_rel.out = s.l.to<u64>();
  out_rel.cost = s.cost;

  for (auto &d : data.x) {
    out_rel.in.emplace_back();
    conv_one_data(s, d, out_rel.in.back());
  }
  for (auto &d : data.x_diff) {
    out_rel.in.emplace_back();
    conv_one_data(s, d.first, out_rel.in.back());
    out_rel.in.emplace_back();
    conv_one_data(s, d.second, out_rel.in.back());
  }
}

FeistelSolver::Result FeistelLinearSolver::add_rel(const RelState &rel_) {
  RelState rel = rel_;

  const Result kDone = Result(false);
  const Result kNotDone = Result(true);
  OPA_DISP("try add rel ", rel.l, rel.r, rel.lin, rel.rin);
  if (rel.l.size() == 0)
    return kNotDone;
  if (seen.count(rel.l))
    return kNotDone;

  int n = m_pending_entries.size();
  if (n > 30 && m_pending_entries[n / 2].rel.cost > rel.cost)
    return kNotDone;
  // m_pending_entries are always sorted by cost

  auto lst = target->get_pre_key(rel.l);
  Basis basis;
  basis.init(lst, N);
  OPA_DISP("BSIS dim for >> ", lst.size(), basis.dim(), rel.l, lst, N);
  if (basis.dim() == 0)
    return kNotDone;
  m_pending_entries.pb(Entry(rel, basis));

  seen.insert(rel.l);

  struct Comparator {
    bool operator()(const Entry &a, const Entry &b) const {
      if (fabs(a.rel.cost - b.rel.cost) > 1e-8)
        return a.rel.cost > b.rel.cost;
      return 0;
    }
  };

  sort(ALL(m_pending_entries), Comparator());

  int c2 = evaluate_cost();
  // every rel, we expect a decrease of cost of v
  // give up when cost < nrel * v
  // TODO: take into account minimum possible cost
  double v = 0.3;
  OPA_DISP0(m_pending_entries.size(), c2, v);
  if (m_pending_entries.size() >= 3 * target->key_size() &&
      c2 < m_pending_entries.size() * v) {
    return kDone;
  }
  return kNotDone;
}

void FeistelLinearSolver::setup_solver_rels() {
  fill_rels();

  vi order;
  evaluate_cost(&order);
  vector<Entry> final_entries;
  set<int> used;
  REP (i, order.size()) {
    used.insert(order[i]);
    final_entries.pb(m_pending_entries[order[i]]);
  }

  REP (i, m_pending_entries.size()) {
    if (used.count(i))
      continue;
    final_entries.pb(m_pending_entries[i]);
  }

  std::vector<RelState> rels;
  for (auto &e : final_entries)
    rels.pb(e.rel);
  finish_rels(final_entries, order.size(), rels);
}

void FeistelLinearSolver::fill_rels() {

  if (this->nround == 1) {
    REP (i, feistel->blk_size()) {
      RelState cur;
      cur.lin.add(i);
      cur.l.add(i);
      cur.cost = 1;
      Result res = add_rel(cur);
      if (!res.more)
        break;
    }

  } else if (this->nround == 2) {
    REP (i, feistel->blk_size()) {
      RelState cur;
      cur.rin.add(i);
      cur.l.add(i);
      cur.cost = 1;
      Result res = add_rel(cur);
      if (!res.more)
        break;
    }

  } else {
    RelState s;
    s.cost = 1;
    s.round = 0;
    Result x = rec(s, -1);
  }
}

FeistelSolver::Result FeistelLinearSolver::rec(const RelState &s, int parid) {
  OPA_DISP0(s, parid);
  if (m_rec_cache.count(s)) {
    OPA_DISP0("CACHED");
    return Result();
  }

  int r = m_rec_cache[s] = m_rec_par.size();
  m_rec_par.pb(parid);

  if (s.round == feistel->nround() - 1) {
    RelState ns = s;
    std::swap(ns.r, ns.l);
    Result res = add_rel(ns);
    return res;
  }

  const Relations &rels = feistel->round()->get_input_relations(false);
  if (s.r.size() == 0) {
    RelState ns = s;
    ns.l = s.r;
    ns.r = s.l;
    ns.cost = s.cost;
    ns.round = s.round + 1;
    Result res = rec(ns, r);
    if (!res.more)
      return res;
  }

  for (auto &x : rels.tb) {
    RelState ns = s;
    ns.r = s.l;
    ns.l = s.r;

    ++ns.round;
    ns.cost *= x.cost.prob();
    if (s.round < 2)
      ;
    else {
      if (ns.l != x.out)
        continue;
    }
    if (s.round == 0) {
      ns.rin.do_xor(x.in);
      ns.lin.do_xor(x.out);
      ns.l = x.out;
      ns.r = x.in; // cancel next xor
    } else if (s.round == 1) {
      ns.rin.do_xor(x.out);
      ns.l.do_xor(x.out);
    }
    ns.r.do_xor(x.in);
    Result res = rec(ns, r);
    if (!res.more)
      return res;
  }

  return Result();
}

double FeistelLinearSolver::get_score_count(const FindKeysContext *ctx, int pos,
                                            u64 key, int round_id) const {
  const OutRel &orel = ctx->out_rels->at(pos);
  u32 cnt = 0;
  u32 tot = 0;

  for (auto &x : orel.in) {
    u64 ov = target->fast_eval2(x.first, key);
    int c = opa::utils::dot(ov, orel.out) ^ x.second;
    cnt += c;
    tot += 1;
  }

  OPA_CHECK0(tot > 0);
  double res = 1. * cnt / tot;
  res = fabs(res - 0.5) * 2;
  return res;
}

struct Scorer {
  double best_prob;
  Scorer(double best_prob) { this->best_prob = best_prob; }

  // higher the better
  double get_score(double prob, int sz) {
    return pow(10., prob / best_prob) / pow(1.5, sz);
  }
};

int FeistelLinearSolver::evaluate_cost(std::vector<int> *order) {
  /*
     Build incremental basis from each relation
     For a new relation,
     */
  std::vector<Basis> incremental_basis;

  incremental_basis.emplace_back(target->key_size());
  Scorer scorer(m_pending_entries[0].rel.cost);

  if (order) {
    for (auto &x : m_pending_entries)
      OPA_DISP0(x.rel);
  }

  int cost = 0;
  vi used(m_pending_entries.size(), 0);
  while (1) {
    int best = -1;
    int bestv = 0;
    int first = -1;
    REP (i, m_pending_entries.size()) {
      if (used[i])
        continue;

      Basis tmp = m_pending_entries[i].basis;
      tmp.add(incremental_basis.back());
      int cnd = tmp.dim() - incremental_basis.back().dim();
      if (cnd == 0)
        continue;
      if (first == -1)
        first = i;
      if (order)
        OPA_DISP("Cost for ", m_pending_entries[i].rel.cost, cnd,
                 scorer.get_score(m_pending_entries[i].rel.cost, cnd));

      if (best == -1 ||
          scorer.get_score(m_pending_entries[first].rel.cost, bestv) <
            scorer.get_score(m_pending_entries[i].rel.cost, cnd))
        best = i, bestv = cnd;
    }
    if (best == -1)
      break;
    if (order)
      OPA_DISP("best codst ", bestv);
    cost = std::max(cost, bestv);

    used[best] = 1;
    Basis next = incremental_basis.back();
    next.add(m_pending_entries[best].basis);
    incremental_basis.pb(next);
    if (order)
      order->pb(best);
  }
  cost =
    std::max(cost, int(target->key_size() - incremental_basis.back().dim()));
  if (0 && order) {
    REP (i, m_pending_entries.size())
      OPA_DISP("GOT ENTRY ", i, m_pending_entries[i].basis.dim(),
               m_pending_entries[i].rel.l, m_pending_entries[i].basis.mat());
    OPA_DISP("COSTING >> ", cost);
  }
  return cost;
}

OPA_NAMESPACE_END(opa, crypto, la)

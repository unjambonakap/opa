#include <opa/crypto/la/feistel_diff_driver.h>

using namespace std;
using namespace opa::utils;
using namespace opa::math::common;

OPA_NAMESPACE(opa, crypto, la)
RequiredDiffs FeistelDiffSolver::required_data_to_diffs(
  const RequiredData &required_data) const {
  map<BitVec, u64> diff_to_count;

  for (auto &per_round : required_data.data_per_round) {
    for (auto &rel : per_round.rels) {
      BitVec diff = feistel->range_lr_to_bitvec(rel.lin, rel.rin);
      // TODO: I need something real
      u64 need = 10;
      diff_to_count[diff] = need;
    }
  }
  RequiredDiffs res;
  for (auto &e : diff_to_count) {
    res.entries.emplace_back(e.first, e.second);
  }
  return res;
}

void FeistelDiffSolver::data_to_outrel(const RelState &s,
                                       const CipherData &data,
                                       OutRel &out_rel) const {
  out_rel.cost = s.cost;
  out_rel.out = s.l.to<u64>();

  for (auto &e : data.x_diff) {
    BitVec diff = e.ST.ST.xorz(e.ND.ST);
    if (!diff.xorz(feistel->range_lr_to_one(s.lin, s.rin)).isz())
      continue;
    OutRel::DiffEntry entry;
    entry.a1 = e.ST.ST.to<u64>();
    entry.a2 = e.ND.ST.to<u64>();
    entry.b1 = e.ST.ND.to<u64>(feistel->blk_size());
    entry.b2 = e.ND.ND.to<u64>(feistel->blk_size());
    entry.diff_res =
      e.ST.ND.xorz(e.ND.ND).to<u64>(0, feistel->blk_size()) ^ out_rel.out;
    out_rel.in_diff.push_back(entry);
  }
  OPA_CHECK(out_rel.in_diff.size() > 0, s.lin, s.rin);
}

double FeistelDiffSolver::eval_one_outrel(const OutRel &out_rel, u64 key,
                                          u64 check_mask) const {
  u32 cnt = 0;
  OPA_CHECK0(out_rel.in_diff.size() > 0);
  for (auto &e : out_rel.in_diff) {
    u64 v1 = target->fast_eval2(e.b1, key);
    u64 v2 = target->fast_eval2(e.b2, key);
    u64 r = v1 ^ v2;
    r ^= e.diff_res;
    r &= check_mask;
    cnt += r == 0;
  }
  return 1. * cnt / out_rel.in_diff.size();
}

double FeistelDiffSolver::get_score_count(const FindKeysContext *ctx, int pos,
                                          u64 key, int round_id) const {
  double res = 0.;

  u64 check_mask = m_buckets[round_id].to<u64>();
  u32 tot = 0;
  for (const auto &out_rel : *ctx->out_rels) {
    double prob = eval_one_outrel(out_rel, key, check_mask);
    res += prob;
    tot += 1;
  }
  OPA_CHECK0(!std::isnan(res));
  res /= tot;
  return res;
}

void FeistelDiffSolver::setup_solver_rels() {
  // TODO: API for that
  m_alloc.reset(Timer::froms(2));
  m_best_rels.init(20);

  fill_rels();
  m_final_rels.clear();
  m_best_rels.get_items(m_final_rels);

  int bucket_size = 4;
  int ks = feistel->round()->key_size();

  m_buckets.clear();
  for (int i = 0; i < ks; i += bucket_size) {
    int ni = min(ks, i + bucket_size);
    RangeCoverage cur;
    cur.add_interval(i, ni - i);
    m_buckets.pb(cur);
  }

  vector<Entry> final_entries;
  for (auto &bucket : m_buckets) {
    Entry cur_entry;
    // TODO: a fuckign real value
    cur_entry.rel.cost = 1.;
    auto lst = target->get_pre_key(bucket);
    cur_entry.basis.init(lst, N);
    final_entries.pb(cur_entry);
  }
  finish_rels(final_entries, final_entries.size(), m_final_rels);
}

FeistelSolver::Result FeistelDiffSolver::add_rel(const RelState &rel_) {
  const Result kDone = Result(false);
  const Result kNotDone = Result(true);
  if (rel_.lin.size() == 0 && rel_.rin.size() == 0)
    return kNotDone;

  double resv = this->test_r2(200, rel_, true);
  OPA_DISP("adding rel >> ", this->nround, resv, rel_);
  OPA_CHECK0(resv > 0.99);

  m_best_rels.maybe_add(rel_.cost, rel_);

  if (m_best_rels.size() < 3)
    return kNotDone;
  if (m_alloc.elapsed())
    return kDone;

  return kNotDone;
}

void FeistelDiffSolver::fill_rels() {
  OPA_CHECK0(this->nround > 2);
  RelState s;
  s.cost = 1;
  s.round = 0;
  Result x = rec(s, -1);
}

FeistelSolver::Result FeistelDiffSolver::rec(const RelState &s, int parid) {
  OPA_DISP0(s, parid);
  if (m_rec_cache.count(s)) {
    OPA_DISP0("CACHED");
    return Result();
  }

  int r = m_rec_cache[s] = m_rec_par.size();
  m_rec_par.pb(parid);

  if (s.round == feistel->nround() - 2) {
    RelState ns = s;
    ns.r.clear();
    Result res = add_rel(ns);
    return res;
  }

  const Relations &rels = feistel->round()->get_input_relations(true);

  for (auto &x : rels.tb) {
    RelState ns = s;
    ns.r = s.l;
    ns.l = s.r;

    ++ns.round;
    ns.cost *= x.cost.prob();
    if (s.round == 0) {
      ns.r = x.in;
      ns.rin = x.in;
    } else if (s.round == 1) {
      ns.lin.do_xor(ns.r ^ x.in);
      ns.r = x.in;
    } else {
      if (ns.r != x.in)
        continue;
    }

    ns.l.do_xor(x.out);
    Result res = rec(ns, r);
    if (!res.more)
      return res;
  }

  return Result();
}
OPA_NAMESPACE_END(opa, crypto, la)

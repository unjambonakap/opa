#include <opa/crypto/la/feistel_linear_driver.h>
#include <opa/crypto/la/feistel.h>
#include <opa/crypto/la/blocks.h>
#include <opa/math/common/Utils.h>

using namespace opa::OR;

OPA_NAMESPACE(opa, crypto, la)

void FeistelSolver::init(FeistelBlock *feistel) {
  this->init(feistel->nround(), feistel, feistel->round().get());
}

void FeistelSolver::init(int nround, FeistelBlock *feistel,
                         const CipherBlock *target) {
  Initable::init();
  N = target->key_size();
  this->nround = nround;
  this->target = target;
  this->feistel = feistel;

  if (nround > 1) {
    FeistelBlock::Params params = feistel->params();
    params.nround -= 1;
    if (params.nround <= 2) {
      m_next_solver.reset(new FeistelLinearSolver());
    } else {
      m_next_solver.reset(this->duplicate());
    }

    m_next_blk.reset(new FeistelBlock);
    m_next_blk->init(params);
    m_next_solver->init(nround - 1, m_next_blk.get(), target);
  }

  m_incremental_basis.clear();
}

const FeistelSolver::RequiredData &FeistelSolver::first_phase_find_rels() {
  if (!m_required_data) {
    m_required_data.reset(new RequiredData);
    if (m_next_solver) {
      *m_required_data = m_next_solver->first_phase_find_rels();
    }

    setup_solver_rels();

    OPA_CHECK0(m_rels.size() > 0);
    puts("FINISH RELS");
    m_required_data->data_per_round.push_back(RequiredDataPerRound{ m_rels });
    // OPA_DISP0(opa::utils::Join("\n", round_rels));
    // exit(0);
  }
  return *m_required_data;
}

u64 FeistelSolver::key_to_base(u64 k) {
  u64 res = 0;
  REP (i, N) {
    int tmp =
      opa::utils::dot(RangeCoverage::From_binary_vec(
                        m_incremental_basis.back().mat().get_row(i).tovec())
                        .to<u64>(),
                      k);
    res |= (u64)tmp << i;
  }
  return res;
}

u64 FeistelSolver::base_to_key(u64 b) {
  u64 res = 0;
  REP (i, N) {
    if (b >> i & 1) {
      u64 tmp =
        RangeCoverage::From_binary_vec(inv.get_col(i).tovec()).to<u64>();
      res ^= tmp;
    }
  }
  return res;
}

void FeistelSolver::finish_rels(const std::vector<Entry> &entries_,
                                int main_size,
                                const std::vector<RelState> &rels) {
  m_entries = entries_;
  m_rels = rels;

  for (auto &x : m_entries)
    OPA_DISP("BEFORE >> ", x.rel);

  m_steps.resize(main_size);
  OPA_DISP("Main size IS >> ", main_size);

  std::set<int> used;
  REP (i, main_size)
    used.insert(i);

  m_incremental_basis.clear();
  // empty basis
  m_incremental_basis.emplace_back(target->key_size());

  double cum_score = 1;
  int cost2 = 0;
  REP (i, main_size) {
    m_steps[i].main = i;
    auto &entry = m_entries[i];
    if (i <= 1)
      evaluate_entry_real_cost(entry);
    OPA_DISP("MAIN STEP ", i, entry.rel);

    {
      cum_score *= entry.rel.cost;
      entry.cum_score = cum_score;
    }

    Basis last = m_incremental_basis.back();
    Basis cur = entry.basis;
    Basis other(target->key_size());
    int nx = 0;
    entry.bit_pos = last.dim();
    REP (j, cur.dim()) {
      auto vj = cur.get(j);
      if (last.spans(vj)) {
        continue;
      }
      entry.used_in_basis.insert(j);
      last.add(vj);
      OPA_CHECK0(last.dim() == entry.bit_pos + entry.used_in_basis.size());
    }
    cost2 += entry.used_in_basis.size();
    m_incremental_basis.pb(last);
    OPA_CHECK(last.dim() - entry.bit_pos == entry.used_in_basis.size(),
              last.dim(), entry.bit_pos, entry.used_in_basis);

    REP (j, m_entries.size()) {
      if (!used.count(j) && cur.spans(m_entries[j].basis)) {
        cum_score *= m_entries[j].rel.cost;
        m_entries[j].cum_score = cum_score;
        used.insert(j);
        m_steps[i].checks.pb(j);
      }
    }
    // if (m_steps[i].checks.size() > 2)
    //  m_steps[i].checks.resize(2);
    OPA_DISP("CHECKER at ", m_steps[i].checks);
  }

  Basis final_basis = m_incremental_basis.back();

  std::vector<int> last_vecs;
  int prev_dim = final_basis.dim();
  REP (i, N) {
    if (final_basis.dim() == N)
      break;
    auto xx = RangeCoverage().add(i).to_mat(N);
    OPA_DISP(" >> ", xx);
    if (final_basis.spans(xx))
      continue;
    OPA_DISP("Try adding ", i);
    last_vecs.pb(i);
    final_basis.add(RangeCoverage().add(i));
  }
  OPA_CHECK_EQ0(final_basis.dim(), N);

  m_incremental_basis.pb(final_basis);
  bool ok = final_basis.mat().invert(&inv);
  // OPA_DISP0(final_basis.mat(), final_basis.dim());
  OPA_CHECK0(ok);

  // OPA_DISP("last mat >> ", final_basis.mat());
  // OPA_DISP("INV mat >> ", inv);
  // REP (i, N) {
  //  OPA_DISP("FOR ", i, inv.eval(RangeCoverage().add(i).to_binary_vec(N)));
  //}
  // OPA_DISP("last size >> ", last_vecs.size());

  REP (i, m_steps.size()) {
    auto &entry = m_entries[m_steps[i].main];
    REP (j, entry.used_in_basis.size()) {
      entry.bit_to_key.pb(
        RangeCoverage::From_binary_vec(inv.get_col(entry.bit_pos + j).tovec())
          .to<u64>());
      OPA_DISP("wanted bit ", entry.bit_pos + j);
      OPA_CHECK0(key_to_base(entry.bit_to_key[j]) ==
                 1ull << (entry.bit_pos + j));
    }
  }

  // int ncheck = 20;
  // REP (nn, ncheck) {
  //  u32 v = opa::math::common::rng() % (1ull << target->raw_input_size());
  //  u64 t1 = 0;
  //  REP (i, N) {
  //    auto r =
  //      RangeCoverage::From_binary_vec(final_basis.mat().get_row(i).tovec());
  //    u64 x = opa::utils::dot(r.to<u32>(), v);
  //    t1 |= x << i;
  //  }
  //  u32 t2 = 0;
  //  REP (i, N) {
  //    if ((t1 >> i & 1) == 0)
  //      continue;
  //    auto r = RangeCoverage::From_binary_vec(inv.get_col(i).tovec());
  //    t2 ^= r.to<u32>();
  //  }
  //  OPA_DISP0(v, t1, t2);
  //  OPA_CHECK_EQ0(t2, v);
  //}

  REP (i, m_entries.size()) {
    if (!used.count(i))
      continue;
    auto &entry = m_entries[i];

    REP (j, entry.basis.dim()) {
      if (entry.used_in_basis.count(j))
        continue;
      entry.key_to_rel.pb(
        RangeCoverage::From_binary_vec(entry.basis.get(j).tovec()).to<u64>());
    }
    int nx = entry.used_in_basis.size();
    REP (j, nx) {
      u64 mod = 0;
      REP (k, entry.key_to_rel.size()) {
        mod |= opa::utils::dot(entry.key_to_rel[k], entry.bit_to_key[j]) << k;
      }
      entry.bit_to_cachekey.pb(mod << nx | 1 << j);
    }
  }

  Entry last_entry;
  last_entry.bit_pos = prev_dim;
  OPA_CHECK0(last_vecs.size() + prev_dim == feistel->blk_size());
  REP (i, last_vecs.size()) {
    last_entry.bit_to_key.pb(
      RangeCoverage::From_binary_vec(inv.get_col(prev_dim + i).tovec())
        .to<u64>());
  }
  m_entries.pb(last_entry);

  Step last_step;
  last_step.main = m_entries.size() - 1;
  m_steps.pb(last_step);
}

void FeistelSolver::second_phase_entry(const CipherData &cipher_data,
                                       FeistelSolverResult *output_result) {

  OPA_CHECK0(output_result != nullptr);
  FeistelSolverSecondPhaseData second_phase_data;
  FeistelSolverSharedContext shared_context;
  shared_context.result = output_result;
  second_phase_data.key.init(0);
  second_phase_data.shared_context = &shared_context;
  second_phase_find_keys(cipher_data, second_phase_data);
}

void FeistelSolver::second_phase_find_keys(
  const CipherData &cipher_data,
  const FeistelSolverSecondPhaseData &second_phase_data) {
  FindKeysContext context;

  std::unique_ptr<FeistelSolverSecondPhaseData> second_phase_data_owned;

  context.cipher_data = &cipher_data;
  context.second_phase_data = &second_phase_data;
  std::vector<OutRel> out_rels;
  for (auto &rel : m_required_data->data_per_round.back().rels) {
    out_rels.emplace_back();
    data_to_outrel(rel, cipher_data, out_rels.back());
  }

  BeamSearch<RecState>::Params params;
  params.node_lim = 1000;
  params.debug = false;
  params.search_params.func = [&](Search<RecState> *self, const RecState &s) {
    this->find_keys_rec(s);
  };
  context.search.init(params);

  context.out_rels = &out_rels;
  find_keys_start(&context);
}

void FeistelSolver::find_keys_start(FindKeysContext *context) {
  context->cached_scores.clear();
  REP (i, m_steps.size()) {
    auto &entry = m_entries[m_steps[i].main];
    // OPA_DISP("step ", i, entry.cum_score, entry.rel, entry.basis.mat());
    // for (auto &x : m_steps[i].checks) {
    //  auto &e2 = entries[x];
    //  OPA_DISP("check ", e2.cum_score, e2.rel);
    //}
    // puts("");
  }
  RecState s(context, 0, 0, 0);
  context->search.add(s);
  context->search.start();
}

void FeistelSolver::find_keys_end(const RecState &s) {
  CipherData ndata;
  BitVec found_key = BitVec::From(s.vkey, feistel->round()->key_size());
  feistel->undo_last_round(*s.context->cipher_data, found_key, &ndata);
  // OPA_DISP0(found_key.to<u64>());
  // OPA_DISP0(target->context()->target_key.to<u64>());

  FeistelSolverSecondPhaseData new_second_phase_data =
    *s.context->second_phase_data;
  new_second_phase_data.key =
    found_key.concat(s.context->second_phase_data->key);

  if (this->m_next_solver) {
    this->m_next_solver->second_phase_find_keys(ndata, new_second_phase_data);
    if (new_second_phase_data.shared_context->should_stop) {
      s.context->search.stop();
    }

  } else {
    BitVec xor_key(feistel->blk_size() * 2);
    BitVec final_key = new_second_phase_data.key;

    for (auto &e : ndata.x) {
      // swap back
      e.second = e.second.extract(feistel->blk_size())
                   .concat(e.second.extract(0, feistel->blk_size()));
    }

    for (auto &x : ndata.x_diff) {
      for (auto &e : { &x.ST, &x.ND }) {
        // swap back
        e->second = e->second.extract(feistel->blk_size())
                     .concat(e->second.extract(0, feistel->blk_size()));
      }
    }

    if (feistel->mid_xor()) {
      if (ndata.x.size() > 0) {
        xor_key = ndata.x[0].first.xorz(ndata.x[0].second);
      } else {
        OPA_CHECK0(ndata.x_diff.size() > 0);
        xor_key =
          ndata.x_diff[0].first.first.xorz(ndata.x_diff[0].first.second);
      }
      final_key = final_key.concat(xor_key);
    }

    bool ok = true;
    for (auto &e : ndata.x) {
      BitVec tmp = e.first.xorz(xor_key);
      if (tmp != e.second) {
        ok = false;
        break;
      }
    }

    if (ok) {
      for (auto &x : ndata.x_diff) {
        for (auto &e : { x.ST, x.ND }) {
          BitVec tmp = e.first.xorz(xor_key);
          OPA_DISP0(tmp, e.second);
          if (tmp != e.second) {
            ok = false;
            break;
          }
        }
      }
    }

    if (ok) {
      puts("OK KAPPA");
      OPA_DISP("TRY FINAL 1 here", final_key, ok, final_key.to<u64>(), s.cost);
      OPA_DISP0(target->context()->target_key.to<u64>());
      new_second_phase_data.shared_context->result->keys.pb(final_key);
      if (new_second_phase_data.shared_context->result->keys.size() >= 2) {
        new_second_phase_data.shared_context->should_stop = true;
        s.context->search.stop();
      }
    }
  }
}

void FeistelSolver::find_keys_rec(const RecState &s) {
  if (s.round == m_steps.size()) {
    find_keys_end(s);
    return;
  }
  bool last = s.round == m_steps.size() - 1;

  auto &step = m_steps[s.round];
  int basis_pos = m_incremental_basis[s.round].dim();
  auto entry = m_entries[m_steps[s.round].main];
  int sz = entry.bit_to_key.size();
  OPA_CHECK_EQ0(entry.bit_pos, basis_pos);

  if (sz == 0) {
    RecState ns = s;
    ++ns.round;
    return find_keys_rec(ns);
  }

  std::vector<std::pair<double, RecState> > scores;

  u64 cache_key = 0;

  if (!last) {
    REP (i, entry.key_to_rel.size())
      cache_key |= (u64)opa::utils::dot(entry.key_to_rel[i], s.vkey) << i;
    cache_key <<= sz;
  }

  bool CHECKING = 0;
  u64 want = 0;
  if (CHECKING) {
    u64 ans_key = 0;
    ans_key =
      feistel->get_round_key(target->context()->target_key, this->nround - 1)
        .to<u64>();
    OPA_DISP0("TRY GO FOR KEY ", ans_key, this->nround);
    REP (j, sz) {
      auto tmp = m_incremental_basis.back().mat().get_row(basis_pos + j);
      int res = RangeCoverage()
                  .From_binary_vec(tmp.tovec())
                  .dot(RangeCoverage::From(ans_key, N));
      OPA_DISP0(res, tmp, ans_key);
      want |= (u64)res << j;
    }
    OPA_DISP0("WANT >> ", want, ans_key, key_to_base(ans_key));
  }
  // OPA_DISP("rec step ", cache_key, sz, s, basis_pos, entry.rel.cost);
  std::set<u64> lim;

  REP (i, 1 << sz) {
    bool is_want = i == want && CHECKING;
    if (CHECKING && !is_want)
      continue;
    u64 cache_key2 = cache_key;
    RecState ns = s;
    ns.vbase = s.vbase | ((u64)i << basis_pos);
    ns.vkey = s.vkey;
    ns.cost = s.cost;
    ns.round = s.round + 1;

    REP (j, sz)
      if (i >> j & 1) {
        ns.vkey ^= entry.bit_to_key[j];
        if (!last)
          cache_key2 ^= entry.bit_to_cachekey[j];
      }

    if (lim.size() && !lim.count(ns.vkey))
      continue;

    if (0 && !last && is_want) {
      u64 check_key = 0;
      int pos1 = 0;
      int pos2 = 0;
      REP (j, entry.basis.dim()) {
        u64 res = RangeCoverage()
                    .From_binary_vec(entry.basis.get(j).tovec())
                    .dot(RangeCoverage::From(ns.vkey, N));

        if (entry.used_in_basis.count(j)) {
          check_key |= res << pos1++;
        } else
          check_key |= res << (entry.used_in_basis.size() + pos2++);
      }
      OPA_CHECK_EQ0(check_key, cache_key2);
    }

    if (is_want) {
      REP (j, sz) {
        auto tmp = m_incremental_basis.back().mat().get_row(basis_pos + j);
        u64 res = RangeCoverage()
                    .From_binary_vec(tmp.tovec())
                    .dot(RangeCoverage::From(ns.vkey, N));
        // OPA_DISP0(j, tmp, ns.vkey, res, want);

        OPA_CHECK0((res) == (want >> j & 1));
      }
    }

    if (last) {
      find_keys_rec(ns);

    } else {
      double score =
        get_score(s.context, m_steps[s.round].main, ns, cache_key2, s.round);
      ns.update_score(score, entry.rel.cost);
      OPA_DISP("NS >> ", ns, score, want, i, entry.rel.cost);
      scores.pb(MP(score, ns));
    }
  }

  // OPA_DISP("filtering ", scores.size());
  if (!last) {

    if (filter_res(s.context, scores, entry)) {

      for (auto check : m_steps[s.round].checks) {
        // OPA_DISP("Filter status >>> ", scores.size(), 1 << sz);
        const auto &entry = m_entries[check];
        for (auto &e : scores) {
          u64 cache_key = 0;
          REP (j, entry.key_to_rel.size())
            cache_key |= opa::utils::dot(entry.key_to_rel[j], e.second.vkey)
                         << j;
          double score =
            get_score(s.context, check, e.second, cache_key, s.round);
          e.second.update_score(score, entry.rel.cost);
          e.first *= score;
          // OPA_DISP0(score, e.second.cost, entry.cum_score);
        }
        if (!filter_res(s.context, scores, entry))
          break;
      }
    }
    OPA_DISP("End status >>> ", scores.size(), 1 << sz);
  }

  for (auto e : scores) {
    // puts("PUsHING heRe");
    s.context->search.add(e.second);
  }
}

bool FeistelSolver::filter_res(const FindKeysContext *context,
                               std::vector<std::pair<double, RecState> > &res,
                               const Entry &entry) const {
  sort(ALL(res), std::greater<std::pair<double, RecState> >());
  // puts("RES OF FILTER");
  // for (auto &e : res)
  //  OPA_DISP0(e);
  // TODO: comment this + api to set parameters
  double tlow = -3.090232;
  double thigh = 1.644854;
  int n = context->out_rels->at(0).in.size();
  double p0 = entry.rel.cost / 2 + 0.5;
  double tmp = sqrt(std::max<double>(0., p0 * (1 - p0) / n));
  double vlow = tlow * tmp + p0;
  double vhigh = thigh * tmp + p0;
  OPA_DISP0(tlow, p0, tmp, vlow, vhigh, n);
  // vlow *= 0.9;

  for (int i = 0; i < res.size(); ++i) {
    double v0 = res[i].first;
    OPA_DISP0(v0, vlow, res[i].first, p0, res[i].second.vkey);
    if (i != -1 && v0 + 1e-7 < vlow) {
      res.resize(i);
      break;
    }
  }
  return res.size() > 2;
}

double FeistelSolver::get_score(FindKeysContext *context, int data_pos,
                                const RecState &s, u64 cache_key,
                                int round_id) const {
  bool has = true;

  if (!context->cached_scores.count(MP(data_pos, cache_key)))
    has = false;
  double &res = context->cached_scores[MP(data_pos, cache_key)];

  if (!has) {
    res = this->get_score_count(context, data_pos, s.vkey, round_id);
    OPA_CHECK0(!std::isnan(res));
  }
  return res;
}

void FeistelSolver::evaluate_entry_real_cost(Entry &entry) {
  double res = test_r2(200, entry.rel, m_diff_mode);
  // OPA_DISP0("evaluate real antry cost", res);
  entry.rel.cost = res;
}

// OPA_CLOUDY_REGISTER_BASE_TMPL(testr2mapper, opa::threading::MapJob,
//                              opa::utils::DoubleRes, int,
//                              SPTR(const opa::crypto::la::FeistelBlock),
//                              opa::crypto::la::RelState);

double FeistelSolver::test_r2(int nr, const RelState &rel, bool diff) const {
  int sz = 3;
  std::vector<opa::utils::DoubleRes> res;
  REP (i, sz) { res.emplace_back(test_r2_internal(nr, rel, diff)); }
  double mean = 0;
  for (auto x : res)
    mean += x.v;
  mean /= res.size();
  // OPA_DISP("TEST_R2 >> ", rel.cost, mean);
  return mean;
}

double FeistelSolver::test_r2_internal(int nr, const RelState &rel,
                                       bool diff) const {
  check_init();
  BitVec kv = BitVec::rand(feistel->key_size());
  CipherData data;

  if (diff) {
    feistel->gen_rand_diff(feistel->range_lr_to_bitvec(rel.lin, rel.rin), nr,
                           kv, &data);

  } else {
    data = feistel->gen_rand(nr, kv);
  }
  OutRel d2;
  data_to_outrel(rel, data, d2);
  double res = feistel->round()->test_outrel(
    feistel->get_round_key(kv, this->nround - 1), d2);
  return res;
}

OPA_NAMESPACE_END(opa, crypto, la)

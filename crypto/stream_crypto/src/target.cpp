#include <opa/crypto/stream/target.h>

using namespace opa::crypto::la;
using namespace opa::math::common;
using namespace std;

OPA_NAMESPACE(opa, crypto, stream)

void SboxLfsr::init(const Params &params) {
  opa::utils::Initable::init();
  m_params = params;
  m_n = sbox().input_size();
}

u32 SboxLfsr::get() {
  if (!m_bs.has_more()) {
    int v = 0;
    REP (i, m_n) {
      int u = 0;
      for (auto &lfsr : lfsrs()) {
        u ^= lfsr.get_next();
      }
      v |= u << i;
    }
    m_bs.push_multiple(sbox().get(v), sbox().output_size());
  }
  return m_bs.get();
}
void SboxLfsr::set_rand() {
  reset();
  for (auto &lfsr : lfsrs()) {
    lfsr.set_rand();
  }
}

opa::math::common::BitVec SboxLfsr::get_state() const {
  int m_input_len = 0;
  for (auto &lfsr : lfsrs())
    m_input_len += lfsr.size();

  opa::math::common::BitVec res;
  res.init(m_input_len);
  int pos = 0;
  for (auto &lfsr : lfsrs()) {
    REP (j, lfsr.size()) { res.set(pos++, lfsr.state >> j & 1); }
  }
  return res;
}
void SboxLfsr::reset() { m_bs.reset(); }

void SboxLfsrSolver::test_rel(const la::Relation &rel) {
  /*
  int ntry = 100000;
  int cnt = 0;
  int pre = 100 / no * no;
  REP (tt, ntry) {
    int x = 0;
    m_params.lfsr_spec.set_rand();

    REP (lfsrid, m_params.lfsr_spec.lfsrs().size()) {
      auto &lfsr = m_params.lfsr_spec.lfsrs()[lfsrid];
      for (auto &i : rel.in.all()) {
        REP (j, lfsr.size())
          x ^= m_tb[lfsrid].get(j, pre / no * ni + i);
      }
    }

    REP (i, pre)
      m_params.lfsr_spec.get();
    REP (k, no) {
      int v = m_params.lfsr_spec.get();
      if (rel.out.in(k))
        x ^= v;
    }
    cnt += x == 1;
  }
  OPA_DISP(">>> %d\n", abs(cnt - ntry / 2), m_input_len, m_output_len);
  */
}

MulMatrixF2<false> analyse_lfsr(const LFSR_GF2_small &lfsr, int N) {
  MulMatrixF2<false> res;
  int n = lfsr.size();
  res.init(N, n);
  REP (ib, n) {
    LFSR_GF2_small tmp = lfsr;
    tmp.set_state(1ull << ib);
    REP (j, N) { res.set(j, ib, tmp.get_next()); }
  }

  return res;
}

void SboxLfsrSolver::init(const Params &params) {
  m_params = params;
  m_nlfsr = params.lfsr.lfsrs().size();
  m_ni = params.lfsr.sbox().input_size();
  m_no = params.lfsr.sbox().output_size();

  m_input_pos.pb(0);
  REP (i, m_nlfsr)
    m_input_pos.pb(m_input_pos.back() + params.lfsr.lfsrs()[i].size());

  m_input_len = m_input_pos.back();
  OPA_DISP("input len >> ", m_input_len, m_params.ans_input.size(), m_nlfsr);
}

void SboxLfsrSolver::load_plan(SolverPlan *plan) {
  m_plan = plan;
  configure_output_bits(m_plan->m_known_data);
}

std::vector<la::Relation> SboxLfsrSolver::get_used_relations() const {
  const auto &sbox_desc = m_params.lfsr.sbox();
  SboxBlock blk;
  blk.init(sbox_desc);
  auto relations = blk.get_relations().tb;

  sort(ALL(relations), la::Relation::CmpVecByCost());
  reverse(ALL(relations));

  int n = sbox_desc.input_size() + sbox_desc.output_size();
  OPA_DISP("got basis ", n);
  Basis basis(n);
  vector<la::Relation> used_rels;
  for (auto &rel : relations) {
    if (rel.cost.prob() < relations[0].cost.prob() / 2)
      break;
    OPA_DISP("trying to add rel", rel.in, rel.out, rel.cost.prob());
    RangeCoverage vec(n);
    vec.add(rel.in);
    vec.add(rel.out.shift(sbox_desc.input_size()));
    if (basis.add(vec)) {
      puts("pushing");
      used_rels.pb(rel);
    }
  }
  return used_rels;
}

void SboxLfsrSolver::configure_output_bits(const std::set<int> &data) {
  OPA_CHECK0(m_output_bit_list.size() == 0);

  for (auto &e : data)
    m_output_bit_list.pb(e);
  sort(ALL(m_output_bit_list));

  m_output_bit_rmp.clear();
  REP (i, m_output_bit_list.size())
    m_output_bit_rmp[m_output_bit_list[i]] = i;
}

void SboxLfsrSolver::load_plan_from_data(const std::set<int> &known_data) {
  m_own_plan.reset(new SolverPlan);
  int output_len = known_data.size();
  m_plan = m_own_plan.get();
  auto used_rels = get_used_relations();

  m_plan->init(m_input_len, known_data);
  RelsStore &store = m_plan->m_rels_store;

  configure_output_bits(known_data);

  int end = m_output_bit_list.back() + 1;
  REP (i, m_nlfsr) {
    OPA_DISP0("ANALYSE", i, end);
    m_tb.pb(analyse_lfsr(m_params.lfsr.lfsrs()[i], end));
  }

  int n2 = (m_output_bit_list.size() + 1 + m_no - 1) / m_no;
  REP (i, used_rels.size()) {
    auto &rel = used_rels[i];
    store.add_new_typ(rel.cost.bias());

    for (int i = 0; i < n2 / m_no; i += 1) {
      BitVec cur(m_input_len + output_len);
      bool ok = true;
      for (auto &out : rel.out.all()) {
        int oid = i * m_no + out;
        if (!m_output_bit_rmp.count(oid)) {
          ok = false;
          break;
        }
        cur.toggle(m_input_len + m_output_bit_rmp[oid]);
      }
      if (!ok)
        continue;

      for (auto &in : rel.in.all()) {
        REP (j, m_nlfsr) {
          REP (k, m_input_pos[j + 1] - m_input_pos[j]) {
            cur.xorv(m_input_pos[j] + k, m_tb[j].get(i * m_ni + in, k));
          }
        }
      }
      store.add_new_rel(SolverRel(cur, rel.c));
    }
    OPA_DISP("pushed type >> has size ", i, store.get_num_rels_by_type(i));
  }
  m_plan->build_plan();
}

BitVec SboxLfsrSolver::solve(const ObservedData &obs_data) {
  vector<int> obits;
  for (auto &e : m_output_bit_list) {
    OPA_CHECK0(obs_data.count(e));
    obits.pb(obs_data.find(e)->ND);
  }
  m_plan->m_rels_store.update_rels_store(obits);

  if (m_params.check) {
    OPA_DISP0(m_plan->m_rels_store.full_rels.size());
    int cnt = 0;
    REP (i, m_plan->m_rels_store.get_num_rels_by_type(0)) {
      auto &r = m_plan->m_rels_store.full_rels[i];
      int v = r.c;
      REP (j, obits.size())
        v ^= obits[j] & r.v.get(m_input_len + j);
      REP (j, m_input_len)
        v ^= m_params.ans_input.get(j) & r.v.get(j);
      cnt += v == 0;
    }

    int cnt2 = 0;
    REP (i, m_plan->m_rels_store.get_num_rels_by_type(0)) {
      auto &r = m_plan->m_rels_store.rels[i];
      int v = r.c;
      REP (j, m_input_len)
        v ^= m_params.ans_input.get(j) & r.v.get(j);
      cnt2 += v == 0;
    }
    OPA_DISP0(cnt, cnt2);
  }

  Solver solver;
  Solver::Params params;
  params.check = m_params.check;
  params.ans_input = m_params.ans_input;
  params.input_len = m_input_len;
  params.plan = m_plan;

  if (m_params.check) {
    check_key(m_params.ans_input, obs_data);
  }

  solver.init(params);
  BitVec key = solver.solve();

  check_key(key, obs_data);
  return key;
}

void SboxLfsrSolver::check_key(const BitVec &key,
                               const ObservedData &obs_data) const {
  int pos = 0;
  SboxLfsr target = m_params.lfsr;
  target.reset();
  for (auto &lfsr : target.lfsrs()) {
    u64 state = 0;
    REP (j, lfsr.n) { state |= u64(key.get(pos++)) << j; }
    lfsr.set_state(state);
  }

  int end = m_output_bit_list.back() + 1;
  REP (j, end) {
    int cur = target.get();
    if (obs_data.count(j)) {
      OPA_CHECK(obs_data.find(j)->ND == cur, "fail at ", j);
    }
  }
}

void get_sbox_lfsr_plan(const SboxLfsr &lfsr, const std::set<int> &known_data,
                        SolverPlan *plan) {}
void solve_sbox_lfsr(const SboxLfsr &lfsr, const ObservedData &obs_data,
                     const SolverPlan &plan) {}

OPA_NAMESPACE_END(opa, crypto, stream)

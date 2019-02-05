#include <opa/crypto/stream/utils.h>
#include <opa/math/common/bignum.h>
#include <opa/math/common/rng.h>
#include <opa/math/common/stats.h>
#include <opa/threading/auto_job.h>
#include <opa/threading/runner.h>

using namespace std;
using namespace opa::math::common;
DEFINE_bool(quadrisection_local, false, "");

OPA_NAMESPACE(opa, crypto, stream)

u64 get_count_for_sizes(const std::vector<pii> &take,
                        const std::vector<int> &sizes) {
  bignum score = 1;
  for (auto &e : take) {
    score *= nchoosek(sizes[e.ST], e.ND);
  }
  return score.fitsu64() ? score.getu64() : -1;
}

bool OutputRel::add(u64 nkey, MitmRel &res) {
  // sort(ALL(res));
  // REP (u, res.size() - 1) {
  //  if (res[u] == res[u + 1]) {
  //    return;
  //  }
  //}
  if (nwant == -1 || out->size() < nwant) {
    out->pb(MP(nkey, res));
    return true;
  }
  return false;
}

class RelFinder {
public:
  RelFinder(const RelsStore &rels_store, const vector<int> &ids,
            const KeyRelsList &key_rels, MidRelDesc &out_rels)
      : m_rels_store(rels_store), ids(ids), key_rels(key_rels),
        out_rels(out_rels) {
    state.resize(ids.size());
    state_key.resize(ids.size() + 1, 0);
  }

  void go() {
    rec(0);
    sort(ALL(out_rels));
  }

  void rec(int pos, int last = -1) {
    if (pos == 0 || ids[pos] != ids[pos - 1]) last = -1;

    if (pos == ids.size()) {
      out_rels.pb(MP(state_key[pos], state));
      return;
    }

    int nx = m_rels_store.get_num_rels_by_type(ids[pos]);
    FOR (i, last + 1, nx) {
      u16 id = m_rels_store.get_pos_by_type(ids[pos], i);
      state_key[pos + 1] = state_key[pos] ^ key_rels[id];
      state[pos] = id;
      rec(pos + 1, i);
    }
  }

  FinalRel state;
  vector<KeyRel> state_key;

  const vector<int> &ids;
  MidRelDesc &out_rels;
  const KeyRelsList &key_rels;
  const RelsStore &m_rels_store;
};

class CutPlanHelper {
public:
  CutPlan compute_cut(const vector<int> &sizes,
                      const vector<pii> &take_per_bucket) {
    m_sizes = sizes;
    m_take_per_bucket = take_per_bucket;

    best_plan = CutPlan();
    CutPlan plan;
    plan.take_left.resize(take_per_bucket.size());
    plan.take_right.resize(take_per_bucket.size());

    find_cut_plan_rec(0, plan);
    best_plan.can_left = get_cut_score_one_side(best_plan.take_left);
    best_plan.can_right = get_cut_score_one_side(best_plan.take_right);
    best_plan.to_ids();
    return best_plan;
  }

  void find_cut_plan_rec(int pos, CutPlan &plan) {
    if (pos == m_take_per_bucket.size()) {
      u64 best_score = get_cut_score(best_plan);
      u64 cur_score = get_cut_score(plan);
      if (best_score > cur_score) {
        best_plan = plan;
      }
      return;
    }

    int w = m_take_per_bucket[pos].ST;
    REP (i, m_take_per_bucket[pos].ND + 1) {
      plan.take_left[pos] = MP(w, i);
      plan.take_right[pos] = MP(w, m_take_per_bucket[pos].ND - i);
      find_cut_plan_rec(pos + 1, plan);
    }
  }

  u64 get_cut_score_one_side(const vector<pii> &take) {
    return get_count_for_sizes(take, m_sizes);
  }

  u64 get_cut_score(const CutPlan &plan) {
    if (plan.take_left.size() == 0) return -1;
    return max(get_cut_score_one_side(plan.take_left),
               get_cut_score_one_side(plan.take_right));
  }

  CutPlan best_plan;
  vector<int> m_sizes;
  vector<pii> m_take_per_bucket;
};

void StepHelper::solve(u64 lim, int type) {
  OPA_CHECK0(m_params.zero_ub <= 64);

  CutPlanHelper helper;
  CutPlan cut_plan =
    helper.compute_cut(m_params.rels_store->get_sizes(), m_params.target->tb);
  OPA_DISP("Procing plan ", cut_plan.left, cut_plan.right, cut_plan.can_left,
           cut_plan.can_right, m_params.target->tb);

  const u64 mitm_lim = 1e6;
  MidRelDesc res;
  bool want_quad = max(cut_plan.can_left, cut_plan.can_right) >= mitm_lim;
  want_quad |= type == 0;
  want_quad &= type != 1;

  if (want_quad) {
    solve_quadrisection(cut_plan, res);
  } else {
    solve_mitm(cut_plan, lim, res);
  }

  for (auto &rel : res) {
    m_params.out_desc->rels.pb(rel.ND);
  }
  OPA_DISP0("End of plan, got ", m_params.out_desc->rels.size());
}

void mitm_proc(const MitmEntryList &tleft, const MitmEntryList &tright,
               OutputRel &out, u64 mask, bool is_sym) {
  OPA_DISP("mitm proc ", tleft.size(), tright.size());
  for (int il = 0, ir = 0; il < tleft.size() && ir < tright.size();) {
    for (; il < tleft.size() && (tleft[il].ST & mask) < (tright[ir].ST & mask);
         ++il) {
    }
    for (; ir < tright.size() && (tright[ir].ST & mask) < (tleft[il].ST & mask);
         ++ir) {
    }

    int nl, nr;
    // no need to mask for second
    for (nl = il;
         nl < tleft.size() && (tleft[nl].ST & mask) == (tright[ir].ST & mask);
         ++nl)
      ;
    for (nr = ir;
         nr < tright.size() && (tright[nr].ST & mask) == (tleft[il].ST & mask);
         ++nr)
      ;

    if (nl != il && nr != ir) {
      FOR (i, il, nl) {
        FOR (j, ir, is_sym ? i : nr) {
          MitmRel res;

          if (is_sym && tleft[i].ND == tright[j].ND) continue;
          res.ST = tleft[i].ND.ST | tright[j].ND.ST;
          res.ND = tleft[i].ND.ND | tright[j].ND.ND;
          if (!out.add(tleft[i].ST ^ tright[j].ST, res)) return;
        }
      }
    }
    il = nl;
    ir = nr;
  }
}

inline u32 get_mitm_rel_elem(const MitmRel &x, u32 type) {
  if (type == 0)
    return u32(x.ST);
  else if (type == 1)
    return u32(x.ST >> 32);
  else if (type == 2)
    return u32(x.ND);
  else if (type == 3)
    return u32(x.ND >> 32);
  return 0;
}

inline void set_mitm_rel_elem(MitmRel &x, u32 type, u32 v) {
  if (type == 0)
    x.ST |= v;
  else if (type == 1)
    x.ST |= u64(v) << 32;
  else if (type == 2)
    x.ND |= v;
  else if (type == 3)
    x.ND |= u64(v) << 32;
}

struct FrelToMitmRelMap {
  std::map<FinalRel, int> mp;
  vector<FinalRel> rmp;

  // first elem is empty vec
  FrelToMitmRelMap() { get({}); }
  int get(const FinalRel &rel) {
    if (!mp.count(rel)) {
      mp[rel] = rmp.size();
      rmp.pb(rel);
    }
    return mp[rel];
  }

  const FinalRel &rget(int id) const { return rmp[id]; }
};

void midrel_to_mitmrel(const MidRelDesc &mid_rels, MitmEntryList &mitm_rels,
                       int type, FrelToMitmRelMap &mp) {
  for (auto &ie : mid_rels) {
    MitmRel rel;
    set_mitm_rel_elem(rel, type, mp.get(ie.ND));
    mitm_rels.emplace_back(ie.ST, rel);
  }
}

void mitmrels_to_midrels(const MitmEntryList &mitm_rels, MidRelDesc &mid_rels,
                         const vector<FinalRel> &mp, int nentries) {
  OPA_CHECK0(nentries == 2 || nentries == 4);
  for (auto &e : mitm_rels) {
    FinalRel rel;
    rel.insert(rel.end(), ALL(mp[get_mitm_rel_elem(e.ND, 0)]));
    rel.insert(rel.end(), ALL(mp[get_mitm_rel_elem(e.ND, 1)]));
    if (nentries == 4) {
      rel.insert(rel.end(), ALL(mp[get_mitm_rel_elem(e.ND, 2)]));
      rel.insert(rel.end(), ALL(mp[get_mitm_rel_elem(e.ND, 3)]));
    }

    sort(ALL(rel));
    bool ok = 1;
    REP (i, rel.size() - 1) {
      if (rel[i] == rel[i + 1]) {
        ok = 0;
        break;
      }
    }
    if (!ok) continue;

    mid_rels.emplace_back(e.ST, rel);
  }

  sort(ALL(mid_rels));
  mid_rels.resize(unique(ALL(mid_rels)) - mid_rels.begin());
}

struct QuadrisectionTaskParams : public opa::utils::ProtobufParams {
  u64 mask;
  OPA_TGEN_IMPL(mask);
};

struct QuadrisectionTask
    : public opa::threading::AutoCollectJob<QuadrisectionTaskParams,
                                            MidRelEntry> {
  struct Params : public opa::utils::ProtobufParams {
    MitmEntryList a, b, c, d;
    std::vector<FinalRel> rels;

    struct Data : public opa::utils::ProtobufParams {
      int t;
      int nbits;
      bool left_sym;
      bool right_sym;
      bool sym;
      OPA_TGEN_IMPL(t, nbits, left_sym, right_sym, sym);
    } data;
    OPA_TGEN_IMPL(a, b, c, d, rels, data);
  };

  void init(const Params &params) { m_params = params; }

  virtual void auto_worker_do_work(const QuadrisectionTaskParams &data,
                                   MidRelDesc &out_res) override {
    int nshift = m_params.data.nbits - m_params.data.t;
    u64 M = data.mask << nshift;
    u64 mask_t = ((1ull << m_params.data.t) - 1) << nshift;

    MitmEntryList b2, d2;
    for (auto &e : m_params.b) b2.pb(MP(e.ST ^ M, e.ND));
    for (auto &e : m_params.d) d2.pb(MP(e.ST ^ M, e.ND));
    sort(ALL(b2));
    sort(ALL(d2));
    MitmEntryList l1, l2;
    OutputRel r1(&l1), r2(&l2);

    mitm_proc(m_params.a, b2, r1, mask_t,
              m_params.data.left_sym && M == 0); // 20 %

    mitm_proc(m_params.c, d2, r2, mask_t, m_params.data.right_sym && M == 0);
    // taking 40-% at this point
    if (m_params.data.left_sym) {
      for (int i = 0; i < l1.size(); ++i) {
        if (get_mitm_rel_elem(l1[i].ND, 0) == get_mitm_rel_elem(l1[i].ND, 1)) {
          swap(l1[i], l1.back());
          l1.pop_back();
          --i;
        }
      }
    }

    if (m_params.data.right_sym) {
      for (int i = 0; i < l2.size(); ++i) {
        if (get_mitm_rel_elem(l2[i].ND, 2) == get_mitm_rel_elem(l2[i].ND, 3)) {
          swap(l2[i], l2.back());
          l2.pop_back();
          --i;
        }
      }
    }

    sort(ALL(l1));
    sort(ALL(l2));
    l1.resize(unique(ALL(l1)) - l1.begin());
    l2.resize(unique(ALL(l2)) - l2.begin());

    MitmEntryList mitm_list;
    OutputRel out_rel(&mitm_list);
    mitm_proc(l1, l2, out_rel, -1, m_params.data.sym);

    if (m_params.data.sym) {
      for (int i = 0; i < mitm_list.size(); ++i) {
        bool ok = 1;
        REP (i1, 4) {
          int x1 = get_mitm_rel_elem(mitm_list[i].ND, i1);
          // skip empty vec from check
          if (x1 == 0) continue;
          REP (i2, i1) {
            if (x1 == get_mitm_rel_elem(mitm_list[i].ND, i2)) {
              ok = 0;
              break;
            }
          }
          if (!ok) break;
        }
        if (!ok) {
          swap(mitm_list[i], mitm_list.back());
          mitm_list.pop_back();
          --i;
        }
      }
    }
    mitmrels_to_midrels(mitm_list, out_res, m_params.rels, 4);
    OPA_DISP("on mask ", M, 1 << m_params.data.t, l1.size(), l2.size(),
             m_params.data.nbits, out_res.size(), b2.size(), d2.size());
  }

  virtual void auto_server_get_work() override { run(false); }

  void run(bool local = true) {
    REP (M2, 1 << m_params.data.t) {
      QuadrisectionTaskParams task;
      task.mask = M2;
      if (local) {
        auto_worker_do_work(task, m_res_list);
      } else {
        bool more;
        cb()(task, more);
        if (!more) break;
      }
    }
  }

  OPA_TGEN_IMPL(m_params);
  OPA_CLOUDY_JOB_DECL;

  Params m_params;
};

OPA_CLOUDY_REGISTER_BASE(QuadrisectionTask);
OPA_CLOUDY_JOB_IMPL(QuadrisectionTask);
void StepHelper::solve_quadrisection(const CutPlan &cut_plan, MidRelDesc &res) {
  CutPlanHelper helper;
  CutPlan left_plan =
    helper.compute_cut(m_params.rels_store->get_sizes(), cut_plan.take_left);
  CutPlan right_plan =
    helper.compute_cut(m_params.rels_store->get_sizes(), cut_plan.take_right);
  OPA_CHECK0(cut_plan.left.size() > 0 && cut_plan.right.size() > 0);

  const u64 max_sz = 1e5;
  const int t_pw = max(cut_plan.can_left, cut_plan.can_right) / max_sz;
  int t = log2_high_bit(max(1, t_pw - 1));
  if (m_params.dispatcher) {
    t = max<int>(t, log2_high_bit(m_params.dispatcher->nthread()));
  }

  OPA_DISP("gogo quadrisection ", cut_plan.left, cut_plan.right, left_plan.left,
           left_plan.right, right_plan.left, right_plan.right, t,
           m_params.zero_ub);

  QuadrisectionTask::Params params;
  MidRelDesc a, b, c, d;
  FrelToMitmRelMap frel_to_mitmrel;

  get_all_rels(left_plan.left, a);
  get_all_rels(left_plan.right, b);
  get_all_rels(right_plan.left, c);
  get_all_rels(right_plan.right, d);

  int pos = 0;
  for (auto &e : { MP(&a, &params.a), MP(&b, &params.b), MP(&c, &params.c),
                   MP(&d, &params.d) }) {
    midrel_to_mitmrel(*e.ST, *e.ND, pos, frel_to_mitmrel);
    sort(ALL(*e.ND));
    ++pos;
  }

  if (params.a.size() < params.b.size()) swap(params.a, params.b);
  if (params.c.size() < params.d.size()) swap(params.c, params.d);

  params.rels = frel_to_mitmrel.rmp;
  params.data.nbits = m_params.zero_ub;
  params.data.t = t;
  params.data.left_sym = left_plan.is_sym();
  params.data.right_sym = right_plan.is_sym();
  params.data.sym = cut_plan.is_sym();

  OPA_DISP("buckets", params.a.size(), params.b.size(), params.c.size(),
           params.d.size(), t);
  UPTR(QuadrisectionTask)
  task(opa::threading::Runner::GetJob<QuadrisectionTask>(
    QuadrisectionTask::JobName));
  task->init(params);

  if (!FLAGS_quadrisection_local && m_params.dispatcher) {
    m_params.dispatcher->process_job(*task);
  } else {
    task->run();
  }
  puts("dispatcher done");

  res = task->res_list();
  sort(ALL(res));
  res.resize(unique(ALL(res)) - res.begin());

  OPA_DISP("got rels", res.size(), params.data.left_sym, params.data.right_sym,
           params.data.sym);
  OPA_CHECK0(res.size() > 0);
}

void StepHelper::solve_mitm(const CutPlan &cut_plan, const u64 lim,
                            MidRelDesc &res) {
  OPA_DISP("solve mitm ", cut_plan.left, cut_plan.right, cut_plan.can_left,
           cut_plan.can_right);

  MidRelDesc tleft;
  MidRelDesc tright;
  u64 want_rels = lim * (1ull << m_params.zero_ub);
  u64 nleft = min(cut_plan.can_left,
                  max(want_rels / cut_plan.can_right, u64(sqrt(want_rels))));
  u64 nright = min(cut_plan.can_right,
                   max(want_rels / cut_plan.can_left, u64(sqrt(want_rels))));

  puts("");
  OPA_DISP0(nleft, nright, want_rels, cut_plan.can_left, cut_plan.can_right,
            lim, m_params.zero_ub, cut_plan.left, cut_plan.right);

  if (lim == 0 || nleft >= cut_plan.can_left / 2) {
    get_all_rels(cut_plan.left, tleft);
  } else {
    get_lim_rels(cut_plan.left, tleft, nleft);
  }

  if (lim == 0 || nright >= cut_plan.can_right / 2) {
    get_all_rels(cut_plan.right, tright);
  } else {
    get_lim_rels(cut_plan.right, tright, nright);
  }

  MitmEntryList l1, l2, mitm_res;
  OutputRel rel_res(&mitm_res, lim);
  FrelToMitmRelMap frel_to_mitmrel;
  midrel_to_mitmrel(tleft, l1, 0, frel_to_mitmrel);
  midrel_to_mitmrel(tright, l2, 1, frel_to_mitmrel);

  mitm_proc(l1, l2, rel_res, -1, cut_plan.is_sym());
  OPA_DISP0("DOOONE");
  mitmrels_to_midrels(mitm_res, res, frel_to_mitmrel.rmp, 2);
  OPA_DISP0("done conv1");
}

void StepHelper::get_lim_rels(const vector<int> &ids, MidRelDesc &out_rels,
                              u64 nrels) {
  for (u64 relid = 0; relid < nrels; ++relid) {
    FinalRel state(ids.size());
    u64 key = 0;
    REP (pos, ids.size()) {
      int w = rng() % m_params.rels_store->get_num_rels_by_type(ids[pos]);
      state[pos] = w;
      key ^=
        (*m_params.key_rels)[m_params.rels_store->get_pos_by_type(ids[pos], w)];
    }
    out_rels.pb(MP(key, state));
  }
  sort(ALL(out_rels));
}

void StepHelper::get_all_rels(const std::vector<int> &ids,
                              MidRelDesc &out_rels) {
  RelFinder finder(*m_params.rels_store, ids, *m_params.key_rels, out_rels);
  finder.go();
}

OPA_NAMESPACE_END(opa, crypto, stream)

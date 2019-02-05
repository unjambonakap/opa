#pragma once
#include <opa/crypto/stream/common.h>
#include <opa/math/common/fast_gf2.h>
#include <opa/threading/dispatcher.h>

OPA_NAMESPACE(opa, crypto, stream)

typedef u64 KeyRel;
typedef std::vector<KeyRel> KeyRelsList;
typedef std::vector<u32> FinalRel;
typedef std::pair<u64, u64> MitmRel;
typedef std::pair<u64, MitmRel> MitmEntry;
typedef std::vector<MitmEntry> MitmEntryList;

typedef std::pair<u64, FinalRel> MidRelEntry;
typedef std::vector<MidRelEntry> MidRelDesc;
u64 get_count_for_sizes(const std::vector<pii> &take,
                        const std::vector<int> &sizes);

inline double bias_to_prob(double bias) { return (1 + bias) / 2; }
inline double prob_to_bias(double prob) { return 2 * prob - 1; }

inline static u64 get_key(const opa::math::common::BitVec &bv, int rem) {
  u64 res = 0;
  // OPA_CHECK0(rem <= 32);
  REP (i, rem)
    res |= u64(bv.get(i)) << i;
  return res;
}

struct SolverRel;
struct SolverRelParams : public opa::utils::ProtobufParams {
  opa::math::common::BitVecRepr v;
  bool c;
  SolverRelParams() {}
  SolverRelParams(const SolverRel &rel);
  OPA_TGEN_IMPL(v, c);
};

struct SolverRel {
  opa::math::common::BitVec v;
  bool c;
  SolverRel() {}
  SolverRel(const SolverRelParams &params) {
    v = params.v.to_bitvec();
    c = params.c;
  }

  SolverRel(int len) { init(len); }
  SolverRel(const opa::math::common::BitVec &v, bool c) : v(v), c(c) {}

  void init(int len) {
    v.init(len);
    c = 0;
  }

  bool check(const opa::math::common::BitVec &key) const {
    return key.dot(v) == c;
  }

  SolverRel &sadd(const SolverRel &other) {
    v.sxorz(other.v);
    c ^= other.c;
    return *this;
  }

  SolverRel add(const SolverRel &other) const {
    SolverRel res;
    res.v = v.xorz(other.v);
    res.c = c ^ other.c;
    return res;
  }
};

inline SolverRelParams::SolverRelParams(const SolverRel &rel) {
  c = rel.c;
  v = rel.v.to_repr();
}

struct SolverState {
  std::vector<pii> tb;
  double cost;
  SolverState(const std::vector<pii> &tb, double cost) {
    this->tb = tb;
    this->cost = cost;
  }
  SolverState() {}
  bool operator<(const SolverState &s) const { return cost < s.cost; }
};

struct StepDescription : public opa::utils::ProtobufParams {
  int n_bruteforce = -1;
  std::vector<int> fixed_input;
  std::vector<FinalRel> rels;
  u32 get_fixed_key_from_bitvec(const opa::math::common::BitVec &key) const;

  void update_bitvec_from_fixed_key(u32 fixed_key,
                                    opa::math::common::BitVec &out_key) const;
  OPA_TGEN_IMPL(n_bruteforce, fixed_input, rels);
};

class RelsStore : public opa::utils::ProtobufParams {
public:
  void init(int input_len) { this->input_len = input_len; }

  inline const SolverRel &get_rel_by_type(int typ, int pos) const {
    return rels[rels_type_pos[typ] + pos];
  }

  inline int get_num_rels_by_type(int typ) const {
    return rels_type_length[typ];
  }

  inline int get_pos_by_type(int typ, int pos) const {
    return rels_type_pos[typ] + pos;
  }

  int get_num_rels() const {
    return rels_type_pos.back() + rels_type_length.back();
  }
  int get_num_types() const { return type_bias.size(); }

  double get_type_bias(int typ) const { return type_bias[typ]; }

  KeyRelsList to_key_rels(int rem) const {
    KeyRelsList key_rels(get_num_rels());
    REP (i, key_rels.size()) { key_rels[i] = get_key(rels[i].v, rem); }
    return key_rels;
  }

  void add_new_typ(double bias) {
    type_bias.pb(bias);
    rels_type_pos.pb(rels_type_length.size() == 0 ? 0
                                                  : rels_type_length.back());
    rels_type_length.pb(0);
  }

  int get_rel_type(int relid) const {
    REP (typ, rels_type_length.size()) {
      relid -= rels_type_length[typ];
      if (relid < 0) return typ;
    }
    OPA_CHECK0(false);
  }

  void add_new_rel(const SolverRel &rel) {
    SolverRel rel_only_input;
    rel_only_input.init(input_len);
    REP (i, input_len)
      rel_only_input.v.set(i, rel.v.get(i));
    rel_only_input.c = rel.c;

    rels.pb(rel_only_input);
    OPA_CHECK(rels.size() <= 1ull << 32, "Too many rels, id has to fit on u32");
    full_rels.pb(rel);
    rels_type_length.back() += 1;
  }

  void update_rels_store(const std::vector<int> &obs_data) {
    REP (i, full_rels.size()) {
      int c = full_rels[i].c;
      REP (j, obs_data.size()) {
        c ^= obs_data[j] & full_rels[i].v.get(input_len + j);
      }
      rels[i].c = c;
    }
  }

  double get_score(const opa::math::common::BitVec &bv) const {
    double score = 0;
    for (auto &x : rels) {
      score += x.check(bv);
    }
    return score;
  }

  std::vector<int> get_sizes() const {
    std::vector<int> res;
    REP (i, get_num_types())
      res.pb(get_num_rels_by_type(i));
    return res;
  }

  std::vector<SolverRel> rels;
  std::vector<SolverRel> full_rels;

  std::vector<int> rels_type_pos;
  std::vector<int> rels_type_length;
  std::vector<double> type_bias;
  int input_len = -1;

  mutable std::vector<SolverRelParams> rels_params;
  mutable std::vector<SolverRelParams> full_rels_params;
  virtual void after_load() {
    for (auto &rel : rels_params) {
      rels.emplace_back(rel);
    }
    for (auto &rel : full_rels_params) {
      full_rels.emplace_back(rel);
    }
    rels_params.clear();
    full_rels_params.clear();
  }

  virtual void before_store() const {
    for (auto &rel : rels) {
      rels_params.emplace_back(rel);
    }
    for (auto &rel : full_rels) {
      full_rels_params.emplace_back(rel);
    }
  }

  OPA_TGEN_IMPL(input_len, rels_params, full_rels_params, rels_type_pos,
                type_bias, rels_type_length);
};

struct CutPlan {
  std::vector<int> left;
  std::vector<int> right;
  std::vector<pii> take_left;
  std::vector<pii> take_right;
  u64 can_left;
  u64 can_right;
  bool is_sym() const { return left.size() > 0 && left == right; }

  void to_ids(const std::vector<pii> &tb, std::vector<int> &res) {
    res.clear();
    for (auto &e : tb) REP (j, e.ND)
        res.pb(e.ST);
  }

  /*
      Once the final cutplan is selected, format of left and right changes.
      Instead of being pos[i] -> number of time ith bucket is selected
      We have pos[i] -> which bucket for the ith selection (so a list of
      buckets (with possible repetitions))
      */
  void to_ids() {
    to_ids(take_left, left);
    to_ids(take_right, right);
  }
};

struct OutputRel {
  OutputRel(MitmEntryList *out, u64 nwant = -1) : out(out), nwant(nwant) {}
  MitmEntryList *out;
  u64 nwant;

  bool add(u64 nkey, MitmRel &res);
};

class StepHelper {
public:
  struct Params {
    const RelsStore *rels_store;
    const KeyRelsList *key_rels;
    const SolverState *target;
    int zero_ub;
    opa::threading::Dispatcher *dispatcher;
    StepDescription *out_desc;
    Params() {}
    Params(const RelsStore *rels_store, const KeyRelsList *key_rels,
           const SolverState *target, int zero_ub, StepDescription *out_desc,
           opa::threading::Dispatcher *dispatcher = nullptr)
        : rels_store(rels_store), key_rels(key_rels), target(target),
          zero_ub(zero_ub), out_desc(out_desc), dispatcher(dispatcher) {}
  };

  StepHelper(const Params &params) { m_params = params; }

  void get_all_rels(const std::vector<int> &ids, MidRelDesc &out_rels);

  void get_lim_rels(const std::vector<int> &ids, MidRelDesc &out_rels,
                    u64 nrels);

  void solve(u64 lim, int type = -1);

  void solve_quadrisection(const CutPlan &cut_plan, MidRelDesc &res);

  void solve_mitm(const CutPlan &cut_plan, const u64 lim, MidRelDesc &res);

  Params m_params;
};

OPA_NAMESPACE_END(opa, crypto, stream)

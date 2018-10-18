#pragma once
#include <opa/crypto/la/context.h>
#include <opa/math/common/Matrix.h>
#include <opa/math/common/algebra/basis.h>
#include <opa/math/common/base/range_coverage.h>
#include <opa/math/common/base/value.h>
#include <opa/math/common/fast_gf2.h>
#include <opa/utils/misc.h>
#include <opa_common.h>

OPA_NAMESPACE_DECL3(opa, crypto, la)
typedef opa::math::common::BitVec BitVec;
typedef opa::math::common::RangeCoverage RangeCoverage;
typedef opa::math::common::Basis Basis;

class CipherContext;
class JitContext;
class JitBuilder;
class RelationDriver;
class Tester;

struct OutRel {
  struct DiffEntry {
    // An input pair with the correct differential.
    u64 a1, a2;
    u64 b1, b2;
    u64 diff_res;
  };

  // Bias of the rel
  double cost;
  // Output bits that belongs to the relations
  u64 out;
  // Data pairs <input, constant>.
  // Constant is for linear relations, taking into account input + ignored
  // part of output.
  typedef std::pair<u64, int> InputLinearRel;
  std::vector<InputLinearRel> in;
  std::vector<DiffEntry> in_diff;
};

struct CipherData {
  typedef std::pair<BitVec, BitVec> IOPair;
  std::vector<IOPair> x;
  std::vector<std::pair<IOPair, IOPair> > x_diff;
};

struct CallInfo {
  CallInfo(int input_base_size = 0, int key_base_size = 0, bool inlined = true,
           bool force_array_key = false) {
    this->input_base_size = input_base_size;
    this->key_base_size = key_base_size;
    this->inlined = inlined;
    this->force_array_key = force_array_key;
  }

  void setup(int input_size, int key_size, int output_size) {
    this->input_size = input_size;
    this->key_size = key_size;
    this->output_size = output_size;
    n_key_args = 0;
    n_input_args = 0;

    if (input_size) {
      if (input_base_size == 0)
        input_base_size = input_size;
      OPA_CHECK_EQ0(input_size % input_base_size, 0);
      n_input_args = input_size / input_base_size;
    }

    if (key_size) {
      if (key_base_size == 0)
        key_base_size = key_size;
      OPA_CHECK_EQ0(input_size % input_base_size, 0);
      n_key_args = key_size / key_base_size;
    }
    input_array = false;
    key_array = n_key_args > 1 || force_array_key;
  }

  int n_input_args;
  int n_key_args;
  bool input_array;
  bool key_array;

  int input_size;
  int key_size;
  int output_size;
  int input_base_size;
  int key_base_size;
  bool inlined;
  bool force_array_key;
};

struct Relation {
  struct RelationCost {
    static constexpr double kEps = 1e-8;
    OPA_DECL_COUT_OPERATOR3(RelationCost, m_bias);

    bool operator<(const RelationCost &x) const { OPA_LT_OP(x, m_bias); }

    bool is_exact() const { return m_bias >= 1 - kEps; }
    void from_count(int cnt, int sz) {
      this->m_bias = 2 * fabs(1. * cnt / sz - 0.5);
    }
    void from_walsh_count(int cnt, int sz) {
      this->m_bias = fabs(1. * cnt / sz);
    }
    double bias() const { return m_bias; }

    void set_bias(double p) { this->m_bias = p; }
    OPA_ACCESSOR_R(double, m_bias, prob);

  private:
    // TODO: rename to bias
    double m_bias = 0.;
  };

  bool operator==(const Relation &a) const {
    return a.in == in && a.out == out;
  }

  bool operator<(const Relation &a) const { OPA_LT_OP(a, in, out, cost); }
  OPA_DECL_COUT_OPERATOR3(Relation, in, out, c, cost);


  struct CmpVecByCost {
    bool operator()(const Relation &a, const Relation &b) const {
      return a.cost < b.cost;
    }
  };

  struct CmpVec {
    bool operator()(const Relation &a, const Relation &b) const {
      if (a.in != b.in)
        return a.in < b.in;
      return a.out < b.out;
    }
  };

  struct EqVec {
    bool operator()(const Relation &a, const Relation &b) const {
      return a == b;
    }
  };

  RangeCoverage in, out;
  bool c = 0;
  RelationCost cost;
};

struct ExactRelations : public opa::utils::Initable {
public:
  virtual void init(int input_size, int output_size,
                    const std::vector<Relation> &rels);
  bool try_rel(const Basis &out, Basis &res) const;

private:
  opa::math::common::Matrix<u32> m_mat;
  int input_size;
  int output_size;
};

struct Relations : public opa::utils::Initable {
  virtual void init(int size) {
    opa::utils::Initable::init();
    diff_basis.init(size);
  }

  void add_rel(const Relation &r) { tb.push_back(r); }
  void add_basis_vec(const BitVec &r) { diff_basis.add(r.to_mat()); }

  std::vector<Relation> tb;
  Basis diff_basis;

  friend std::ostream &operator<<(std::ostream &os, const Relations &v) {
    os << "Relations diff=" << v.diff << "\n";
    for (auto &e : v.tb)
      os << e << "\n";
    os << v.diff_basis;
    return os;
  }
  OPA_DECL_STR_FROM_COUT();

  bool diff;
};

class CipherBlock : public opa::utils::Initable,
                    public opa::utils::BaseStorable {
public:
  struct Params {
    CallInfo m_call;
    u32 m_input_size = 0;
    u32 m_output_size = 0;
    u32 m_key_size = 0;
    bool m_fast_eval = false;
    mutable SPTR(RelationDriver) m_rel_driver;

    Params();

    Params &call(const CallInfo &call) {
      m_call = call;
      return *this;
    }
    Params &input_size(u32 input_size) {
      m_input_size = input_size;
      return *this;
    }
    Params &output_size(u32 output_size) {
      m_output_size = output_size;
      return *this;
    }
    Params &key_size(u32 key_size) {
      m_key_size = key_size;
      return *this;
    }
    Params &fast_eval(bool fast_eval) {
      m_fast_eval = fast_eval;
      return *this;
    }
    // ownership transfered
    Params &rel_driver(RelationDriver *driver);
  };

  void create_jit_builder(JitBuilder &builder, JitContext &jc) const;

  typedef bool (*SolveCb)(const BitVec &v);
  virtual void set_context(const CipherContext *context) {
    m_context = context;
  }

  virtual ~CipherBlock() {}
  virtual void init(const Params &params) {
    opa::utils::Initable::init();
    m_params = params;
    m_params.m_input_size += m_params.m_key_size;
    m_func = 0;
    call().setup(raw_input_size(), key_size(), output_size());
  }
  virtual void fini() override;

  virtual void verify() const {}
  virtual void build() {}

  virtual BitVec evaluate(const BitVec &iv_) const;

  virtual u64 fast_eval1(u64 a) const { OPA_CHECK(0, desc()); }
  virtual u64 fast_eval2(u64 a, u64 b) const { OPA_CHECK(0, desc()); }
  virtual u64 fast_eval_tb(u64 a, const u64 *b) const { OPA_CHECK(0, desc()); }
  BitVec fast_eval_fe(const BitVec &iv) const;
  BitVec fast_eval_fe(const BitVec &input, const BitVec &kv) const;

  virtual void do_evaluate(const BitVec &iv, const BitVec &kv,
                           BitVec &ov) const = 0;

  const Relations &get_relations(bool diff = false) const {
    Relations &cur_rels = diff ? m_relations_diff : m_relations;
    if (!cur_rels.is_init()) {
      cur_rels.init(input_size() + output_size());
      if (diff)
        do_get_relations_diff(cur_rels);
      else
        do_get_relations(cur_rels);
      // OPA_CHECK0(m_relations.size() > 0);
      cur_rels.diff = diff;
    }
    return cur_rels;
  }

  const ExactRelations &get_exact_relations() const {
    if (!m_exact_relations.is_init()) {
      const auto &relations = get_relations();
      m_exact_relations.init(input_size(), output_size(), relations.tb);
    }
    return m_exact_relations;
  }

  const Relations &get_input_relations(bool diff = false) const {
    Relations &cur_rels = diff ? m_input_relations_diff : m_input_relations;
    if (!cur_rels.is_init()) {
      cur_rels = get_relations(diff);
      for (auto &rel : cur_rels.tb)
        rel.in.filter_greater(raw_input_size() - 1);
      // OPA_CHECK0(m_relations.size() > 0);
      sort(ALL(cur_rels.tb), Relation::CmpVec());
      opa::utils::make_unique(cur_rels.tb, Relation::EqVec());
      cur_rels.diff = diff;
    }
    return cur_rels;
  }

  void test_rel(int nr, const RangeCoverage &irel, const RangeCoverage &orel,
                bool diff = false) const;

  void get_relations_dumb(Relations &rels) const;
  void get_relations_walsh(Relations &rels) const;
  void get_relations_walsh_diff(Relations &rels) const;

  virtual void do_get_relations(Relations &rels) const = 0;
  virtual void do_get_relations_diff(Relations &rels) const = 0;

  virtual BitVec evaluate2(const BitVec &iv, const BitVec &kv) const {
    BitVec iv2;
    OPA_CHECK(iv.size() + kv.size() == input_size(), iv.size(), kv.size(),
              input_size());
    iv2.init(input_size());
    iv2.copy(iv);
    iv2.copy(kv, iv.size());
    return this->evaluate(iv2);
  }

  virtual void setup_jit(const JitBuilder &builder) const = 0;
  uintptr_t get_jit_func() const;

  virtual std::vector<RangeCoverage> get_pre(const RangeCoverage &out) const;

  std::vector<RangeCoverage> get_pre_key(const RangeCoverage &x) const {
    std::vector<RangeCoverage> tb, res;
    tb = get_pre(x);
    for (auto &a : tb) {
      res.pb(a.sshift(-raw_input_size()).filter_less(0));
    }
    return res;
  }

  virtual void do_get_pre(const Basis &out, Basis &res) const;

  virtual void do_get_pre_fail(const Basis &out, Basis &res) const {
    OPA_DISP("PRE FAIL DEFAULT ", desc());
    REP (i, input_size())
      res.add(RangeCoverage().add(i));
  }

  double test_outrel(const BitVec &kv, const OutRel &orel) const;

  CipherData gen_rand(int nr, const BitVec &kv) const;
  void gen_rand_diff(const BitVec &diff, int nr, const BitVec &kv,
                     CipherData *out_cipher_data) const;
  bool verify_data(const CipherData &data, const BitVec &kv) const;

  void do_evals(u32 *tb) const;

  void test_pre(int nr, const RangeCoverage &r);
  RelationDriver *driver() const { return m_params.m_rel_driver.get(); }

  BitVec evaluate_jit(const BitVec &iv, const BitVec &kv) const;

  OPA_ACCESSOR_R2(u32, m_params.m_input_size - m_params.m_key_size,
                  raw_input_size);
  OPA_ACCESSOR_R(u32, m_params.m_input_size, input_size);
  OPA_ACCESSOR_R(u32, m_params.m_output_size, output_size);
  OPA_ACCESSOR_R(u32, m_params.m_key_size, key_size);
  OPA_ACCESSOR_R(bool, m_params.m_fast_eval, fast_eval);
  OPA_ACCESSOR(CallInfo, m_params.m_call, call);
  OPA_ACCESSOR_PTR_R(CipherContext, m_context, context);

  OPA_ACCESSOR(std::string, m_desc, desc);

private:
  mutable Relations m_relations;
  mutable Relations m_relations_diff;
  mutable ExactRelations m_exact_relations;
  mutable Relations m_input_relations;
  mutable Relations m_input_relations_diff;
  const CipherContext *m_context = nullptr;
  std::string m_desc;
  Params m_params;
  mutable uintptr_t m_func;
};

#define OPA_TGEN_CIPHERBLOCK0() OPA_TGEN_CODE({ init(); }, { ; })

#define OPA_TGEN_CIPHERBLOCK(params_typ)                                       \
  OPA_TGEN_CODE(                                                               \
    {                                                                          \
      params_typ params;                                                       \
      opa::utils::tgen_load(params, any.any(0));                               \
      this->init(params);                                                      \
    },                                                                         \
    {                                                                          \
      any.add_any();                                                           \
      opa::utils::tgen_store(this->params(), *any.mutable_any(0));             \
    })

#define OPA_REG_BLOCK(cl)                                                      \
  _OPA_CLASS_STORE_REGISTER_BY_T(cl, []() {                                    \
    return opa::crypto::la::CipherContext::singleton.get()                     \
      ->instanciate<cl>(true)                                                  \
      .release();                                                              \
  })
OPA_NAMESPACE_DECL3_END

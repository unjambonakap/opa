#pragma once

#include <opa_common.h>
#include <opa/math/common/UtilsPoly.h>
#include <opa/math/common/Types.h>
#include <opa/math/common/FractionField.h>
#include <opa/math/common/Matrix.h>
#include <opa/math/common/GF_pBN.h>
#include <opa/math/common/GF_q.h>
#include <opa/math/common/Zn_BG.h>

OPA_NAMESPACE(opa, math, adv)
typedef OPA_MATH::bignum bignum;

template <class T> class RingExtension {
public:
  typedef OPA_MATH::Fraction<T> QE;
  void init(OPA_MATH::Ring<T> *ring, const OPA_MATH::Poly<T> &f) {
    m_ring = ring;
    m_fraction_field.reset(new OPA_MATH::FractionField<T>(ring));
    m_f = f;
    m_pr.init(m_ring);
    m_n = m_f.deg() + 1;
  }

  T get_norm(const T &a, const T &b) const {
    int d = m_f.deg();
    T bd = m_ring->getE();
    T n_b = m_ring->neg(b);
    OPA_MATH::Poly<T> rf = m_f;

    REP (i, d + 1) {
      rf[d - i] = m_ring->mul(rf[d - i], bd);
      bd = m_ring->mul(bd, n_b);
    }
    return m_pr.eval(rf, a);
  }

  T get_norm(const OPA_MATH::Poly<T> &a) const {
    OPA_MATH::Matrix<QE> norm_mat;
    OPA_MATH::PolyRing<QE> q_pr(m_fraction_field.get());

    std::vector<QE> q_a_vec;
    for (auto &v : a.toVector()) {
      q_a_vec.emplace_back(m_fraction_field->import_integer(v));
    }
    OPA_MATH::Poly<QE> q_a = q_pr.import(q_a_vec);

    norm_mat.init(m_n, m_n, &q_pr);
    REP (i, m_n) {
      OPA_MATH::Poly<QE> e = q_pr.xpw(m_fraction_field->getE());
      e = q_pr.mul(e, q_a);
      norm_mat.setRow(i, e.toVector());
    }
    QE dx = norm_mat.det();
    OPA_CHECK0(m_fraction_field.is_integer(dx));
    return dx.p;
  }

  int m_n;
  OPA_MATH::Ring<T> *m_ring;
  OPA_MATH::PolyRing<T> m_pr;
  std::unique_ptr<OPA_MATH::FractionField<T> > m_fraction_field;
  OPA_MATH::Poly<T> m_f;
};

struct NfsParams {
  NfsParams() {}
  NfsParams(const OPA_MATH::Poly<bignum> &f, const bignum &n, const bignum &m)
      : f(f), n(n), m(m) {
    f_deriv = f.get_poly_ring()->derivate(f);
  }

  OPA_MATH::Poly<bignum> f;
  OPA_MATH::Poly<bignum> f_deriv;
  bignum n, m;
};

struct NfsAlgebraicIdeal {
  u64 p;
  u64 r;
  bool operator<(const NfsAlgebraicIdeal &other) const {
    OPA_LT_OP(other, p, r);
  }
  OPA_DECL_COUT_OPERATOR2(NfsAlgebraicIdeal, a.p, a.r);
};

struct NfsRationalIdeal {
  u64 p;
  u64 m;
  bool operator<(const NfsRationalIdeal &other) const {
    OPA_LT_OP(other, p, m);
  }
  OPA_DECL_COUT_OPERATOR2(NfsRationalIdeal, a.p, a.m);
};

struct NfsIdealParams {
  std::vector<NfsRationalIdeal> rational_ideals;
  std::vector<NfsAlgebraicIdeal> algebraic_ideals;
  std::vector<NfsAlgebraicIdeal> quadratic_ideals;

  std::string str() const {
    std::ostringstream ss;
    ss << "Rational ideals:\n" << opa::utils::Join("\n", rational_ideals)
       << "\nAlgebraic ideals:\n" << opa::utils::Join("\n", algebraic_ideals)
       << "\nQuadratic ideals:\n" << opa::utils::Join("\n", quadratic_ideals)
       << "\n";

    return ss.str();
  }

  int get_dim() const {
    return rational_ideals.size() + algebraic_ideals.size() +
           quadratic_ideals.size();
  }

  OPA_DECL_COUT_OPERATOR(NfsIdealParams);
};

class NfsIdealHelper {
public:
  void setup(const NfsParams &params);
  void find_algebraic_ideals(int maxp);
  void find_rational_ideals(int maxp);
  void find_quadratic_ideals(int lowp, int maxp);

  void get_algebraic_ideals_for_prime(u32 p, std::vector<u32> *rlist);
  void get_algebraic_ideals(int lowp, int maxp,
                            std::set<NfsAlgebraicIdeal> *algebraic_ideals);

  void set_nfs_ideal_params(NfsIdealParams *ideal_params) const;

  std::set<NfsRationalIdeal> rational_ideals;
  std::set<NfsAlgebraicIdeal> algebraic_ideals;
  std::set<NfsAlgebraicIdeal> quadratic_ideals;

  NfsParams params;
};

struct SieveRationalEntry {
  s32 logv;
};

struct SieveAlgebraicEntry {
  s32 logv;
};

struct SieveResultEntry {
  bignum a = 0;
  bignum b = 0;
  SieveResultEntry(const bignum &a, const bignum &b) : a(a), b(b) {}
  SieveResultEntry() {}
};

class SieveComputationHelper {
public:
  void init(const RingExtension<bignum> *re, const NfsParams *params) {
    m_re = re;
    m_params = params;
  }

  s32 get_algebraic_log_bound(const SieveResultEntry &entry) const {
    bignum normv = m_re->get_norm(entry.a, entry.b);
    s32 logv = normv.get_size(2);
    return get_bound_from_log(logv);
  }

  s32 get_rational_log_bound(const SieveResultEntry &entry) const {
    bignum v = m_params->m * entry.b + entry.a;
    s32 logv = v.get_size(2);
    return get_bound_from_log(logv);
  }

  s32 get_bound_from_log(s32 logv) const {
    return std::max<s32>(logv * 9. / 10 - 5, logv * 2. / 3);
  }

  const NfsParams *m_params;
  const RingExtension<bignum> *m_re;
};

struct SieveLineStructure {
  std::vector<SieveAlgebraicEntry> algebraic_entries;
  std::vector<SieveRationalEntry> rational_entries;
};

class NfsSieveResult {
public:
  std::vector<SieveResultEntry> maybe_entries;
  std::vector<SieveResultEntry> checked_entries;
};

class NfsSquare {
public:
  std::vector<SieveResultEntry> pairs;
};

struct NfsRelation {
  std::vector<u32> rel;
  SieveResultEntry entry;
};

class NfsRelationCollector {
public:
  void init(int m) { this->m = m; }
  void get_kernel(OPA_MATH::Matrix<u32> *mat) const;
  void add_relation(const NfsRelation &relation) {
    relations.emplace_back(relation);
  }

  NfsSquare vector_to_square(const std::vector<u32> &vec) const;

  mutable OPA_MATH::Matrix<u32> orig_rel_mat;
  std::vector<NfsRelation> relations;
  int m = -1;
};

struct NfsSquareDataEntry {
  OPA_MATH::GF_pBN gfp;
  OPA_MATH::GF_q<bignum> gfq;
};

struct NfsSquareData {
  // pointer remains valid if only push_back
  std::deque<NfsSquareDataEntry> entries;
};

class NfsSieveHelper {
public:
  void setup(const NfsParams &params, const NfsIdealParams &ideal_params) {
    m_params = params;
    m_ideal_params = ideal_params;
    m_re.init(&OPA_MATH::Ring_Z, params.f);
    m_comp.init(&m_re, &m_params);
    m_collector.init(m_ideal_params.get_dim());
    zn.init(m_params.n);
    pr_zn.init(&zn);
    setup_square_data();
  }
  void setup_square_data();

  void do_sieve_stupid(const bignum &al, const bignum &ah, const bignum &bl,
                       const bignum &bh);
  void do_sieve(const bignum &al, const bignum &ah, const bignum &bl,
                const bignum &bh);
  void sieve_b(const bignum &b, const bignum &al, const bignum &ah);

  void ensure_valid_results();
  bool
  get_algebraic_decomposition(const SieveResultEntry &entry,
                              std::vector<u32> *decomposition = nullptr) const;
  bool check_algebraic_entry(const SieveResultEntry &entry,
                             NfsRelation *relation = nullptr) const;
  bool check_rational_entry(const SieveResultEntry &entry,
                            NfsRelation *relation = nullptr) const;
  bool check_entry(const SieveResultEntry &entry) const;

  bignum get_rational_val(const SieveResultEntry &entry) const;
  bignum get_algebraic_val(const SieveResultEntry &entry) const;

  void setup_quadratic_rel(const SieveResultEntry &entry,
                           NfsRelation *relation) const;

  void entry_to_relation(const SieveResultEntry &entry,
                         NfsRelation *relation) const;
  void setup_relations(const NfsSieveResult &result,
                       NfsRelationCollector *collector) const;

  int get_rational_offset(int rational_pos) const { return rational_pos; }
  int get_algebraic_offset(int algebraic_pos) const {
    return get_rational_offset(m_ideal_params.rational_ideals.size()) +
           algebraic_pos;
  }

  int get_quadratic_offset(int quadratic_pos) const {
    return get_algebraic_offset(m_ideal_params.algebraic_ideals.size()) +
           quadratic_pos;
  }

  bignum get_algebraic_root(const NfsSquare &square) const;
  bignum get_rational_root(const NfsSquare &square) const;
  bool maybe_add_square_data(const bignum &p);
  OPA_MATH::Poly<bignum>
  compute_square_root(const NfsSquareDataEntry &entry, const NfsSquare &square,
                      const std::vector<u32> &square_decomposition) const;

  std::vector<u32> compute_square_decomposition(const NfsSquare &square) const;

  bool verify_square(const NfsSquare &square) const;
  bignum try_square(const NfsSquare &square) const;

  NfsParams m_params;
  NfsIdealParams m_ideal_params;
  RingExtension<bignum> m_re;
  SieveComputationHelper m_comp;

  NfsSieveResult m_result;
  NfsRelationCollector m_collector;
  OPA_MATH::Zn_BG zn;
  OPA_MATH::PolyRing<bignum> pr_zn;

  NfsSquareData m_square_data;
};

OPA_NAMESPACE_END(opa, math, adv)

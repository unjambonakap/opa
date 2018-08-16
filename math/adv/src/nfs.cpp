#include <opa/math/adv/nfs.h>
#include <opa/math/common/Utils.h>
#include <opa/math/common/UtilsGFq.h>
#include <opa/math/common/algo.h>

using namespace OPA_MATH;
using namespace std;

OPA_NAMESPACE(opa, math, adv)

void NfsIdealHelper::setup(const NfsParams &params) { this->params = params; }

void NfsIdealHelper::find_rational_ideals(int maxp) {
  for (auto &p : OPA_MATH::pl) {
    if (p >= maxp)
      break;
    rational_ideals.insert(NfsRationalIdeal({ p, (params.m % p).getu32() }));
  }
}

void NfsIdealHelper::get_algebraic_ideals_for_prime(u32 p,
                                                    std::vector<u32> *rlist) {

  OPA_MATH::GF_p gfp(p);
  OPA_MATH::PolyRing<u32> pr;
  pr.init(&gfp);

  std::vector<u32> tb;
  for (auto &e : params.f.toVector()) {
    tb.emplace_back((e % p).getu32());
  }

  Poly<u32> f2 = pr.import(tb);
  vector<Poly<u32> > factors = pr.factor(f2, 1);
  for (auto &factor : factors) {
    if (factor.deg() != 1)
      continue;
    u32 proot = pr.linear_root(factor);
    rlist->emplace_back(proot);
  }
}

void NfsIdealHelper::get_algebraic_ideals(
  int lowp, int maxp, std::set<NfsAlgebraicIdeal> *algebraic_ideals) {
  for (auto &p : OPA_MATH::pl) {
    if (p >= maxp)
      break;
    if (p < lowp)
      continue;
    std::vector<u32> rlist;
    get_algebraic_ideals_for_prime(p, &rlist);
    for (auto &r : rlist) {
      algebraic_ideals->insert(NfsAlgebraicIdeal{ p, r });
    }
  }
}

void NfsIdealHelper::find_algebraic_ideals(int maxp) {
  get_algebraic_ideals(0, maxp, &algebraic_ideals);
}

void NfsIdealHelper::find_quadratic_ideals(int lowp, int maxp) {
  std::set<NfsAlgebraicIdeal> quad_dirty;
  get_algebraic_ideals(lowp, maxp, &quad_dirty);
  auto *poly_ring = params.f_deriv.get_poly_ring();

  for (auto &ideal : quad_dirty) {
    u32 resv = (poly_ring->eval(params.f_deriv, ideal.r) % ideal.p).getu32();
    if (resv == 0)
      continue;
    quadratic_ideals.insert(ideal);
  }
}

void NfsIdealHelper::set_nfs_ideal_params(NfsIdealParams *ideal_params) const {
  ideal_params->rational_ideals =
    std::vector<NfsRationalIdeal>(ALL(this->rational_ideals));
  ideal_params->algebraic_ideals =
    std::vector<NfsAlgebraicIdeal>(ALL(this->algebraic_ideals));
  ideal_params->quadratic_ideals =
    std::vector<NfsAlgebraicIdeal>(ALL(this->quadratic_ideals));
}

void NfsSieveHelper::do_sieve_stupid(const bignum &al, const bignum &ah,
                                     const bignum &bl, const bignum &bh) {
  bignum bpos = bl;
  for (; bpos < bh; bpos += 1) {
    bignum apos = al;
    for (; apos < ah; apos += 1) {
      SieveResultEntry entry(apos, bpos);
      if (this->check_entry(entry)) {

        m_result.maybe_entries.pb(entry);
      }
    }
  }
}

void NfsSieveHelper::do_sieve(const bignum &al, const bignum &ah,
                              const bignum &bl, const bignum &bh) {
  bignum bpos = bl;
  for (; bpos < bh; bpos += 1) {
    sieve_b(bpos, al, ah);
  }
}

void NfsSieveHelper::sieve_b(const bignum &b, const bignum &al,
                             const bignum &ah) {
  SieveLineStructure line_struct;
  u64 len = (ah - al).getu32();
  line_struct.algebraic_entries.resize(len);
  line_struct.rational_entries.resize(len);
  s32 rational_logv = m_comp.get_rational_log_bound(SieveResultEntry(al, b));

  for (int i = 0; i < len; ++i) {
    line_struct.rational_entries[i].logv = rational_logv;
    line_struct.algebraic_entries[i].logv =
      m_comp.get_algebraic_log_bound(SieveResultEntry(al + i, b));
  }

  for (const NfsRationalIdeal &rational_ideal :
       m_ideal_params.rational_ideals) {
    u32 p = rational_ideal.p;
    const auto &gfp = OPA_MATH::get_gfp(p);
    u64 b2 = (b % p).getu32();
    u64 a0 = gfp.mul(rational_ideal.m, gfp.inv(gfp.neg(b2)));
    u32 cur_log = log(p) / log(2.);

    s64 cur_a = a0 - (al % p).getu32();
    if (cur_a < 0)
      cur_a += p;

    for (; cur_a < len; cur_a += p) {
      line_struct.rational_entries[cur_a].logv -= cur_log;
    }
  }

  for (const NfsAlgebraicIdeal &algebraic_ideal :
       m_ideal_params.algebraic_ideals) {
    u32 p = algebraic_ideal.p;
    u32 r = algebraic_ideal.r;
    const auto &gfp = OPA_MATH::get_gfp(p);
    u32 b2 = (b % p).getu32();
    u32 a0 = (p - 1ll * b2 * r % p) % p;
    u32 cur_log = log(p) / log(2.);

    s64 cur_a = a0 - (al % p).getu32();
    if (cur_a < 0)
      cur_a += p;

    for (; cur_a < len; cur_a += p) {
      line_struct.algebraic_entries[cur_a].logv -= cur_log;
    }
  }

  REP (i, len) {
    if (line_struct.algebraic_entries[i].logv < 0 &&
        line_struct.rational_entries[i].logv < 0) {
      m_result.maybe_entries.pb(SieveResultEntry(al + i, b));
    }
  }
}

bignum NfsSieveHelper::get_rational_val(const SieveResultEntry &entry) const {
  return entry.a + entry.b * m_params.m;
}
bignum NfsSieveHelper::get_algebraic_val(const SieveResultEntry &entry) const {
  return m_re.get_norm(entry.a, entry.b);
}

bool NfsSieveHelper::check_rational_entry(const SieveResultEntry &entry,
                                          NfsRelation *relation) const {
  bignum v = get_rational_val(entry);
  int pos = 0;
  for (auto &ideal : m_ideal_params.rational_ideals) {
    int cnt = 0;
    while (v % ideal.p == 0)
      v /= ideal.p, cnt ^= 1;
    if (relation)
      relation->rel[this->get_rational_offset(pos)] = cnt;
    ++pos;
  }
  return v == 1;
}

// if decomposition is empty, initialized to zero vector
// otherwise values will be ADDED to the vector
bool NfsSieveHelper::get_algebraic_decomposition(
  const SieveResultEntry &entry, std::vector<u32> *decomposition) const {
  bignum norm = get_algebraic_val(entry);
  if (decomposition && decomposition->size() == 0) {
    *decomposition =
      std::vector<u32>(m_ideal_params.algebraic_ideals.size(), 0);
  }

  int pos = 0;
  for (auto &ideal : m_ideal_params.algebraic_ideals) {
    if (norm % ideal.p == 0) {
      if ((entry.a + entry.b * ideal.r) % ideal.p == 0) {
        u32 cnt = 0;
        while (norm % ideal.p == 0)
          norm /= ideal.p, cnt += 1;
        if (decomposition) // adding here, for easier updates
          (*decomposition)[pos] += cnt;
      }
    }
    ++pos;
  }
  return norm.abs() == 1;
}

bool NfsSieveHelper::check_algebraic_entry(const SieveResultEntry &entry,
                                           NfsRelation *relation) const {
  std::vector<u32> tb;
  bool res = get_algebraic_decomposition(entry, relation ? &tb : nullptr);
  if (relation) {
    REP (i, tb.size())
      relation->rel[this->get_algebraic_offset(i)] = tb[i] & 1;
  }
  return res;
}

bool NfsSieveHelper::check_entry(const SieveResultEntry &entry) const {
  return check_rational_entry(entry) && check_algebraic_entry(entry);
}

void NfsSieveHelper::setup_quadratic_rel(const SieveResultEntry &entry,
                                         NfsRelation *relation) const {
  int pos = 0;
  for (auto &ideal : m_ideal_params.quadratic_ideals) {
    bignum v = entry.a + entry.b * ideal.r;
    v %= ideal.p;
    bignum legendre = v.legendre(ideal.p);
    OPA_CHECK(legendre.abs() == 1, v, ideal.p, legendre);
    relation->rel[this->get_quadratic_offset(pos)] = (legendre != 1);
    ++pos;
  }
}

void NfsSieveHelper::ensure_valid_results() {
  for (auto &entry : m_result.maybe_entries) {
    if (entry.a == 0 && entry.b == 0)
      continue;
    if (check_entry(entry)) {
      m_result.checked_entries.emplace_back(entry);
    }
  }
}
void NfsSieveHelper::entry_to_relation(const SieveResultEntry &entry,
                                       NfsRelation *relation) const {
  relation->rel.resize(m_ideal_params.get_dim());
  relation->entry = entry;
  check_rational_entry(entry, relation);
  check_algebraic_entry(entry, relation);
  setup_quadratic_rel(entry, relation);
}

void NfsSieveHelper::setup_relations(const NfsSieveResult &result,
                                     NfsRelationCollector *collector) const {
  for (auto &entry : result.checked_entries) {
    NfsRelation relation;
    entry_to_relation(entry, &relation);
    collector->add_relation(relation);
  }
}

void NfsRelationCollector::get_kernel(OPA_MATH::Matrix<u32> *mat) const {
  OPA_MATH::Matrix<u32> rel_mat;
  rel_mat.initialize(&GF2, relations.size(), this->m + relations.size());
  REP (i, rel_mat.getN()) {
    REP (j, this->m) { rel_mat.get(i, j) = relations[i].rel[j]; }
    rel_mat.get(i, this->m + i) = 1;
  }

  this->orig_rel_mat =
    rel_mat.get_submatrix(0, 0, relations.size(), this->m).clone();
  int pos = rel_mat.row_echelon(-1, this->m);
  OPA_DISP0(rel_mat.getN(), rel_mat.getM(), this->m, pos);
  // exit(0);
  *mat = rel_mat.get_submatrix(pos, this->m).clone();
}

NfsSquare
NfsRelationCollector::vector_to_square(const std::vector<u32> &vec) const {
  NfsSquare res;
  REP (i, vec.size())
    if (vec[i])
      res.pairs.emplace_back(this->relations[i].entry);
  return res;
}

void NfsSieveHelper::setup_square_data() {
  bignum tot = 1;
  if (1) {
    bignum start = std::min(m_params.n + 10, bignum(1e9));
    while (tot < m_params.n.pow(20) * 2) {
      start = OPA_MATH::next_prime(start);
      if (maybe_add_square_data(start)) {
        tot *= start;
      }
    }
  } else {
    // OPA_CHECK0(maybe_add_square_data(int(1000037)));
    // OPA_CHECK0(maybe_add_square_data(int(1000039)));
    // OPA_CHECK0(maybe_add_square_data(int(1000081)));
    // OPA_CHECK0(maybe_add_square_data(int(1e9+7)));
    OPA_CHECK0(maybe_add_square_data(9851));
    OPA_CHECK0(maybe_add_square_data(9907));
    OPA_CHECK0(maybe_add_square_data(9929));
  }
}

bool NfsSieveHelper::maybe_add_square_data(const bignum &p) {
  GF_pBN gfp(p);
  PolyRing<bignum> pr(&gfp);

  Poly<bignum> pf = pr.import(m_params.f);
  if (pf.deg() != m_params.f.deg())
    return false;
  LOG(INFO) << "Trying for " << p << " " << pf;
  if (!pr.isIrred(pf))
    return false;

  m_square_data.entries.emplace_back();
  NfsSquareDataEntry *entry = &m_square_data.entries.back();
  entry->gfp = gfp;
  entry->gfq.init(&entry->gfp, pf);
  LOG(INFO) << "Ok for " << p << " " << pr.factor(pf);

  return true;
}

bignum NfsSieveHelper::get_algebraic_root(const NfsSquare &square) const {
  std::vector<u32> decomposition;
  for (auto &square_elem : square.pairs) {
    OPA_CHECK0(this->get_algebraic_decomposition(square_elem, &decomposition));
  }

  // valid square, check holds
  for (auto &e : decomposition) {
    OPA_CHECK0(e % 2 == 0);
    e >>= 1;
  }

  std::vector<CRT> crt_tb(m_params.f.deg());

  for (int i = 0; i < this->m_square_data.entries.size(); ++i) {
    printf("\nprocessing square %d %lu\n", i,
           this->m_square_data.entries.size());
    const auto &split_field = this->m_square_data.entries[i];
    Poly<bignum> ev = compute_square_root(split_field, square, decomposition);

    REP (i, crt_tb.size()) {
      crt_tb[i].add(split_field.gfp.getSize(), 1, ev.get_safe(i));
    }
    // OPA_DISP0(ev, ev % m_params.n, ev * ev % m_params.n);
  }

  std::vector<bignum> lifted_coeffs;
  for (auto &crt : crt_tb) {
    bignum bound = crt.get_bound();
    bignum v = crt.solve();
    if (v > bound / 2)
      v -= bound;
    lifted_coeffs.emplace_back(v);
  }
  Zn_BG zn(m_params.n);
  PolyRing<bignum> pr_zn(&zn);
  bignum res = pr_zn.import(lifted_coeffs)(m_params.m);
  return res;
}

OPA_MATH::Poly<bignum> NfsSieveHelper::compute_square_root(
  const NfsSquareDataEntry &entry, const NfsSquare &square,
  const std::vector<u32> &square_decomposition) const {
  typedef Poly<bignum> FieldT;
  const GF_q<bignum> &gfq = entry.gfq;
  FieldT p_square = gfq.getE();

  bignum target_norm = 1;
  REP (i, square_decomposition.size()) {
    REP (j, square_decomposition[i])
      target_norm =
        target_norm * m_ideal_params.algebraic_ideals[i].p % gfq.getCar();
  }

  for (auto &entry : square.pairs) {
    FieldT imported_entry = gfq.import_vec({ entry.a, entry.b });
    p_square = gfq.mul(p_square, imported_entry);
  }

  PolyRing<FieldT> pr(&gfq);
  vector<FieldT> fpoly_vec;
  for (auto &e : m_params.f_deriv.toVector()) {
    fpoly_vec.emplace_back(gfq.from_base_field(e));
  }

  Poly<FieldT> fpoly = pr.import(fpoly_vec);
  FieldT ex = pr.eval(fpoly, gfq.theta());
  p_square = gfq.mul(p_square, gfq.mul(ex, ex));

  OPA_CHECK0(has_square_root(gfq, p_square));

  FieldT p_root;
  find_square_root(gfq, p_square, &p_root);
  target_norm = target_norm * gfq.get_norm(ex) % gfq.getCar();
  // OPA_DISP0(gfq.get_norm(p_root), gfq.get_norm(gfq.neg(p_root)), target_norm,
  // gfq.getCar());
  if (gfq.get_norm(p_root) != target_norm) {
    p_root = gfq.neg(p_root);
    OPA_CHECK(gfq.get_norm(p_root) == target_norm, target_norm,
              gfq.get_norm(p_root));
  }
  return p_root;
}

bignum NfsSieveHelper::get_rational_root(const NfsSquare &square) const {
  // TODO: use smooth basis to compute this hist
  bignum res = 1;
  for (const auto &entry : square.pairs) {
    bignum v = entry.a + entry.b * m_params.m;
    res = res * v;
  }

  OPA_DISP0(m_params.f(m_params.m) % m_params.n);

  bignum fdm = m_params.f_deriv(m_params.m);
  OPA_DISP0(res);

  bignum root = res.sqrt();
  OPA_DISP0(res * fdm * fdm);
  OPA_DISP0(root * fdm % m_params.n);
  OPA_CHECK(root * root == res, root, res);
  root = root * fdm % m_params.n;
  return root;
}

bool NfsSieveHelper::verify_square(const NfsSquare &square) const {
  Matrix<u32> v =
    Matrix<u32>::fromzerovec(&GF2, m_ideal_params.get_dim(), VecType::Row);
  for (auto &entry : square.pairs) {
    NfsRelation relation;
    this->entry_to_relation(entry, &relation);
    v = v.add(v.fromvec(relation.rel, VecType::Row));
  }
  return v.is_null();
}

bignum NfsSieveHelper::try_square(const NfsSquare &square) const {
  bignum algebraic_root = get_algebraic_root(square);
  bignum rational_root = get_rational_root(square);

  bignum v1 = algebraic_root * algebraic_root % m_params.n;
  bignum v2 = rational_root * rational_root % m_params.n;
  OPA_CHECK(v1 == v2, v1, v2);
  bignum d1 = (algebraic_root - rational_root) % m_params.n;
  bignum d2 = (algebraic_root + rational_root) % m_params.n;
  bignum d = m_params.n.gcd(d1);
  return d;
}

OPA_NAMESPACE_END(opa, math, adv)

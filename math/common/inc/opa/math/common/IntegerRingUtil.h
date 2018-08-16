#pragma once

#include <opa/math/common/GF_p.h>
#include <opa/math/common/Poly.h>
#include <opa/math/common/PolyRing.h>
#include <opa/math/common/Ring.h>
#include <opa/math/common/Types.h>
#include <opa/math/common/Utils.h>
#include <opa/math/common/UtilsPoly.h>
#include <opa/math/common/UtilsRing.h>
#include <opa/math/common/Zn_BG.h>
#include <opa/math/common/bignum.h>
#include <opa/math/common/stats.h>

OPA_NAMESPACE_DECL3(opa, math, common)

class MignotteBoundHelper {
public:
  static bignum ComputeBound(const P_Z &poly) {
    MignotteBoundHelper helper(&poly);
    return helper.get_bound();
  }

  MignotteBoundHelper(const P_Z *poly) { init(poly); }

  void init(const P_Z *poly) {
    this->poly = poly;
    poly_norm = compute_norm(poly);
  }

  bignum compute_norm(const P_Z *poly) const {
    bignum sum = 0;
    for (auto &e : poly->to_vec()) {
      sum += e * e;
    }
    sum = sum.sqrt();
    return sum;
  }

  // bignum compute_coeff_bound(int pos) const {
  //  OPA_CHECK0(false);
  //  bignum v = nchoosek(n, pos) * poly_norm + nchoosek(n - 1, pos - 1);
  //}

  bignum get_bound() const {
    int d = poly->deg();
    return poly_norm * nchoosek(d - 1, d / 2) + nchoosek(d - 1, d / 2 - 1);
  }

  const P_Z *poly;
  bignum poly_norm;
};

class HenselLifting {
public:
  typedef P_Z P_t;
  struct Params {
    const bignum *bound = nullptr;
    u32 p = 0;
    const P_t *z_poly = nullptr;
    const std::vector<Poly<u32> > *orig_factors = nullptr;

    Params() {}
    Params(const bignum *bound, u32 p, const P_t *z_poly,
           const std::vector<Poly<u32> > *orig_factors)
        : bound(bound), p(p), z_poly(z_poly), orig_factors(orig_factors) {}
  };

  struct LiftNode {
    P_t p, gcd_coeff;
    LiftNode *left = nullptr;
    LiftNode *right = nullptr;
  };

  void init(const Params &params) { this->params = params; }
  void set_cur_ring(const bignum &pe) {
    zn.init(pe * params.p);
    pr_zn.init(&zn);
    this->cur_pe = pe;
  }

  void set_gcd_coeffs(const P_t &a, const P_t &b, P_t *u, P_t *v) const {
    pr_zn.egcd(a, b, *u, *v);
  }

  LiftNode *init_lift_tree(int T, int H) {
    nodes.emplace_back();
    LiftNode *cur = &nodes.back();

    if (T + 1 == H) {
      cur->p = factors[T];
    } else {
      int M = (T + H) / 2;
      cur->left = init_lift_tree(T, M);
      cur->right = init_lift_tree(M, H);
      cur->p = pr_zn.mul(cur->left->p, cur->right->p);
      // OPA_DISP("INIT >> ", cur->left->p, cur->right->p, cur->p);
      set_gcd_coeffs(cur->left->p, cur->right->p, &cur->left->gcd_coeff,
                     &cur->right->gcd_coeff);
    }
    return cur;
  }

  void lift_one(const P_t &c, P_Z &a, P_t &b, P_Z &u, P_Z &v) const {
    a.unsafe_change_ring(&pr_zn);
    b.unsafe_change_ring(&pr_zn);
    u.unsafe_change_ring(&pr_zn);
    v.unsafe_change_ring(&pr_zn);
    OPA_CHECK0(u.deg() < b.deg());
    OPA_CHECK0(v.deg() < a.deg());
    OPA_CHECK0(a.lc().gcd(params.p) == 1);
    // OPA_DISP0(c, a, b, cur_pe, c - a * b);
    auto nc = a * b;
    auto tmp = c - nc;
    for (auto &x : tmp) OPA_CHECK(x % cur_pe == 0, nc, c, a, b, tmp, cur_pe);

    P_t f = (c - a * b) / cur_pe;
    P_t t = v * f / a;
    P_t a0 = v * f - a * t;
    P_t b0 = u * f + b * t;

    a = a + a0 * cur_pe;
    b = b + b0 * cur_pe;
  }

  void process_lift_tree(LiftNode *node) {
    if (node->left == nullptr) return;
    lift_one(node->p, node->left->p, node->right->p, node->left->gcd_coeff,
             node->right->gcd_coeff);

    process_lift_tree(node->left);
    process_lift_tree(node->right);
  }

  // factors are ovewritten
  void do_lift() {
    bignum pe = params.p;
    while (true) {
      if (pe > *params.bound * 2 * params.z_poly->lc()) break;
      set_cur_ring(pe);
      root->p = pr_zn.import(*params.z_poly).monic();
      // OPA_DISP("LIFT STEP >> ", root->p, *params.z_poly);
      process_lift_tree(root);
      pe *= params.p;
    }
    cur_pe *= params.p;

    factors.clear();
    for (auto &node : nodes) {
      if (node.left != nullptr) continue;
      factors.emplace_back(node.p);
    }
  }

  void lift(std::vector<P_t> *out_lifted_factors) {
    this->out_lifted_factors = out_lifted_factors;
    set_cur_ring(1);

    for (auto &factor : *params.orig_factors) {
      factors.emplace_back(
        pr_zn.import(conv_vector_t_to_bignum<u32>(factor.to_vec())));
    }
    root = init_lift_tree(0, factors.size());
    do_lift();

    find_orig_factors(out_lifted_factors);
  }

  struct RecState {
    std::vector<int> used;
    int rem;
    std::vector<int> cur_cnd;
    P_t target_rem;
    std::vector<P_t> result;
  };

  int do_rec(int pos, int cnd_size, const P_t &cnd) {
    if (cnd_size > 0) {
      P_t z_cnd;
      // OPA_DISP("Candidate ", cnd * rec_state.target_rem.lc(), cnd, cur_pe);
      if (back_to_z(cnd * rec_state.target_rem.lc(), &z_cnd)) {
        z_cnd = z_cnd.pp();

        if (rec_state.target_rem[0] == 0 ||
            (z_cnd[0] != 0 && rec_state.target_rem[0] % z_cnd[0] == 0)) {
          // OPA_DISP0("TRYing ", z_cnd, rec_state.rem, rec_state.target_rem,
          // cnd_size);
          bool take_cnd = false;
          take_cnd = (rec_state.target_rem % z_cnd).deg() == -1;

          if (take_cnd) {
            rec_state.target_rem = rec_state.target_rem / z_cnd;
            rec_state.rem -= rec_state.cur_cnd.size();
            for (auto &e : rec_state.cur_cnd) {
              rec_state.used[e] = 1;
            }
            rec_state.result.emplace_back(z_cnd);
            return rec_state.cur_cnd[0];
          }
        }
      }
    }

    if (cnd_size >= rec_state.rem / 2) return pos - 1;
    if (pos == rec_state.used.size()) return pos - 1;
    int resv = do_rec(pos + 1, cnd_size, cnd);
    if (resv < pos) return resv;
    if (rec_state.used[pos]) return pos;

    P_t cnd_take = cnd * factors[pos];
    rec_state.cur_cnd.pb(pos);
    resv = do_rec(pos + 1, cnd_size + 1, cnd_take);
    rec_state.cur_cnd.pop_back();
    if (resv < pos) return resv;
    return pos;
  }

  void find_orig_factors(std::vector<P_t> *out_lifted_factors) {
    // OPA_DISP("BOUND >> ", *params.bound, factors, *params.z_poly, cur_pe);
    rec_state.used.resize(factors.size(), 0);
    rec_state.rem = factors.size();
    rec_state.target_rem = *params.z_poly;
    do_rec(0, 0, pr_zn.getE());

    if (rec_state.target_rem.deg() >= 1) {
      P_t last_cnd = pr_zn.getE();
      REP (i, rec_state.used.size()) {
        if (!rec_state.used[i]) last_cnd = last_cnd * factors[i];
      }

      OPA_CHECK0(pr_zn.is_monic(last_cnd)); // only handling monic
      last_cnd = last_cnd * rec_state.target_rem.lc();
      // OPA_DISP0(last_cnd);
      OPA_CHECK0(back_to_z(last_cnd, &last_cnd));
      // OPA_DISP0(last_cnd);
      OPA_CHECK(last_cnd % rec_state.target_rem == (u32)0, last_cnd,
                rec_state.target_rem, cur_pe);
      OPA_CHECK(PR_Z.isZ(last_cnd % rec_state.target_rem), last_cnd,
                rec_state.target_rem, cur_pe);
      out_lifted_factors->emplace_back(last_cnd);
    }
    // OPA_DISP0(*out_lifted_factors, rec_state.target_rem);

    for (const P_t &factor : rec_state.result) {
      out_lifted_factors->emplace_back(factor);
    }
  }

  bool norm_coeff(const bignum &coeff, bignum *res) const {
    *res = coeff <= *params.bound ? coeff : coeff - cur_pe;
    return (-*res <= *params.bound);
  }

  bool back_to_z(const P_t &poly, P_t *res) const {
    std::vector<bignum> remapped_coeffs(poly.size());

    REP (i, poly.size()) {
      if (!norm_coeff(poly[i], &remapped_coeffs[i])) return false;
    }

    *res = PR_Z.import(remapped_coeffs);
    if (res->lc() < 0) *res = -*res;
    // OPA_DISP("OK convert ", *res);
    return true;
  }

  RecState rec_state;
  LiftNode *root;
  bignum cur_pe;

  Zn_BG zn;
  PolyRing<bignum> pr_zn;
  std::deque<LiftNode> nodes;
  std::vector<P_t> factors;
  std::vector<P_t> *out_lifted_factors;
  Params params;
};

void factor_zpoly_squarefree(const P_Z &poly, std::vector<P_Z> *result) {
  if (poly.deg() <= 0) return;
  OPA_CHECK(poly.lc() > 0, poly);
  bignum bound = MignotteBoundHelper::ComputeBound(poly);
  int p = 2;
  while (true) {
    p = nextPrimeSmall(p + 1);
    GF_p gfp(p);
    PolyRing<u32> gfp_pr(&gfp);
    std::vector<u32> coeffs;
    for (auto &v : poly.to_vec()) coeffs.emplace_back((v % p).getu32());

    Poly<u32> gfp_poly = gfp_pr.import(coeffs);
    if (gfp_poly.deg() != poly.deg()) continue;
    OPA_TRACE("Selected prime ", p);
    gfp_poly = gfp_poly.monic();

    std::vector<Poly<u32> > factors = gfp_pr.factor(gfp_poly);
    std::set<bignum> seen;
    for (auto &factor : factors) {
      gfp_pr.smonic(factor);
      bignum v = gfp_pr.to_base(gfp_pr.monic(factor));
      seen.insert(v);
    }
    // not squarefree for this prime, nop
    // OPA_DISP("Factor for ", p, factors, bound, gfp_poly,
    //         seen.size() == factors.size());
    if (seen.size() != factors.size()) continue;
    if (factors.size() == 1) {
      result->push_back(poly);
      return;
    }

    // Found good prime: do lifting and return
    // OPA_DISP("Doing lifting for ", p);
    HenselLifting lifting;
    lifting.init(HenselLifting::Params(&bound, p, &poly, &factors));
    lifting.lift(result);
    for (auto &factor : *result) {
      if (factor.lc() < 0) factor = -factor;
    }
    return;
  }
}

bignum factor_zpoly(const P_Z &poly,
                    std::vector<std::pair<P_Z, int> > *result) {
  bignum pc = poly.cont();
  P_Z rem;
  P_Z u = make_squarefree(poly.pp(), &rem);
  OPA_TRACE("After squarefree", poly, u);
  rem = rem * pc;
  if (u.lc() < 0) u = -u, rem = -rem;
  std::vector<P_Z> fact_squarefree;
  factor_zpoly_squarefree(u, &fact_squarefree);
  for (auto fact : fact_squarefree) result->emplace_back(fact, 1);

  for (auto &e : *result) {
    while (rem % e.first == (u32)0) rem = rem / e.first, ++e.second;
  }
  std::sort(ALL(*result));
  return rem.lc();
}

bool pz_is_irred(const P_Z &poly) {
  // TODO: so inefficient
  std::vector<std::pair<P_Z, int> > tmp;
  factor_zpoly(poly, &tmp);
  return tmp.size() == 1 && tmp[0].second == 1;
}

template <typename Q>
bignum factor_qpoly(const Poly<Q> &poly,
                    std::vector<std::pair<P_Z, int> > *result) {
  P_Z pz = PR_Z.get_poly();
  bignum d = force_poly_to_frac_base(poly, &pz);
  factor_zpoly(pz, result);
  return d;
}

template <typename Q>
bignum factor_qpoly_squarefree(const Poly<Q> &poly, std::vector<P_Z> *result) {
  P_Z pz = PR_Z.get_poly();
  bignum d = force_poly_to_frac_base(poly, &pz);
  if (pz.lc() < 0) d = -d, pz = -pz;
  // OPA_DISP0(pz);
  factor_zpoly_squarefree(pz, result);
  return d;
}

OPA_NAMESPACE_DECL3_END

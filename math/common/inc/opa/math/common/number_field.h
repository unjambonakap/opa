#pragma once

#include <opa/math/common/Field.h>
#include <opa/math/common/IntegerRingUtil.h>
#include <opa/math/common/PolyModRing.h>
#include <opa/math/common/Utils.h>
#include <opa/math/common/UtilsPoly.h>
#include <opa/math/common/Z.h>
#include <opa/math/common/algebra/basis.h>

OPA_NAMESPACE_DECL3(opa, math, common)

template <typename FType> class RootFinder {
  typedef Poly<FType> PolyT;
  typedef std::complex<FType> CType;
  typedef Poly<CType> PolyCT;
  typedef CType ComplexRoot;
  typedef FType RealRoot;

  enum RootType {
    FAIL = 0,
    REAL = 1,
    COMPLEX = 2,
  };

public:
  RootFinder(const FType &precision) {
    m_precision = precision;
    m_pr.init(&m_rfield);
    m_cpr.init(&m_cfield);
  }

  RootType find_one_root(const PolyT &rpoly, RealRoot *real_root,
                         ComplexRoot *complex_root) const {
    PolyCT poly = import_poly(m_cpr, rpoly);
    CType x(1.3, 0.314159);
    PolyCT dpoly = poly.derivate();

    CType v = poly(x);
    CType dv = dpoly(x);
    FType norm = std::norm(v);
    CType dx = v / dv;

    for (int nstep = 0;; ++nstep) {
      if (nstep > this->max_step) {
        OPA_DISP("Failed to converge", rpoly, v, x, dx, v, dv, dpoly);
        return RootType::FAIL;
      }

      if (std::abs(dx) < m_precision) {
        break;
      }

      CType nx = x - dx;
      CType nv = poly(nx);
      FType nnorm = std::norm(nv);
      if (nnorm < norm) {
        norm = nnorm;
        x = nx;
        v = nv;
        dx = v / dpoly(x);
        --nstep;
      } else {
        dx = dx / Float(2);
      }
    }

    if (std::abs(std::imag(x)) * std::abs(std::imag(x)) < m_precision) {
      *real_root = std::real(x);
      return RootType::REAL;
    }
    *complex_root = x;
    return RootType::COMPLEX;
  }

  void find_poly_roots(const P_Z &poly) {
    real_roots.clear();
    complex_roots.clear();
    all_roots.clear();
    PolyT cur_poly = import_poly(m_pr, poly);

    while (cur_poly.deg() > 0) {
      ComplexRoot complex_root;
      RealRoot real_root;
      RootType root_type = find_one_root(cur_poly, &real_root, &complex_root);
      OPA_CHECK0(root_type != RootType::FAIL);
      std::vector<FType> root_poly;

      if (root_type == RootType::REAL) {
        real_roots.push_back(real_root);
        root_poly = { -real_root, 1 };
      } else {
        complex_roots.push_back(complex_root);
        root_poly = { std::norm(complex_root), -(std::real(complex_root) * 2),
                      1 };
      }
      PolyT npoly = m_pr.get_poly();
      m_pr.ediv(cur_poly, m_pr.import(root_poly), &npoly, nullptr);
      cur_poly = npoly;
    }
    all_roots = complex_roots;
    for (auto &x : real_roots) all_roots.emplace_back(x);
    std::sort(ALL(all_roots));
  }

  std::vector<RealRoot> real_roots;
  std::vector<ComplexRoot> complex_roots;
  std::vector<ComplexRoot> all_roots;
  FType m_precision;
  RealField<FType> m_rfield;
  ComplexField<FType> m_cfield;

  PolyRingBaseField<FType> m_pr;
  PolyRingBaseField<CType> m_cpr;
  int max_step = 40;
};

class NumberField : public PolyModField<Q> {
public:
  typedef P_Q K;
  typedef Poly<K> P_K;
  typedef Poly<P_K> P_KK;

  struct NumberFieldElem {
    bignum d;
    P_Z a;
  };

  void init(const Field<Q> *field, const P_Z &root) {
    m_field = field;
    m_root = root;
    m_dim = root.deg();
    pq_root = to_pq(root).monic();
    pqq_root = import_change_var(Q_xy, pq_root);
    PolyModField<Q>::init(field, pq_root);
    OPA_CHECK0(this->get_underlying_ring() != nullptr);

    K_x.init(this);
    OPA_CHECK0(K_x.get_underlying_ring() != nullptr);

    K_xx.init(&K_x);
    pk_root = import_change_var(K_x, pq_root);
    OPA_TRACE(pq_root, pk_root, root);
  }

  P_Q get_char_poly(const NumberFieldElem &e) const {
    Poly<P_Z> dx_ay = Z_xy.import({ PR_Z.x() * e.d - e.a });
    Poly<P_Z> ty = import_change_var(Z_xy, m_root);
    P_Z res = resultant2(ty, dx_ay);

    std::vector<Q> tb_q;
    bignum dn = e.d.pow(m_root.size());
    for (auto &x : res) tb_q.push_back(QF(x, dn));
    return Q_x.import(tb_q);
  }

  Q get_norm(const NumberFieldElem &e) const {
    bignum dn = e.d.pow(m_root.size());
    bignum res = resultant2(m_root, e.a);
    return QF.import(res, dn);
  }

  P_Q get_norm_p_qq(const P_QQ &poly) const {
    OPA_DISP0(pqq_root, poly);
    P_Q res = resultant2(pqq_root, poly);
    return res;
  }

  // A(X - kY, O=Y)
  P_QQ get_axy(const P_K &poly, int k) const {
    P_QQ poly_xy = Q_xy.import(poly);
    // P_QQ poly_yx = poly_switch_vars(poly_xy);
    P_QQQ poly_xzy = import_change_var(Q_xyz, poly_xy);
    P_QQ eval_at = Q_xy.x() - Q_xy.constant(Q_x.x() * Q_x.constant(QF(k)));
    P_QQ a_rmp_xy = poly_xzy(eval_at);
    P_QQ a_rmp_yx = poly_switch_vars(a_rmp_xy);
    OPA_DISP0(eval_at, poly_xy, a_rmp_xy);
    return a_rmp_yx;
  }

  std::pair<P_Q, int> find_nonsquarefree_norm(const P_K &poly) const {

    for (int k = 0;; ++k) {
      P_QQ poly_xy = get_axy(poly, k);
      P_Q norm = get_norm_p_qq(poly_xy);
      P_K norm_pk = import_change_var(K_x, norm);
      if (!is_squarefree(norm_pk)) {
        OPA_DISP("Norm not squarefree, bad ", k);
        continue;
      }
      return { norm, k };
    }
  }

  void factor_squarefree(const P_K &p, std::vector<P_K> *factors) const {

    OPA_DISP("factoring squarefree >> ", p);
    P_Q norm;
    int k;
    std::tie(norm, k) = find_nonsquarefree_norm(p);

    OPA_DISP("FACTORQ >> ", norm, k);
    OPA_CHECK(norm.deg() >= p.deg(), norm, p);
    P_K norm_pk = import_change_var(K_x, norm);
    P_K eval_for_remap = K_x.x() + K_x.constant(this->x() * QF(k));
    P_KK norm_pkk = import_change_var(K_xx, norm_pk);
    P_K norm_pk_remap = norm_pkk(eval_for_remap);
    //OPA_DISP0(norm_pk_remap, p, norm_pk_remap % p,
    //          is_squarefree(norm_pk_remap));
    //OPA_CHECK_EQ(norm_pk_remap % p, K_x.getZ(), norm_pk_remap, p);

    std::vector<P_Z> q_factors;
    OPA_DISP("Factor Q >> ", norm);
    factor_qpoly_squarefree(norm, &q_factors);
    OPA_DISP0(q_factors);

    for (const auto &factor : q_factors) {
      P_Q fact_pq = to_pq(factor);
      P_K fact_pk = import_change_var(K_x, fact_pq);
      P_KK fact_pkk = import_change_var(K_xx, fact_pk);

      P_K fact_pk_remap = fact_pkk(eval_for_remap);
      P_K cnd = K_x.gcd(fact_pk_remap, p);
      OPA_DISP0(fact_pk_remap, p, cnd);

      if (cnd.deg() > 0) factors->push_back(cnd.monic());
    }
  }

  void factor_k(const P_K &_p, std::vector<std::pair<P_K, int> > *res) const {
    P_K rem;
    P_K p = _p.monic();
    P_K u = make_squarefree(p, &rem);
    std::vector<P_K> fact_squarefree;
    factor_squarefree(u, &fact_squarefree);
    for (auto fact : fact_squarefree) res->emplace_back(fact, 1);

    for (auto &e : *res) {
      OPA_DISP("Kappa ", rem, e.first);
      while (rem % e.first == K_x.getZ()) rem = rem / e.first, ++e.second;
    }
  }

  std::vector<Q> pk_to_vec(int nmod_deg, const P_K &v) const {
    std::vector<Q> x;
    REP (i, nmod_deg) {
      REP (j, pk_root.deg()) { x.push_back(v.get_safe(i).get_safe(j)); }
    }
    return x;
  }

  void extend(const P_K &np, NumberField *res, K *out_r1, K *out_r2) const {

    PolyModField<K> ef;
    ef.init(this, np);

    int k = 1;

    for (;; ++k) {
      int m = np.deg() * pk_root.deg();
      BasisT<Q> basis(m, m_field);
      P_K a = ef.x() + ef.constant(this->x() * QF(k));
      P_K cur = ef.getE();
      std::vector<Q> x;
      std::vector<P_K> cur_tb;
      while (true) {
        x = pk_to_vec(np.deg(), cur);
        cur_tb.push_back(cur);
        if (!basis.add(x)) break;
        OPA_DISP0(cur);
        cur = cur * a;
      }
      OPA_DISP0(k, basis.dim(), basis.full_span());

      if (basis.full_span()) {
        std::vector<Q> new_poly = basis.get_coord(x);
        P_Q new_poly_q = Q_x.xpw(m) - Q_x.import(new_poly);
        P_Z new_poly_z = PR_Z.get_poly();
        force_poly_to_frac_base(new_poly_q, &new_poly_z);
        OPA_DISP0(new_poly, new_poly_z);

        res->init(m_field, new_poly_z);

        P_K r1 = ef.constant(this->x());
        P_K r2 = ef.x();
        std::vector<Q> t1 = pk_to_vec(np.deg(), r1);
        std::vector<Q> t2 = pk_to_vec(np.deg(), r2);

        for (auto &v : cur_tb) v = res->K_x.import(v);
        P_K tmp = res->K_x.getZ();
        REP (i, t1.size())
          tmp = tmp + cur_tb[i] * res->constant(t1[i]);
        OPA_DISP0(tmp);

        std::vector<Q> x1 = basis.get_coord(t1);
        std::vector<Q> x2 = basis.get_coord(t2);
        K r1k, r2k;
        r1k = res->import(x1);
        r2k = res->import(x2);
        OPA_DISP0(x1, x2, r1k, r2k);
        OPA_DISP0(ef.import(pk_root)(this->import(r1k)));
        OPA_DISP0(ef.import(pk_root)(this->x()));
        OPA_DISP0(ef.import(np)(this->import(r2k)));

        OPA_DISP0(res->K_x.import(pk_root) %
                  (res->K_x.x() - res->K_x.constant(r1k)));
        OPA_DISP0(res->K_x.import(np) %
                  (res->K_x.x() - res->K_x.constant(r2k)));
        OPA_DISP0(res->K_x.import(pk_root)(r1k));
        OPA_DISP0(res->K_x.import(np)(r2k));
        *out_r1 = r1k;
        *out_r2 = r2k;

        // std::vector<std::pair<P_K, int> > factor_and_pw;
        // res->factor_k(import_change_var(res->K_x, pq_root), &factor_and_pw);
        // OPA_DISP0(factor_and_pw);
        return;
      }
    }
  }

  P_Q minimal_poly(K e) const {
    BasisT<Q> basis(m_dim, m_field);
    K x = getE();
    basis.add(x.to_vec(m_dim));
    int deg = 1;
    for (;; ++deg) {
      x = x * e;
      if (!basis.add(x.to_vec(m_dim))) break;
    }
    std::vector<Q> vpoly = basis.get_coord(x.to_vec(m_dim));
    P_Q p = Q_x.import(vpoly);
    OPA_DISP0(p, deg);
    p.get_force(deg) = QF(-1);
    OPA_DISP0(basis.mat().eval(vpoly), p, x.to_vec(m_dim));
    OPA_TRACES(import_change_var(K_x, p)(e));
    return -p;
  }

  const Field<Q> *m_field;
  PolyRing<P_Q> K_x;
  PolyRing<P_K> K_xx;
  P_K pk_root;
  P_Q pq_root;
  P_QQ pqq_root;
  P_Z m_root;
  int m_dim;
};

OPA_NAMESPACE_DECL3_END

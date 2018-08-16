#pragma once

#include <opa_common.h>
#include <opa/math/common/bignum.h>
#include <opa/math/common/Matrix.h>
#include <opa/math/common/Poly.h>
#include <opa/math/common/PolyRing.h>
#include <opa/math/common/GF_q.h>
#include <opa/math/common/UtilsGFq.h>

OPA_NAMESPACE_DECL3(opa, math, common)

template <class T> class LinearRecurrenceSolver {
public:
  typedef Poly<T> TPol;
  typedef Poly<T> TExt;
  typedef GF_q<T> SplitField;

private:
  class RootElem {
  public:
    RootElem(std::shared_ptr<SplitField> split_field, const Poly<T> &poly,
             const TExt &root)
        : m_poly(poly), m_root(root) {
      m_split_field = split_field;
      n = m_poly.deg();
      setup();
    }

    void setup() {
      m_mat.initialize(m_split_field->getBaseField(), n, n);

      PolyRing<T> pr(m_split_field->getBaseField());

      TExt cur = m_root;
      REP (i, n) {
        m_mat.setCol(i, pr.toVector(cur, n));
        cur = m_split_field->mul(cur, m_split_field->import(pr.x()));
      }
    }

    TExt get(const bignum &p) const {
      Matrix<T> tmp;
      tmp.affect_const(m_mat.faste(p));
      PolyRing<T> pr(m_split_field->getBaseField());
      return m_split_field->import(pr.import(tmp.getCol(0).tovec()));
      // return m_split_field->import(pr.import({}));
    }

  private:
    Poly<T> m_poly;
    std::shared_ptr<SplitField> m_split_field;
    Matrix<T> m_mat;
    TExt m_root;
    int n;
  };

  struct DegInternal {
    std::vector<TPol> polys;
    std::shared_ptr<SplitField> big_split_field;
    std::shared_ptr<SplitField> split_field;
    Matrix<T> toRootBase;
    Matrix<T> toBigField;
    const Field<T> *base_field;
    int cur_deg;
    int big_deg;

    std::vector<RootElem> roots;
    std::shared_ptr<GF_q<T> > field;

    void setup(const Field<T> *_base_field,
               std::shared_ptr<SplitField> big_split_field) {
      base_field = _base_field;
      this->big_split_field = big_split_field;

      const TPol &basePol = polys[0];
      split_field.reset(new SplitField(base_field, basePol));
      PolyRing<T> pr(base_field);
      PolyRing<TExt> split_pr(split_field.get());
      PolyRing<TExt> big_pr(big_split_field.get());

      big_deg = big_split_field->getModPoly().deg();
      cur_deg = basePol.deg();

      for (auto &poly : polys) {
        Poly<TExt> ext_poly = toExtField(split_field.get(), pr.import(poly));
        std::vector<Poly<TExt> > factors = split_pr.factor(ext_poly);
        assert(factors.size() > 0);
        assert(factors[0].deg() == 1);
        split_pr.smonic(factors[0]);

        TExt root = split_field->neg(factors[0].get(0));
        REP (i, poly.deg()) {
          roots.emplace_back(split_field, basePol, root);
          root =
            split_field->faste(root, split_field->getBaseField()->getSize());
        }
      }

      Matrix<T> tmp(base_field, cur_deg, cur_deg);

      TExt x = split_field->import(pr.x());
      TExt cur = x;

      REP (i, cur_deg) {
        tmp.setCol(i, pr.toVector(cur, cur_deg));
        cur = split_field->mul(cur, x);
      }

      bool ok = tmp.invert(&toRootBase);
      assert(ok);

      std::vector<Poly<TExt> > factors =
        big_pr.factor(toExtField(big_split_field.get(), basePol));
      assert(factors.size() > 0);
      assert(factors[0].deg() == 1);
      big_pr.smonic(factors[0]);

      x = big_split_field->neg(factors[0].get(0));
      cur = x;

      toBigField.initialize(base_field, big_deg, cur_deg);
      REP (i, cur_deg) {
        toBigField.setCol(i, pr.toVector(cur, big_deg));
        cur = big_split_field->mul(cur, x);
      }
    }

    int get_size() const { return roots.size(); }

    void get(TExt *tb, const bignum &p) const {
      PolyRing<T> pr(base_field);
      REP (i, roots.size()) {
        TExt cur = roots[i].get(p);
        std::vector<T> vec = pr.toVector(cur, cur_deg);
        vec = toRootBase.eval(vec);
        tb[i] = big_split_field->import(pr.import(toBigField.eval(vec)));
        // import just for form
      }
    }
  };

public:
  void initialize(const Matrix<T> &m, const Field<T> &field) {
    m_matrix.affect(m.clone());
    m_field = &field;
    m_pr.reset(new PolyRing<T>(m_field));
  }

  void setup2(const std::vector<Poly<T> > &factors) {
    m_factors = factors;
    u32 big_field_deg = 1;
    m_n_coeff = 0;
    for (const auto &x : m_factors) {
      m_n_coeff += x.deg();
      big_field_deg = u32_lcm(big_field_deg, x.deg());
    }

    m_big_split_field.reset(new GF_q<T>(m_field, big_field_deg));

    for (auto &x : m_factors) {
      m_pr->smonic(x);
      m_deg_data[x.deg()].polys.pb(x);
    }

    for (auto &x : m_deg_data)
      x.ND.setup(m_field, m_big_split_field);

    m_tmp_vals.resize(m_n_coeff);

    m_coeff_mat.initialize(m_big_split_field.get(), m_n_coeff, m_n_coeff);
    m_cur_rank = 0;
  }

  void setup() {
    auto char_poly = gfx_char_poly(m_matrix);
    auto factors = m_pr->factor(char_poly);
    setup2(factors);
  }

  bool add_init_cond(T _val, u32 p) {
    TExt val = m_pr->constant(_val);

    int pos = 0;
    for (auto &x : m_deg_data) {
      x.ND.get(m_tmp_vals.data() + pos, p);
      pos += x.ND.get_size();
    }

    m_coeff_mat.setRow(m_cur_rank, m_tmp_vals);
    int old = m_cur_rank;
    m_cur_rank = m_coeff_mat.rank(m_cur_rank + 1, m_n_coeff);
    assert(old == m_cur_rank || old + 1 == m_cur_rank);

    if (m_cur_rank == old + 1)
      m_val_tb.pb(val);

    if (m_cur_rank != m_n_coeff)
      return false;

    m_coeff_vec.initialize(m_big_split_field.get(), 1, m_n_coeff);
    std::vector<TExt> coeffs = m_coeff_mat.solve(m_val_tb);
    m_coeff_vec.setRow(0, coeffs);

    return true;
  }

  TExt ext_get(const bignum &p) const {

    int pos = 0;
    for (auto &x : m_deg_data) {
      x.ND.get((TExt *)m_tmp_vals.data() + pos, p);
      pos += x.ND.get_size();
    }

    TExt res = m_coeff_vec.eval(m_tmp_vals)[0];
    return res;
  }

  T get(const bignum &p) const {
    TExt tmp = ext_get(p);
    assert(tmp.deg() <= 0);
    if (tmp.deg() < 0)
      return 0;
    return tmp.get(0);
  }

private:
  Matrix<T> m_matrix;
  const Field<T> *m_field;
  std::shared_ptr<PolyRing<T> > m_pr;
  std::vector<Poly<T> > m_factors;
  std::map<int, DegInternal> m_deg_data;

  std::shared_ptr<SplitField> m_big_split_field;
  std::vector<TExt> m_tmp_vals;
  std::vector<TExt> m_val_tb;
  Matrix<TExt> m_coeff_vec;
  Matrix<TExt> m_coeff_mat;
  int m_cur_rank;
  int m_n_coeff;
};

template <class T>
bool checkGenSeq(const Ring<T> &ring, const std::vector<T> &gen,
                 const std::vector<T> &tb, int sz) {
  // printf("CHECKING >>> \n");
  // utils::out(gen);
  // utils::out(tb,sz);
  for (int i = gen.size() - 1; i < sz; ++i) {
    T tmp = ring.getZ();
    for (int j = 0; j < gen.size(); ++j)
      tmp = ring.add(tmp, ring.mul(gen[j], tb[i - j]));

    if (!ring.isZ(tmp))
      return false;
  }
  return true;
}

template <class T>
std::vector<T> findMinLinearRecursion_Slow(const Ring<T> &ring,
                                           const std::vector<T> &tb, int maxd) {
  maxd = std::min<int>(tb.size() / 2, maxd);
  // TODO: fix loop order
  for (int i = maxd; i >= 1; --i) {

    Matrix<T> matrix(&ring, i, i);
    std::vector<T> b(i, ring.getZ());
    for (int j = 0; j < i; ++j) {
      for (int k = 0; k < i; ++k)
        matrix(j, k) = tb[j + i - 1 - k];
      b[j] = tb[j + i];
    }

    std::vector<T> x = matrix.solve(b);

    if (x.size() != 0) {
      auto check = matrix.eval(x);
      REP(j,b.size()){
        OPA_CHECK(b[j] == check[j], j, b[j], check[j], b, check);
      }
      std::vector<T> res(i + 1);
      res[0] = ring.getE();
      for (int j = 0; j < i; ++j)
        res[j + 1] = ring.neg(x[j]);
      while(res.size()>0 && ring.isZ(res.back())) res.pop_back();
      if (checkGenSeq(ring, res, tb, tb.size())) return res;
    }
  }
  return std::vector<T>();
}

template <class T>
std::vector<T> findMinLinearRecursion_Massey(const Field<T> &field,
                                             const std::vector<T> &tb,
                                             int maxd = -1) {
  if (maxd == -1)
    maxd = tb.size();

  std::vector<T> a, b;

  int i = 0;
  for (; tb[i] == field.getZ() && i < tb.size(); ++i)
    ;

  if (i == tb.size())
    return std::vector<T>(1, field.getE());

  int lastPos = 0;
  T lastDiff;

  a.resize(i + 2, field.getZ());
  a[0] = field.getE();
  lastDiff = tb[i];
  lastPos = i;
  b.push_back(field.getE());

  for (++i; i < tb.size(); ++i) {

    T tmp = field.getZ();
    for (int j = 0; j < a.size(); ++j)
      tmp = field.add(tmp, field.mul(a[j], tb[i - j]));

    if (tmp != field.getZ()) {

      std::vector<T> na = a;
      int ns = i - lastPos + b.size();
      na.resize(std::max((int)a.size(), ns), field.getZ());
      T coeff = field.neg(field.div(tmp, lastDiff));
      for (int j = 0; j < b.size(); ++j)
        na[i - lastPos + j] =
          field.add(na[i - lastPos + j], field.mul(coeff, b[j]));

      if (a.size() != na.size()) {
        b = a;
        lastDiff = tmp;
        lastPos = i;
      }
      a = na;
    }

    // assert(checkGenSeq(field, a, tb, i+1));

    if (a.size() > maxd + 1) {
      return std::vector<T>();
    }
  }

  return a;
}

template <class T> bool has_square_root(const Field<T> &field, const T &a) {
  bignum order = field.getSize() - 1;
  OPA_CHECK0(order % 2 == 0);
  T tmp1 = field.faste(a, order / 2);
  return field.isE(tmp1);
}

template <class T>
bool find_square_root(const Field<T> &field, const T &a, T *answer) {
  if (!has_square_root(field, a))
    return false;

  bignum order = field.getSize() - 1;
  int s = 0;
  {
    bignum cur = order;
    for (; cur.get_bit(s) == 0; ++s)
      ;
  }

  bignum q = order.rshift(s);
  if (s == 1) {
    *answer = field.faste(a, (q + 1) / 2);
    return true;
  }

  T z;

  while (true) {
    z = field.getRandRaw();
    if (!field.isZ(z) && !has_square_root(field, z))
      break;
  }

  T c = field.faste(z, q);
  T res = field.faste(a, (q + 1) / 2);
  T t = field.faste(a, q);
  int m = s;
  while (!field.isE(t)) {
    // invariant res^2=t*a
    T tmp = t;
    int s2 = 0;
    for (; !field.isE(tmp); ++s2, tmp = field.mul(tmp, tmp))
      ;
    OPA_CHECK0(s2 < m);
    T b = field.faste(c, bignum(2).pow(m - s2 - 1)); // b order 2^(s2+1)
    res = field.mul(res, b);
    c = field.faste(b, 2); // c order 2^s2
    t = field.mul(t, c);   // t order <= s2-1
    m = s2;
  }
  *answer = res;
  return true;
}
OPA_NAMESPACE_DECL3_END

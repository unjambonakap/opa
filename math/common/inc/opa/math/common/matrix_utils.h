#pragma once

#include <opa/math/common/Matrix.h>
#include <opa/math/common/Types.h>
#include <opa/math/common/UtilsRing.h>

OPA_NM_MATH_COMMON

void gcd_simplify(Matrix<Z> &m) {
  GcdBuilder<Z> gcd(Ring_Z);
  REP (i, m.n)
    REP (j, m.m)
      gcd.update(m(i, j));
  if (Ring_Z.isInv(gcd.v)) return;
  m.scdiv(gcd.v);
}

void hnf(const Matrix<Z> &a, Matrix<Z> *b, Matrix<Z> *u_res = nullptr,
         Z *det_res = nullptr, bool nofail = false) {
  int n = a.getN();
  int m = a.getM();
  // working on rows operation, should be better for memory
  *b = a.transpose();
  Matrix<Z> u = Matrix<Z>::identity(a.get_ring(), m);
  Z det = a.get_ring()->getE();

  REP (i, std::min(n, m)) {
    std::vector<Z> vals;
    FOR (j, i, m)
      vals.push_back(b->get(j, i));

    // OPA_DISP0(*b);
    auto plan = get_gcd_plan(*a.get_ring(), vals);
    if (plan.gcd == 0 && nofail) {
      OPA_CHECK0(det_res != nullptr);
      *det_res = 0;
      return;
    }
    // OPA_DISP0(plan.steps, plan.pos, plan.gcd, vals);
    OPA_CHECK0(plan.gcd != 0);
    for (auto mat : { b, &u }) {
      for (const auto &step : plan.steps) {
        mat->row_op(i + step.r1, i + step.r2, step.a11, step.a12, step.a21,
                    step.a22);
      }
      if (plan.pos != 0) {
        mat->swap_rows(i, plan.pos + i);
        det = -det;
      }
    }

    OPA_CHECK(plan.gcd == b->get(i, i), plan, *b, vals);
    if (plan.gcd < 0) {
      for (auto mat : { b, &u }) {
        mat->get_mutable_row(i).scmul(-1);
      }
      plan.gcd = -plan.gcd;
      det = -det;
    }
    det = det * b->get(i, i);

    REP (j, m) {
      if (i != j) {
        Z coeff = -b->get(j, i) / plan.gcd;
        if (coeff == 0) continue;
        for (auto mat : { b, &u }) {
          mat->get_mutable_row(j).self_elem_addmul(mat->get_row(i), coeff);
        }
      }
    }
  }

  u = u.transpose();
  *b = b->transpose();

  Matrix<Z> tmp = a * u;
  OPA_CHECK(*b == tmp, a, tmp);
  if (u_res) *u_res = std::move(u);
  if (det_res) *det_res = det;
}

Z hnf_det(const Matrix<Z> &mat) {
  Matrix<Z> b;
  Z det;
  hnf(mat, &b, nullptr, &det, true /*nofail*/);
  return det;
}

bool reduce_hnf_vec(const Matrix<Z> &mat, std::vector<Z> &invec,
                    std::vector<Z> *expl_res = nullptr, bool force = false) {
  int n = mat.n;
  int m = mat.m;
  std::vector<Z> expl(n, 0);
  OPA_CHECK_EQ0(mat.n, invec.size());

  int c = 0;
  REP (i, n) {
    if (invec[i] == 0) continue;
    if (!force && (c == m || mat(i, c) == 0 || invec[i] % mat(i, c) != 0))
      return false;
    if (c == m || mat(i, c) == 0) continue;
    Z v = invec[i] / mat(i, c);
    expl[c] = v;
    REP (j, n)
      invec[j] -= v * mat(j, c);
    c++;
  }

  if (expl_res != nullptr) *expl_res = expl;
  REP (i, n)
    if (invec[i] != 0) return false;
  return true;
}

void snf(const Matrix<Z> &a, Matrix<Z> *d_res, Matrix<Z> *p_res,
         Matrix<Z> *q_res) {
  Matrix<Z> d, p, q;
  int n = a.n;
  int m = a.m;
  // q * a * p = d
  d = a.clone();
  q = a.identity(n);
  p = a.identity(m);

  int c = 0;
  for (int i = 0; i < std::min(n, m); ++i) {
    int pivot = -1;
    FOR (j, i, n)
      if (d(j, c) != 0) {
        pivot = i;
        break;
      }
    if (pivot == -1) {
      continue;
    }
    q.swap_rows(pivot, i);
    d.swap_rows(pivot, i);

    d.swap_cols(c, i);
    p.swap_cols(c, i);

    bool row_step = true;
    while (true) {
      std::vector<Z> vals;
      FOR (j, i + 1, m)
        if (d.get(i, j) != 0) goto NOTDONE;
      FOR (j, i + 1, n)
        if (d.get(j, i) != 0) goto NOTDONE;
      break;

    NOTDONE:
      if (row_step) {
        FOR (j, i, n)
          vals.push_back(d.get(j, i));
      } else {
        FOR (j, i, m)
          vals.push_back(d.get(i, j));
      }

      auto plan = get_gcd_plan(*a.get_ring(), vals);
      // OPA_DISP0(vals, plan, d);
      for (const auto &step : plan.steps) {
        int i1 = i + step.r1;
        int i2 = i + step.r2;
        if (row_step) {
          d.row_op(i1, i2, step.a11, step.a12, step.a21, step.a22);
          q.row_op(i1, i2, step.a11, step.a12, step.a21, step.a22);
        } else {
          d.col_op(i1, i2, step.a11, step.a12, step.a21, step.a22);
          p.col_op(i1, i2, step.a11, step.a12, step.a21, step.a22);
        }
      }
      if (plan.pos != 0) {
        if (row_step) {
          d.swap_rows(i, plan.pos + i);
          q.swap_rows(i, plan.pos + i);
        } else {
          d.swap_cols(i, plan.pos + i);
          p.swap_cols(i, plan.pos + i);
        }
      }
      OPA_CHECK(d(i, i) == plan.gcd, i, d, plan);

      if (row_step) {
        FOR (j, i + 1, n) {
          Z coeff = -d.get(j, i) / plan.gcd;
          if (coeff == 0) continue;
          d.get_mutable_row(j).self_elem_addmul(d.get_row(i), coeff);
          q.get_mutable_row(j).self_elem_addmul(q.get_row(i), coeff);
        }
      } else {
        FOR (j, i + 1, m) {
          Z coeff = -d.get(i, j) / plan.gcd;
          if (coeff == 0) continue;
          d.get_mutable_col(j).self_elem_addmul(d.get_col(i), coeff);
          p.get_mutable_col(j).self_elem_addmul(p.get_col(i), coeff);
        }
      }

      row_step = !row_step;
    }
    ++c;
  }

  if (d_res != nullptr) *d_res = std::move(d);
  if (q_res != nullptr) *q_res = std::move(q);
  if (p_res != nullptr) *p_res = std::move(p);
}

template <class T>
Matrix<T> expand_unimodular(const Ring<T> *ring, const std::vector<T> &tb,
                            int target_pos) {
  int n = tb.size();
  OPA_CHECK(target_pos >= 0 && target_pos < tb.size(), tb.size(), target_pos);
  std::vector<T> res = tb;
  Matrix<T> m = Matrix<T>::identity(ring, n);
  auto plan = get_gcd_plan(*ring, tb);
  OPA_CHECK0(ring->isInv(plan.gcd));

  for (auto step : plan.steps) {
    T va = ring->add(ring->mul(step.a11, tb[step.r1]),
                     ring->mul(step.a12, tb[step.r2]));
    T vb = ring->add(ring->mul(step.a21, tb[step.r1]),
                     ring->mul(step.a22, tb[step.r2]));
    res[step.r1] = va;
    res[step.r2] = vb;
    m.col_op(step.r1, step.r2, step.a11, step.a12, step.a21, step.a22);
  }

  T xx = ring->inv(res[plan.pos]);
  m.get_mutable_col(plan.pos).scmul(xx);
  OPA_DISP0(m, res, plan.pos);

  REP (i, n) {
    if (plan.pos == i) continue;
    T v = ring->neg(res[i]);
    m.get_mutable_col(i).self_elem_addmul(m.get_col(plan.pos), v);
  }
  m.swap_cols(plan.pos, target_pos);
  return m;
}

template <class T> void gram_schmidt(Matrix<T> &a) {
  int n = a.n;
  auto ring = a.get_ring();
  REP (i, n) {
    REP (j, i) {
      T dt = a.get_row(i).dot(a.get_row(j));
      T tt = a.get_row(j).dot(a.get_row(j));
      T d = ring->gcd(dt, tt);
      dt = ring->neg(ring->div(dt, d));
      tt = ring->div(tt, d);
      a.row_op(i, j, tt, dt, ring->getZ(), ring->getE());
      OPA_CHECK0(ring->isZ(a.get_row(i).dot(a.get_row(j))));
    }
  }
}

template <class T> const Ring<T> *GuessRing();
template <> const Ring<Z> *GuessRing<Z>() { return &Ring_Z; }
template <> const Ring<Q> *GuessRing<Q>() { return &QF; }

template <class T> struct Module_t {
  typedef Module_t<T> Type;
  typedef T value_type;

  std::vector<T> elems;
  static Type zero(int n) { return Type(n); }

  Type zero_like() const { return Type::zero(elems.size()); }

  Module_t(int n) {
    ring = GuessRing<T>();
    elems = std::vector<T>(n, ring->getZ());
  }

  Module_t() { ring = GuessRing<T>(); }
  Module_t(const Matrix<T> &m) {
    ring = GuessRing<T>();
    OPA_CHECK(m.is_vec(), m);
    elems = m.tovec();
  }

  typename std::vector<T>::const_iterator begin() const {
    return elems.begin();
  }
  typename std::vector<T>::const_iterator end() const { return elems.end(); }
  typename std::vector<T>::iterator begin() { return elems.begin(); }
  typename std::vector<T>::iterator end() { return elems.end(); }

  template <class Container> Module_t(const Container &elems_) {
    ring = GuessRing<T>();
    for (const auto &e : elems_) elems.push_back(e);
  }

  bool operator==(const Type &peer) const { return peer.elems == elems; }
  int size() const { return elems.size(); }

  const T &operator[](int a) const { return elems[a]; }
  T &operator[](int a) { return elems[a]; }

  Type operator-() const {
    Type res;
    for (auto &e : elems) res.elems.push_back(-e);
    return res;
  }

#define MODULE_FE(res, a, b, op)                                               \
  REP (i, (a).elems.size()) {                                                  \
    T x = (a).elems[i];                                                        \
    T y = (b).elems[i];                                                        \
    res.elems[i] = op;                                                         \
  }
  Type operator*(const T &x) const {
    Type res = *this;
    for (auto &e : res.elems) e = e * x;
    return res;
  }

  Type operator+(const T &x) const {
    Type res = *this;
    for (auto &e : res.elems) e = e + x;
    return res;
  }

  Type operator*(const Type &b) const {
    Type res = zero_like();
    MODULE_FE(res, *this, b, x * y);
    return res;
  }

  Type operator-(const Type &b) const {
    Type res = zero_like();
    MODULE_FE(res, *this, b, x - y);
    return res;
  }

  Type operator+(const Type &b) const {
    Type res = zero_like();
    MODULE_FE(res, *this, b, x + y);
    return res;
  }

  T dot(const Type &b) const {
    OPA_CHECK_EQ(elems.size(), b.elems.size(), *this, b);
    T res = ring->getZ();
    REP (i, elems.size()) { res = res + elems[i] * b.elems[i]; }
    return res;
  }

  OPA_DECL_COUT_OPERATOR2(Type, a.elems);
  const Ring<T> *ring;
};

using ZModule_t = Module_t<Z>;
using QModule_t = Module_t<Q>;
QModule_t lift_zmodule(const ZModule_t &x) {
  QModule_t res = QModule_t::zero(x.size());
  REP (i, x.size()) { res[i] = QF.import(x[i]); }
  return res;
}
ZModule_t force_zmodule(const QModule_t &x) {
  return ZModule_t(make_integers(Ring_Z, x.elems));
}

std::vector<QModule_t> lift_zmodules(const std::vector<ZModule_t> &x) {
  std::vector<QModule_t> res;
  for (auto &e : x) res.push_back(lift_zmodule(e));
  return res;
}

Matrix<Z> generate_lattice(int n, int m, Z bound, s64 maxstep = -1) {
  Matrix<Z> mat(&Ring_Z, n, m);
  REP (i, std::min(n, m))
    mat(i, i) = 1;
  int M = n + m;
  std::vector<Z> possibilities = { 0, 0, 1, -1 };
  for (; maxstep != 0; --maxstep) {
    int randv = rng() % M;
    if (randv < n) {
      int r = randv;
      ZModule_t cur = ZModule_t::zero(m);
      REP (i, n - 1) {
        Z mul = possibilities[rng() % possibilities.size()];
        cur = cur + ZModule_t(mat.get_row(i + (i >= r)).tovec()) * mul;
      }
      mat.get_mutable_row(r).addv(cur);
      if (mat.get_row(r).l_inf() >= bound) break;
    } else {
      int c = randv - n;
      ZModule_t cur = ZModule_t::zero(n);
      REP (i, m - 1) {
        Z mul = possibilities[rng() % possibilities.size()];
        cur = cur + ZModule_t(mat.get_col(i + (i >= c)).tovec()) * mul;
      }
      mat.get_mutable_col(c).addv(cur);
      if (mat.get_col(c).l_inf() >= bound) break;
    }
  }

  return mat;
}

OPA_NM_MATH_COMMON_END

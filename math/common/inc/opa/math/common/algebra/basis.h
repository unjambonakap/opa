#pragma once
#include <opa/math/common/Field.h>
#include <opa/math/common/GF_p.h>
#include <opa/math/common/Matrix.h>
#include <opa/math/common/Types.h>
#include <opa/math/common/base/range_coverage.h>
#include <opa/utils/misc.h>
#include <opa_common.h>

OPA_NAMESPACE(opa, math, common)
class RangeCoverage;

class Basis : public opa::utils::Initable {
public:
  Basis() {}
  Basis(int sz_) { init(sz_); }

  void copy(const Basis &a) {
    m_mat.affect(a.m_mat.clone());
    m_reduced.affect(a.m_reduced.clone());
    m_dim = a.m_dim;
    m_sz = a.m_sz;
  }

  Basis(const Basis &a) : Initable(a) {
    Initable::operator=(a);
    copy(a);
  }

  Basis &operator=(const Basis &a) {
    Initable::operator=(a);
    copy(a);
    return *this;
  }

  virtual void init(int sz_) {
    Initable::init();
    m_sz = sz_;
    m_mat.initialize(&GF2, m_sz, m_sz);
    m_reduced.initialize(&GF2, m_sz, m_sz);
    m_dim = 0;
  }

  virtual void init(const std::vector<RangeCoverage> &lst, int sz) {
    init(sz);
    for (auto &x : lst) add(x);
  }

  void add(const std::vector<Matrix<u32> > &tb) {
    for (auto &x : tb) {
      if (!add(x)) break;
    }
  }

  bool add(const RangeCoverage &x) {
    return add(
      mat().fromvec(x.to_binary_vec(sz()), opa::math::common::VecType::Row));
  }

  void add(const Basis &a) {
    REP (i, a.dim())
      add(a.mat().get_row(i));
  }

  bool add(const Matrix<u32> &x) {
    check_init();
    if (m_dim == m_sz) return false;
    reduced().set_submatrix(x, m_dim, 0);
    mat().set_submatrix(x, m_dim, 0);
    m_dim = reduced().row_echelon(m_dim + 1);
    return m_dim < m_sz;
  }

  Matrix<u32> get(int i) const { return m_mat.get_row(i); }

  // vector
  bool spans(const Matrix<u32> &v, Matrix<u32> *output = nullptr) const {
    Matrix<u32> tmp = v.clone();
    bool res = reduced().reduce(tmp, -1, true, output);
    return res;
  }

  bool spans(const Basis &b) const {
    REP (i, dim())
      if (!spans(b.get(i))) return false;
    return true;
  }

  friend std::ostream &operator<<(std::ostream &os, const Basis &v) {
    os << "Basis sz=" << v.sz() << " dim=" << v.dim() << "\n";
    os << v.mat().str(v.dim());
    return os;
  }

  void reduce() { m_dim = mat().row_echelon(); }

  OPA_ACCESSOR_NOSETTER(Matrix<u32>, m_mat, mat);
  OPA_ACCESSOR_NOSETTER(Matrix<u32>, m_reduced, reduced);
  OPA_ACCESSOR_R(int, m_dim, dim);
  OPA_ACCESSOR_R(int, m_sz, sz);

private:
  Matrix<u32> m_mat, m_reduced;
  int m_dim;
  int m_sz;
};

template <class T> class BasisT : public opa::utils::Initable {
public:
  typedef BasisT<T> SelfType;
  BasisT() {}
  BasisT(int sz_, const Field<T> *ring) { init(sz_, ring); }

  void copy(const SelfType &a) {
    m_mat.affect(a.m_mat.clone());
    m_reduced.affect(a.m_reduced.clone());
    m_dim = a.m_dim;
    m_sz = a.m_sz;
  }

  BasisT(const SelfType &a) : Initable(a) {
    Initable::operator=(a);
    copy(a);
  }

  SelfType &operator=(const SelfType &a) {
    Initable::operator=(a);
    copy(a);
    return *this;
  }

  virtual void init(int sz_, const Field<T> *field) {
    Initable::init();
    m_sz = sz_;
    m_field = field;
    m_mat.initialize(field, m_sz, m_sz);
    m_reduced.initialize(field, m_sz, 2 * m_sz);
    clear();
  }
  void clear() {
    m_dim = 0;
    m_mat.clear_zeroes();
    m_reduced.clear_zeroes();
    REP (i, m_sz)
      m_reduced(i, m_sz + i) = m_field->getE();
  }

  bool add(const std::vector<T> &x) {
    check_init();
    if (m_dim == m_sz) return false;
    reduced().setRow(m_dim, x);
    mat().setRow(m_dim, x);
    int old_dim = m_dim;
    m_dim = reduced().row_echelon(m_dim + 1, m_sz);
    return old_dim != m_dim;
  }

  Matrix<T> get(int i) const { return m_mat.get_col(i); }

  // vector
  bool spans(const Matrix<T> &v, Matrix<T> *output = nullptr) const {
    Matrix<T> tmp = v.clone();
    bool res = reduced().reduce(tmp, -1, true, output);
    return res;
  }

  bool spans(const SelfType &b) const {
    REP (i, dim())
      if (!spans(b.get(i))) return false;
    return true;
  }

  friend std::ostream &operator<<(std::ostream &os, const SelfType &v) {
    os << "Basis sz=" << v.sz() << " dim=" << v.dim() << "\n";
    os << v.mat().str(v.dim());
    return os;
  }

  std::vector<T> get_coord(const std::vector<T> &v) const {
    Matrix<T> imat = m_reduced.get_submatrix(0, m_sz);
    std::vector<T> res = imat.evalT(v);
    OPA_CHECK_EQ(v, m_mat.evalT(res), res, m_mat, m_reduced);
    return res;
  }

  void reduce() { m_dim = mat().row_echelon(); }
  bool full_span() const { return m_dim == m_sz; }

  OPA_ACCESSOR_NOSETTER(Matrix<T>, m_mat, mat);
  OPA_ACCESSOR_NOSETTER(Matrix<T>, m_reduced, reduced);
  OPA_ACCESSOR_R(int, m_dim, dim);
  OPA_ACCESSOR_R(int, m_sz, sz);

private:
  Matrix<T> m_mat, m_reduced;
  int m_dim = 0;
  int m_sz;
  const Field<T> *m_field;
};

OPA_NAMESPACE_END(opa, math, common)

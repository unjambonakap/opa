#pragma once

#include <opa/math/common/FractionField.h>
#include <opa/math/common/Poly.h>
#include <opa/math/common/PolyRing.h>
#include <opa/math/common/Ring.h>
#include <opa/math/common/Utils.h>
#include <opa/math/common/UtilsRing.h>
#include <opa/math/common/base.h>

#include <opa/utils/base.h>
#include <opa/utils/misc.h>
#include <opa/utils/string.h>

OPA_NM_MATH_COMMON

enum class VecType : int { Row, Col };

template <class T> class Matrix {

public:
  Matrix() {
    ring = 0;
    n = m = 0;
  }

  void init_objs();
  void init_data();
  void clear_zeroes() {
    REP (i, n)
      REP (j, m)
        get(i, j) = ring->getZ();
  }

  Matrix<T> clone() const {
    Matrix<T> res;
    res.copy(*this);
    return res;
  }

  Matrix<T>(Matrix<T> &&a) { affect(a); }
  Matrix<T> &operator=(Matrix<T> &&a) {
    affect(a);
    return *this;
  }
  Matrix<T>(const Matrix<T> &a) { copy(a); }

  void copy(const Matrix<T> &a);
  void prepare_set(const Matrix<T> &a, int r0, int c0, int nr, int nc);
  void set_from(const Matrix<T> &a, int r0 = 0, int c0 = 0, int nr = -1,
                int nc = -1);
  void set_mutable_from(Matrix<T> &a, int r0 = 0, int c0 = 0, int nr = -1,
                        int nc = -1);
  void set_mutable_from(Matrix<T> &&a, int r0 = 0, int c0 = 0, int nr = -1,
                        int nc = -1);

  std::vector<std::vector<T> > as_rows(int r0 = 0, int nr = -1) const {
    if (nr == -1) nr = n;
    std::vector<std::vector<T> > res;
    FOR (i, r0, nr)
      res.push_back(this->get_row(i).tovec());
    return res;
  }

  std::vector<std::vector<T> > as_cols(int c0 = 0, int nc = -1) const {
    if (nc == -1) nc = m;
    std::vector<std::vector<T> > res;
    FOR (i, c0, nc)
      res.push_back(this->get_col(i).tovec());
    return res;
  }

  Matrix<T> get_submatrix(int r0, int c0, int nr = -1, int nc = -1) const;
  Matrix<T> get_mutable(int r0, int c0, int nr = -1, int nc = -1);
  Matrix<T> get_mutable_col(int c);
  Matrix<T> get_mutable_row(int r);
  Matrix(const Ring<T> *ring, int n, int m);
  Matrix<T> cmul(const T &a) const {
    Matrix<T> res = this->clone();
    res.scmul(a);
    return res;
  }

  template <class Container> void set_cols(const Container &container) {
    OPA_CHECK(container.size() <= m, container.size(), m);
    REP (i, container.size())
      this->set_col(i, container[i]);
  }

  template <class Container> void set_rows(const Container &container) {
    OPA_CHECK(container.size() <= n, container.size(), n);
    REP (i, container.size())
      this->set_row(i, container[i]);
  }

  void row_op(int i1, int i2, const T &a11, const T &a12, const T &a21,
              const T &a22) {
    REP (i, m) {
      T v1 = get(i1, i);
      T v2 = get(i2, i);
      get(i1, i) = ring->add(ring->mul(v1, a11), ring->mul(v2, a12));
      get(i2, i) = ring->add(ring->mul(v1, a21), ring->mul(v2, a22));
    }
  }

  void col_op(int i1, int i2, const T &a11, const T &a12, const T &a21,
              const T &a22) {
    REP (i, n) {
      T v1 = get(i, i1);
      T v2 = get(i, i2);
      get(i, i1) = ring->add(ring->mul(v1, a11), ring->mul(v2, a12));
      get(i, i2) = ring->add(ring->mul(v1, a21), ring->mul(v2, a22));
    }
  }

  Matrix<T> faste(bignum p) const {
    Matrix<T> x, a;
    x.affect(identity());
    a.affect(this->clone());
    for (; p != 0; p.srshift(1), a.affect(a * a))
      if (p.get_bit(0)) x.affect(x * a);
    return x;
  }

  void set_submatrix(const Matrix<T> &a, int r0 = 0, int c0 = 0);

  template <class Container> void addv(const Container &a) {
    OPA_CHECK0(this->is_vec());
    REP (i, size())
      get(i) = ring->add(get(i), a[i]);
  }

  Matrix<T> cadd(const T &a) const {
    Matrix<T> res = this->clone();
    res.scadd(a);
    return res;
  }
  Matrix<T> &scmul(const T &a) {
    REP (i, n)
      REP (j, m)
        get(i, j) = ring->mul(get(i, j), a);
    return *this;
  }

  Matrix<T> &scadd(const T &a) {
    REP (i, n)
      REP (j, m)
        get(i, j) = ring->add(get(i, j), a);
    return *this;
  }

  Matrix<T> &scdiv(const T &a) {
    REP (i, n)
      REP (j, m)
        get(i, j) = ring->div(get(i, j), a);
    return *this;
  }

  Matrix<T> &ssub(const Matrix<T> &b) {
    REP (i, n)
      REP (j, m)
        get(i, j) = ring->sub(get(i, j), b(i, j));
    return *this;
  }

  Matrix<T> &sadd(const Matrix<T> &b) {
    REP (i, n)
      REP (j, m)
        get(i, j) = ring->add(get(i, j), b(i, j));
    return *this;
  }

  Matrix<T> &sdiv(const Matrix<T> &b) {
    REP (i, n)
      REP (j, m)
        get(i, j) = ring->div(get(i, j), b(i, j));
    return *this;
  }

  Matrix<T> &scmul(const Matrix<T> &b) {
    REP (i, n)
      REP (j, m)
        get(i, j) = ring->mul(get(i, j), b(i, j));
    return *this;
  }

  void initialize(const Ring<T> *ring, int n, int m);
  void finalize();
  ~Matrix();
  const T &get(int a) const {
    OPA_CHECK0(n == 1 || m == 1);
    return n == 1 ? get(0, a) : get(a, 0);
  }

  int size() const {
    OPA_CHECK0(this->is_vec());
    return std::max(n, m);
  }

  T &get(int a) {
    OPA_CHECK0(n == 1 || m == 1);
    if (n == 1) return get(0, a);
    return get(a, 0);
  }

  const T &get(int a, int b) const { return (*m_data.get())[a][b]; }
  void set(int a, int b, const T &v) { get(a, b) = v; }
  T dot(const Matrix<T> &a) const {
    T res = ring->getZ();
    REP (i, n)
      REP (j, m)
        res = ring->add(res, ring->mul(get(i, j), a(i, j)));
    return res;
  }

  T &get(int a, int b) {
    OPA_CHECK0(!m_mem.is_ro());
    return (*m_data.get())[a][b];
  }

  T &operator[](int a) { return get(a); }
  const T &operator[](int a) const { return get(a); }
  T &operator()(int a, int b) { return get(a, b); }
  const T &operator()(int a, int b) const { return get(a, b); }

  T &operator()(int a) { return get(a); }
  const T &operator()(int a) const { return get(a); }

  int getN() const;
  int getM() const;

  static Matrix<T> identity(const Ring<T> *ring, int n);
  Matrix<T> identity(int n = -1) const;
  Matrix<T> zeroes(int nr = -1, int nc = -1) const {
    return constant(ring->getZ(), nr, nc);
  }
  Matrix<T> ones(int nr = -1, int nc = -1) const {
    return constant(ring->getE(), nr, nc);
  }
  Matrix<T> constant(const T &cv, int nr = -1, int nc = -1) const;
  static Matrix<T> rand(const Ring<T> *ring, int n, int m);

  template <typename Func>
  static Matrix<T> rand(const Ring<T> *ring, int n, int m, const Func &func) {
    Matrix<T> res(ring, n, m);

    REP (i, n)
      REP (j, m)
        res(i, j) = func();
    return res;
  }

  Matrix<T> &self_elem_addmul(const Matrix<T> &a, const T &c);
  Matrix<T> &self_elem_mul(const Matrix<T> &a);
  Matrix<T> mul(const Matrix<T> &a) const;

  Matrix<T> &self_mul(const Matrix<T> &a) {
    Matrix<T> res = *this * a;
    copy_coeffs(res);
    return *this;
  }

  Matrix<T> add(const Matrix<T> &a) const;
  Matrix<T> sub(const Matrix<T> &a) const;

  void copy_coeffs(const Matrix<T> &a) {
    OPA_CHECK_EQ0(a.n, n);
    OPA_CHECK_EQ0(a.m, m);
    REP (i, n)
      REP (j, m)
        get(i, j) = a(i, j);
  }

  void swap_cols(int a, int b) {
    OPA_CHECK(a < m && b < m, a, b, m);
    if (a == b) return;
    REP (i, m)
      std::swap(get(i, a), get(i, b));
  }

  void swap_rows(int a, int b) {
    OPA_CHECK(a < n && b < n, a, b, n);
    if (a == b) return;
    REP (i, m)
      std::swap(get(a, i), get(b, i));
  }

  void setCol(int colId, const std::vector<T> &tb, int start = 0);
  void setRow(int rowId, const std::vector<T> &tb, int start = 0);

  // DEPRECATED
  Matrix<T> getCol(int c) const { return get_col(c); }
  Matrix<T> getRow(int r) const { return get_row(r); }

  Matrix<T> get_col(int c) const;
  Matrix<T> get_row(int r) const;

  T l_inf() const {
    T res = ring->getZ();
    REP (i, n)
      REP (j, m)
        if (ring->compareRank(res, get(i, j))) res = get(i, j);
    return ring->abs(res);
  }

  template <class Container> void set_col(int c, const Container &container) {
    REP (i, n)
      get(i, c) = container[i];
  }

  template <class Container> void set_row(int r, const Container &container) {
    REP (i, m)
      get(r, i) = container[i];
  }

  Matrix<T> transpose() const;

  Matrix<T> operator*(const Matrix<T> &b) const { return mul(b); }
  Matrix<T> operator+(const Matrix<T> &b) const { return add(b); }
  Matrix<T> operator-(const Matrix<T> &b) const { return sub(b); }
  bool operator==(const Matrix<T> &b) const {
    if (n != b.n || m != b.m) return false;
    REP (i, n)
      REP (j, m)
        if (get(i, j) != b(i, j)) return false;
    return true;
  }

  Matrix<T> kernel_basis() const {
    Matrix<T> a(ring, n, m + n);
    REP (i, n)
      REP (j, m)
        a(i, j) = get(i, j);
    REP (i, n)
      a(i, m + i) = ring->getE();
    int d = a.row_echelon(n, m);
    if (d == n) return Matrix<T>();
    Matrix<T> res = a.get_submatrix(d, m).clone();

    return res;
  }

  Matrix<T> null_basis() const {
    Matrix<T> a = this->clone();
    int d = a.row_echelon();
    int null_dim = m - d;

    Matrix<T> res(ring, null_dim, m);
    std::vector<int> pos;
    int count = 0;
    int c = 0;
    OPA_DISP0(a);

    REP (i, d + 1) {
      int lastc = c;
      for (; c < m && ring->isZ(a(i, c)); ++c)
        ;
      FOR (j, lastc, c) {
        res(count, j) = ring->getME();
        REP (k, i)
          res(count, pos[k]) = ring->div(a(k, j), a(k, pos[k]));
        ++count;
      }
      pos.push_back(c);
      ++c;
    }
    return res;
  }

  bool reduce(Matrix<T> &other, int last_col = -1, bool early_exit = true,
              Matrix<T> *output = nullptr) const;

  bool left_inverse(Matrix<T> *res) const {
    Matrix<T> cur;
    if (!(this->transpose() * *this).invert(&cur)) return false;
    *res = cur * this->transpose();
    return true;
  }

  Matrix<T> inverse() const {
    Matrix<T> res;
    OPA_CHECK0(this->invert(&res));
    return res;
  }

  bool invert(Matrix<T> *res, T *det = nullptr) const;
  bool is_null() const {
    REP (i, n)
      REP (j, m)
        if (!ring->isZ(this->get(i, j))) return false;
    return true;
  }

  template <typename Container> std::vector<T> eval(const Container &x) const {
    OPA_CHECK0(x.size() == m);

    std::vector<T> res(n);
    for (int i = 0; i < n; ++i) {
      T u = ring->getZ();
      for (int j = 0; j < m; ++j) u = ring->add(u, ring->mul(x[j], get(i, j)));
      res[i] = u;
    }
    return res;
  }

  template <typename Container> std::vector<T> evalT(const Container &x) const {
    OPA_CHECK0(x.size() == n);

    std::vector<T> res(m);
    REP (i, m) {
      T u = ring->getZ();
      REP (j, n)
        u = ring->add(u, ring->mul(x[j], get(j, i)));
      res[i] = u;
    }
    return res;
  }

  std::vector<T> solve(std::vector<T> y) const;
  T get_det() const;
  T get_det_slow() const;
  T get_det_row_echelon();
  bool is_square() const { return n == m; }
  bool is_vec() const { return n == 1 || m == 1; }
  int rank(int n2, int m2) const;
  int row_echelon(int n2 = -1, int m2 = -1);
  int row_echelon_gcd(T *gcd_mul = nullptr, int n2 = -1, int m2 = -1);

  const Ring<T> *get_ring() const { return ring; }

  T trace() const;

  std::string str(int nrows = -1, bool header = true) const;
  void disp() const;

  Matrix<T> fromvec(const std::vector<T> &vec, VecType type = VecType::Col) {
    return Matrix<T>::fromvec(ring, vec, type);
  }

  Matrix<T> from(const Ring<T> *ring, const std::vector<std::vector<T> > &vec) {
    Matrix<T> res;
    res.initialize(ring, vec.size(), vec.size() > 0 ? vec[0].size() : 0);
    REP (i, res.n)
      REP (j, res.m)
        res(i, j) = vec[i][j];
    return res;
  }

  static Matrix<T> fromzerovec(const Ring<T> *ring, int n,
                               VecType type = VecType::Col) {
    Matrix<T> res;
    if (type == VecType::Row)
      res.initialize(ring, 1, n);
    else
      res.initialize(ring, n, 1);
    REP (i, n)
      res(i) = ring->getZ();
    return res;
  }

  static Matrix<T> fromvec(const Ring<T> *ring, const std::vector<T> &vec,
                           VecType type = VecType::Col) {
    Matrix<T> res = Matrix<T>::fromzerovec(ring, vec.size(), type);
    REP (i, vec.size())
      res(i) = vec[i];
    return res;
  }

  std::vector<T> tovec() const {
    OPA_CHECK0(this->is_vec());
    std::vector<T> res;
    REP (i, std::max(n, m))
      res.pb(get(i));
    return res;
  }
  OPA_DECL_COUT_OPERATOR(Matrix<T>);

  Matrix<T> &affect(Matrix<T> &a) {
    set_mutable_from(a);
    return *this;
  }

  Matrix<T> &affect(Matrix<T> &&a) {
    set_mutable_from(a);
    return *this;
  }

  Matrix<T> &affect_const(const Matrix<T> &a) {
    set_from(a);
    return *this;
  }

  void add_row() {
    OPA_CHECK0(!m_mem.is_ro());
    ++n;
    m_mem.objs->resize(n * m);
    // wont work with mutable stuff
    m_data->push_back(nullptr);
    reset_ptrs();
  }

  template <class RingType>
  Matrix<typename RingType::Type> lift(const RingType &r) const {
    Matrix<typename RingType::Type> res(&r, n, m);
    REP (i, n)
      REP (j, m)
        res(i, j) = r.import(get(i, j));
    return res;
  }

  template <class RingType>
  Matrix<typename RingType::BaseType> project(const RingType &r) const {
    Matrix<typename RingType::BaseType> res(r.get_base_ring(), n, m);
    REP (i, n)
      REP (j, m)
        res(i, j) = r.project(get(i, j));
    return res;
  }

  std::map<int, int> row_echelon_map;

  const Ring<T> *ring;
  int n, m; // n rows, m cols
private:
  void reset_ptrs() {
    REP (i, n)
      (*m_data.get())[i] = m_mem.get_ptr(m_off.r() + i, m_off.c());
  }

  Matrix<T> &operator=(const Matrix<T> &a) = delete;

  SPTR(std::vector<T *>) m_data;

  /*
  struct Offset {
    int offr;
    int offc;

    Offset operator+(const Offset &peer) const {
      Offet res;
      res.offr = offr + peer.offr;
      res.offc = offc + peer.offc;
      return res;
    }
  };
  */
  OPA_NAMED_PAIR(Offset, r, c, int, int);

  struct MatrixMem {
    SPTR(std::vector<T>) objs;
    int stride;
    bool ro;

    void reset() { objs.reset(); }
    bool is_ro() const { return ro; }

    void set_const(const MatrixMem &mem) {
      stride = mem.stride;
      objs = mem.objs;
      ro = true;
    }

    void set(const MatrixMem &mem) {
      set_const(mem);
      ro = false;
    }

    void set_own(std::vector<T> *tb, int stride) {
      objs.reset(tb);
      this->stride = stride;
      ro = false;
    }

    T *get_ptr(int row, int col) {
      return (T *)objs->data() + stride * row + col;
    }

    const T *get_ptr(int row, int col) const {
      return (const T *)objs->data() + stride * row + col;
    }
  };

  Offset m_off;

  MatrixMem m_mem;
};

template <class T> void Matrix<T>::init_data() {
  if (n && m) {
    m_data.reset(new std::vector<T *>(n));
    reset_ptrs();
  }
}

template <class T> void Matrix<T>::init_objs() {
  if (n && m) {
    m_mem.set_own(new std::vector<T>(n * m, ring->getZ()), m);
  }
  m_off = Offset(0, 0);
}

template <class T>
void Matrix<T>::prepare_set(const Matrix<T> &a, int r0, int c0, int nr,
                            int nc) {
  if (nr == -1) nr = a.getN() - r0;
  if (nc == -1) nc = a.getM() - c0;
  OPA_CHECK0(r0 + nr <= a.n);
  OPA_CHECK0(c0 + nc <= a.m);
  n = nr;
  m = nc;
  ring = a.ring;
  row_echelon_map = a.row_echelon_map;

  m_off = a.m_off + Offset(r0, c0);
}

template <class T> void Matrix<T>::copy(const Matrix<T> &a) {
  finalize();
  prepare_set(a, 0, 0, -1, -1);
  m_off = Offset(0, 0);
  init_objs();
  init_data();
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < m; ++j) get(i, j) = a(i, j);
}

template <class T>
void Matrix<T>::set_from(const Matrix<T> &a, int r0, int c0, int nr, int nc) {
  finalize();
  m_mem.set_const(a.m_mem);
  prepare_set(a, r0, c0, nr, nc);
  init_data();
}

template <class T>
void Matrix<T>::set_mutable_from(Matrix<T> &&a, int r0, int c0, int nr,
                                 int nc) {
  finalize();
  m_mem.set(a.m_mem);
  prepare_set(a, r0, c0, nr, nc);
  init_data();
}

template <class T>
void Matrix<T>::set_mutable_from(Matrix<T> &a, int r0, int c0, int nr, int nc) {
  finalize();
  m_mem.set(a.m_mem);
  prepare_set(a, r0, c0, nr, nc);
  init_data();
}

template <class T>
Matrix<T> Matrix<T>::get_mutable(int r0, int c0, int nr, int nc) {
  Matrix<T> res;
  res.set_mutable_from(*this, r0, c0, nr, nc);
  return res;
}
template <class T>
Matrix<T> Matrix<T>::get_submatrix(int r0, int c0, int nr, int nc) const {
  Matrix<T> res;
  res.set_from(*this, r0, c0, nr, nc);
  return res;
}

template <class T> Matrix<T> Matrix<T>::get_mutable_col(int c) {
  Matrix<T> res;
  res.set_mutable_from(*this, 0, c, -1, 1);
  return res;
}

template <class T> Matrix<T> Matrix<T>::get_mutable_row(int r) {
  Matrix<T> res;
  res.set_mutable_from(*this, r, 0, 1, -1);
  return res;
}

template <class T> Matrix<T> Matrix<T>::get_col(int c) const {
  Matrix<T> res;
  res.set_from(*this, 0, c, -1, 1);
  return res;
}

template <class T> Matrix<T> Matrix<T>::get_row(int r) const {
  Matrix<T> res;
  res.set_from(*this, r, 0, 1, -1);
  return res;
}

template <class T> Matrix<T>::Matrix(const Ring<T> *ring, int n, int m) {
  initialize(ring, n, m);
}

template <class T>
void Matrix<T>::initialize(const Ring<T> *ring, int n, int m) {
  OPA_CHECK0(ring != nullptr);
  this->ring = ring;
  this->n = n;
  this->m = m;
  m_off = Offset(0, 0);
  init_objs();
  init_data();
}

template <class T> void Matrix<T>::finalize() {
  m_data.reset();
  m_mem.reset();
}

template <class T> Matrix<T>::~Matrix() { finalize(); }

template <class T> int Matrix<T>::getN() const { return n; }

template <class T> int Matrix<T>::getM() const { return m; }

template <class T> Matrix<T> Matrix<T>::identity(const Ring<T> *ring, int n) {
  Matrix<T> res(ring, n, n);
  for (int i = 0; i < n; ++i) res.get(i, i) = ring->getE();
  return res;
}

template <class T>
Matrix<T> Matrix<T>::rand(const Ring<T> *ring, int n, int m) {
  Matrix<T> res(ring, n, m);

  for (int i = 0; i < n; ++i)
    for (int j = 0; j < m; ++j) res.get(i, j) = ring->getRand();
  return res;
}

template <class T> bool Matrix<T>::invert(Matrix<T> *res, T *res_det) const {
  OPA_CHECK0(res || res_det);
  T det = ring->getE();

  OPA_CHECK0(n == m);
  if (res) *res = Matrix<T>::identity(ring, n);

  Matrix<T> cur = this->clone();

  for (int i = 0; i < n; ++i) {
    int u = -1;

    T curGcd = ring->getZ();
    if (ring->isField()) {
      for (int j = i; j < n; ++j) {
        T &x = cur.get(j, i);
        if (ring->isZ(x)) continue;
        curGcd = ring->gcd(curGcd, x);
        if (ring->isInv(curGcd)) {
          u = j;
          break;
        }
      }
      if (u == -1) {
        return false;
      }

      if (i != u) {
        det = ring->neg(det);
        for (int j = i; j < n; ++j) std::swap(cur.get(i, j), cur.get(u, j));
        if (res)
          for (int j = 0; j < n; ++j) std::swap(res->get(i, j), res->get(u, j));
      }
      det = ring->mul(det, cur.get(i, i));
      OPA_CHECK0(ring->isInv(cur.get(i, i)));
      T iv = ring->inv(cur.get(i, i));

      // cur.get(i,i) = ring->mul(cur.get(i,i), iv);
      // std::cout<<iv<<std::endl;
      // std::cout<<cur.get(i,i)<<std::endl;
      // exit(0);
      for (int j = 0; j < n; ++j) cur.get(i, j) = ring->mul(cur.get(i, j), iv);
      if (res)
        for (int j = 0; j < n; ++j)
          res->get(i, j) = ring->mul(res->get(i, j), iv);
    } else {
      std::vector<T> vals;
      FOR (j, i, n)
        vals.push_back(cur.get(j, i));
      auto plan = get_gcd_plan(*ring, vals);
      // OPA_DISP0(plan.steps, plan.pos, plan.gcd, vals);
      OPA_CHECK0(!ring->isZ(plan.gcd));
      for (auto mat : { &cur, res }) {
        if (mat == nullptr) continue;
        for (const auto &step : plan.steps) {
          mat->row_op(i + step.r1, i + step.r2, step.a11, step.a12, step.a21,
                      step.a22);
        }
        if (plan.pos != 0) {
          mat->swap_rows(i, plan.pos + i);
          det = ring->neg(det);
        }
      }

      if (!ring->isInv(plan.gcd)) return false;

      OPA_CHECK(plan.gcd == cur.get(i, i), plan.gcd, cur, vals);
      if (!ring->isE(plan.gcd)) {
        T mul = ring->inv(plan.gcd);
        for (auto mat : { &cur, res }) {
          if (mat == nullptr) continue;
          mat->get_mutable_row(i).scmul(mul);
        }
        det = ring->mul(det, mul);
      }
    }

    OPA_CHECK(ring->isE(cur.get(i, i)), cur.get(i, i));

    for (int j = 0; j < n; ++j) {
      if (i != j) {
        T coeff = ring->neg(cur.get(j, i));
        if (!ring->isZ(coeff)) {
          for (int k = i; k < n; ++k)
            cur.get(j, k) =
              ring->add(cur.get(j, k), ring->mul(coeff, cur.get(i, k)));
          if (res)
            for (int k = 0; k < n; ++k)
              res->get(j, k) =
                ring->add(res->get(j, k), ring->mul(coeff, res->get(i, k)));
          OPA_CHECK0(ring->isZ(cur.get(j, i)));
        }
      }
    }
  }

  if (res_det) *res_det = det;

  return true;
}

template <class T> Matrix<T> &Matrix<T>::self_elem_mul(const Matrix<T> &a) {
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < m; ++j) get(i, j) = ring->mul(get(i, j), a.get(i, j));
  return *this;
}

template <class T>
Matrix<T> &Matrix<T>::self_elem_addmul(const Matrix<T> &a, const T &c) {
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < m; ++j)
      get(i, j) = ring->add(get(i, j), ring->mul(c, a.get(i, j)));
  return *this;
}

template <class T> Matrix<T> Matrix<T>::mul(const Matrix<T> &a) const {
  OPA_CHECK(m == a.n, m, n, a.n, a.m);

  Matrix<T> res(ring, n, a.m);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < a.m; ++j) {
      T v = ring->getZ();
      for (int k = 0; k < m; ++k)
        v = ring->add(v, ring->mul(get(i, k), a.get(k, j)));
      res.get(i, j) = v;
    }
  return res;
}

template <class T> Matrix<T> Matrix<T>::add(const Matrix &a) const {
  OPA_CHECK(m == a.m, m, n, a.n, a.m);
  OPA_CHECK(n == a.n, m, n, a.n, a.m);

  Matrix<T> res(ring, n, m);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < a.m; ++j)
      res.get(i, j) = ring->add(get(i, j), a.get(i, j));
  return res;
}

template <class T> Matrix<T> Matrix<T>::sub(const Matrix &a) const {
  OPA_CHECK0(m == a.m);
  OPA_CHECK0(n == a.n);

  Matrix<T> res(ring, n, m);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < a.m; ++j)
      res.get(i, j) = ring->sub(get(i, j), a.get(i, j));
  return res;
}

template <class T> std::string Matrix<T>::str(int nrows, bool header) const {
  std::ostringstream ss;
  ss << std::hex << std::showbase;
  if (nrows == -1) nrows = n;
  if (header)
    ss << utils::stdsprintf("# Matrix(%d %d, shorten=%d)\n", n, m, nrows);
  ss << '[';
  for (int i = 0; i < nrows; ++i) {
    ss << '[';
    for (int j = 0; j < m; ++j) {
      ss << get(i, j);
      if (j != m - 1) ss << ", ";
    }
    ss << ']';
    if (i != nrows - 1)
      ss << ',';
    else
      ss << ']';
    ss << std::endl;
  }
  return ss.str();
}

template <class T> void Matrix<T>::disp() const { std::cout << str(); }

template <class T> std::vector<T> Matrix<T>::solve(std::vector<T> y) const {
  OPA_CHECK0(y.size() == n);
  OPA_CHECK0(m <= n);
  std::vector<T> x(m, ring->getZ());

  Matrix tmp = this->clone();
  for (int i = 0; i < m; ++i) {
    int u = -1;
    for (int j = i; j < n; ++j)
      if (!ring->isZ(tmp.get(j, i))) {
        if (u != -1 && ring->compareRank(tmp(u, i), tmp(j, i))) continue;
        u = j;
      }

    // OPA_DISP0(tmp, u, i, y);
    if (u == -1) {
      break;
    }

    std::swap(y[i], y[u]);
    for (int j = i; j < m; ++j) std::swap(tmp.get(i, j), tmp.get(u, j));

    if (1 && ring->isInv(tmp.get(i, i))) {
      T iv = ring->inv(tmp.get(i, i));
      for (int j = i; j < m; ++j) tmp.get(i, j) = ring->mul(iv, tmp.get(i, j));
      y[i] = ring->mul(iv, y[i]);
      for (int j = i + 1; j < n; ++j) {
        if (j != i && !ring->isZ(tmp.get(j, i))) {
          T c = ring->neg(tmp.get(j, i));
          for (int k = i; k < m; ++k)
            tmp.get(j, k) =
              ring->add(tmp.get(j, k), ring->mul(tmp.get(i, k), c));
          y[j] = ring->add(y[j], ring->mul(y[i], c));
        }
      }
    } else {

      for (int j = i + 1; j < n; ++j) {
        if (0) {
          while (!ring->isZ(tmp(j, i))) {
            T q;
            ring->ediv(tmp(i, i), tmp(j, i), &q, nullptr);
            y[i] = ring->sub(y[i], ring->mul(y[j], q));
            std::swap(y[i], y[j]);
            for (int k = i; k < m; ++k) {
              tmp(i, k) = ring->sub(tmp(i, k), ring->mul(tmp(j, k), q));
              std::swap(tmp(i, k), tmp(j, k));
            }
          }
        } else {

          T u, v, d;
          d = ring->egcd(tmp(i, i), tmp(j, i), u, v);
          OPA_CHECK(ring->isInv(ring->gcd(u, v)), u, v);
          if (!ring->isInv(u)) {
            for (int k = i; k < m; ++k) std::swap(tmp(i, k), tmp(j, k));
            std::swap(u, v);
            std::swap(y[i], y[j]);
          }
          OPA_CHECK(ring->isInv(u), u, v, d, tmp(i, i), tmp(j, i));
          OPA_DISP0(ring->gcd(u, v), u, v, d);

          OPA_CHECK(!ring->isZ(u), tmp(i, i), tmp(j, i), u, v, d);
          T c = ring->neg(ring->div(tmp(j, i), d));
          for (int k = i; k < m; ++k) {
            tmp(i, k) =
              ring->add(ring->mul(tmp(i, k), u), ring->mul(tmp(j, k), v));
            tmp(j, k) = ring->add(tmp(j, k), ring->mul(c, tmp(i, k)));
          }
          OPA_CHECK(ring->isZ(tmp(j, i)), i, j, u, v, tmp(i, i), tmp(j, i), d,
                    c);
          y[i] = ring->add(ring->mul(y[i], u), ring->mul(y[j], v));
          y[j] = ring->add(y[j], ring->mul(y[i], c));
        }
      }
    }
  }

  OPA_DISP0(tmp, y);
  // solve failed
  for (int i = m; i < n; ++i)
    if (!ring->isZ(y[i])) {
      return {};
    }

  REPV (i, m) {
    T sum = ring->getZ();
    FOR (j, i + 1, m)
      sum = ring->add(sum, ring->mul(tmp(i, j), x[j]));
    sum = ring->add(y[i], ring->neg(sum));
    if (ring->isZ(tmp(i, i))) {
      if (!ring->isZ(sum)) return {};
      x[i] = ring->getZ();
      continue;
    }

    if (!ring->isZ(ring->mod(sum, tmp(i, i)))) return {};
    x[i] = ring->div(sum, tmp(i, i));
    OPA_CHECK(ring->mul(tmp(i, i), x[i]) == sum, i, x[i], tmp(i, i), y[i], sum);
  }
  OPA_DISP0(tmp.eval(x), y);

  return x;
}

template <class T> int Matrix<T>::row_echelon(int n2, int m2) {
  if (n2 == -1) n2 = n;
  if (m2 == -1) m2 = m;
  row_echelon_map.clear();

  OPA_CHECK0(ring->isField());
  int col = 0;
  int nc = 0;
  for (int i = 0; i < n2;) {
    int cur = col++;
    if (cur == m2) return i;

    int u = -1;
    FOR (j, i, n2)
      if (ring->isInv(get(j, cur))) {
        u = j;
        break;
      }
    if (u == -1) continue;
    ++nc;
    FOR (j, cur, getM())
      std::swap(get(u, j), get(i, j));
    row_echelon_map[cur] = i;

    T inv = ring->inv(get(i, cur));
    FOR (j, cur, getM())
      get(i, j) = ring->mul(inv, get(i, j));
    REP (j, n2)
      if (j != i) {
        if (ring->isZ(get(j, cur))) continue;
        T coeff = ring->neg(get(j, cur));
        FOR (k, cur, getM())
          get(j, k) = ring->add(get(j, k), ring->mul(get(i, k), coeff));
      }

    ++i;
  }
  return nc;
}

template <class T> int Matrix<T>::row_echelon_gcd(T *det_res, int n2, int m2) {
  // Need to make use of step.coeff_mul and step.final_mul. Invalid in UFD

  if (n2 == -1) n2 = n;
  if (m2 == -1) m2 = m;
  row_echelon_map.clear();
  T det;

  int col = 0;
  int nc = 0;
  for (int i = 0; i < n2;) {
    int cur = col++;
    if (cur == m2) return i;

    std::vector<T> vals;
    FOR (j, i, n2) { vals.push_back(get(j, cur)); }
    auto plan = get_gcd_plan(*ring, vals);
    if (plan.pos == -1) {
      det = ring->getZ();
      continue;
    }

    ++nc;
    for (const auto &step : plan.steps) {
      row_op(i + step.r1, i + step.r2, step.a11, step.a12, step.a21, step.a22);
    }
    if (plan.pos != 0) {
      swap_rows(i, plan.pos + i);
      det = ring->neg(-det);
    }
    T cv = get(i, i);
    det = ring->mul(det, cv);
    row_echelon_map[cur] = i;
    FOR (j, i + 1, n2) {
      T coeff = -get(j, i) / cv;
      if (ring->isZ(coeff)) continue;
      get_mutable_row(j).self_elem_addmul(get_row(i), coeff);
    }
    ++i;
  }

  if (det_res) *det_res = det;
  return nc;
}

template <class T> T Matrix<T>::get_det_row_echelon() {
  OPA_CHECK0(is_square());
  T det_mul;
  this->row_echelon_gcd(&det_mul);
  OPA_DISP0(*this);
  T det = ring->getE();
  REP (i, getN()) {
    det = ring->mul(det, get(i, i));
    OPA_DISP0(i, get(i, i), det);
  }
  return det;
}

template <class T> T Matrix<T>::get_det() const {
  OPA_CHECK0(is_square());
  T res;
  bool ok = invert(0, &res);
  if (!ok) return ring->getZ();
  return res;
}

template <class T>
void Matrix<T>::setCol(int colId, const std::vector<T> &tb, int start) {
  OPA_CHECK0(start + tb.size() <= n);
  REP (i, tb.size())
    get(start + i, colId) = tb[i];
}

template <class T>
void Matrix<T>::setRow(int rowId, const std::vector<T> &tb, int start) {
  OPA_CHECK0(start + tb.size() <= m);
  REP (i, tb.size())
    get(rowId, start + i) = tb[i];
}

template <class T> Matrix<T> Matrix<T>::transpose() const {
  Matrix<T> res(ring, m, n);
  REP (i, n)
    REP (j, m)
      res.get(j, i) = get(i, j);
  return res;
}

template <class T> int Matrix<T>::rank(int n2, int m2) const {
  Matrix<T> x = this->clone();
  return x.row_echelon(n2, m2);
}

template <class T> T Matrix<T>::trace() const {
  OPA_CHECK0(n == m);
  T res = ring->getZ();
  REP (i, n)
    res = ring->add(res, get(i, i));
  return res;
}

template <class T> Matrix<T> Matrix<T>::identity(int n) const {
  if (n == -1) n = this->n;
  return Matrix<T>::identity(ring, n);
}

template <class T>
Matrix<T> Matrix<T>::constant(const T &cv, int nr, int nc) const {
  if (nr == -1) nr = getN();
  if (nc == -1) nc = getM();
  Matrix<T> res;
  res.initialize(ring, nr, nc);
  REP (i, nr)
    REP (j, nc)
      res(i, j) = cv;
  return res;
}

template <class T>
bool Matrix<T>::reduce(Matrix<T> &other, int last_col, bool early_exit,
                       Matrix<T> *output) const {
  OPA_CHECK0(other.getM() <= getM());
  // should be a reduced row echelon matrix
  // That is, row_echelon must have been called before or row_echelon_map must
  // be set.
  if (last_col == -1) last_col = getM();
  bool ok = true;
  if (output != nullptr) {
    output->initialize(this->get_ring(), other.getN(), this->getN());
  }

  REP (i, other.getN()) {
    REP (j, last_col) {
      if (ring->isZ(other(i, j))) continue;
      if (!row_echelon_map.count(j)) {
        if (early_exit) {
          // OPA_DISP("Fail reduce at ", i, j, last_col);
          return false;
        }
        ok = false;
        break;
      }

      int pos = row_echelon_map.find(j)->ND;
      T coeff = ring->neg(ring->div(other.get(i, j), get(pos, j)));
      other.get_mutable_row(i).self_elem_addmul(get_row(pos), coeff);

      if (output) {
        output->get(i, pos) = ring->neg(coeff);
      }
    }
    REP (j, last_col)
      if (!ring->isZ(other(i, j))) {
        if (early_exit) return false;
        ok = false;
      }
  }
  return ok;
}

template <class T> T Matrix<T>::get_det_slow() const {
  OPA_CHECK0(is_square());
  T res = ring->getZ();
  auto perm = utils::Range<int>::StepRange(0, n, 1).tb();

  do {
    T entry = ring->getE();
    REP (i, n)
      entry = ring->mul(entry, get(i, perm[i]));
    if (get_permutation_signature(perm) == -1) entry = ring->neg(entry);
    res = ring->add(res, entry);
  } while (std::next_permutation(ALL(perm)));
  return res;
}

template <class T>
void Matrix<T>::set_submatrix(const Matrix<T> &a, int r0, int c0) {
  OPA_CHECK0(r0 + a.getN() <= getN());
  OPA_CHECK0(c0 + a.getM() <= getM());

  REP (i, a.getN())
    REP (j, a.getM()) { get(r0 + i, c0 + j) = a(i, j); }
}

OPA_NM_MATH_COMMON_END

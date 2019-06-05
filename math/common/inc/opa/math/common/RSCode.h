#pragma once

#include <opa/math/common/CyclicCode.h>
#include <opa/math/common/FFT.h>
#include <opa/math/common/GF_q.h>
#include <opa/math/common/UtilsGFq.h>
#include <opa/math/common/UtilsPoly.h>
#include <opa/math/common/algo.h>

OPA_NAMESPACE_DECL3(opa, math, common)

// reed solomon code
template <class U> class RSCode {
public:
  typedef Poly<U> T;
  typedef Poly<T> PT;
  RSCode() {}

  RSCode(const GF_q<U> *field, int n, int k) { init(field, n, k); }

  void init(const GF_q<U> *field, int n, int k) {
    m_field = field;
    this->m_n = n;
    this->m_k = k;

    m_t = (m_n - m_k) / 2;
    m_pr.init(m_field);
    advanced_init(m_field->getPrimitiveElem(), 1);
  }

  void advanced_init(const T &w, int pw_chk) {
    m_g = m_pr.getE();
    m_pw_chk = pw_chk;
    m_w = w;
    T v = m_field->faste(m_w, pw_chk);
    REP (i, 2 * m_t) {
      m_g = m_g * m_pr.import_vec({ -v, m_field->getE() }); // x - a^i
      v = v * m_w;
    }
  }

public:
  virtual std::vector<T> encode(const std::vector<T> &msg) {
    std::vector<T> cv(m_n, m_field->getZ());
    OPA_CHECK0(msg.size() == m_k);
    std::copy(ALL(msg), cv.begin() + (m_n - m_k));
    PT c = m_pr.import_vec(cv);
    c = c - (c % m_g);
    return m_field->import_vecs(c.to_vec(m_n));
  }

  bool decode_correct(std::vector<T> &rv) const {
    PT r = m_pr.import_vec(rv);
    m_nerrs = 0;

    std::vector<T> syndroms(2 * m_t);
    bool nz = false;
    std::vector<T> pts;
    {
      T v = m_field->faste(m_w, m_pw_chk);
      REP (i, 2 * m_t) {
        pts.pb(v);
        syndroms[i] = r(v);
        nz |= !m_field->isZ(syndroms.back());
        v = v * m_w;
      }
    }
    std::vector<T> locator_poly_tb;

    locator_poly_tb = findMinLinearRecursion_Massey(*m_field, syndroms, m_t);

    PT locator_poly = m_pr.import_vec(locator_poly_tb);

    if (locator_poly.deg() <= 0) return !nz;

    std::vector<int> err_pos;
    std::vector<T> wlist;
    {
      T v = m_field->faste(m_w, m_pw_chk);
      REP (i, m_n) {
        if (locator_poly(v) == m_field->getZ()) {
          err_pos.push_back((2 * m_n - m_pw_chk - i) % m_n);
          wlist.push_back(m_field->inv(v));
        }
        v = v * m_w;
      }
    }

    if (locator_poly.deg() != err_pos.size()) return false;
    int nerr = err_pos.size();
    m_nerrs = nerr;

    Matrix<T> err_mat(m_field, nerr, nerr);
    REP (i, nerr) {
      REP (j, nerr)
        err_mat(i, j) = m_field->faste(pts[i], err_pos[j]);
    }

    std::vector<T> b_vec(syndroms.begin(), syndroms.begin() + nerr);
    std::vector<T> err_tb = err_mat.solve(b_vec);
    REP (i, nerr)
      rv[err_pos[i]] = rv[err_pos[i]] - err_tb[i];
    return true;
  }

  std::vector<T> decode(const std::vector<T> &rv) const {
    std::vector<T> c = rv;
    if (!decode_correct(c)) return {};
    return std::vector<T>(c.begin() + m_n - m_k, c.end());
  }
  OPA_ACCESSOR_R(int, m_n, n);
  OPA_ACCESSOR_R(int, m_k, k);
  OPA_ACCESSOR_R(int, m_t, t);
  OPA_ACCESSOR_R(PT, m_g, g);
  OPA_ACCESSOR_R(T, m_w, w);
  OPA_ACCESSOR_PTR_R(GF_q<U>, m_field, field);
  OPA_ACCESSOR_R(int, m_field->p(), p);
  OPA_ACCESSOR_R(int, m_field->q(), q);
  OPA_ACCESSOR_R(int, m_nerrs, nerrs);

private:
  PT m_g;
  T m_w;
  int m_n;
  int m_k;
  int m_t;
  int m_pw_chk;
  mutable int m_nerrs;

  const GF_q<U> *m_field;
  PolyRing<T> m_pr;
};

typedef RSCode<u32> RSCode_u32;
typedef RSCode<Poly<u32> > RSCode_P_u32;

OPA_NAMESPACE_DECL3_END

#pragma once

#include <opa/crypto/base.h>
#include <opa/math/common/Field.h>
#include <opa/math/common/Poly.h>
#include <opa/math/common/PolyRing.h>
#include <opa/math/common/Ring_col.h>
#include <opa/math/common/Types.h>
#include <opa/math/common/UtilsPoly.h>
#include <opa/math/common/algo.h>
#include <opa_common.h>

OPA_NM_CRYPTO
template <class T> class LFSR_Base {
public:
  T get_next() {
    OPA_STATIC_ASSERT_IMPL(T);
    return T();
  }
  int size() const {
    OPA_STATIC_ASSERT_IMPL(T);
    return -1;
  }
};

template <class T> class LFSR : public LFSR_Base<T> {
public:
  typedef OPA_MATH::Poly<T> TExt;
  LFSR() {}

  void init(const OPA_MATH::Field<T> *field, const TExt &state,
            const OPA_MATH::Poly<T> &poly) {
    m_field = field;
    m_pr.init(m_field);
    m_poly = m_pr.import(poly);
    m_ring.init(m_field, m_poly);
    m_state = m_ring.import(state);
    m_x = m_pr.x();
    OPA_CHECK0(m_poly.deg() >= 1);
    m_n = m_poly.deg();
  }

  LFSR(const OPA_MATH::Field<T> *field, const TExt &state,
       const OPA_MATH::Poly<T> &poly) {
    init(field, state, poly);
  }

  const TExt &get_xinv() const {
    if (!m_xinv_init) {
      // if (!m_pr.isZ(m_pr.mod(m_x, m_poly))) {
      m_xinv = m_ring.inv(m_x);
      m_xinv_init = true;
    }
    return m_xinv;
  }

  void set_prev() { m_state = m_ring.mul(get_xinv(), m_state); }
  void advance(opa::math::common::bignum v) {
    TExt y = m_x;
    if (v < 0) {
      y = get_xinv();
      v = -v;
    }
    m_state = m_ring.mul(m_state, m_ring.faste(y, v));
  }

  T get_next() {
    T res = m_field->neg(m_state.get_safe(m_poly.deg() - 1));
    m_state = m_ring.mul(m_x, m_state);
    return res;
  }
  T get_next_v2() {
    std::vector<T> tmp = m_state.to_vec(m_n);
    REP(i, m_n) tmp[i] = m_state.get_safe(i+1);
    T output = m_state.get(0);
    FOR(i,1, m_n) tmp[i] = m_field->add(tmp[i], m_field->mul(m_poly.get_safe(i+1), output));
    m_state = m_pr.import(tmp);
    return output;
  }

  T get_next_non_galois() {
    T res = m_field->getZ();
    m_state = m_pr.mul(m_x, m_state);
    REP (i, m_poly.deg() + 1)
      res = m_field->add(res, m_field->mul(m_poly[i], m_state.get_safe(i)));
    m_pr.set1(m_state, m_poly.deg(), m_field->getZ());
    m_pr.set1(m_state, 0, res);
    return res;
  }

  T get_next_non_galois2() {
    T res = m_field->getZ();
    T output = m_state.get_safe(0);

    std::vector<T> tmp = m_state.to_vec(m_n);
    REP (i, m_poly.deg()) {
      res = m_field->add(res, m_field->mul(m_poly[i], m_state.get_safe(i)));
      tmp[i] = m_state.get_safe(i+1);
    }
    tmp[m_poly.deg()-1] = res;
    m_state = m_pr.import(tmp);
    return output;
  }


  T get_next_non_galois_gps_style() {
    T res = m_state.get_safe(m_poly.deg() - 1);
    get_next_non_galois();
    return res;
  }

  OPA_MATH::bignum max_period() const { return m_ring.order; }

  TExt get_initial_state(const std::vector<T> &obs) const {
    OPA_CHECK0(obs.size() >= m_n);
    TExt base_state = m_pr.getZ();
    TExt cur_toggle = m_pr.getZ();
    REP (i, m_n) {
      cur_toggle = m_pr.mul(cur_toggle, m_x);
      T diff = m_field->sub(m_field->neg(obs[i]), cur_toggle.get_safe(m_n));
      m_pr.set1(base_state, m_n - i - 1, diff);
      m_pr.set1(cur_toggle, m_n, m_field->neg(obs[i]));
      cur_toggle = m_pr.mod(cur_toggle, m_poly);
    }
    return base_state;
  }

  int size() const { return m_poly.deg(); }
  bool init_from_seq(const OPA_MATH::Field<T> *field,
                     const std::vector<T> &seq) {

    OPA_MATH::PolyRing<T> pr(field);
    auto vec =
      OPA_MATH::findMinLinearRecursion_Massey(*field, seq, seq.size() / 2);
    if (vec.size() == 0) {
      OPA_DISP0("failed to find recursion");
      return false;
    }
    if (vec.size() == 1) {
      // all zero sequence
      vec.pb(field->getE());
    }
      OPA_DISP0("KAA  ", vec);
    TExt poly = pr.import(vec, true);
    TExt rseq = pr.import(seq, true);

    TExt state = pr.mul(rseq, poly);
    pr.sresize(state, poly.size());
    OPA_DISP0(poly, state);
    this->init(field, state, poly);
    return true;
  }

  const OPA_MATH::Poly<T> get_poly() const { return m_poly; }
  const TExt get_state() const { return m_state; }

  void init_rand(int n, const OPA_MATH::Field<T> *field) {
    OPA_MATH::Poly<T> poly;
    poly = OPA_MATH::find_primitive_poly(*field, n);
    OPA_MATH::PolyRing<T> pr(field);
    TExt init_state = pr.rand(n - 1);
    this->init(field, init_state, poly);
  }
  void set_state(const TExt &state) { m_state = state; }

  OPA_ACCESSOR_R(OPA_MATH::PolyRing<T>, m_pr, pr);

private:
  const OPA_MATH::Field<T> *m_field;
  TExt m_state;
  OPA_MATH::PolyRing<T> m_pr;
  OPA_MATH::Ring_col<T> m_ring;
  OPA_MATH::Poly<T> m_poly;
  mutable OPA_MATH::Poly<T> m_xinv;
  OPA_MATH::Poly<T> m_x;
  int m_n;
  mutable bool m_xinv_init = false;
};

typedef LFSR<u32> LFSR_u32;

OPA_NM_CRYPTO_END

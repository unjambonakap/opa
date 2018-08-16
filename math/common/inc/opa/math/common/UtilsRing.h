#pragma once

#include <opa/math/common/FractionField.h>
#include <opa/math/common/Poly.h>
#include <opa/math/common/Ring.h>
#include <opa/utils/misc.h>

OPA_NM_MATH_COMMON

template <class T> class GcdBuilder {
public:
  const Ring<T> &m_ring;
  T v;

  GcdBuilder(const Ring<T> &ring) : m_ring(ring) { v = ring.getZ(); }
  const T &update(const T &x) {
    v = m_ring.gcd(x, v);
    return v;
  }

  template <class U> const T &update_iter(const U &u) {
    std::for_each(
      ALL(u), std::bind(&GcdBuilder<T>::update, this, std::placeholders::_1));
    return v;
  }
};

template <class T, class U> T gcd_list(const Ring<T> &ring, const U &u) {
  return GcdBuilder<T>(ring).update_iter(u);
}

template <class T, class U> T gcd_simplify(const Ring<T> &ring, U &u) {
  T gcd = gcd_list(ring, u);
  OPA_CHECK(!ring.isZ(gcd), gcd, u);
  for (auto &x : u) x = ring.div(x, gcd);
  return gcd;
}

template <class T, class U> T lcm_list(const Ring<T> &ring, const U &u) {
  T res = ring.getE();
  for (const T &a : u) {
    T d = ring.gcd(a, res);
    res = ring.mul(res, ring.div(a, d));
  }
  return res;
}

template <class T> struct GcdComputationStep {
  int r1, r2;
  T a11, a12;
  T a21, a22;
  T mul;
  bool gcd_on_1;
  OPA_DECL_COUT_OPERATOR2(GcdComputationStep<T>, a.r1, a.r2, a.a11, a.a12,
                          a.a21, a.a22, a.mul, a.gcd_on_1);
};

template <class T> struct GcdComputationPlan {
  std::vector<GcdComputationStep<T> > steps;
  T gcd;
  int pos;
  T final_mul;
  OPA_DECL_COUT_OPERATOR2(GcdComputationPlan<T>, a.steps, a.gcd, a.pos,
                          a.final_mul);
};

template <class T>
Poly<T> compute_gcd_step(const PolyRing<T> &ring, const Poly<T> &a,
                         const Poly<T> &b, GcdComputationStep<Poly<T> > *step) {
  T u, v, u2, v2;
  auto ctx = ring.egcd2_ctx(a, b);

  step->mul = ring.constant(ctx.coeff_mul);
  step->a11 = ctx.u;
  step->a12 = ctx.v;
  step->a21 = ctx.u2;
  step->a22 = ctx.v2;
  return ctx.res * ctx.final_mul;
}

template <class RingType, class T>
T compute_gcd_step(const RingType &ring, const T &a, const T &b,
                   GcdComputationStep<T> *step) {
  T u, v, u2, v2;
  bool gcd_on_1 = true;
  T nxt = ring.egcd2(a, b, u, v, u2, v2, gcd_on_1);
  if (!gcd_on_1) std::swap(u, u2), std::swap(v, v2);

  step->mul = ring.getE();
  step->a11 = u;
  step->a12 = v;
  step->a21 = u2;
  step->a22 = v2;
  step->gcd_on_1 = gcd_on_1;
  if (gcd_on_1) {
    OPA_CHECK(ring.add(ring.mul(step->a11, a), ring.mul(step->a12, b)) == nxt,
              step, a, b);
  } else {

    OPA_CHECK(ring.add(ring.mul(step->a21, a), ring.mul(step->a22, b)) == nxt,
              step, a, b);
  }
  return nxt;
}

template <class RingType, class T>
GcdComputationPlan<T> get_gcd_plan(const RingType &ring,
                                   const std::vector<T> &vals) {
  // TODO: something less costly
  GcdComputationPlan<T> plan;
  T curgcd = ring.getZ();
  int gcd_pos = -1;
  REP (i, vals.size()) {
    GcdComputationStep<T> step;
    T nxt = compute_gcd_step(ring, curgcd, vals[i], &step);
    if (nxt == curgcd) continue;
    curgcd = nxt;
    if (gcd_pos == -1) {
      gcd_pos = i;
      continue;
    }

    step.r1 = gcd_pos;
    step.r2 = i;
    plan.steps.push_back(step);
    if (!step.gcd_on_1) gcd_pos = i;
  }
  plan.gcd = curgcd;
  plan.pos = gcd_pos;
  return plan;
}

template <class T> void reduce_gcd(const Ring<T> &ring, std::vector<T> &vals) {
  auto plan = get_gcd_plan(ring, vals);
  if (plan.pos == -1) return;

  for (const auto &step : plan.steps) {
    auto &v1 = vals[step.r1];
    auto &v2 = vals[step.r2];
    T t1 = ring.add(ring.mul(v1, step.a11), ring.mul(v2, step.a12));
    T t2 = ring.add(ring.mul(v1, step.a21), ring.mul(v2, step.a22));
    v1 = t1;
    v2 = t2;
  }
  OPA_CHECK_EQ(vals[plan.pos], plan.gcd, plan.pos, vals);
}

template <class T>
void reduce_gcd2(const PolyRing<T> &ring, Poly<T> &a, Poly<T> &b) {
  GcdComputationStep<Poly<T> > step;
  Poly<T> nxt = compute_gcd_step(ring, a, b, &step);

  Poly<T> t1 = ring.add(ring.mul(a, step.a11), ring.mul(b, step.a12));
  Poly<T> t2 = ring.add(ring.mul(a, step.a21), ring.mul(b, step.a22));
  a = t1;
  b = t2;
  OPA_CHECK_EQ(a, nxt, a, b, step);
}

template <class T>
std::vector<T> make_integers(const Ring<T> &base_ring,
                             const std::vector<Fraction<T> > &container) {
  std::vector<T> qs;
  for (auto &x : container) qs.push_back(x.q);
  T lcm = lcm_list(base_ring, qs);
  REP (i, container.size()) {
    qs[i] = base_ring.mul(container[i].p, base_ring.div(lcm, container[i].q));
  }
  return qs;
}

OPA_NM_MATH_COMMON_END

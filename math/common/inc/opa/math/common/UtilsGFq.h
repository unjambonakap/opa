#ifndef _H_OPA_MATH_COMMON_UTILS_GFQ
#define _H_OPA_MATH_COMMON_UTILS_GFQ

#include <opa/math/common/Field.h>
#include <opa/math/common/GF_q.h>
#include <opa/math/common/Poly.h>
#include <opa/math/common/PolyRing.h>

OPA_NAMESPACE_DECL3(opa, math, common)

template <class T>
bool isInBaseField(const GF_q<T> *extensionField, const Poly<Poly<T> > &poly) {
  for (int i = 0; i < poly.size(); ++i)
    if (poly.get(i).size() > 1) return false;
  return true;
}

template <class T>
T toBaseField1(const GF_q<T> *extensionField, const Poly<T> &u) {
  typedef T TYPE_BASE;
  typedef Poly<T> TYPE_EXT;

  const Field<TYPE_BASE> *baseField = extensionField->getBaseField();

  OPA_CHECK0(u.size() <= 1);
  if (u.size() == 0)
    return baseField->getZ();
  else
    return u.get(0);
}

template <class T>
Poly<T> toBaseField(const GF_q<T> *extensionField, const Poly<Poly<T> > &poly) {
  typedef T TYPE_BASE;
  typedef Poly<T> TYPE_EXT;

  const Field<TYPE_BASE> *base_field = extensionField->getBaseField();
  Poly<TYPE_BASE> res = base_field->get_poly_ring()->xpw(poly.deg());
  for (int i = 0; i < poly.size(); ++i)
    res.get(i) = toBaseField1(extensionField, poly.get(i));
  return res;
}

template <class T>
Poly<Poly<T> > toExtField(const GF_q<T> *extensionField, const Poly<T> &poly) {
  // DO NOT USE
  OPA_CHECK0(false);
  typedef T TYPE_BASE;
  typedef Poly<T> TYPE_EXT;

  // lifetime scope issues
  PolyRing<TYPE_EXT> pr(extensionField);
  PolyRing<TYPE_BASE> pr2(extensionField->getBaseField());
  Poly<TYPE_EXT> res = pr.xpw(poly.deg());

  for (int i = 0; i < poly.size(); ++i) res.get(i) = pr2.constant(poly.get(i));
  return res;
}

template <typename T>
Poly<T> findMinimalPoly(const GF_q<T> *extensionField, const Poly<T> &wElem) {
  typedef Poly<T> TYPE_EXT;
  typedef T TYPE_BASE;

  const Field<T> *baseField = extensionField->getBaseField();

  int qm1 = extensionField->getSizeU32() - 1;
  int q = baseField->getSizeU32();

  int r = 1;
  TYPE_EXT cur = wElem;

  PolyRing<TYPE_EXT> pr(extensionField);
  Poly<TYPE_EXT> p = pr.getE();

  for (;; ++r) {
    if (r != 1 && cur == wElem) break;
    p = pr.mul(p, pr.add(pr.x(), pr.neg(pr.constant(cur))));
    cur = extensionField->faste(cur, q);
  }

  for (int i = 0; i < p.size(); ++i) assert(p.get(i).size() <= 1);
  Poly<TYPE_EXT> xn1 =
    pr.add(pr.xpw(qm1), pr.neg(pr.constant(extensionField->getE())));

  OPA_CHECK0(pr.mod(xn1, p) == pr.getZ());
  OPA_CHECK0(pr.eval(p, wElem) == extensionField->getZ());

  return toBaseField(extensionField, p);
}

template<typename T>
Poly<T> gfx_char_poly(const Matrix<T> &mat) {
  // Poly ring is not valid anymore after return of function.
  OPA_CHECK0(mat.ring->isField());
  OPA_CHECK0(mat.n == mat.m);
  int n = mat.n;
  Field<T> *field = (Field<T> *)mat.ring;
  typedef GF_q<T> EF;
  EF ef(field, n + 1);
  Matrix<Poly<T> > em(&ef, n, n);
  REP (i, n) {
    REP (j, n) {
      std::vector<T> ex;
      ex.push_back(mat.ring->neg(mat.get(i, j)));
      if (i == j) ex.push_back(field->getE());
      em(i, j) = ef.import_vec(ex);
    }
  }

  Poly<T> res = em.get_det();
  return res;
}


OPA_NAMESPACE_DECL3_END

#endif

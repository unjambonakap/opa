#ifndef _H_OPA_MATH_COMMON_CYCLICCODE
#define _H_OPA_MATH_COMMON_CYCLICCODE

#include <opa/math/common/Field.h>
#include <opa/math/common/PolyRing.h>
#include <opa/math/common/Matrix.h>

OPA_NAMESPACE_DECL3(opa, math, common)

template <class T> class CyclicCode {
private:
  const Field<T> &field;
  int n, nk;
  Poly<T> g, h;
  PolyRing<T> pr;
  Poly<T> xn1;

public:
  int getN() const { return n; }
  int getNK() const { return nk; }
  int getK() const { return n - nk; }
  CyclicCode(const Field<T> &f, const Poly<T> &g, int n);

  std::vector<T> encode(const std::vector<T> &msg) const;
  std::vector<T> decode(const std::vector<T> &msg) const;

  bool isCodeword(const std::vector<T> &msg) const;
};
template <class T>
CyclicCode<T>::CyclicCode(const Field<T> &f, const Poly<T> &g, int n)
    : field(f), pr(f) {
  this->g = g;
  this->n = n;
  nk = n - g.deg();
  xn1 = pr.add(pr.xpw(n), pr.neg(pr.getE()));

  Poly<T> r;
  pr.ediv(xn1, g, &h, &r);
  assert(pr.isZ(r));
}

template <class T>
std::vector<T> CyclicCode<T>::encode(const std::vector<T> &msg) const {
  Poly<T> in = pr.import(msg);
  Poly<T> res = pr.mulmod(in, g, xn1);

  Poly<T> tmp = pr.mulmod(res, h, xn1);

  std::vector<T> ans = res.toVector();
  ans.resize(n, field.getZ());
  return ans;
}

template <class T>
std::vector<T> CyclicCode<T>::decode(const std::vector<T> &msg) const {
  Matrix<T> m(&field, n, nk);
  std::vector<T> tb = g.toVector();

  for (int j = 0; j < nk; ++j)
    for (int i = 0; i < n - nk + 1; ++i)
      m.get(i + j,j) = g.get(i);
  std::vector<T> res = m.solve(msg);
  return res;
}

template <class T>
bool CyclicCode<T>::isCodeword(const std::vector<T> &msg) const {
  Poly<T> res = pr.mulmod(pr.import(msg), h, xn1);
  return pr.isZ(res);
}

OPA_NAMESPACE_DECL3_END

#endif

#pragma once

#include <opa/math/common/UtilsPoly.h>
#include <opa/math/common/GF_q.h>
#include <opa/math/common/UtilsGFq.h>
#include <opa/math/common/CyclicCode.h>
#include <opa/math/common/FFT.h>
#include <opa/math/common/algo.h>

OPA_NAMESPACE_DECL3(opa, math, common)

std::vector<std::vector<int> > findCosets(int q, int n) {
  std::set<int> seen;
  std::vector<std::vector<int> > res;
  for (int i = 0; i < n; ++i) {
    if (seen.count(i))
      continue;
    std::vector<int> cur;
    cur.push_back(i);
    for (int x = (i * q) % n; x != i; x = (x * q) % n)
      cur.push_back(x);
    res.push_back(cur);
    for (int j = 0; j < cur.size(); ++j)
      seen.insert(cur[j]);
  }
  return res;
}

template <class T>
Poly<T> findBCHPoly(const Field<T> *baseField, int n, int j0, int t,
                    GF_q<T> *&extField, Poly<T> &w) {
  typedef T BASE_TYPE;
  typedef Poly<T> EXT_TYPE;

  assert(j0 + 2 * t < n);

  int q = baseField->getSizeU32();
  int m = 1;
  int tmp = q;
  for (; (tmp - 1) % n != 0 && m < 100; ++m, tmp = tmp * q)
    ;
  assert(m != 100);

  printf(">>> %d\n", m);
  extField = new GF_q<BASE_TYPE>(baseField, m);
  w = extField->getNthRoot(n);

  std::vector<std::vector<int> > cosets = findCosets(q, n);
  PolyRing<EXT_TYPE> pr(extField);

  Poly<EXT_TYPE> res = pr.getE();

  for (int i = 0; i < cosets.size(); ++i) {
    for (int j = 0; j < cosets[i].size(); ++j) {
      int e = cosets[i][j];
      if (e >= j0 && e <= j0 + 2 * t) {
        for (int k = 0; k < cosets[i].size(); ++k) {
          EXT_TYPE root = extField->faste(w, cosets[i][k]);
          Poly<EXT_TYPE> linearPoly = pr.add(pr.x(), pr.neg(pr.constant(root)));
          res = pr.mul(res, linearPoly);
        }
        break;
      }
    }
  }
  OPA_CHECK0(isInBaseField(extField, res));
  return toBaseField(extField, res);
}

template <class T> class BCHCode : public CyclicCode<T> {
  Poly<T> m_poly;
  Poly<T> m_w;
  int n, j0, t;

  const Field<T> *m_baseField;
  GF_q<T> *m_extField;

  BCHCode(const Field<T> *baseField, int n, int j0, int t, GF_q<T> *extField,
          const Poly<T> &w, const Poly<T> &poly)
      : CyclicCode<T>(baseField, poly, n), m_baseField(baseField), n(n), j0(j0),
        t(t), m_poly(poly), m_w(w), m_extField(extField) {

    // extField->getNthRoot(n).disp();

    // Poly<Poly<int> > extPoly=utils::toExtField(*extField, poly);
    // std::vector<Poly<int> > itb=fft.doifft(tb);
  }

public:
  ~BCHCode() { delete m_extField; }

  static BCHCode *getNew(const Field<T> *baseField, int n, int j0, int t) {
    GF_q<T> *extField;
    Poly<T> w;
    Poly<T> poly = findBCHPoly(baseField, n, j0, t, extField, w);
    return new BCHCode(baseField, n, j0, t, extField, w, poly);
  }

  virtual std::vector<T> decode(const std::vector<T> &msg) {
    typedef T BASE_TYPE;
    typedef Poly<T> EXT_TYPE;

    std::vector<BASE_TYPE> res;
    PolyRing<BASE_TYPE> pr(m_baseField);
    PolyRing<EXT_TYPE> extPr(m_extField);
    Poly<EXT_TYPE> m2 = toExtField(m_extField, pr.import(msg));

    std::vector<EXT_TYPE> tb(2 * t + 1);
    EXT_TYPE ww = m_extField->faste(m_w, j0);
    for (int i = j0; i < j0 + 2 * t + 1; ++i) {
      tb[i - j0] = extPr.eval(m2, ww);
      ww = m_extField->mul(ww, m_w);
    }

    out(tb);

    FFT<EXT_TYPE> fft;
    fft.init(m_extField, n, m_w);

    std::vector<T> correctMsg;
    std::vector<EXT_TYPE> x;

    if (0)
      x = findMinLinearRecursion_Slow(*m_extField, tb, t);
    else {
      x = findMinLinearRecursion_Massey(*m_extField, tb, t);
      printf("RETURN >> %d\n", x.size());
    }
    printf("GOT >> %d %d\n", x.size(), t);
    out(x);

    // std::reverse(x.begin(), x.end());
    out(x);

    int sz = x.size() - 1;
    std::vector<int> errorPos;
    Poly<EXT_TYPE> locatorPoly = extPr.import(x);
    Poly<BASE_TYPE> msgPoly;
    Poly<BASE_TYPE> errPoly;

    if (sz <= 0) {
      correctMsg = msg;
      goto out;
    }

    ww = m_extField->faste(m_w, j0);

    for (int j = 0; j < n; ++j) {
      if (extPr.eval(locatorPoly, ww) == m_extField->getZ())
        errorPos.push_back(n - 1 - j);
      ww = m_extField->mul(ww, m_w);
    }
    printf("ERROR POS >>> \n");
    out(errorPos);
    printf(">>> GOT %d for %d\n", sz, errorPos.size());
    if (sz != errorPos.size())
      goto fail;

    if (m_baseField->getSize() == 2) {
      for (int j = 0; j < errorPos.size(); ++j)
        pr.set1(errPoly, errorPos[j], m_baseField->getE());
    } else {
      Matrix<EXT_TYPE> errMatrix(m_extField, sz, sz);
      for (int j = 0; j < sz; ++j)
        errMatrix(0, j) = m_extField->faste(m_w, errorPos[j]);

      for (int j = 1; j < sz; ++j)
        for (int k = 0; k < sz; ++k)
          errMatrix.get(j, k) =
            m_extField->mul(errMatrix.get(j - 1, k), errMatrix.get(0, k));

      tb.resize(sz);
      x = errMatrix.solve(tb);
      assert(x.size() != 0);
      out(x);
      for (int j = 0; j < sz; ++j)
        pr.set1(errPoly, errorPos[j], toBaseField1(m_extField, x[j]));
    }

    msgPoly = pr.import(msg);
    msgPoly = pr.add(msgPoly, pr.neg(errPoly));
    correctMsg = pr.toVector(msgPoly, n);
    goto out;

    correctMsg = msg;
  out:;
    return CyclicCode<T>::decode(correctMsg);
  fail:
    return std::vector<T>();
  }
};
typedef BCHCode<u32> BCHCode_u32;
typedef BCHCode<Poly<u32>> BCHCode_P_u32;

OPA_NAMESPACE_DECL3_END

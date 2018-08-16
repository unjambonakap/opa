#include <stdio.h>

#include <opa/math/common/GF_q.h>
#include <opa/math/common/GF_p.h>
#include <opa/math/common/PolyRing.h>
#include <opa/math/common/Utils.h>
#include <opa/math/common/Matrix.h>
#include <opa/math/common/CyclicCode.h>
#include <opa/math/common/BCHCode.h>
#include <opa/math/common/RSCode.h>
#include <opa/math/common/UtilsPoly.h>
#include <opa/math/common/UtilsGFq.h>
#include <opa/math/common/BCHCode.h>
#include <opa/math/common/FFT.h>
#include <opa/math/common/algo.h>

using namespace opa::math::common;
using namespace std;

void testMatrix() {
  GF_p a(7);
  PolyRing<u32> pr(a);

  std::vector<u32> tb = pr.rand(4).toVector();
  out(tb);
  rng.seed(14);
  Matrix<u32> m;
  m.affect(Matrix<u32>::rand(&a, 5, 5));
  // m=Matrix::identity(&a,5);
  m.scmul(3);

  Matrix<u32> res;
  bool ok;

  ok = m.invert(&res);

  puts("");
  printf("=============\n");
  printf("INIT MATRIX\n");
  m.disp();

  puts("");
  printf("RES OF INVERT >> %d\n", ok);
  if (ok) {
    printf("RESULT MATRIX >>> \n");
    res.disp();
    puts("++++");

    Matrix<u32> prod;
    prod.affect(m.mul(res));
    printf("PROD WITH INV >>> \n");
    prod.disp();
  }
}

void testCyclicCode() {
  GF_p baseField(3);
  int d = 2;
  GF_q<u32> extensionField(baseField, d);

  Poly<u32> tmp =
    findMinimalPoly(extensionField, extensionField.getPrimitiveElem());
  printf("MINIMAL POLY >>> \n");
  tmp.disp();

  int n = u32_faste(baseField.getSize().getu32(), d) - 1;
  int k = tmp.size();

  CyclicCode<u32> c(baseField, tmp, n);

  std::vector<u32> tb = genRand(baseField.getSizeU32(), c.getNK());
  std::vector<u32> res;

  printf("ENCODE INPUT\n");
  out(tb);
  res = c.encode(tb);
  // res[0]=baseField.add(res[0], baseField.getE());
  printf("ENCODED WORD\n");
  out(res);
  printf(">>> %d\n", c.isCodeword(res));

  res = c.decode(res);
  printf("ENCODE OUTPUT\n");
  out(res);
  out(c.encode(res));
}

template <class T>
std::vector<T> genSeq(const Field<T> &field, std::vector<T> start,
                      const std::vector<T> &gen, int sz) {
  assert(start.size() == gen.size());

  for (int i = 0; i < sz; ++i) {
    T tmp = field.getZ();
    for (int j = 0; j < gen.size(); ++j)
      tmp = field.add(tmp, field.mul(gen[j], start[start.size() - 1 - j]));

    start.push_back(tmp);
  }
  return start;
}

void test() {

  if (0) {

    for (int ntry = 0; ntry < 100; ++ntry) {
      GF_p a(5);
      int n = 50;
      int sz = 300;
      int maxd = (n + sz) / 2;
      std::vector<u32> ta = genRandVector(n, a);

      std::vector<u32> st = genRandVector(n, a);

      std::vector<u32> res = genSeq(a, st, ta, sz);

      ta.insert(ta.begin(), a.neg(a.getE()));

      std::vector<u32> r1, r2;
      assert(checkGenSeq(a, ta, res, res.size()));

      // r1=findMinLinearRecursion_Slow(a, res, maxd);
      r2 = findMinLinearRecursion_Massey(a, res, n + sz);

      puts("ORIG WAS:");
      out(ta);
      // out(r1);
      out(r2);
      puts("ND");
      assert(ta.size() >= r2.size());
    }

    exit(0);
  }

  if (1) {
    GF_p a(11);
    int n = 120;
    int j0 = 1, t = 20;
    BCHCode<u32> *bch = BCHCode<u32>::getNew(a, n, 1, t);
    printf("BCH CODE >>> n=%d, k=%d\n", bch->getN(), bch->getK());

    for (int tr = 0; tr < 100; ++tr) {
      std::vector<u32> tb = genRand(a.getSizeU32(), bch->getNK());
      std::vector<u32> res;

      res = bch->encode(tb);
      PolyRing<u32> pr(a);

      int weight = t;
      std::vector<u32> err = genRandWeight(n, weight, a);
      int cnt = 0;
      for (int i = 0; i < err.size(); ++i)
        if (err[i] != a.getZ())
          ++cnt;
      printf("REAL WEIGHT >> %d (wanted %d\n", cnt, weight);

      std::vector<u32> msg = pr.add(pr.import(res), pr.import(err)).toVector();
      printf("RES >>> ");
      out(res);
      msg.resize(bch->getN(), a.getZ());

      printf("ERR POLY >> ");
      out(err);
      std::vector<u32> sol = bch->decode(msg);
      if (sol.size() == 0) {
        printf("DECODING FALIED\n");
        exit(0);

      } else {
        printf("RESULT OF DECODING>>\n");
        out(sol);
        out(tb);
        assert(sol == tb);
      }
    }

    delete bch;
    exit(0);
  }

  //{
  //    GF_p a(3);
  //    GF_q<int> b(a,3);

  //    testCyclicCode();
  //    exit(0);
  //}
}

void test_codechef_stepnum() {
#define CODECHEF_MOD 4294967143
  const u32 mod = 4294967143;
  const int N = 11;
  GF_p field(mod);
  Matrix<u32> m(&field, N, N);
  REP(i, N) REP(j, N) m(i, j) = (i + 1 == j || i - 1 == j) ? 1 : 0;
  Poly<u32> charPol = gfx_char_poly(m);

  PolyRing<u32> pr(field);
  charPol = pr.import(charPol);
  auto x = pr.factor(charPol);
  charPol.disp();
  puts("\n====FACTORS");

  printf(">> %u\n", pr.eval(charPol, 0));
  printf(">> %u\n", pr.eval(charPol, 1));
  printf(">> %u\n", pr.eval(charPol, field.importu32(mod - 1)));
  printf(">> %u\n", pr.eval(charPol, field.importu32(mod - 3761656586)));
  printf(">> %u\n", pr.eval(charPol, field.importu32(mod - 533310557)));

  REP(i, x.size())
  pr.smonic(x[i]);

  sort(ALL(x));
  puts("################");
  for (auto y : x)
    std::cout << y << std::endl;
  puts("################");

  {
    bignum u = mod;
    u.smul(u);
    u.sadd(bignum::froms32(-1));
    printf("try factor %s\n", u.str(10).c_str());
    BGFactors factors = factor_large(u);
    printf("%s has %lu factors\n", u.str(10).c_str(), factors.size());
    for (auto e : factors) {
      printf(">> %s %d\n", e.ST.str(10).c_str(), e.ND);
    }
  }

  auto poly = x.back();
  GF_q<u32> gfq(field, poly);
  typedef Poly<u32> T2;
  puts("PRIMITIVE ELEM");
  cout << gfq.getPrimitiveElem() << endl;

  PolyRing<T2> r2(gfq);
  vector<T2> vecPoly;
  for (auto a : poly.toVector())
    vecPoly.pb(pr.constant(a));
  puts("EVAL");
  auto polyGFQ = r2.import(vecPoly);
  auto xGFQ = gfq.import(pr.x());
  cout << r2.eval(polyGFQ, xGFQ) << endl;
  auto factor2 = r2.factor(polyGFQ);
  cout << polyGFQ << endl;
  for (auto a : factor2) {
    puts("NEW FACTOR");
    r2.smonic(a);
    cout << a << endl;
    cout << r2.eval(a, xGFQ) << endl;
  }

  puts("\n\nGO SOLVER");
  LinearRecurrenceSolver<u32> solver;
  solver.initialize(m, field);
  solver.setup();

  int pos = 0;
  vector<u32> check_vals;
  while (1) {
    Matrix<u32> m2;
    m2.affect(m.faste(pos));
    u32 cnt = 0;

    REP(i, N) {
      if (i - 1 > 0)
        cnt += m2.get(i, i - 1);
      if (i + 1 < N)
        cnt += m2.get(i, i + 1);
    }
    check_vals.pb(cnt);

    bool ok = solver.add_init_cond(cnt, pos);
    if (ok)
      break;

    ++pos;
  }

  bignum maxv = bignum(1 << 16).pow(4);
  int ntry = 20000;
  REP(i, ntry) {
    bignum test = maxv.rand();
    u32 res = solver.get(test);
    printf(">> got %u\n", res);
  }
  out(check_vals);
}

void testMinimal() {
  const int MOD = 13;

  typedef u32 T;
  typedef Poly<T> TExt;
  GF_p field(MOD);
  PolyRing<T> pr(field);

  Poly<T> irred = findIrred(field, 3);
  GF_q<T> extField(field, irred);
  PolyRing<TExt> extPr(extField);

  Poly<TExt> irredEmbed = toExtField(extField, irred);
  Matrix<T> mat1(&field, 3, 3);
  std::vector<T> eVec;
  {
    auto factors = extPr.factor(irredEmbed);
    REP(i, factors.size()) {
      auto &x = factors[i];
      assert(x.deg() == 1);
      extPr.smonic(x);
      mat1.setCol(i, pr.toVector(extField.neg(x.get(0)), 3));
    }
    std::vector<T> res = mat1.solve(pr.toVector(pr.E(), 3));
    std::cout << mat1.get_det() << std::endl;
    eVec = res;
    out(eVec);
  }

  std::vector<vector<T> > mpIrred2;
  Poly<T> irred2 = findIrred(field, 3);
  {
    Poly<TExt> embed = toExtField(extField, irred2);
    auto factors = extPr.factor(embed);
    for (auto x : factors) {
      extPr.smonic(x);
      auto y = mat1.solve(pr.toVector(extField.neg(x.get(0)), 3));
      out(y);
      mpIrred2.pb(y);
    }
  }

  Poly<T> irred4 = findIrred(field, 4);

  Poly<T> irred12 = findIrred(field, 12);

  {
    GF_q<T> bigField(field, irred12);
    PolyRing<TExt> bigPr(bigField);
    Poly<TExt> a = toExtField(bigField, irred);
    auto factors = bigPr.factor(a);

    Matrix<T> mat(&field, 12, 3);
    REP(i, factors.size()) {
      auto &x = factors[i % factors.size()];
      std::cout << x << std::endl;
      assert(x.deg() == 1);
      bigPr.smonic(x);
      mat.setCol(i, pr.toVector(bigField.neg(x.get(0)), 12));
    }

    std::vector<T> res = mat.eval(eVec);
    out(res);

    puts("FOR IRRED2");
    factors = bigPr.factor(toExtField(bigField, irred2));
    GF_q<T> gf2(field, irred2);
    Poly<TExt> poly2 = toExtField(gf2, irred2);
    PolyRing<TExt> pr2(gf2);

    for (auto x : factors) {
      bigPr.smonic(x);
      TExt p = bigField.neg(x.get(0));
      cout << bigPr.eval(toExtField(bigField, irred2), p) << endl;
      cout << p << endl;
      TExt q = pr.mod(p, irred2);
      cout << q << endl;
      cout << pr2.eval(poly2, pr.x()) << endl;
      puts("");
    }

    puts("eval");
    for (auto x : mpIrred2) {
      std::vector<T> res = mat.eval(x);
      out(res);
    }
  }
}

int main() {
  initMathCommon();

  int mode = 0;
  if (mode == 0) {
    test_codechef_stepnum();
  } else if (mode == 1) {
    testMinimal();
  } else {

    test();

    GF_p a(7);
    int b = a.mul(5, 6);

    PolyRing<u32> pr(a);
    Poly<u32> pmod = findIrred(a, 4);
    GF_q<u32> gfq(a, pmod);
    Poly<u32> w = gfq.getPrimitiveElem();

    printf("IRRED: ");
    pmod.disp();
    printf("PRIMITIVE: ");
    w.disp();

    // testMatrix();
    testCyclicCode();
  }

  return 0;
}

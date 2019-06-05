#include <gtest/gtest.h>

#include <opa/math/common/BCHCode.h>
#include <opa/math/common/CyclicCode.h>
#include <opa/math/common/FFT.h>
#include <opa/math/common/FractionField.h>
#include <opa/math/common/GF_p.h>
#include <opa/math/common/GF_pBN.h>
#include <opa/math/common/GF_q.h>
#include <opa/math/common/IntegerRingUtil.h>
#include <opa/math/common/Matrix.h>
#include <opa/math/common/PolyRing.h>
#include <opa/math/common/RSCode.h>
#include <opa/math/common/Utils.h>
#include <opa/math/common/UtilsGFq.h>
#include <opa/math/common/UtilsPoly.h>
#include <opa/math/common/UtilsRing.h>
#include <opa/math/common/Z.h>
#include <opa/math/common/Zn_BG.h>
#include <opa/math/common/algo.h>
#include <opa/math/common/float.h>
#include <opa/math/common/matrix_utils.h>
#include <opa/math/common/number_field.h>
#include <opa/math/common/utils_num.h>
#include <opa/or/grid_search.h>
#include <opa/utils/misc.h>

using namespace opa;

using namespace opa::math::common;
using namespace std;
typedef NumberField::P_K P_K;
typedef NumberField::P_KK P_KK;
typedef NumberField::K K;

DEFINE_int32(opa_roots_prec, -500, "");
DEFINE_int32(opa_roots_steps, 2000, "");

template <class T> bool has_linear_term(const std::vector<T> &tb) {
  for (auto &e : tb)
    if (e.deg() == 1) return true;
  return false;
}

template <class T>
std::vector<T> get_linear_roots(const std::vector<Poly<T> > &tb) {
  std::vector<T> res;
  for (auto &e : tb)
    if (e.deg() == 1) res.push_back(e.linear_root());
  return res;
}

void test_init() { initMathCommon(0); }

template <class T> T extract(const Matrix<T> &m, const Matrix<T> &extractor) {
  const Ring<T> *ring = m.get_ring();
  T res = ring->getZ();
  REP (i, m.getN())
    REP (j, m.getM())
      res = ring->add(res, ring->mul(m.get(i, j), extractor.get(i, j)));

  return res;
}

template <class T>
T get_simple(const Matrix<T> &m, const Matrix<T> &extractor, const bignum &p) {
  Matrix<T> m2 = m.faste(p);
  T val = extract(m2, extractor);
  return val;
}

template <class T> void doTest(int n, const Field<T> &field) {
  int ntry = 10;

  REP (tt, ntry) {
    Matrix<T> m = Matrix<T>::rand(&field, n, n);
    Matrix<T> extractor = Matrix<T>::rand(&field, n, n);

    LinearRecurrenceSolver<T> solver;
    solver.initialize(m, field);
    solver.setup();

    int pos = 0;
    while (1) {

      T res = get_simple(m, extractor, pos);
      bool ok = solver.add_init_cond(res, pos);
      if (ok) break;

      ++pos;
    }

    REP (ii, ntry) {
      u32 pw = rng();
      T correct = get_simple(m, extractor, pw);
      T fast = solver.get(pw);
      ASSERT_EQ(fast, correct);
    }
  }
}

TEST(LinearRecurrenceTest, RandTest_Small) {
  test_init();

  int min_size = 2;
  int max_size = 7;
  int ntest = 10;

  REP (tt, ntest) {
    int n = rng() % (max_size - min_size + 1) + min_size;
    u32 mod = gen_prime(bignum(1e3)).getu32();
    printf(">> ON FIELD >> %u\n", mod);
    GF_p field(mod);
    ASSERT_NO_FATAL_FAILURE(doTest(n, field));
  }
}

TEST(PolyTest1, PolyTest1) {
  typedef FractionField<P_Q> QPZ;
  typedef PolyRing<QPZ::Type> PZ2;

  QPZ qpz1(&Q_x);
  PZ2 pr2(&qpz1);

  auto f1 = pr2.get_poly();
  auto f2 = pr2.get_poly();

  {
    Q x3[] = { QF(0), QF(1) };
    Q x2[] = { QF(0), QF(1), QF(2) };
    Q x1[] = { QF(1), QF(0), QF(1), QF(1) };
    Q x0[] = { QF(1), QF(1) };
    pr2.set1(f1, 0, qpz1.import(Q_x.import(LIST2VEC(x0))));
    pr2.set1(f1, 1, qpz1.import(Q_x.import(LIST2VEC(x1))));
    pr2.set1(f1, 2, qpz1.import(Q_x.import(LIST2VEC(x2))));
    pr2.set1(f1, 3, qpz1.import(Q_x.import(LIST2VEC(x3))));
  }

  {
    Q x4[] = { QF(0), QF(1) };
    Q x3[] = { QF(0), QF(0), QF(2) };
    Q x2[] = { QF(1), QF(0), QF(0), QF(1) };
    Q x1[] = { QF(0), QF(1) };
    Q x0[] = { QF(0) };
    pr2.set1(f2, 0, qpz1.import(Q_x.import(LIST2VEC(x0))));
    pr2.set1(f2, 1, qpz1.import(Q_x.import(LIST2VEC(x1))));
    pr2.set1(f2, 2, qpz1.import(Q_x.import(LIST2VEC(x2))));
    pr2.set1(f2, 3, qpz1.import(Q_x.import(LIST2VEC(x3))));
    pr2.set1(f2, 4, qpz1.import(Q_x.import(LIST2VEC(x4))));
  }
  printf("pr2 >> %d\n", pr2.get_poly_pos());
  cout << f1 << endl;
  cout << f2 << endl;

  auto x = pr2.gcd(f1, f2);
  cout << pr2.div(f1, x) << endl;
}

TEST(PolyTest2, PolyTest) {
  typedef Zn_BG R;
  typedef PolyRing<R::Type> PR;
  typedef PR::Type Poly;
  //    ('N',
  // 740765548979273098467598803958212385151570053921334237430171491357308450305938925395058048571558613364002948004291135518240329572789525487495147870779619379982865011328775565850048248526863374376024296921937798169737860584047065593928295857417452372744936947544816804233701992919611488140593397159150152160920639)
  //        \n('E', 3)
  //        \n('c1',
  // 321451913057900142348436621563079153898436032837412854031246697790410602040147332179869901737501439750726664592096795391349241878705910327861241059454661619432324907631836944052325945666694131891489395959762339277410381692975197784562565409741727780043108506458930540459411390110527535022558117942647833465024816)
  //        \n('c2',
  // 245544492996888727164841815357590445824184017819212225646984254796592976347385430003123536033004742032759416589790020081404988144261423483934413815011827391460021382138980857256740580080876659438050231270325521568176877577234140954433459361598359898768186238850013895811016147701584167992547871319386741583536303)
  //        \n('sol', None)\n"
  OPA_BG n = OPA_BG::fromstr(
    "740765548979273098467598803958212385151570053921334237430171491357308"
    "4503059389253950580485715586133640029480042911355182403295727895254874"
    "9514787077961937998286501132877556585004824852686337437602429692193779"
    "8169737860584047065593928295857417452372744936947544816804233701992919"
    "611488140593397159150152160920639",
    10);

  OPA_BG e = 3;

  OPA_BG c1 = OPA_BG::fromstr(
    "3214519130579001423484366215630791538984360328374128540312466977904106"
    "0204014733217986990173750143975072666459209679539134924187870591032786"
    "1241059454661619432324907631836944052325945666694131891489395959762339"
    "2774103816929751977845625654097417277800431085064589305404594113901105"
    "27535022558117942647833465024816",
    10);

  OPA_BG c2 = OPA_BG::fromstr(
    "2455444929968887271648418153575904458241840178192122256469842547965929"
    "7634738543000312353603300474203275941658979002008140498814426142348393"
    "4413815011827391460021382138980857256740580080876659438050231270325521"
    "5681768775772341409544334593615983598987681862388500138958110161477015"
    "84167992547871319386741583536303",
    10);
  OPA_BG sol = OPA_BG::fromu64(2449967717ll);

  OPA_BG expected_m = OPA_BG::fromstr(
    "2782365860388557286310809887683786504528738655865781703464278179636746"
    "2742441257825190950207363331440302148015785208096779553247828949275151"
    "5506113530367705961747287052814906081099788439953915040340986457201263"
    "045837650841031752009184680515477175397898656180618898643429",
    10);

  R r(n);
  PR pr(&r);

  Poly p1 = pr.sub(pr.faste(pr.x(), e), pr.constant(c1));
  cout << pr.constant(sol) << endl;
  Poly p2 =
    pr.sub(pr.faste(pr.add(pr.x(), pr.constant(sol)), e), pr.constant(c2));
  cout << pr.eval(p1, expected_m) << endl;
  cout << pr.eval(p2, expected_m) << endl;

  auto res = pr.gcd(p1, p2);
  pr.smonic(res);
  cout << (r.neg(res.get(0)) >> 32) << endl;
  cout << expected_m << endl;
}

TEST(FloatTest, Test1) {
  puts("HERE");
  Float a;
  Float b("-0.12");
  cout << a.add(b).to_str() << endl;
}

template <class T> void test_tonelli_field(const Field<T> &field) {
  T x;
  while (true) {
    x = field.getRandRaw();
    if (has_square_root(field, x)) break;
  }
  T sq;
  ASSERT_TRUE(find_square_root(field, x, &sq));
  OPA_DISP0(x, sq, field.mul(sq, sq));
  ASSERT_EQ(field.mul(sq, sq), x);
}

void test_tonelli(const bignum &p) {
  GF_pBN gfp(p);
  test_tonelli_field(gfp);
}

TEST(TonelliShanks, Tiny) {
  GF_p gfp(7);
  u32 sq;
  u32 x = gfp.mul(3, 3);
  ASSERT_TRUE(find_square_root(gfp, x, &sq));
  ASSERT_EQ(gfp.mul(sq, sq), x);
}

TEST(TonelliShanks, Test1) {
  bignum p = gen_prime(1 << 30);
  REP (i, 10) { test_tonelli(p); }
}

TEST(TonelliShanks, TestLarge) {
  REP (prime_test, 10) {
    bignum p = gen_prime(bignum(2).lshift(500));
    REP (i, 10) { test_tonelli(p); }
  }
}

TEST(TonelliShanks, TestGfq) {
  std::uniform_int_distribution<> deg_rng(3, 8);
  REP (i, 10) {
    u32 p = 2;
    while (p == 2) p = gen_prime(1 << 7).getu32();
    int deg = deg_rng(rng);
    GF_p gfp(p);
    GF_q<u32> gfq(&gfp, deg);
    OPA_DISP0(p, deg, gfq.getSize());
    test_tonelli_field(gfq);
  }
}

TEST(ZPolyFact, Hensel_LiftOne) {
  typedef Poly<bignum> P_t;
  P_t a, b, c;

  GF_pBN p_big(nextPrimeSmall(1e6 + 3));
  a = findIrred(p_big, 2);
  b = findIrred(p_big, 3);

  a.unsafe_change_ring(&PR_Z);
  b.unsafe_change_ring(&PR_Z);
  c = a * b;
  OPA_DISP0("Try factor ", c, a, b);

  bignum p = 13;

  HenselLifting lift;
  lift.params.p = p.getu32();
  lift.set_cur_ring(1);

  {
    P_t pa, pb, u, v;
    pa = lift.pr_zn.import(a);
    pb = lift.pr_zn.import(b);
    lift.set_gcd_coeffs(pa, pb, &u, &v);
    P_t tmp = pa * u + pb * v;
    OPA_CHECK(lift.pr_zn.isE(tmp), tmp, pa, pb, u, v);

    lift.set_cur_ring(p);
    P_t pc = lift.pr_zn.import(c);
    lift.lift_one(pc, pa, pb, u, v);
    OPA_CHECK(pc == pa * pb, pc, pa * pb, pa, pb);
    // OPA_CHECK(lift.pr_zn.isE(pa * u + pb * v), pa * u + pb * v, pa, u, pb,
    // v);
  }

  {
    std::vector<Poly<bignum> > factors;
    factor_zpoly_squarefree(c, &factors);
    OPA_DISP0(factors);
    sort(ALL(factors));
    std::vector<Poly<bignum> > expected = { a, b };
    sort(ALL(expected));
    ASSERT_EQ(factors, expected);
  }
}

TEST(PolyFFT, TestFromReal) {
  GF_p gfp(107);
  PolyRing<u32> pr(&gfp);
  typedef Poly<u32> PT;
  PT a = pr.import({ 57, 1, 2, 105 });
  PT b = pr.import({ 84, 33, 8 });

  PT c1 = pr.mul(a, b);
  PolyFFTHelper_WithReal<u32, double> fft_helper;
  fft_helper.init(&gfp, a.deg() + b.deg() + 1);
  std::vector<u32> t_c2 = fft_helper.multiply(a.to_vec(), b.to_vec());
  PT c2 = pr.import(t_c2);

  OPA_DISP0(c1, c2, a, b, c1 == c2);
}

TEST(PolyFFT, TestReal) {
  FFT_Real<double> fft;
  fft.init(10);
  std::vector<double> a = { 1., 3., 5., 7. };
  auto tmp = fft.fft(a);
  auto res = fft.ifft(tmp);
  OPA_DISP0(a, res);
}

TEST(PolyFFT, Test2) {
  GF_p gfp(13);
  PolyRing<u32> pr(&gfp);
  typedef Poly<u32> PT;
  PT a = pr.import({ 1, 2, 4 });
  PT b = pr.import({ 2, 3 });

  int wanted_order = a.deg() + b.deg() + 1;
  wanted_order = 12;
  PolyFFTPlanner<u32> planner;
  planner.init(&gfp, wanted_order);
  OPA_DISP0(planner.fft(a.to_vec()));

  PolyFFTPlanner<u32> planner2;
  planner2.init(&gfp, wanted_order, false, false);
  OPA_DISP0(planner2.fft(a.to_vec()));

  std::vector<u32> t_c2 =
    planner.main_executor->multiply(a.to_vec(), b.to_vec());

  FFT<u32> fft2;
  fft2.init(&gfp, wanted_order, gfp.getNthRoot(wanted_order));
  auto a1 = fft2.fft(a.to_vec());
  auto a2 = fft2.fft(b.to_vec());
  OPA_DISP0(a1, a2);
  REP (i, a1.size())
    a1[i] = gfp.mul(a1[i], a2[i]);
  OPA_DISP0(a1);
  OPA_DISP0(fft2.ifft(a1));

  PT c1 = pr.mul(a, b);
  PT c2 = pr.import(t_c2);

  OPA_DISP0(c1, c2, a, b, c1 == c2);
}

TEST(PolyFFT, TestMulMultiLarge) {
  int prime = 1e9 + 2;
  REP (i, 1) {
    prime = next_prime(prime + 1).getu32();
    GF_p gfp(prime);
    PolyRing<u32> pr(&gfp);
    REP (j, 2) {
      int pa = rng() % int(1e5) + 1;
      int pb = rng() % int(1e5) + 1;
      OPA_DISP0(i, j, prime, pa, pb);
      int wanted_order = pa + pb + 1;

      PolyFFTPlanner<u32> planner;
      planner.init(&gfp, wanted_order, false, false);

      REP (k, 5) {
        typedef Poly<u32> PT;
        PT a = pr.rand(pa);
        PT b = pr.rand(pb);

        std::vector<u32> t_c2 =
          planner.main_executor->multiply(a.to_vec(), b.to_vec());

        if (0) {
          PT c1 = pr.mul(a, b);
          PT c2 = pr.import(t_c2);

          OPA_CHECK(c1 == c2, c1, c2, a, b, c1 == c2);
        }
      }
    }
  }
}

TEST(PolyFFT, TestMulMulti) {
  int prime = 1e6 + 2;
  REP (i, 5) {
    prime = nextPrimeSmall(prime + 1);
    GF_p gfp(prime);
    PolyRing<u32> pr(&gfp);
    REP (j, 5) {
      OPA_DISP0(i, j);
      int pa = rng() % int(1e3) + 1;
      int pb = rng() % int(1e3) + 1;
      int wanted_order = pa + pb + 1;

      PolyFFTPlanner<u32> planner;
      planner.init(&gfp, wanted_order, false, false);

      REP (k, 10) {
        typedef Poly<u32> PT;
        PT a = pr.rand(pa);
        PT b = pr.rand(pb);

        std::vector<u32> t_c2 =
          planner.main_executor->multiply(a.to_vec(), b.to_vec());

        PT c1 = pr.mul(a, b);
        PT c2 = pr.import(t_c2);

        OPA_CHECK(c1 == c2, c1, c2, a, b, c1 == c2);
      }
    }
  }
}

TEST(PolyFFT, TestMulSimple) {
  int p = 0x2d000001;
  int order = 30;
  int root = 0x1e8fa7f5;

  if (0) {
    p = 257;
    order = p - 1;
    root = -1;
  }

  GF_p gfp(p);
  if (root == -1) root = gfp.getPrimitiveElem();
  OPA_DISP0(gfp.getOrderOf(root));
  OPA_DISP0(gfp.faste(root, 5));

  PolyRing<u32> pr(&gfp);
  typedef Poly<u32> PT;
  PT a = pr.import({ 1, 2, 4 });
  PT b = pr.import({ 2, 3 });

  int wanted_order = order;
  PolyFFTPlanner<u32> planner;
  planner.init(&gfp, wanted_order, !PolyFFTPlanner<u32>::kForceSimple,
               PolyFFTPlanner<u32>::kWantDft);
  PolyFFTPlanner<u32> simple_planner;
  simple_planner.init(&gfp, wanted_order, PolyFFTPlanner<u32>::kForceSimple,
                      PolyFFTPlanner<u32>::kWantDft);

  OPA_DISP0(planner.fft(a.to_vec()));
  OPA_DISP0(simple_planner.fft(a.to_vec()));
  OPA_DISP0(planner.ifft(planner.fft(a.to_vec())));
  OPA_DISP0(simple_planner.ifft(simple_planner.fft(a.to_vec())));
  OPA_DISP0(planner.ifft(planner.fft(a.to_vec())));
  std::vector<u32> t_c2 =
    planner.main_executor->multiply(a.to_vec(), b.to_vec());

  PT c1 = pr.mul(a, b);
  PT c2 = pr.import(t_c2);

  OPA_DISP0(c1, c2, a, b, c1 == c2);
}

TEST(PolyFFT, TestMulSplit) {
  GF_p gfp(107);
  PolyRing<u32> pr(&gfp);
  typedef Poly<u32> PT;
  PT a = pr.import({ 1, 2, 4 });
  PT b = pr.import({ 2, 3 });

  int wanted_order = a.deg() + b.deg() + 1;
  wanted_order = gfp.getSizeU32() - 1;
  PolyFFTPlanner<u32> planner;
  planner.init(&gfp, wanted_order, !PolyFFTPlanner<u32>::kForceSimple,
               PolyFFTPlanner<u32>::kWantDft);

  std::vector<u32> t_c2 =
    planner.main_executor->multiply(a.to_vec(), b.to_vec());

  PT c1 = pr.mul(a, b);
  PT c2 = pr.import(t_c2);

  OPA_DISP0(c1, c2, a, b, c1 == c2);
}

TEST(PolyFFT, TestMul) {
  GF_p gfp(1e6 + 3);
  PolyRing<u32> pr(&gfp);
  typedef Poly<u32> PT;
  PT a = pr.import({ 1, 2, 4 });
  PT b = pr.import({ 2, 3 });

  int wanted_order = a.deg() + b.deg() + 1;
  wanted_order = 30;
  PolyFFTPlanner<u32> planner;
  planner.init(&gfp, wanted_order, false, false);

  std::vector<u32> t_c2 =
    planner.main_executor->multiply(a.to_vec(), b.to_vec());

  PT c1 = pr.mul(a, b);
  PT c2 = pr.import(t_c2);

  OPA_DISP0(c1, c2, a, b, c1 == c2);
}

TEST(PolyFFT, Test1) {
  PolyFFTPlanner<u32> planner;
  GF_p gfp(107);
  PolyRing<u32> pr(&gfp);
  typedef Poly<u32> PT;
  PT a = pr.import({ 57, 1, 2, 105 });
  PT b = pr.import({ 84, 33, 8 });

  int wanted_order = a.deg() + b.deg() + 1;
  wanted_order = 106;
  planner.init(&gfp, wanted_order);
  std::vector<u32> t_c2 =
    planner.main_executor->multiply(a.to_vec(), b.to_vec());

  FFT<u32> fft2;
  fft2.init(&gfp, wanted_order, gfp.getNthRoot(wanted_order));
  auto a1 = fft2.fft(a.to_vec());
  auto a2 = fft2.fft(b.to_vec());
  OPA_DISP0(a1, a2);
  REP (i, a1.size())
    a1[i] = gfp.mul(a1[i], a2[i]);
  OPA_DISP0(a1);
  OPA_DISP0(fft2.ifft(a1));

  PT c1 = pr.mul(a, b);
  PT c2 = pr.import(t_c2);

  OPA_DISP0(c1, c2, a, b, c1 == c2);
}

TEST(SmoothPrimeOrder, Test1) {
  puts("LA");
  OPA_DISP0(spoh_for_gfp_split.data().size());
  for (auto &e : spoh_for_gfp_split.data()) {
    OPA_DISP0(e);
  }
}

TEST(ZX, Test1) {
  auto poly = PR_Z.import({ 5, 2 });
  OPA_DISP0(PR_Z.gcd(poly, poly));
}

TEST(ZX, GCD1) {
  zrand_bound = 20;
  FOR (maxdeg, 1, 7) {
    REP (nstep, 100) {
      int da = 1 + rng() % maxdeg;
      int db = 1 + rng() % maxdeg;
      int dc = 1 + rng() % maxdeg;
      P_Z a = PR_Z.rand(da);
      P_Z b = PR_Z.rand(db);
      P_Z c = PR_Z.rand(dc);
      P_Z p = a * b;
      P_Z p2 = a * c;
      auto ctx = PR_Z.egcd2_ctx(p, p2);
      OPA_DISP0(ctx.res, ctx.gcd_mul, ctx.coeff_mul, p, p2, a);
      P_Z d = ctx.res;
      OPA_CHECK0(PR_Z.mod(p, d) == PR_Z.getZ());
      OPA_CHECK0(PR_Z.mod(p2, d) == PR_Z.getZ());
    }
  }
}

TEST(ZX, MignotteBound) {
  bignum bound = 1e5;
  REP (nstep, 100) {
    int da = 1 + rng() % 10;
    int db = 1 + rng() % 10;
    P_Z a = PR_Z.rand(da);
    P_Z b = PR_Z.rand(db);
    P_Z p = a * b;
    bignum bound = MignotteBoundHelper::ComputeBound(p);
    bignum diff = bignum::fromFloat(1e100);
    for (auto &x : a) diff = std::min(diff, bound - x.abs());
    for (auto &x : b) diff = std::min(diff, bound - x.abs());
    OPA_DISP0(bound, a, b, diff);
    OPA_CHECK(diff >= 0, diff, a, b, bound);
  }
}

TEST(ZX, Factor1) {
  zrand_bound = 20;
  FOR (nfacts, 1, 10) {
    REP (nstep, 10) {
      puts("\n\n");
      std::vector<P_Z> factors;
      REP (i, nfacts) {
        int factor_max_deg = 4;
        int deg = rng() % factor_max_deg + 1;
        auto tmp = PR_Z.rand(deg).pp();
        if (tmp.lc() < 0) tmp = -tmp;
        factors.push_back(tmp);
      }

      P_Z poly = PR_Z.getE();
      for (auto &fact : factors) poly = poly * fact;
      std::vector<std::pair<P_Z, int> > found_factors;
      OPA_DISP("Trying ", poly, factors);
      bignum unit = factor_zpoly(poly, &found_factors);
      OPA_CHECK(unit == 1, unit); // took primitive part + positive
      auto factors_with_count = utils::list_to_count(factors);

      P_Z res_poly = PR_Z.getE() * unit;
      for (auto &factor_and_pw : found_factors) {
        REP (i, factor_and_pw.second)
          res_poly = res_poly * factor_and_pw.first;
      }
      OPA_CHECK(res_poly == poly, res_poly, poly, found_factors,
                factors_with_count);

      // Cannot check like this, some factors maybe not be irreducible
      // OPA_DISP0(found_factors, factors_with_count);
      // OPA_CHECK_EQ0(found_factors, factors_with_count);
    }
  }
}

TEST(GcdPlan, ZX1) {
  auto poly = PR_Z.import({ 5, 2 });
  auto plan = get_gcd_plan(PR_Z, std::vector<P_Z>{ poly, PR_Z.getZ() });
  OPA_DISP0(plan.gcd);
  OPA_DISP0(plan.pos);
}

TEST(GcdPlan, QX1) {
  P_Q a = Q_x.xpw(2) * QF(-1, 2) + Q_x.x() * QF(5) + Q_x.constant(QF(-19, 2));
  P_Q b = Q_x.x() - Q_x.constant(QF(5));
  //(((-19 / 2)*A^0 5*A^1 (1 / -2)*A^2),(-5*A^0 1*A^1)),
  reduce_gcd2(Q_x, a, b);
  OPA_DISP0(a, b);
}

TEST(PolyBasic, SwitchVars) {
  PolyRing<P_Z> pr_zz(&PR_Z);
  typedef Poly<P_Z> P_ZZ;
  P_ZZ pxy = pr_zz.import_vec(
    { PR_Z.import_vec({ 1, 4, 5 }), PR_Z.import_vec({ 2, 8 }) });
  P_ZZ pyx =
    pr_zz.import_vec({ PR_Z.import_vec({ 1, 2 }), PR_Z.import_vec({ 4, 8 }),
                       PR_Z.import_vec({ 5 }) });
  OPA_CHECK_EQ0(pyx, poly_switch_vars(pxy));
}

TEST(PolyBasic, Resultant) {
  P_QQ p1 =
    Q_xy.import({ Q_x.getE(), Q_x.getZ(), Q_x.constant(QF(-2)) }, kPolyRev);
  P_QQ p2 =
    Q_xy.import({ Q_x.constant(QF(-1, 2)), Q_x.x() - Q_x.constant(QF(5)),
                  Q_x.xpw(2) * QF(-1, 2) + Q_x.x() * QF(5) - QF(0x17, 2) },
                kPolyRev);
  OPA_DISP0(resultant2(p1, p2));
  //  (1*A^0) () (-2*A^0) ()
  //    () (1*A^0) () (-2*A^0)
  //    ((1 / -2)*A^0) (-5*A^0 1*A^1) ((-17 / 2)*A^0 5*A^1 (1 / -2)*A^2) ()
  //    () ((1 / -2)*A^0) (-5*A^0 1*A^1) ((-17 / 2)*A^0 5*A^1 (1 / -2)*A^2)
  //    ,
  //
}

TEST(PolyBasic, Resultant2) {
  zrand_bound = 10;
  REP (nstep, 100) {
    int d1 = rng() % 5 + 1;
    int d2 = rng() % 5 + 1;
    P_Q p1 = Q_x.rand(d1);
    P_Q p2 = Q_x.rand(d2);
    OPA_DISP0(nstep, p1, p2);
    Q r1 = resultant2(p1, p2);
    Q r2 = resultant_slow(p1, p2);
    OPA_CHECK_EQ(r1, r2, p1, p2);
  }
}

TEST(PolyBasic, ResultantQQ) {
  zrand_bound = 10;
  Q_x.poly_rand_deg = 2;
  REP (nstep, 100) {
    int d1 = rng() % 3 + 1;
    int d2 = rng() % 3 + 1;
    P_QQ p1 = Q_xy.rand(d1);
    P_QQ p2 = Q_xy.rand(d2);
    OPA_DISP0(nstep, p1, p2);
    P_Q r1 = resultant2(p1, p2);
    P_Q r2 = resultant_slow(p1, p2);
    REP (i, r1.size())
      OPA_CHECK_EQ(r1[i], r2[i], i);
    OPA_CHECK_EQ(r1, r2, p1, p2);
  }
}

TEST(Base, PermutationSignature) {
  OPA_CHECK_EQ0(1, get_permutation_signature({ 0, 1 }));
  OPA_CHECK_EQ0(-1, get_permutation_signature({ 1, 0 }));
}

TEST(MatrixBasic, DeterminantSimple) {
  Matrix<bignum> m1 = Matrix<bignum>::rand(&Ring_Z, 2, 2);
  Matrix<Q> m1_q = m1.lift(QF);
  bignum r3 = QF.to_base_or_fail(m1_q.get_det());
  bignum r1 = m1(0, 0) * m1(1, 1) - m1(1, 0) * m1(0, 1);
  bignum r2 = m1.get_det_slow();
  OPA_CHECK_EQ(r1, r2, m1);
  OPA_CHECK_EQ(r3, r2, m1);
}

TEST(MatrixBasic, DeterminantSimple3x3) {
  REP (step, 100) {
    const Matrix<bignum> m1 = Matrix<bignum>::rand(&Ring_Z, 3, 3);
    Matrix<Q> m1_q = m1.lift(QF);
    bignum r3 = QF.to_base_or_fail(m1_q.get_det());
    bignum r1 = m1(0, 0) * (m1(1, 1) * m1(2, 2) - m1(1, 2) * m1(2, 1)) -
                m1(1, 0) * (m1(0, 1) * m1(2, 2) - m1(0, 2) * m1(2, 1)) +
                m1(2, 0) * (m1(0, 1) * m1(1, 2) - m1(0, 2) * m1(1, 1));
    bignum r2 = m1.get_det_slow();
    OPA_CHECK_EQ(r1, r2, m1, r1, r2, r3);
    OPA_CHECK_EQ(r3, r2, m1, r1, r2, r3);
  }
}

TEST(MatrixBasic, Determinant) {
  zrand_bound = 10;
  REP (nstep, 100) {
    OPA_DISP0(nstep);
    int n = rng() % 8 + 1;
    Matrix<bignum> m1 = Matrix<bignum>::rand(&Ring_Z, n, n);
    ;
    Matrix<Q> m1_q = m1.lift(QF);
    bignum r1 = QF.to_base_or_fail(m1_q.get_det());
    bignum r2 = m1.get_det_slow();
    OPA_CHECK_EQ(r1, r2, m1, "\n\n", m1_q);
  }
}

TEST(NumberField, Test1) {
  typedef NumberField NF;
  typedef NumberField::NumberFieldElem NFE;
  PolyRing<Q> pr;
  pr.init(&QF);

  P_Z p1 = PR_Z.import({ 1, 0, -2 }, kPolyRev);
  NF nf;
  nf.init(&QF, p1);

  NFE e;
  e.a = PR_Z.import({ 12, 5 });
  e.d = 1;
  auto cp_q = nf.get_char_poly(e);
  OPA_DISP0(cp_q);

  P_K pk = nf.K_x.import(
    { nf.import({ QF(1) }), nf.import({ QF(-10) }), nf.import({ QF(23) }) },
    kPolyRev);

  K pk_sol_root1 = nf.import({ QF(5), QF(1) }); // 5 + sqrt(2)
  OPA_CHECK_EQ(pk(pk_sol_root1), Q_x.getZ(), pk_sol_root1, pk);
  P_K pk_sol_factor1 = nf.K_x.x() - pk_sol_root1;
  OPA_CHECK_EQ(pk % pk_sol_factor1, nf.K_x.getZ(), pk, pk_sol_factor1);
  P_K pk_sol_factor2 = pk / pk_sol_factor1;
  K pk_sol_root2 = pk_sol_factor2.monic()[0];
  OPA_DISP0(pk_sol_root1, pk_sol_root2);

  int k = 2;
  P_QQ r1_rmp = nf.get_axy(pk_sol_factor1, k);
  P_QQ r2_rmp = nf.get_axy(pk_sol_factor2, k);
  OPA_DISP0(r1_rmp, r2_rmp);
  P_Q n1 = nf.get_norm_p_qq(r1_rmp);
  P_Q n2 = nf.get_norm_p_qq(r2_rmp);
  OPA_DISP0(n1, n2);
  OPA_CHECK0(is_squarefree(n1));
  OPA_CHECK0(is_squarefree(n2));

  std::vector<P_Z> fx1, fx2;
  factor_qpoly_squarefree(n1, &fx1);
  factor_qpoly_squarefree(n2, &fx2);
  OPA_DISP0(fx1, fx2);

  P_K eval_for_remap = nf.K_x.x() + nf.K_x.constant(nf.x() * QF(k));

  P_Q fq1 = to_pq(fx1[0]);
  P_K fk1 = import_change_var(nf.K_x, fq1);
  P_KK fpk1 = import_change_var(nf.K_xx, fk1);
  P_K fk1_rmp = fpk1(eval_for_remap);
  OPA_DISP0(fk1_rmp % pk_sol_factor1);
  OPA_DISP0(fk1_rmp % pk_sol_factor2);

  P_Q fq2 = to_pq(fx2[0]);
  P_K fk2 = import_change_var(nf.K_x, fq2);
  P_KK fpk2 = import_change_var(nf.K_xx, fk2);
  P_K fk2_rmp = fpk2(eval_for_remap);
  OPA_DISP0(fk2_rmp % pk_sol_factor1);
  OPA_DISP0(fk2_rmp % pk_sol_factor2);

  P_QQ r1r2_rmp = nf.get_axy(pk, k);
  P_Q n1n2 = nf.get_norm_p_qq(r1r2_rmp);
  OPA_DISP0(r1r2_rmp, n1n2);
  std::vector<P_Z> fx1x2;
  factor_qpoly_squarefree(n1n2, &fx1x2);
  OPA_DISP0(fx1x2, n1n2, n1 * n2);

  std::vector<std::pair<P_K, int> > factors;
  nf.factor_k(pk, &factors);
  OPA_DISP0(factors);
  for (auto &factor_and_pw : factors) {
    const auto &factor = factor_and_pw.first;
    if (factor.deg() != 1) continue;
    auto ifx = nf.K_x.import(factor);
    K e = nf.neg(ifx.monic()[0]);
    OPA_DISP0(e);
  }
}

TEST(NumberField, Test2) {
  typedef NumberField NF;
  typedef NumberField::NumberFieldElem NFE;
  typedef NumberField::P_K P_K;
  typedef NumberField::P_KK P_KK;
  typedef NumberField::K K;
  PolyRing<Q> pr;
  pr.init(&QF);

  P_Z p1 = PR_Z.import({ 1, 0, -3 }, kPolyRev);
  NF nf;
  nf.init(&QF, p1);

  P_K pk = import_change_var(nf.K_x, to_pq(p1));
  std::vector<std::pair<P_K, int> > factors;
  nf.factor_k(pk, &factors);
  OPA_DISP0(factors);
  P_K rem;
  for (auto &factor_and_pw : factors) {
    const auto &factor = factor_and_pw.first;
    if (factor.deg() != 2) continue;
    rem = factor.monic();
  }
  // OPA_CHECK0(rem.deg() == 2);

  P_Z p2 = PR_Z.import({ 1, 0, -2 });
  P_K rem2 = import_change_var(nf.K_x, to_pq(p2));

  K r1, r2;
  NF nf2;
  nf.extend(rem2, &nf2, &r1, &r2);
  OPA_DISP0(r1, r2);
  // P_K np1 = nf2.K_x.x() - nf2.K_x.constant(r1);
  // P_K mod_lift = import_change_var(nf2.K_x, nf2.mod);
  // OPA_DISP0(mod_lift % np1);
}

TEST(NumberField, FindRoot) {
  P_Z p1 = Z_x.import({ 1, 0, -2 }, kPolyRev);
  RootFinder<Float> rf(Float::Float_10pw(FLAGS_opa_roots_prec));
  rf.max_step = FLAGS_opa_roots_steps;
  rf.find_poly_roots(p1);
  OPA_DISP0(rf.real_roots);
  OPA_DISP0(rf.complex_roots);
  OPA_CHECK0(rf.real_roots.size() == 2);

  p1 = Z_x.import({ 1, 0, 0, -3 }, kPolyRev);
  rf.find_poly_roots(p1);
  OPA_DISP0(rf.real_roots);
  OPA_DISP0(rf.complex_roots);
  OPA_CHECK0(rf.real_roots.size() == 1);
}

TEST(NumberField, RootComposition) {
  P_Z p1 = Z_x.import({ 1, 0, -2 }, kPolyRev);
  RootFinder<Float> rf(Float::Float_10pw(FLAGS_opa_roots_prec));
  rf.max_step = FLAGS_opa_roots_steps;
  rf.find_poly_roots(p1);
  Float r1 = rf.real_roots[0];
  auto t1 = rf.all_roots;

  P_Z p2 = Z_x.import({ 1, 0, 0, -3 }, kPolyRev);
  rf.find_poly_roots(p2);
  Float r2 = rf.real_roots[0];
  auto t2 = rf.all_roots;

  P_Q pq1 = to_pq(p1);
  P_QQ pqq1 = import_change_var(Q_xy, pq1);
  P_Q pq2 = to_pq(p2);
  P_QQ pqq2 = import_change_var(Q_xy, pq2);
  P_QQQ pqqq2 = import_change_var(Q_xyz, pqq2);
  P_QQ eval_at = Q_xy.constant(Q_x.x()) - Q_xy.x();
  P_QQ tx = pqqq2(eval_at); // p2(X-Y)

  P_Q res = resultant2(pqq1, tx);
  std::vector<std::pair<P_Z, int> > lst;
  factor_qpoly(res, &lst);
  OPA_DISP0(r1, r2);
  OPA_DISP0(lst);

  rf.find_poly_roots(lst[0].first);
  OPA_DISP0(rf.real_roots);

  std::vector<std::complex<Float> > tb;
  for (auto &a1 : t1)
    for (auto &a2 : t2) tb.push_back(a1 + a2);
  std::sort(ALL(tb));
  auto cmp = FloatEQCompare<Float>(1e-10);
  opa::utils::make_unique(tb, cmp);
  OPA_CHECK_EQ(tb.size(), rf.all_roots.size(), tb, rf.all_roots);

  REP (i, tb.size()) {
    OPA_CHECK(cmp(tb[i], rf.all_roots[i]), i, tb[i], rf.all_roots[i]);
  }
  OPA_DISP0(tb);
  OPA_DISP0(rf.all_roots);
}

TEST(NumberField, PrimElem) {
  P_Z p1 = Z_x.import({ 1, 0, -2 }, kPolyRev);
  P_Z p2 = Z_x.import({ 1, 0, 0, -3 }, kPolyRev);
  P_Q pq1 = to_pq(p1);
  P_QQ pqq1 = import_change_var(Q_xy, pq1);
  P_Q pq2 = to_pq(p2);
  P_QQ pqq2 = import_change_var(Q_xy, pq2);
  P_QQQ pqqq2 = import_change_var(Q_xyz, pqq2);

  RootFinder<Float> rf(Float::Float_10pw(FLAGS_opa_roots_prec));
  rf.max_step = FLAGS_opa_roots_steps;
  rf.find_poly_roots(p1);
  OPA_CHECK(rf.real_roots.size() > 0, rf.all_roots);
  Float r1v = rf.real_roots[0];
  rf.find_poly_roots(p2);
  OPA_CHECK(rf.real_roots.size() > 0, rf.all_roots);
  Float r2v = rf.real_roots[0];

  for (int k = 1;; ++k) {
    P_QQ eval_at = Q_xy.constant(Q_x.x()) - Q_xy.xpwv(1, Q_x.constant(QF(k)));
    P_QQ tx = pqqq2(eval_at); // p2(X-Y)
    P_Q res = resultant2(pqq1, tx);
    std::vector<std::pair<P_Z, int> > lst;
    factor_qpoly(res, &lst);
    OPA_DISP0(k, lst);

    auto px = lst[0].first;
    NumberField nf;
    nf.init(&QF, px);
    std::vector<P_K> f1, f2;
    nf.factor_squarefree(import_change_var(nf.K_x, pq1), &f1);
    nf.factor_squarefree(import_change_var(nf.K_x, pq2), &f2);
    OPA_DISP0(f1, f2);

    if (has_linear_term(f1) && has_linear_term(f2)) {
      auto r1 = get_linear_roots(f1)[0];
      auto r2 = get_linear_roots(f2)[0];
      P_Qf pr1 = import_poly(Qf_x, r1);
      P_Qf pr2 = import_poly(Qf_x, r2);

      Float rxv = r2v + r1v * k;
      rf.find_poly_roots(px);
      OPA_DISP0(rxv, rf.real_roots);
      OPA_DISP0(r1v, pr1(rxv));
      OPA_DISP0(r2v, pr2(rxv));
      break;
    }
    puts("\n\n\n");
  }
}

P_Z pz_find_irred(int deg, const bignum &bound) {
  utils::ScopedPush<Z> s1(zrand_bound, bound);
  while (true) {
    P_Z cnd = Z_x.rand(deg);
    OPA_DISP("Trying ", cnd);
    if (pz_is_irred(cnd)) return cnd;
  }
}

TEST(NumberField, Prob1) {
  int p_deg = 4;
  int n_vec = 5;
  int subspace_dim = 3;
  int kernel_dim = n_vec - subspace_dim;
  OPA_CHECK0(subspace_dim <= p_deg);
  int maxv = 1000000;
  int lim2 = 1000000;
  auto mat = generate_lattice(p_deg, subspace_dim, maxv);
  P_Z p = pz_find_irred(p_deg, bignum(lim2));
  auto vecs_lat = generate_lattice(subspace_dim, n_vec, maxv);
  auto vecs = mat * vecs_lat;

  OPA_DISP0(mat, vecs_lat, vecs);
  NumberField nf;
  nf.init(&QF, p);
  std::vector<std::vector<Q> > lst;
  REP (i, n_vec) {
    std::vector<Q> vv = import_vec(QF, vecs.get_col(i).tovec());
    auto v = Q_x.import(vv);

    auto mv = nf.minimal_poly(v);
    OPA_DISP0(v, mv);
    std::vector<std::pair<P_K, int> > factors;
    auto px = nf.K_x.x() - nf.K_x.constant(v);
    auto m_lifted = import_change_var(nf.K_x, mv);
    OPA_TRACES(nf.K_x.mod(m_lifted, px));
    nf.factor_k(m_lifted, &factors);
    OPA_DISP0(v, factors);

    lst.push_back(vv);
  }

  Matrix<Q> mm(&QF, lst.size(), p_deg);
  mm.set_rows(lst);
  auto kernel = mm.kernel_basis();
  OPA_DISP0(kernel);
  OPA_DISP0(kernel * mm);
}

TEST(BasicMat, HNF) {
  int n = 3, m = 6;
  Matrix<bignum> m1(&Ring_Z, n, m);
  m1.setRow(0, { 2, 1, 1, 2, 8, 12 });
  m1.setRow(1, { -3, 5, 1, 10, -2, 88 });
  m1.setRow(2, { 12, 1, 77, 2, 74, 33 });
  Matrix<bignum> h, u;
  bignum det;
  hnf(m1, &h, &u, &det);
  OPA_DISP0(h);
  OPA_DISP0(u);
  OPA_DISP0(u.get_det_slow());
  OPA_DISP0(det);

  auto null_space = u.get_submatrix(0, n);
  OPA_DISP0(m1 * null_space);
  null_space = null_space.transpose();

  Matrix<bignum> h2, u2;
  hnf(null_space, &h2, &u2, &det);
  OPA_DISP0(null_space, h2, u2, null_space.get_det_slow());
}

TEST(BasicMat, SNF) {
  int n = 3, m = 2;
  Matrix<bignum> m1(&Ring_Z, n, m);
  m1.setCol(0, { 8, 4, 8 });
  m1.setCol(1, { 4, 8, 4 });
  Matrix<bignum> p, d, q, ip, iq;
  snf(m1, &d, &p, &q);
  auto res = q * m1 * p;
  OPA_DISP0(p);
  OPA_DISP0(d);
  OPA_DISP0(q);
  OPA_DISP0(res);

  Matrix<bignum> lat_hnf;
  hnf(m1, &lat_hnf);

  OPA_CHECK(p.invert(&ip), p);
  OPA_CHECK(q.invert(&iq), q);
  OPA_DISP0(ip, iq, lat_hnf);
  std::vector<std::vector<bignum> > ranges;
  REP (i, n) {
    if (i < m) {
      ranges.push_back(
        utils::Range<bignum>::StepRange(0, d.get(i, i).abs(), 1).tb());
    } else
      ranges.push_back({ 0 });
  }

  auto v1 = m1.get_col(0).tovec();
  auto v2 = m1.get_col(1).tovec();
  auto mq = m1.lift(QF);
  Matrix<Q> imq;
  mq.invert(&imq);

  OPA_DISP0(Ring_Z.dot(v1, v1), Ring_Z.dot(v2, v2));
  OR::CrossProdGen<bignum>(ranges, [&](const std::vector<bignum> &v) {
    auto tmp = iq.eval(v);
    std::vector<bignum> tmp_reduced = tmp;
    reduce_hnf_vec(lat_hnf, tmp_reduced, nullptr, true);
    OPA_DISP0(v, tmp, tmp_reduced, Ring_Z.dot(v1, tmp_reduced),
              Ring_Z.dot(v2, tmp_reduced));
  });
}

TEST(BasicOp, ComplexOp) {
  std::complex<Float> v(-0.4086958773e0, 0.8168134000e0);
  std::complex<Float> dv(0.2600000000e1, 0.6283180000e0);
  OPA_DISP0(v / dv, CF->div(v, dv));
}

TEST(T1, RSCode) {
  RSCode_u32 rscode;
  GF_q_u32 gf256(&GF2, 8);
  OPA_DISP0(gf256.getModPoly(), &gf256);

  int n = 255;
  int k = 223;

  std::vector<u32> data;
  srand(0);
  REP (i, k * 8)
    data.push_back(rand() % 2);
  rscode.init(&gf256, n, k);
  std::vector<Poly_u32> c = rscode.encode(pack_vector_gfq(&gf256, data));
  c[0] = c[0] + gf256.getE();
  c[5] = c[5] + gf256.getE();
  c[n-1] = c[n-1] + gf256.getE();
  std::vector<Poly_u32> res_enc = rscode.decode(c);
  std::vector<u32> res =  unpack_vector_gfq(&gf256, res_enc);
  OPA_CHECK0(data == res);
}

GTEST_API_ int main(int argc, char **argv) {
  printf("Running main() from gtest_main.cc\n");
  testing::InitGoogleTest(&argc, argv);
  puts("1");
  opa::init::opa_init(argc, argv);
  puts("2");
  return RUN_ALL_TESTS();
}

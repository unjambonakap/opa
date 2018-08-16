#include <gtest/gtest.h>

#include <opa/crypto/dlp.h>
#include <opa/crypto/dlp_job.h>
#include <opa/crypto/lfsr.h>
#include <opa/math/common/GF_p.h>
#include <opa/math/common/GF_pBN.h>
#include <opa/math/common/Poly.h>
#include <opa/math/common/PolyRing.h>
#include <opa/math/common/Types.h>
#include <opa/math/common/Utils.h>
#include <opa/math/common/UtilsPoly.h>
#include <opa/threading/dispatcher.h>
#include <opa/threading/runner.h>
#include <opa_common.h>

using namespace opa::threading;
using namespace opa::crypto;
using namespace std;
using namespace std::chrono;
using namespace opa::math::common;

Dispatcher *g_dispatcher;

void init() { initMathCommon(0); }

/*
void do_test(const bignum &mod, const bignum &g, const bignum &x,
             const bignum &y, bool distributed) {

  cout << "mod=" << mod << endl;
  cout << "g=" << g << endl;
  cout << "x=" << x << endl;
  cout << "y=" << y << endl;

  DlpJob dlp_job;
  Dlp *dlp = dlp_job.dlp.get();
  dlp->setup_main(mod, g, y);

  auto t1 = high_resolution_clock::now();
  if (distributed) {
    g_dispatcher->process_job(dlp_job, Runner::GetJobId(DlpJob::JOB_NAME));
  } else {
    dlp->solve();
  }
  auto diff =
    duration_cast<milliseconds>(high_resolution_clock::now() - t1).count();
  cout << "TOOK: " << diff << endl;

  cout << "RESULT=" << dlp->ans << endl;

  ASSERT_EQ(x, dlp->ans);
  ASSERT_TRUE(dlp->check(dlp->ans));
}

void test_batch(const bignum &prime_size, int ntest, bool distributed) {
  init();

  REP (tt, ntest) {
    bignum mod = gen_prime(prime_size);
    bignum g = (mod - 1).rand() + 1;
    GF_pBN f(mod);
    bignum g_order = f.compute_order(g, f.get_factors());

    bignum x = g_order.rand();
    bignum y = f.faste(g, x);
    ASSERT_NO_FATAL_FAILURE(do_test(mod, g, x, y, distributed));
  }
}

void test_smallfixed(bool distributed) {
  init();
  u32 mod = 1e9 + 7;
  u32 x = sqrt(mod) + 2;
  GF_p f(mod);
  u32 g = f.getPrimitiveElem();
  u32 y = f.faste(g, x);
  do_test(mod, g, x, y, distributed);
}

void gen_hard(int nprimes, u64 psize, bignum &mod, bignum &g, bignum &x,
              bignum &y) {
  bignum ub = bignum::fromu64(psize);

  while (1) {
    mod = 2;
    REP (i, nprimes)
      mod *= gen_prime(ub.rand());
    mod += 1;
    if (testPrime(mod))
      break;
  }

  GF_pBN f(mod);
  g = f.getPrimitiveElem();
  x = (mod - 1).rand();
  y = f.faste(g, x);
}

TEST(Dlp, RandTest_Small) { test_batch(1e4, 100, false); }
TEST(Dlp, RandTest_Large) { test_batch(bignum::fromu64(1e14), 10, false); }
TEST(Dlp, RandTest_LargeDistributed) {
  test_batch(bignum::fromu64(1e16), 10, true);
}

TEST(Dlp, RandTest_Hard) {
  init();
  bignum mod, g, x, y;
  gen_hard(20, 1e8, mod, g, x, y);
  do_test(mod, g, x, y, false);
}

void test_hardfixed(bool distributed) {
  init();
  bignum mod, g, x, y;
  mod = bignum("2d1de2a090ffa9582a8cb662599348d2f3838b55ca49812c14127168ee54b"
               "ad3c5927955ac2db126c9d5caf6d6081fd73e4bfe473ef625b6f7a421b3",
               16);
  GF_pBN f(mod);
  g = f.getPrimitiveElem();
  x = (mod - 1).rand();
  y = f.faste(g, x);

  do_test(mod, g, x, y, distributed);
}

TEST(Dlp, RandTest_HardFixed) { test_hardfixed(false); }

TEST(Dlp, RandTest_HardFixedDistributed) { test_hardfixed(true); }

TEST(Dlp, RandTest_SmallFixed) { test_smallfixed(false); }

TEST(Dlp, RandTest_SmallFixedDistributed) { test_smallfixed(true); }
*/

TEST(LFSR, TestLfsr1) {
  typedef u32 T;
  GF_p field(5);
  PolyRing<T> pr(&field);
  T tb_poly[] = { 1, 0, 1, 0, 2, 1, 3, 4, 1, 2, 3 };
  vector<T> vec_poly(tb_poly, tb_poly + sizeof(tb_poly) / sizeof(tb_poly[0]));
  Poly<T> poly = pr.import(vec_poly);
  Poly<T> state = pr.constant(1);
  cout << poly << endl;

  LFSR<T> lfsr(&field, state, poly);
  std::vector<T> seq;
  REP (i, 30)
    seq.pb(lfsr.get_next());

  auto x = LFSR<T>::fromSequence(&field, seq);
  ASSERT_TRUE(x.get() != nullptr);

  cout << x->get_poly() << endl;
  cout << x->get_state() << endl;
  cout << "ORIG >> " << lfsr.get_state() << endl;
}

TEST(LFSR, TestLfsr2) {
  Poly<u32> state = PR_GF2.import({ 1 });
  int deg = 6;
  Poly<u32> lfsr_poly = PR_GF2.import(find_primitive_poly(GF2, deg));
  OPA_DISP0(lfsr_poly);
  vector<u32> t1, t2;
  int n = (1 << lfsr_poly.deg()) - 1;

  {
    LFSR<u32> l1;
    l1.init(&GF2, state, lfsr_poly);
    REP (i, n)
      t1.pb(l1.get_next());

    l1.init(&GF2, state, lfsr_poly);
    REP (i, n)
      t2.pb(l1.get_next_non_galois());
  }

  OPA_DISP0(t1);
  OPA_DISP0(t2);
  REP (i, n) {
    bool ok = 1;
    REP (j, n)
      if (t1[j] != t2[(i + j) % n]) {
        ok = 0;
        break;
      }
    if (!ok)
      continue;
    OPA_DISP("ok at ", i);
  }
}

TEST(LFSR, TestLfsrGF2) {
  REP (i, 10) {
    auto l1 = LFSR<u32>::rand(10 + i, &GF2);
    auto l2 = LFSR_GF2_small::from_lfsr(*l1);
    REP (j, 1000) {
      ASSERT_EQ(l1->get_next(), l2.get_next()) << "Fail at " << j;
    }
  }
}

TEST(LFSR, TestLfsrSeq) {
  REP (ntest, 10) {
    int deg = 6 + ntest;
    Poly<u32> state = PR_GF2.rand(deg - 1);
    Poly<u32> lfsr_poly = PR_GF2.import(find_primitive_poly(GF2, deg));
    OPA_DISP0(lfsr_poly);
    LFSR<u32> l1;
    l1.init(&GF2, state, lfsr_poly);

    vector<u32> obs;
    REP (i, deg + 10) { obs.push_back(l1.get_next()); }

    Poly<u32> istate = l1.get_initial_state(obs);
    LFSR<u32> l2;
    l2.init(&GF2, istate, lfsr_poly);
    vector<u32> nobs;
    REP (i, deg + 10) { nobs.push_back(l2.get_next()); }
    OPA_DISP0(state, istate);
    OPA_CHECK(state == istate, "Wanted ", state, istate);
  }
}

GTEST_API_ int main(int argc, char **argv) {
  // Runner runner;
  testing::InitGoogleTest(&argc, argv);
  opa::init::opa_init(argc, argv);

  // DlpJob::Register();
  // Runner::Build();
  // runner.run_both();
  // Dispatcher *dispatcher = runner.dispatcher();
  // g_dispatcher = dispatcher;
  return RUN_ALL_TESTS();
}

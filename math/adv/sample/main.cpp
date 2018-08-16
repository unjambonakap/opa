#include <opa_common.h>
#include <opa/math/adv/nfs.h>

using namespace std;
using namespace opa::math::adv;
using namespace opa::math::common;

void test_nfs1() {
  bignum n;

  Poly<bignum> f;
  bignum m;

  if (0) {
    n = 45113;
    m = 31;
    f = PR_Z.import({ 1, 15, 29, 8 }, true);
  } else {
    bignum p1 = 1e9 + 9;
    bignum p2 = 1e9 + 7;
    if (0) {
      p1 = bignum::fromstr("100000000000000000039", 10);
      p2 = bignum::fromstr("1000000000000000000117", 10);
    }
    if (1) {
      p1 = 1e4 + 7;
      p2 = 1e4 + 9;
      p1 = 1e6 + 37;
      p2 = 1e6 + 33;
    }

    n = p1 *p2;
    n=4486873;
    int dx = 3;
    m = n.root(dx) - 1;
    OPA_CHECK(m.pow(dx) < n, m, n, dx);
    vector<bignum> f_coeffs;
    n.base_repr(m, &f_coeffs);

    f = PR_Z.import(f_coeffs);
    OPA_CHECK(PR_Z.is_monic(f), f, f_coeffs);
  }

  OPA_DISP0(n, m, f);
  NfsIdealHelper helper;
  NfsParams params(f, n, m);
  helper.setup(params);

  int px = 3000;
  helper.find_algebraic_ideals(px);
  helper.find_quadratic_ideals(px, px + 6000);
  helper.find_rational_ideals(px);

  NfsIdealParams ideal_params;
  helper.set_nfs_ideal_params(&ideal_params);
  OPA_DISP0(ideal_params);

  NfsSieveHelper sieve_helper;
  sieve_helper.setup(params, ideal_params);
  if (1) {
    SieveResultEntry test_entry(61, 1);
    OPA_DISP0(sieve_helper.check_entry(test_entry),
              sieve_helper.get_rational_val(test_entry),
              sieve_helper.get_algebraic_val(test_entry));

    bignum blow = 1;
    bignum alow = 1;
    sieve_helper.do_sieve(alow, alow + 100000, blow, blow + 10);
    sieve_helper.ensure_valid_results();
    for (auto &e : sieve_helper.m_result.checked_entries) {
      OPA_DISP0(e.a, e.b, sieve_helper.get_rational_val(e),
                sieve_helper.get_algebraic_val(e));
    }

    sieve_helper.setup_relations(sieve_helper.m_result,
                                 &sieve_helper.m_collector);

    OPA_MATH::Matrix<u32> kernel_mat;
    sieve_helper.m_collector.get_kernel(&kernel_mat);
    OPA_CHECK0(kernel_mat.getN() > 0);
    int nx = min<int>(20, kernel_mat.getN());
    REP (i, 1 << nx) {
      if (i == 0)
        continue;
      Matrix<u32> select_vec =
        Matrix<u32>::fromzerovec(&GF2, kernel_mat.getN(), VecType::Row);
      REP (j, nx)
        select_vec(j) = i >> j & 1;

      Matrix<u32> row0 = select_vec * kernel_mat;

      Matrix<u32> tmp = row0 * sieve_helper.m_collector.orig_rel_mat;
      NfsSquare square =
        sieve_helper.m_collector.vector_to_square(row0.tovec());

      printf("on vec %d\n", i);
      OPA_CHECK0(sieve_helper.verify_square(square));
      bignum res = sieve_helper.try_square(square);
      if (res != 1 && res != n) {
        OPA_DISP0("Found factor ", res, i);
        break;
      }
    }
    return;
  } else {
    // TODO: compute smooth decomposition for these pairs
    // + evaluate matrix for such pairs

    NfsSquare sq;
    sq.pairs = { { -1, 1 }, { 3, 1 },  { 13, 1 },  { 104, 1 }, { 3, 2 },
                 { 25, 2 }, { -8, 3 }, { 48, 5 },  { 54, 5 },  { -43, 6 },
                 { -8, 7 }, { 11, 7 }, { 856, 11 } };
    OPA_CHECK0(sieve_helper.verify_square(sq));

    bignum algebraic_root = sieve_helper.get_algebraic_root(sq);
    bignum rational_root = sieve_helper.get_rational_root(sq);

    bignum v1 = algebraic_root * algebraic_root % n;
    bignum v2 = rational_root * rational_root % n;
    bignum d1 = (algebraic_root - rational_root) % n;
    bignum d2 = (algebraic_root + rational_root) % n;
    bignum d = n.gcd(d1);
    OPA_DISP0(algebraic_root, rational_root, v1, v2, d1, d2, d);
  }
}

int main(int argc, char **argv) {
  opa::init::opa_init(argc, argv);
  test_nfs1();

  return 0;
}

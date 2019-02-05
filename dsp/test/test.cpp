#include <gtest/gtest.h>
#include <opa_common.h>
#include <opa/dsp/gps/l1ca.h>

using namespace std;
using namespace opa::dsp::gps;

TEST(GPS, L1CA) {
  auto x = Satellites::GetSingleton();
  REP (i, x->nsats()) {
    auto s = x->get(i);
    printf("%02d >> ", i);
    REP (j, 20)
      printf("%d", s.get_l1(j));
    puts("");
  }
}

GTEST_API_ int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  opa::init::opa_init(argc, argv);
  return RUN_ALL_TESTS();
}

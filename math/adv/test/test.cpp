#include <random>
#define OPA_HEX 0
#include "opa/utils/buffer_reader.h"
#include <gtest/gtest.h>

#include <fplll/fplll.h>
#include <opa/math/adv/lie.h>
#include <opa/math/common/matrix_utils.h>
#include <opa/math/common/utils_num.h>
#include <opa/or/grid_search.h>
#include <opa/utils/misc.h>

using namespace opa;
DEFINE_int32(maxd, 2, "");
DEFINE_string(fname_poly3d, "", "");
std::vector<u32> dim_list;

using namespace opa::math::common;
using namespace opa::math::adv;
using namespace std;

TEST(Lie, SO3_1) {}

GTEST_API_ int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  opa::init::opa_init(argc, argv);
  return RUN_ALL_TESTS();
}

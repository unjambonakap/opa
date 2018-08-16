#include <gtest/gtest.h>

#include <opa/utils/files.h>

using namespace opa::utils;
using namespace std;

void init() {}

TEST(FilenameSharder, Test1) {
  OPA_TRACE0;
  FilenameSharder sharder;
  sharder.set_pattern("pattern_{pad: 3}_end").build();
  OPA_TRACE0;
  OPA_CHECK0(sharder.valid());
  OPA_TRACE0;
  REP (i, 5) { OPA_TRACES(sharder.get(i)); }
  REP (i, 5) { OPA_TRACES(sharder.get()); }
  OPA_TRACE0;
}

GTEST_API_ int main(int argc, char **argv) {
  printf("Running main() from gtest_main.cc\n");
  testing::InitGoogleTest(&argc, argv);
  puts("1");
  opa::init::opa_init(argc, argv);
  puts("2");
  return RUN_ALL_TESTS();
}

#include <google/protobuf/io/coded_stream.h>
#include <opa_common.h>

OPA_NAMESPACE_DECL2(opa, init)

struct InitData {
  OpaInitRegister::InitFunc func;
  std::string name;
  bool init = false;
  InitData(OpaInitRegister::InitFunc func, std::string name, bool init)
      : func(func), name(name), init(init) {}
};

static std::vector<InitData> *get_init_funcs() {
  static std::vector<InitData> *tb = new std::vector<InitData>();
  return tb;
}

OpaInitRegister::OpaInitRegister(InitFunc init_func, const std::string &name) {
  get_init_funcs()->emplace_back(init_func, name, false);
}

namespace {
bool is_init = false;
}

void opa_init(int argc, char **argv) {
  if (!is_init) {
    google::ParseCommandLineFlags(&argc, &argv, true);
    google::InitGoogleLogging(argv[0]);
    is_init = true;
  }

  for (auto &x : *get_init_funcs()) {
    if (x.init)
      continue;
    x.func();
    x.init = true;
  }
}

OPA_NAMESPACE_DECL2_END

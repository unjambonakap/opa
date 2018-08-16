#include <opa/common_swig.h>

OPA_NM_INIT

void opa_init_swig(const std::vector<std::string> &tb) {
  std::vector<char *> data(tb.size(), nullptr);
  REP (i, tb.size()) {
    data[i] = new char[tb[i].size() + 1];
    memcpy(data[i], tb[i].data(), tb[i].size() + 1);
  }
  std::vector<char *> data_cp = data;
  opa_init(tb.size(), (char **)data.data());

  for (auto &d : data_cp) {
    delete[] d;
  }
}

OPA_NM_INIT_END

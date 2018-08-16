#include <opa/utils/string_builder.h>

OPA_NAMESPACE(opa, utils)

StringBuilder& StringBuilder::add(opa::stolen::StringRef s) {
  m_res += s;
  return *this;
}

StringBuilder& StringBuilder::add_num(u64 num, int nbytes) {
  return add(std::string((const char*)&num, nbytes));
}

StringBuilder& StringBuilder::add_null(int n) { return add(std::string(n, 0)); }

StringBuilder& StringBuilder::add_repeated(char c, int n) {
  return add(std::string(n, c));
}

StringBuilder& StringBuilder::pad(char c, int n) {
  if (size() > n) {
    m_res.resize(n);
  }
  return add_repeated(c, n - size());
}

StringBuilder& StringBuilder::padmod(char c, int mod) {
  int v = size() % mod;
  if (v == 0) return *this;
  return add_repeated(c, mod - v);
}

OPA_NAMESPACE_END(opa, utils)

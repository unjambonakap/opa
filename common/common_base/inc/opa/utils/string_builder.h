#pragma once

#include <opa_common_base.h>
#include <opa/stolen/StringRef.h>
#include <opa/utils/string_base.h>

OPA_NAMESPACE(opa, utils)

class StringBuilder {
public:
  const std::string &get() const { return m_res; }
  StringBuilder &add(opa::stolen::StringRef s);
  template <class T> StringBuilder &add2(const T &a) {
    return this->add(Conv::to_str(a));
  }

  StringBuilder &add_num(u64 num, int nbytes);
  StringBuilder &add_null(int n);
  StringBuilder &add_repeated(char c, int n);
  StringBuilder &pad(char c, int n);
  StringBuilder &padmod(char c, int mod);
  int size() const { return m_res.size(); }

private:
  std::string m_res;
};
OPA_NAMESPACE_END(opa, utils)

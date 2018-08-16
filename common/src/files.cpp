#include <opa/utils/files.h>
#include <opa/utils/string.h>

OPA_NAMESPACE_DECL2(opa, utils)

FilenameSharder &FilenameSharder::set_defaults() {
  m_valid = false;
  m_pad = 2;
  m_count = 0;
  return *this;
}

FilenameSharder &FilenameSharder::build() {
  m_fmt = "%0" + toStr(m_pad) + "d";
  m_valid = true;
  return *this;
}

FilenameSharder &FilenameSharder::set_number_elem(int ne) {
  int base = 10;
  m_pad = 0;
  while (ne > 0)
    ++m_pad, ne /= base;
  return *this;
}

FilenameSharder &FilenameSharder::set_pattern(const std::string &pattern) {
  if (pattern.size() == 0)
    return *this;
  m_pattern = pattern;
  int pre = pattern.find('{');
  int end = pattern.find('}');
  OPA_CHECK(pre != std::string::npos && end != std::string::npos, pattern);
  m_prefix = pattern.substr(0, pre);
  m_suffix = pattern.substr(end + 1);

  m_config = YAML::Load(pattern.substr(pre + 1, end - pre - 1));
  if (m_config["pad"])
    m_pad = m_config["pad"].as<int>();
  return *this;
}

std::string FilenameSharder::get(int c) {
  OPA_CHECK0(m_valid);
  std::string res = m_prefix;
  res += StringPrintf(StringRef(m_fmt), c);
  res += m_suffix;
  return res;
}

OPA_NAMESPACE_DECL2_END

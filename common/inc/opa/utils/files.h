#pragma once

#include <opa_common.h>

#include <glib/platform/posix/posix_file_system.h>
#include <type_traits>
#include <yaml-cpp/yaml.h>

OPA_NAMESPACE_DECL2(opa, utils)

class FilenameSharder {
public:
  FilenameSharder() { set_defaults(); }
  FilenameSharder &set_defaults();
  FilenameSharder &set_pattern(const std::string &pattern);
  FilenameSharder &set_number_elem(int ne);
  FilenameSharder &build();

  std::string get(int c);
  std::string get() { return get(m_count++); }
  OPA_ACCESSOR_R(bool, m_valid, valid);

private:
  std::string m_pattern;
  std::string m_prefix;
  std::string m_suffix;
  std::string m_fmt;
  int m_pad;
  int m_count;
  YAML::Node m_config;
  bool m_valid = false;
};

OPA_NAMESPACE_DECL2_END

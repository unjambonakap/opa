#pragma once

#include <opa/utils/string.h>
#include <opa/utils/serialize.h>
#include <glib/core/stringpiece.h>
#include <glib/core/raw_coding.h>

OPA_NAMESPACE(opa, utils)

static std::string read_file(glib::StringPiece filename) {
  std::ifstream is{ filename.data(), std::ios::binary | std::ios::ate };
  OPA_CHECK0(is);
  auto size = is.tellg();
  std::string str(size, '\0'); // construct string to stream size
  is.seekg(0);
  is.read(&str[0], size);
  return str;
}

class BufferReader {
public:
  BufferReader(glib::StringPiece str) : m_str(str) {}

  glib::StringPiece consume(int len) { return m_str.substr(m_pos, len); }

  const void *get(int advance) {
    const void *res = (const void *)(m_str.data() + m_pos);
    m_pos += advance;
    return res;
  }

  template <typename T> T read() { return *(T *)get(sizeof(T)); }
  template <typename T> std::vector<T> readn(int count) {
    T *buf = (T *)get(sizeof(T) * count);
    return std::vector<T>(buf, buf + count);
  }
  u8 read_u8() { return read<u8>(); }
  u16 read_u16() { return read<u16>(); }
  u32 read_u32() { return read<u32>(); }
  u64 read_u64() { return read<u64>(); }
  float read_f32() { return read<float>(); }
  double read_f64() { return read<double>(); }

private:
  glib::StringPiece m_str;
  int m_pos=0;
};

OPA_NAMESPACE_END(opa, utils)

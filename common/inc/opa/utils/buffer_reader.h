#pragma once

#include <opa/utils/string.h>
#include <opa/utils/serialize.h>
#include <string_view>

OPA_NAMESPACE(opa, utils)

static std::string read_file(string_view filename) {
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
  BufferReader(string_view str) : m_str(str) {}

  string_view consume(int len) { return m_str.substr(m_pos, len); }

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
  string_view m_str;
  int m_pos=0;
};

class ReaderC8 {
public:
  ReaderC8() {}
  ~ReaderC8();
  ReaderC8(StringRef filename);
  void init(StringRef filename);

  std::complex<s8> read1();

  template <typename T> T read() {
    T res;
    OPA_CHECK0(read(res));
    return res;
  }
  template <typename T> std::vector<T> readn(int count) {
    std::vector<T> res(count);
    size_t got =fread(res.data(), sizeof(T), count, m_f);
    OPA_CHECK(got == count, got, count);
    return res;
  }

  template <typename T> bool read(T &res) {
    return fread(&res, sizeof(res), 1, m_f) == sizeof(T);
  }
  size_t size() const;

  bool read1(std::complex<s8> &e);
  bool is_eof() const;

private:
  FILE *m_f = nullptr;
};

class Writer {
public:
  Writer() {}
  ~Writer();
  Writer(StringRef filename);
  Writer(Writer &&writer){
    m_f = writer.m_f;
    writer.m_f = nullptr;
  }

  Writer(const Writer &writer) = delete;
  Writer& operator=(const Writer & writer) = delete;
  void init(StringRef filename);

  template <typename T> bool write(const T &e) {
    return fwrite(&e, 1, sizeof(e), m_f) == sizeof(e);
  }

private:
  FILE *m_f = nullptr;
};

OPA_NAMESPACE_END(opa, utils)

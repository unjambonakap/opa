#pragma once

#include <absl/strings/substitute.h>
#include "opa/utils/string_base.h"
#include <opa/utils/string.h>
#include <opa/utils/serialize.h>

OPA_NAMESPACE(opa, utils)

static void write_file(std::string_view filename, std::string_view content) {
  std::ofstream os{ filename.data(), std::ios::binary | std::ios::ate };
  OPA_CHECK0(os);
  os << content;
}

class BufferWriter {
public:
  virtual ~BufferWriter() {}
  BufferWriter() {}
  std::string get() { return m_ss.str(); }

  void put(std::string_view content) { m_ss << content; }

  template <typename T> void write(const T &v) { m_ss.write((const char*)&v, sizeof(v)); }
  template <typename T> void write_repeated(const T &v, int count) {
    REP (i, count)
      m_ss.write((const char*)&v, sizeof(v));
  }

  template <typename T> std::vector<T> writen(const std::vector<T> &tb) {
    m_ss.write(tb.data(), tb.size() * sizeof(T));
  }
  void write_u8(u8 v) { return write<u8>(v); }
  void write_u16(u16 v) { return write<u16>(v); }
  void write_u32(u32 v) { return write<u32>(v); }
  void write_u64(u64 v) { return write<u64>(v); }
  void write_f32(float v) { return write<float>(v); }
  void write_f64(double v) { return write<double>(v); }

  template <class T> BufferWriter &operator<<(const T &v) {
    write(v);
    return *this;
  }

private:
  std::ostringstream m_ss;
};

class BufferFileWriter : public BufferWriter {
public:
  BufferFileWriter(std::string_view filename) { m_filename = filename; }

  virtual ~BufferFileWriter() { write_file(m_filename, this->get()); }

  std::string m_filename;
};

template <class T>
void dump_vector_to_files(
  std::string_view file_format, const std::vector<T> &lst,
  std::function<void(BufferWriter &writer, const T &)> func) {
  REP (i, lst.size()) {
    std::string filename = absl::Substitute("$0.$1", file_format, i);
    BufferFileWriter writer(filename);
    func(writer, lst[i]);
  }
}


OPA_NAMESPACE_END(opa, utils)

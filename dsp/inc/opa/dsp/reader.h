#pragma once

#include <opa/dsp/inc.h>

OPA_NAMESPACE(opa, dsp)

class ReaderC8 {
public:
  ReaderC8() {}
  ~ReaderC8();
  ReaderC8(StringRef filename);
  void init(StringRef filename);

  Complex<s8> read1();

  template <typename T> T read() {
    T res;
    OPA_CHECK0(read(res));
    return res;
  }

  template <typename T> bool read(T &res) {
    return fread(&res, 1, sizeof(res), m_f) == sizeof(T);
  }

  bool read1(Complex<s8> &e);
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
OPA_NAMESPACE_END(opa, dsp)

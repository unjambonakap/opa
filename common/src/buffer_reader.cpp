
#include <opa/utils/buffer_reader.h>

OPA_NAMESPACE(opa, utils)
ReaderC8::~ReaderC8() {
  if (m_f) fclose(m_f);
}

ReaderC8::ReaderC8(StringRef filename) { init(filename); }

void ReaderC8::init(StringRef filename) {
  m_f = fopen(filename.data(), "rb");
  OPA_CHECK0(m_f);
}

std::complex<s8> ReaderC8::read1() {
  std::complex<s8> res;
  read1(res);
  return res;
}

bool ReaderC8::read1(std::complex<s8> &e) {
  OPA_CHECK0(!is_eof());
  static_assert(sizeof(e) == 2, "wut");
  return fread(&e, sizeof(e), 1, m_f) == sizeof(e);
}

size_t ReaderC8::size() const {
  fseek(m_f, 0, SEEK_END);
  size_t res = ftell(m_f);
  fseek(m_f, 0, SEEK_SET);
  return res;
}

bool ReaderC8::is_eof() const { return feof(m_f); }

Writer::~Writer() {
  if (m_f) fclose(m_f);
}

Writer::Writer(StringRef filename) { init(filename); }

void Writer::init(StringRef filename) {
  m_f = fopen(filename.data(), "wb");
  OPA_CHECK0(m_f);
}

OPA_NAMESPACE_END(opa, utils)

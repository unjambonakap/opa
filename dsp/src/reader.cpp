#include <opa/dsp/reader.h>

using namespace std;
using namespace opa::utils;

OPA_NAMESPACE(opa, dsp)
ReaderC8::~ReaderC8() {
  if (m_f)
    fclose(m_f);
}

ReaderC8::ReaderC8(StringRef filename) { init(filename); }

void ReaderC8::init(StringRef filename) {
  m_f = fopen(filename.data(), "rb");
  OPA_CHECK0(m_f);
}

Complex<s8> ReaderC8::read1() {
  Complex<s8> res;
  read1(res);
  return res;
}

bool ReaderC8::read1(Complex<s8> &e) {
  OPA_CHECK0(!is_eof());
  static_assert(sizeof(e) == 2, "wut");
  return fread(&e, 1, sizeof(e), m_f) == sizeof(e);
}

bool ReaderC8::is_eof() const { return feof(m_f); }

Writer::~Writer() {
  if (m_f)
    fclose(m_f);
}

Writer::Writer(StringRef filename) { init(filename); }

void Writer::init(StringRef filename) {
  m_f = fopen(filename.data(), "wb");
  OPA_CHECK0(m_f);
}

OPA_NAMESPACE_END(opa, dsp)

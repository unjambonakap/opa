#pragma once

#include <opa/crypto/base.h>

OPA_NM_CRYPTO_CRACKER

class CrackerState {
public:
  virtual void decrypt(const u8 *src, u8 *dest, int n) = 0;
  virtual void encrypt(const u8 *src, u8 *dest, int n) = 0;
  virtual void reset(const u8 *src, int n) = 0;
  virtual ~CrackerState(){}
};

class OutputChecker {
public:
  virtual bool check(const u8 *output, int n) = 0;
  virtual ~OutputChecker(){}
};

class AsciiChecker : public OutputChecker {
public:
  virtual bool check(const u8 *output, int n);
};

class Cracker2 {

public:
  void set_state(CrackerState *state) { m_state = state; }
  void set_checker(OutputChecker *checker) { m_checker = checker; }

  void do_dictionary(const std::string &filename);

  void do_bruteforce(const u8 *charset, int n, int nchar);

  void set_content(const u8 *src, int n);

  Cracker2();
  ~Cracker2();

private:
  void release();
  void go(int pos);
  bool do_check(const u8 *buf, int n);

  OutputChecker *m_checker;
  CrackerState *m_state;
  u8 *m_content;
  u8 *m_encbuf;
  int m_content_len;

  const u8 *m_charset;
  u8 *m_data;
  int m_numchars;
  int m_charset_len;
};

OPA_NM_CRYPTO_CRACKER_END


#include "cracker.h"

OPA_NM_CRYPTO_CRACKER

bool AsciiChecker::check(const u8 *output, int n) {
  int cnt = 0;
  REP (i, n)
    cnt += isprint(output[i]) != 0;

  return cnt >= 9 * n / 10;
}

Cracker2::~Cracker2() { release(); }

Cracker2::Cracker2() {
  m_checker = 0;
  m_state = 0;
  m_content = 0;
  m_encbuf = 0;
}

void Cracker2::release() {
  free(m_encbuf);
  free(m_content);
  m_encbuf = 0;
  m_content = 0;
}

void Cracker2::set_content(const u8 *src, int n) {
  release();
  m_encbuf = (u8 *)malloc(n);
  m_content = (u8 *)malloc(n);
  memcpy(m_content, src, n);
  m_content_len = n;
}

void Cracker2::do_dictionary(const std::string &filename) {
  FILE *f = fopen(filename.c_str(), "r");
  const int buf_len = 256;
  char buf[buf_len];

  int count = 0;
  while (fgets(buf, buf_len, f)) {
    ++count;
    if (count % 1000 == 0)
      printf("step %d\n", count);
    int len = strlen(buf);
    assert(len < buf_len - 1);
    if (buf[len - 1] == '\n')
      buf[--len] = 0;

    if (do_check((u8 *)buf, len))
      printf("found for %s >> %s\n", buf, m_encbuf);
  }

  fclose(f);
}

bool Cracker2::do_check(const u8 *buf, int n) {
  m_state->reset((const u8 *)buf, n);
  m_state->decrypt(m_content, m_encbuf, m_content_len);
  return m_checker->check(m_encbuf, m_content_len);
}

void Cracker2::go(int pos) {

  if (pos == m_numchars) {
    if (do_check(m_data, pos)) {
      printf("ok  for %s\n", m_data);
      printf(">> result %s\n", m_encbuf);
    }
    return;
  }

  REP (i, m_charset_len) {
    m_data[pos] = m_charset[i];
    go(pos + 1);
  }
}

void Cracker2::do_bruteforce(const u8 *charset, int n, int nchar) {
  m_charset = charset;
  m_charset_len = n;
  m_numchars = nchar;
  m_data = (u8 *)alloca(n + 1);
  m_data[nchar] = 0;
  go(0);
}
OPA_NM_CRYPTO_CRACKER_END

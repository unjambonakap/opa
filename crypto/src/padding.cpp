#include <opa/crypto/padding.h>
#include <opa/utils/string.h>

#include <opa/crypto/base.h>

using opa::stolen::StringRef;

OPA_NM_CRYPTO

std::pair<std::string, bool> rpkcs7(StringRef s, int bs) {

  std::string res;
  int n = s.size();
  int x = (u8)s[n - 1];
  if (!x || x > bs || x > n) goto fail;
  for (int i = n - 1; i >= n - x; --i)
    if (s[i] != x) goto fail;
  n -= x;
  res = std::string(s.data(), n);
  return MP(res, true);
fail:
  return MP(std::string(), false);
}

std::string pkcs7(StringRef s, int bs) {
  OPA_CHECK0(bs < 256);
  int u = bs - s.size() % bs;
  std::string res(s.size() + u, 'a');
  memcpy((void*)res.data(), s.data(), s.size());
  REP(i, u) res[s.size() + i] = u;
  return res;
}

std::string pkcs1(StringRef s, int sz) {
  if (sz == -1) sz = s.size();
  OPA_CHECK0(sz >= s.size());
  sz -= s.size();

  std::string res(3 + sz + 4 + s.size(), 'a');
  int pos = 0;
  res[pos++] = 1;
  res[pos++] = 0;

  REP(i, sz) res[pos++] = 0xff;
  res[pos++] = 0;

  *(s32*)(res.data() + pos) = s.size();
  pos += 4;
  memcpy((u8*)res.data() + pos, s.data(), s.size());
  return res;
}

std::pair<std::string, bool> rpkcs1(StringRef s) {
  s32 l;
  int n = s.size();
  int pos = 0;
  if (s[pos++] != 1) goto fail;
  if (s[pos++] != 0) goto fail;
  for (; pos < n && s[pos] == (char)0xff; ++pos)
    ;
  if (pos + 5 >= n || s[pos++] != 0) goto fail;

  l = *(s32*)(s.data() + pos);
  pos += 4;
  // instead of using asn1, just have the length of hash on 4 bytes (little
  // endian)
  if (pos + l != n) goto fail;

  return MP(std::string(s.data() + pos, l), true);

fail:
  return MP(std::string(), false);
}

std::string zeropad(opa::stolen::StringRef s, int sz) {
  OPA_CHECK0(false);
  return "";
}

OPA_NM_CRYPTO_END

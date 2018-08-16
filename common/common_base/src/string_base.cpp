#include <opa/utils/string_base.h>
#include <opa/utils/string_builder.h>

#include <arpa/inet.h>

OPA_NAMESPACE_DECL2(opa, utils)

const char *m64 =
  "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
int rm64[128];
using namespace std;

static INIT_FUNC_ATTR void init_string() {
  REP (i, 64)
    rm64[m64[i]] = i;
}

string b64e(const string &a, int n) {
  if (n == -1)
    n = a.size();

  int pad = (3 - n % 3) % 3;
  int m = 0;
  int sz = (n + 2) / 3 * 4 + pad;
  string res(sz, 0);

  char tmp[5] = { 0 };
  tmp[0] = a[n - 2];
  tmp[1] = a[n - 1];

  int i;
#define DO_ROUND(a, b)                                                         \
  {                                                                            \
    int x = (u8)a[i] << 16 | (u8)a[i + 1] << 8 | (u8)a[i + 2];                 \
    REP (j, 4)                                                                 \
      b[m + 3 - j] = m64[x & 0x3f], x >>= 6;                                   \
    m += 4;                                                                    \
  }
  for (i = 0; i < n - 2; i += 3)
    DO_ROUND(a, res);

  char *tmp2 = tmp - (n - 2);
  for (; i < n; i += 3)
    DO_ROUND(tmp2, res);
  REP (i, pad)
    res[sz - 1 - i] = '=';
  return res;
}

string b64d(const string &a) {
  int pad = 0;
  int n = a.size();
  while (a[n - pad - 1] == '=')
    ++pad;

  int sz = (n - pad) * 3 / 4 - pad;
  string res(sz + pad, 0);

  int m = 0;
  printf(">> %d %d\n", n, pad);
  assert((n - pad) % 4 == 0);
  for (int i = 0; i < n - pad;) {
    int x = 0;
    REP (j, 4)
      x = x << 6 | rm64[a[i++]];
    REP (j, 3)
      res[m + 2 - j] = x & 0xff, x >>= 8;
    m += 3;
  }
  res[sz] = 0;
  res.resize(sz);
  return res;
}

static s8 char2hex(s8 x) { return x <= '9' ? x - '0' : x - 'a' + 10; }

static s8 hex2char(s8 x) { return x <= 9 ? '0' + x : 'a' + x - 10; }

string h2b(const string &a) {
  int n = a.length();
  OPA_CHECK0(!(n & 1));
  string res(n / 2, 0);

  REP (i, n / 2)
    res[i] = char2hex(a[2 * i]) << 4 | char2hex(a[2 * i + 1]);
  return res;
}

string b2h(const u8 *a, int n) {
  string res(n * 2, 0);

  REP (i, n) {
    res[2 * i] = hex2char(a[i] >> 4);
    res[2 * i + 1] = hex2char(a[i] & 0xf);
  }
  return res;
}
string b2h(const string &a, int n) {
  if (n == -1)
    n = a.size();
  return b2h((const u8 *)a.c_str(), n);
}

#define FORMAT_FUNC(FUNC_NAME)                                                 \
  string FUNC_NAME(const char *fmt, ...) {                                     \
    va_list ap;                                                                \
    va_start(ap, fmt);                                                         \
    int n = vsnprintf(NULL, 0, fmt, ap);                                       \
    va_end(ap);                                                                \
    std::string s(n, 0);                                                       \
    va_start(ap, fmt);                                                         \
    vsprintf((char *)s.c_str(), fmt, ap);                                      \
    va_end(ap);                                                                \
    return s;                                                                  \
  }


FORMAT_FUNC(stdsprintf);

OPA_FORMAT_FUNC_STD(StringPrintf, opa::stolen::StringRef);
FORMAT_FUNC(format);

std::string JoinArray(const std::vector<std::string> &tb) {
  StringBuilder builder;

  std::string res;
  REP (i, tb.size()) {
    if (i != 0)
      builder.add2(',');
    builder.add(tb[i]);
  }
  return builder.get();
}

OPA_NAMESPACE_DECL2_END

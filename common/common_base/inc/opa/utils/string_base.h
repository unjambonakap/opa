#pragma once

#include <absl/strings/substitute.h>
#include <opa/stolen/StringRef.h>
#include <opa/utils/stacktrace.h>
#include <opa_common_base.h>
#include <opa_inc.h>
#include <string>

OPA_NAMESPACE_DECL2(opa, utils)

std::string StringPrintf(opa::stolen::StringRef fmt, ...);

#define OPA_FORMAT_FUNC_STD(FUNC_NAME, TYP)                                                        \
  std::string FUNC_NAME(TYP fmt, ...) {                                                            \
    va_list ap;                                                                                    \
    va_start(ap, fmt);                                                                             \
    int n = vsnprintf(NULL, 0, fmt.data(), ap);                                                    \
    va_end(ap);                                                                                    \
    std::string s(n, 0);                                                                           \
    va_start(ap, fmt);                                                                             \
    vsprintf((char *)s.c_str(), fmt.data(), ap);                                                   \
    va_end(ap);                                                                                    \
    return s;                                                                                      \
  }

class Conv {
public:
  template <class T> static void from_str(const opa::stolen::StringRef &a, T &out) { out = T(a); }

  template <class T> static T from_str2(const std::string &a) {
    T out;
    std::istringstream(a) >> out;
    return out;
  }

  template <class T> static std::string to_str(const T &e) {
    std::ostringstream os;
    os << e;
    return os.str();
  }
};

class Conv2 {
public:
  template <class From, class To> static To conv(const From &x) { return To(x); }
};

template <> inline void Conv::from_str(const opa::stolen::StringRef &a, int &out) {
  out = toInt(a);
}
inline u8 char2hex(u8 x) { return x <= '9' ? x - '0' : x - 'a' + 10; }
inline u8 hex2char(u8 x) { return x <= 9 ? '0' + x : 'a' + x - 10; }
inline bool is_hex_char(u8 x) { return (x >= '0' && x <= '9') || (x >= 'a' && x <= 'f'); }

std::string b64e(const std::string &a, int n = -1);
std::string b64d(const std::string &a);

std::string b2h(const std::string &a, int n = -1);
std::string h2b(const std::string &a);
std::string b2h(const u8 *a, int n);

std::string stdsprintf(const char *fmt, ...);
std::string format(const char *fmt, ...);

static const auto &SPrintf = stdsprintf;

template <class T> std::string get_bin_str(const T &v) {
  std::string res;
  REP (i, 8 * sizeof(T)) res += '0' + GETB_BLOCK(&v, i);
  return res;
}

template <class T> std::string dump_hex(const T &v) { return b2h((const u8 *)&v, sizeof(v)); }

std::string JoinArray(const std::vector<std::string> &tb);

static void _Join(const opa::stolen::StringRef &ref, std::ostringstream &os) {}

template <class T, typename... Args>
void _Join(const opa::stolen::StringRef &ref, std::ostringstream &os, const T &value,
           const Args &...args) {
  int size = sizeof...(Args);
  os << Conv::to_str(value);
  if (size > 0) os << ref.data();
  _Join(ref, os, args...);
}

template <typename... Args>
std::string Join(const opa::stolen::StringRef &ref, const Args &...args) {
  std::ostringstream os;
  _Join(ref, os, args...);
  return os.str();
}

template <class T> std::string Join(const opa::stolen::StringRef &ref, const std::vector<T> &tb) {
  std::ostringstream os;
  REP (i, tb.size()) {
    if (i) os << ref.data();
    os << Conv::to_str(tb[i]);
  }
  return os.str();
}

#if OPA_CPP14 == 1
namespace internal {
template <int N, typename... Args> class TemplateIterator {
public:
  template <class T> static void fill1(const opa::stolen::StringRef &ref, T &out) {
    Conv::from_str(ref, out);
  }

  static void fill(const std::vector<opa::stolen::StringRef> &tb, std::tuple<Args...> &out) {
    fill1(tb[N], std::get<N>(out));
    TemplateIterator<N - 1, Args...>::fill(tb, out);
  }

  static void tuple_to_string_vec(const std::tuple<Args...> &a, std::vector<std::string> &out) {
    out[N] = Conv::to_str(std::get<N>(a));
    TemplateIterator<N - 1, Args...>::tuple_to_string_vec(a, out);
  }
};

template <typename... Args> class TemplateIterator<-1, Args...> {
public:
  static void fill(const std::vector<opa::stolen::StringRef> &tb, std::tuple<Args...> &out) {}
  static void tuple_to_string_vec(const std::tuple<Args...> &a, std::vector<std::string> &out) {}
};
} // namespace internal

template <typename... Args>
std::tuple<Args...> string_vec_to_tuple(const std::vector<opa::stolen::StringRef> &tb) {
  std::tuple<Args...> res;
  internal::TemplateIterator<sizeof...(Args) - 1, Args...>::fill(tb, res);
  return res;
}

template <typename... Args>
std::vector<std::string> tuple_to_string_vec(const std::tuple<Args...> &a) {
  std::vector<std::string> res(sizeof...(Args));
  internal::TemplateIterator<sizeof...(Args) - 1, Args...>::tuple_to_string_vec(a, res);
  return res;
}

#endif

OPA_NAMESPACE_DECL2_END

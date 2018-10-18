#pragma once

#include <opa_inc.h>
#include <stdint.h>

#include "config.h"

#ifndef OPA_HAS_STACKTRACE
#define OPA_HAS_STACKTRACE 1
#endif

#define COMMA ,
#define PI std::acos(double(-1.))
#define OPA_MATH opa::math::common
#define OPA_CRYPTO opa::crypto
#define OPA_THREADING opa::threading
#define SPTR(t) std::shared_ptr<t>
#define UPTR(t) std::unique_ptr<t>
#define GETB_BLOCK(a, b)                                                       \
  (BIT_BLOCK(a, b) >> ((b) & (sizeof((a)[0]) * 8 - 1)) & 1)
#define BIT_BLOCK(a, b) ((a)[(b) / sizeof((a)[0]) / 8])
#define BIT_B(a, b) (1ull << ((b) & ((sizeof((a)[0]) * 8) - 1)))

#define OPA_BITSIGN(b) ((b)&1 ? 1 : -1)
#define OPA_SIGN(b) ((b) < 0 ? -1 : 1)

#define BIT0 0x01
#define BIT1 0x02
#define BIT2 0x04
#define BIT3 0x08
#define BIT4 0x10
#define BIT5 0x20
#define BIT6 0x40
#define BIT7 0x80
#define GETB(X, p) ((X) >> (p)&1)
#define SETB(X, p, v) ((X) = (((X) & (~(1 << (p)))) | ((v) << (p))))
#define SET0(X, p) ((X) &= ~MASK(p))
#define SET1(X, p) ((X) |= MASK(p))
#define MASK(p) ((1) << (p))

#if defined(SWIG)
#define OPA_SWIG 1
#else
#define OPA_SWIG 0
#endif

#define LIST2VEC(lst)                                                          \
  (std::vector<__typeof(lst[0])>(lst, lst + LIST_COUNT(lst)))

#define FE_LIST(type, it, action, ...)                                         \
  {                                                                            \
    type OPA_MACRO_TMP_var1[] = __VA_ARGS__;                                   \
    auto OPA_MACRO_TMP_var2 = LIST2VEC(OPA_MACRO_TMP_var1);                    \
    for (auto &it : OPA_MACRO_TMP_var2) {                                      \
      action;                                                                  \
    }                                                                          \
  }

#define OPA_TRACE0 fprintf(stderr, "%s:%s:%d\n", __FILE__, __func__, __LINE__)

#if OPA_HAS_STACKTRACE
#define OPA_PRINT_STACKTRACE()                                                 \
  fprintf(stderr, "OPA_STACKTRACE_BEGIN: \n%s\nOPA_STACKTRACE_END\n",          \
          opa::utils::get_stacktrace().c_str())
#else
#define OPA_PRINT_STACKTRACE()
#endif

#define OPA_NM_UTILS                                                           \
  namespace opa {                                                              \
  namespace utils {
#define OPA_NM_UTILS_END                                                       \
  }                                                                            \
  }
#define OPA_NM_INIT                                                            \
  namespace opa {                                                              \
  namespace init {
#define OPA_NM_INIT_END                                                        \
  }                                                                            \
  }

#define OPA_CASE(typ, action)                                                  \
  case typ: {                                                                  \
    action                                                                     \
  } break;

#define OPA_LINE_INFO                                                          \
  std::string(__FILE__) << ":" << __func__ << ":" << __LINE__

#define OPA_TRACE(msg, ...)                                                    \
  fprintf(stderr, "%s\n", OPA_TRACE_STR(msg, __VA_ARGS__).c_str())

#define OPA_TRACES(...)                                                        \
  fprintf(stderr,                                                              \
          "%s:%s:%d, msg=%s"                                                   \
          "\n",                                                                \
          __FILE__, __func__, __LINE__, OPA_DISP_STR(__VA_ARGS__).c_str())

#define OPA_TRACE_STR(msg, ...)                                                \
  OPA_STREAM_STR(OPA_LINE_INFO << " msg=" << msg << " >> "                     \
                               << OPA_DISP_STR(__VA_ARGS__))

#define OPA_TRACEM(msg)                                                        \
  fprintf(stderr, "%s:%s:%d, msg=" msg "\n", __FILE__, __func__, __LINE__)

#if OPA_DEBUG
#define OPA_ASSERT(cond, fmt, ...)                                             \
  do {                                                                         \
    if (!(cond)) {                                                             \
      fprintf(stderr, "%s:%s:%d, err=%s: " fmt "\n", __FILE__, __func__,       \
              __LINE__, #cond, ##__VA_ARGS__);                                 \
      OPA_PRINT_STACKTRACE();                                                  \
      assert(0);                                                               \
    }                                                                          \
  } while (0)

#define OPA_ASSERTNO(cond, fmt, ...)                                           \
  do {                                                                         \
    if (!(cond)) {                                                             \
      fprintf(stderr, "%s:%s:%d, err=%s, errno=%s: " fmt "\n", __FILE__,       \
              __func__, __LINE__, #cond, strerror(errno), ##__VA_ARGS__);      \
      OPA_PRINT_STACKTRACE();                                                  \
      assert(0);                                                               \
    }                                                                          \
  } while (0)
#else
#define OPA_ASSERT(cond, fmt, ...) assert(cond);
#define OPA_ASSERTNO(cond, fmt, ...) assert(cond);
#endif

#define OPA_ASSERTNO0(cond) OPA_ASSERTNO(cond, "")
#define OPA_ASSERT0(cond) OPA_ASSERT(cond, "")

#define OPA_ABORT(cond, fmt, ...)                                              \
  do {                                                                         \
    if (cond) {                                                                \
      fprintf(stderr, "Abort %s:%s:%d, err=%s: " fmt "\n", __FILE__, __func__, \
              __LINE__, #cond, ##__VA_ARGS__);                                 \
      OPA_PRINT_STACKTRACE();                                                  \
      std::abort();                                                            \
    }                                                                          \
  } while (0)

#define OPA_ABORTS(cond, ...)                                                  \
  do {                                                                         \
    if (cond) {                                                                \
      fprintf(stderr, "Abort %s:%s:%d, cond=%s, errno=%s: \n", __FILE__,       \
              __func__, __LINE__, #cond, strerror(errno));                     \
      fflush(stderr);                                                          \
      std::cerr << RAW_OPA_DISP_VARS(__VA_ARGS__) << std::endl;                \
      OPA_PRINT_STACKTRACE();                                                  \
      std::abort();                                                            \
    }                                                                          \
  } while (0)

#define OPA_ABORTNO(cond, fmt, ...)                                            \
  do {                                                                         \
    if (cond) {                                                                \
      fprintf(stderr, "Abort %s:%s:%d, err=%s, errno=%s: " fmt "\n", __FILE__, \
              __func__, __LINE__, #cond, strerror(errno), ##__VA_ARGS__);      \
      OPA_PRINT_STACKTRACE();                                                  \
      std::abort();                                                            \
    }                                                                          \
  } while (0)
#define OPA_ABORT0(cond) OPA_ABORT(cond, "")
#define OPA_ABORTNO0(cond) OPA_ABORTNO(cond, "")

#define OPA_CHECK0(cond) OPA_ABORT0(!(cond))
#define OPA_CHECK(cond, ...) OPA_ABORTS(!(cond), __VA_ARGS__)
#define OPA_CHECKNO0(cond) OPA_ABORTNO0(!(cond))

#define OPA_CHECK_EQ0(a, b) OPA_ABORTS(!((a) == (b)), "Equality failure ", a, b)
#define OPA_CHECK_EQ(a, b, ...)                                                \
  OPA_ABORTS(!((a) == (b)), "Equality failure ", a, b, __VA_ARGS__)

#define OPA_DECL_EQ_OPERATOR(CL, ...)                                          \
  bool operator==(const CL &other) const { OPA_EQ_OP(other, __VA_ARGS__); }
#define OPA_DECL_LT_OPERATOR(CL, ...)                                          \
  bool operator<(const CL &other) const { OPA_LT_OP(other, __VA_ARGS__); }

#define OPA_DECL_COUT_OPERATOR3(CL, ...)                                       \
  std::string str() const {                                                    \
    return OPA_STREAM_STR(RAW_OPA_DISP_VARS(__VA_ARGS__));                     \
  }                                                                            \
  friend std::ostream &operator<<(std::ostream &os, const CL &a) {             \
    os << a.str();                                                             \
    return os;                                                                 \
  }

#define OPA_DECL_COUT_OPERATOR(CL)                                             \
  friend std::ostream &operator<<(std::ostream &os, const CL &a) {             \
    os << a.str();                                                             \
    return os;                                                                 \
  }

#define OPA_DECL_STR_FROM_COUT()                                             \
  std::string str() const { return OPA_STREAM_STR(*this); }

#define OPA_DECL_COUT_OPERATOR2(CL, ...)                                       \
  friend std::ostream &operator<<(std::ostream &os, const CL &a) {             \
    os << RAW_OPA_DISP_VARS(__VA_ARGS__);                                      \
    return os;                                                                 \
  }

#define OPA_NAMESPACE_DECL1(a) namespace a {

#define OPA_NAMESPACE_DECL1_END }
#define OPA_NAMESPACE_DECL2(a, b)                                              \
  namespace a {                                                                \
  namespace b {

#define OPA_NAMESPACE_DECL2_END                                                \
  }                                                                            \
  }

#define OPA_NAMESPACE_DECL3(a, b, c)                                           \
  namespace a {                                                                \
  namespace b {                                                                \
  namespace c {

#define OPA_NAMESPACE_DECL3_END                                                \
  }                                                                            \
  }                                                                            \
  }

#define _OPA_NAMESPACE1(a) namespace a {
#define _OPA_NAMESPACE2(a, b) _OPA_NAMESPACE1(a) _OPA_NAMESPACE1(b)
#define _OPA_NAMESPACE3(a, b, c) _OPA_NAMESPACE1(a) _OPA_NAMESPACE2(b, c)

#define _OPA_NAMESPACE_END1(a) }
#define _OPA_NAMESPACE_END2(a, b) _OPA_NAMESPACE_END1(a) _OPA_NAMESPACE_END1(b)
#define _OPA_NAMESPACE_END3(a, b, c)                                           \
  _OPA_NAMESPACE_END1(a) _OPA_NAMESPACE_END2(b, c)

#define RAW_OPA_NAMESPACE(...)                                                 \
  OPA_DISPATCH(_OPA_NAMESPACE, __VA_ARGS__)(__VA_ARGS__)

#define RAW_OPA_NAMESPACE_END(...)                                             \
  OPA_DISPATCH(_OPA_NAMESPACE_END, __VA_ARGS__)(__VA_ARGS__)

#if 0 && defined(SWIG)
#define OPA_NAMESPACE(...)
#define OPA_NAMESPACE_END(...)

#else
#define OPA_NAMESPACE(...) RAW_OPA_NAMESPACE(__VA_ARGS__)
#define OPA_NAMESPACE_END(...) RAW_OPA_NAMESPACE_END(__VA_ARGS__)
#endif

#define _OPA_DISP_VAR1(a) #a << "=" << (a) << ","
#define _OPA_DISP_VAR2(a, b) _OPA_DISP_VAR1(a) << _OPA_DISP_VAR1(b)
#define _OPA_DISP_VAR3(a, b, c) _OPA_DISP_VAR2(a, b) << _OPA_DISP_VAR1(c)
#define _OPA_DISP_VAR4(a, b, c, d) _OPA_DISP_VAR3(a, b, c) << _OPA_DISP_VAR1(d)
#define _OPA_DISP_VAR5(a, b, c, d, e)                                          \
  _OPA_DISP_VAR4(a, b, c, d) << _OPA_DISP_VAR1(e)
#define _OPA_DISP_VAR6(a, b, c, d, e, f)                                       \
  _OPA_DISP_VAR5(a, b, c, d, e) << _OPA_DISP_VAR1(f)
#define _OPA_DISP_VAR7(a, b, c, d, e, f, g)                                    \
  _OPA_DISP_VAR6(a, b, c, d, e, f) << _OPA_DISP_VAR1(g)
#define _OPA_DISP_VAR8(a, b, c, d, e, f, g, h)                                 \
  _OPA_DISP_VAR7(a, b, c, d, e, f, g) << _OPA_DISP_VAR1(h)
#define _OPA_DISP_VAR9(a, b, c, d, e, f, g, h, i)                              \
  _OPA_DISP_VAR8(a, b, c, d, e, f, g, h) << _OPA_DISP_VAR1(i)
#define _OPA_DISP_VAR10(a, b, c, d, e, f, g, h, i, j)                          \
  _OPA_DISP_VAR9(a, b, c, d, e, f, g, h, i) << _OPA_DISP_VAR1(j)

#define RAW_OPA_DISP_VARS(...)                                                 \
  std::hex << std::showbase                                                    \
           << OPA_DISPATCH(_OPA_DISP_VAR, __VA_ARGS__)(__VA_ARGS__)

#define OPA_DISP_VARS(...)                                                     \
  (std::cout << RAW_OPA_DISP_VARS(__VA_ARGS__) << std::endl)
#define OPA_DISP(msg, ...)                                                     \
  (std::cout << msg << ": " << RAW_OPA_DISP_VARS(__VA_ARGS__) << std::endl)
#define OPA_DISP0(...) OPA_DISP_VARS(__VA_ARGS__)

#define OPA_DISP_VARSERR(...)                                                  \
  (std::cerr << RAW_OPA_DISP_VARS(__VA_ARGS__) << std::endl)
#define OPA_DISPERR(msg, ...)                                                  \
  (std::cerr << msg << ": " << RAW_OPA_DISP_VARS(__VA_ARGS__) << std::endl)
#define OPA_DISP0ERR(...) OPA_DISP_VARSERR(__VA_ARGS__)

#define OPA_DISPL0(...)                                                        \
  (std::cout << OPA_LINE_INFO << ": " << RAW_OPA_DISP_VARS(__VA_ARGS__)        \
             << std::endl)
#define OPA_DISPL(msg, ...)                                                    \
  (std::cout << OPA_LINE_INFO << "msg:" << msg << ": "                         \
             << RAW_OPA_DISP_VARS(__VA_ARGS__) << std::endl)

#define OPA_STREAM_STR(a)                                                      \
  ((static_cast<std::ostringstream &>(std::ostringstream() << a)).str())

#define OPA_DISP_VARS1(a) (std::cout << _OPA_DISP_VAR1(a) << std::endl)
#define OPA_DISP_VARS2(a, b) (std::cout << _OPA_DISP_VAR2(a, b) << std::endl)
#define OPA_DISP_VARS3(a, b, c)                                                \
  (std::cout << _OPA_DISP_VAR3(a, b, c) << std::endl)
#define OPA_DISP_VARS4(a, b, c, d)                                             \
  (std::cout << _OPA_DISP_VAR4(a, b, c, d) << std::endl)
#define OPA_DISP_VARS5(a, b, c, d, e)                                          \
  (std::cout << _OPA_DISP_VAR5(a, b, c, d, e) << std::endl)
#define OPA_DISP_VARS6(a, b, c, d, e, f)                                       \
  (std::cout << _OPA_DISP_VAR6(a, b, c, d, e, f) << std::endl)
#define OPA_DISP_VARS7(a, b, c, d, e, f, g)                                    \
  (std::cout << _OPA_DISP_VAR7(a, b, c, d, e, f, g) << std::endl)
#define OPA_DISP_VARS8(a, b, c, d, e, f, g, h)                                 \
  (std::cout << _OPA_DISP_VAR8(a, b, c, d, e, f, g, h) << std::endl)
#define OPA_DISP_VARS9(a, b, c, d, e, f, g, h, i)                              \
  (std::cout << _OPA_DISP_VAR9(a, b, c, d, e, f, g, h, i) << std::endl)
#define OPA_DISP_VARS10(a, b, c, d, e, f, g, h, i, j)                          \
  (std::cout << _OPA_DISP_VAR10(a, b, c, d, e, f, g, h, i, j) << std::endl)
#define OPA_DISP_STR(...) (OPA_STREAM_STR(RAW_OPA_DISP_VARS(__VA_ARGS__)))

#define OPA_DECL_SPTR_TMPL(cl, cl2)                                            \
  template <class T> using cl2 = std::shared_ptr<cl<T> >;
#define OPA_DECL_SPTR_TMPL2(cl, cl2)                                           \
  template <class T, class U> using cl2 = std::shared_ptr<cl<T, U> >;
#define OPA_DECL_SPTR_TMPL3(cl, cl2)                                           \
  template <class A, class B, class C>                                         \
  using cl2 = std::shared_ptr<cl<A, B, C> >;

#define OPA_DECL_SPTR(cl, cl2) typedef std::shared_ptr<cl> cl2;
#define OPA_ACCESSOR_PTR_DECL(parCl, cl, field, method)                        \
  cl *method();                                                                \
  const cl *method() const;

#define OPA_ACCESSOR_PTR_IMPL(parCl, cl, field, method)                        \
  cl *parCl::method() { return field; }                                        \
  const cl *parCl::method() const { return field; }

#define OPA_ACCESSOR_PTR(cl, field, method)                                    \
  cl *method() { return field; }                                               \
  OPA_ACCESSOR_PTR_R(cl, field, method)

#define OPA_ACCESSOR_PTR_R(cl, field, method)                                  \
  const cl *method() const { return field; }

#define OPA_ACCESSOR_NOSETTER(cl, field, method)                               \
  OPA_ACCESSOR_R(cl, field, method);                                           \
  OPA_ACCESSOR_W(cl, field, method);

#define OPA_ACCESSOR(cl, field, method)                                        \
  OPA_ACCESSOR_R(cl, field, method);                                           \
  OPA_ACCESSOR_W(cl, field, method);                                           \
  OPA_SETTER(cl, field, method);

#define OPA_ACCESSOR_R(cl, field, method)                                      \
  const cl &method() const { return field; }                                   \
  const cl &method##_const() const { return field; }
#define OPA_ACCESSOR_W(cl, field, method)                                      \
  cl &method() { return field; }

#define OPA_SETTER(cl, field, method)                                          \
  void set_##method(const cl &a) { field = a; }
#define OPA_ACCESSOR_R2(cl, field, method)                                     \
  const cl method() const { return field; }

#define _OPA_FIELDS1(x, a) x.a
#define _OPA_FIELDS2(x, a, b) x.a, _OPA_FIELDS1(x, b)
#define _OPA_FIELDS3(x, a, b, c) x.a, _OPA_FIELDS2(x, b, c)
#define _OPA_FIELDS4(x, a, b, c, d) x.a, _OPA_FIELDS3(x, b, c, d)
#define _OPA_FIELDS5(x, a, b, c, d, e) x.a, _OPA_FIELDS4(x, b, c, d, e)
#define _OPA_FIELDS6(x, a, b, c, d, e, f) x.a, _OPA_FIELDS5(x, b, c, d, e, f)
#define _OPA_FIELDS7(x, a, b, c, d, e, f, g)                                   \
  x.a, _OPA_FIELDS6(x, b, c, d, e, f, g)

#define RAW_OPA_FIELDS(x, ...)                                                 \
  OPA_DISPATCH(_OPA_FIELDS, __VA_ARGS__)(x, __VA_ARGS__)
#define OPA_FIELDS(x, ...) RAW_OPA_FIELDS(x, __VA_ARGS__)

#define _OPA_LT_OP1(x, a)                                                      \
  { return a < x.a; }
#define _OPA_LT_OP2(x, a, b)                                                   \
  {                                                                            \
    if (a != x.a) return a < x.a;                                              \
    _OPA_LT_OP1(x, b);                                                         \
  }
#define _OPA_LT_OP3(x, a, b, c)                                                \
  {                                                                            \
    if (a != x.a) return a < x.a;                                              \
    _OPA_LT_OP2(x, b, c);                                                      \
  }
#define _OPA_LT_OP4(x, a, b, c, d)                                             \
  {                                                                            \
    if (a != x.a) return a < x.a;                                              \
    _OPA_LT_OP3(x, b, c, d);                                                   \
  }
#define _OPA_LT_OP5(x, a, b, c, d, e)                                          \
  {                                                                            \
    if (a != x.a) return a < x.a;                                              \
    _OPA_LT_OP4(x, b, c, d, e);                                                \
  }
#define _OPA_LT_OP6(x, a, b, c, d, e, f)                                       \
  {                                                                            \
    if (a != x.a) return a < x.a;                                              \
    _OPA_LT_OP5(x, b, c, d, e, f);                                             \
  }

#define RAW_OPA_LT_OP(x, ...)                                                  \
  OPA_DISPATCH(_OPA_LT_OP, __VA_ARGS__)(x, __VA_ARGS__)
#define OPA_LT_OP(x, ...) RAW_OPA_LT_OP(x, __VA_ARGS__)

#define _OPA_EQ_OP1(x, a)                                                      \
  { return a == x.a; }
#define _OPA_EQ_OP2(x, a, b)                                                   \
  {                                                                            \
    if (a != x.a) return false;                                                \
    _OPA_EQ_OP1(x, b);                                                         \
  }
#define _OPA_EQ_OP3(x, a, b, c)                                                \
  {                                                                            \
    if (a != x.a) return false;                                                \
    _OPA_EQ_OP2(x, b, c);                                                      \
  }
#define _OPA_EQ_OP4(x, a, b, c, d)                                             \
  {                                                                            \
    if (a != x.a) return false;                                                \
    _OPA_EQ_OP3(x, b, c, d);                                                   \
  }
#define _OPA_EQ_OP5(x, a, b, c, d, e)                                          \
  {                                                                            \
    if (a != x.a) return false;                                                \
    _OPA_EQ_OP4(x, b, c, d, e);                                                \
  }
#define _OPA_EQ_OP6(x, a, b, c, d, e, f)                                       \
  {                                                                            \
    if (a != x.a) return false;                                                \
    _OPA_EQ_OP5(x, b, c, d, e, f);                                             \
  }

#define RAW_OPA_EQ_OP(x, ...)                                                  \
  OPA_DISPATCH(_OPA_EQ_OP, __VA_ARGS__)(x, __VA_ARGS__)
#define OPA_EQ_OP(x, ...) RAW_OPA_EQ_OP(x, __VA_ARGS__)

#define _OPA_COPY_OP1(x, a)                                                    \
  { a = x.a; }
#define _OPA_COPY_OP2(x, a, b)                                                 \
  {                                                                            \
    _OPA_COPY_OP1(x, a);                                                       \
    _OPA_COPY_OP1(x, b);                                                       \
  }
#define _OPA_COPY_OP3(x, a, b, c)                                              \
  {                                                                            \
    _OPA_COPY_OP1(x, a);                                                       \
    _OPA_COPY_OP2(x, b, c);                                                    \
  }
#define _OPA_COPY_OP4(x, a, b, c, d)                                           \
  {                                                                            \
    _OPA_COPY_OP1(x, a);                                                       \
    _OPA_COPY_OP3(x, b, c, d);                                                 \
  }
#define _OPA_COPY_OP5(x, a, b, c, d, e)                                        \
  {                                                                            \
    _OPA_COPY_OP1(x, a);                                                       \
    _OPA_COPY_OP4(x, b, c, d, e);                                              \
  }
#define _OPA_COPY_OP6(x, a, b, c, d, e, f)                                     \
  {                                                                            \
    _OPA_COPY_OP1(x, a);                                                       \
    _OPA_COPY_OP5(x, b, c, d, e, f);                                           \
  }

#define RAW_OPA_COPY_OP(x, ...)                                                \
  OPA_DISPATCH(_OPA_COPY_OP, __VA_ARGS__)(x, __VA_ARGS__)
#define OPA_COPY_OP(x, ...) RAW_OPA_COPY_OP(x, __VA_ARGS__)

#define _TYPE_AND_NAME_TO_ARGS_OP2(t1, n1) t1 n1
#define _TYPE_AND_NAME_TO_ARGS_OP4(t1, n1, t2, n2)                             \
  _TYPE_AND_NAME_TO_ARGS_OP2(t1, n1), _TYPE_AND_NAME_TO_ARGS_OP2(t2, n2)
#define _TYPE_AND_NAME_TO_ARGS_OP6(t1, n1, t2, n2, t3, n3)                     \
  _TYPE_AND_NAME_TO_ARGS_OP2(t1, n1), _TYPE_AND_NAME_TO_ARGS_OP4(t2, n2, t3, n3)
#define _TYPE_AND_NAME_TO_ARGS_OP8(t1, n1, t2, n2, t3, n3, t4, n4)             \
  _TYPE_AND_NAME_TO_ARGS_OP2(t1, n1)                                           \
  , _TYPE_AND_NAME_TO_ARGS_OP4(t2, n2, t3, n3, t4, n4)

#define RAW_TYPE_AND_NAME_TO_ARGS_OP(...)                                      \
  OPA_DISPATCH(_TYPE_AND_NAME_TO_ARGS_OP, __VA_ARGS__)(__VA_ARGS__)
#define TYPE_AND_NAME_TO_ARGS(...) RAW_TYPE_AND_NAME_TO_ARGS_OP(__VA_ARGS__)

#define _OPA_CONSTRUCTOR_OP2(t1, n1) n1(n1)
#define _OPA_CONSTRUCTOR_OP4(t1, n1, t2, n2)                                   \
  _OPA_CONSTRUCTOR_OP2(t1, n1), _OPA_CONSTRUCTOR_OP2(t2, n2)
#define _OPA_CONSTRUCTOR_OP6(t1, n1, t2, n2, t3, n3)                           \
  _OPA_CONSTRUCTOR_OP2(t1, n1), _OPA_CONSTRUCTOR_OP4(t2, n2, t3, n3)
#define _OPA_CONSTRUCTOR_OP8(t1, n1, t2, n2, t3, n3, t4, n4)                   \
  _OPA_CONSTRUCTOR_OP2(t1, n1), _OPA_CONSTRUCTOR_OP4(t2, n2, t3, n3, t4, n4)

#define RAW_OPA_CONSTRUCTOR_OP(...)                                            \
  OPA_DISPATCH(_OPA_CONSTRUCTOR_OP, __VA_ARGS__)(__VA_ARGS__)

#define OPA_CONSTRUCTOR_OP(cl, body, ...)                                      \
  cl(TYPE_AND_NAME_TO_ARGS(__VA_ARGS__))                                       \
      : RAW_OPA_CONSTRUCTOR_OP(__VA_ARGS__) {                                  \
    body;                                                                      \
  }

#define _ARGS_JOIN1(x, a) a
#define _ARGS_JOIN2(x, a, b) a x b
#define _ARGS_JOIN3(x, a, b, c) a x _ARGS_JOIN2(x, b, c)
#define _ARGS_JOIN4(x, a, b, c, d) a x _ARGS_JOIN3(x, b, c, d)
#define _ARGS_JOIN5(x, a, b, c, d, e) a x _ARGS_JOIN4(x, b, c, d, e)
#define ARGS_JOIN(x, ...) OPA_DISPATCH(_ARGS_JOIN, __VA_ARGS__)(x, __VA_ARGS__)

#define OPA_TGEN_DECL_BASE(typ, ...)                                           \
  OPA_NAMESPACE(google, protobuf) class Any;                                   \
  OPA_NAMESPACE_END(google, protobuf) OPA_NAMESPACE(__VA_ARGS__) class typ;    \
  OPA_NAMESPACE_END(__VA_ARGS__);                                              \
  OPA_NAMESPACE(opa, utils);                                                   \
  static void tgen_load(ARGS_JOIN(::, __VA_ARGS__)::typ &,                     \
                        const google::protobuf::Any &);                        \
  static void tgen_store(const ARGS_JOIN(::, __VA_ARGS__)::typ &,              \
                         google::protobuf::Any &);                             \
  OPA_NAMESPACE_END(opa, utils);

#if OPA_CPP
#include <cstdlib>

#include <algorithm>
#include <bitset>
#include <cassert>
#include <climits>
#include <cmath>
#include <complex>
#include <csignal>
#include <cstddef>
#include <cstdio>
#include <ctime>
#include <deque>
#include <errno.h>
#include <fcntl.h>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <queue>
#include <set>
#include <sstream>
#include <stack>
#include <string>
#else
#include <stdlib.h>
#endif

#include <stdarg.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#if OPA_CPP14 == 1
#include <array>
#include <chrono>
#include <condition_variable>
#include <functional>
#include <numeric>
#include <random>
#include <thread>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#endif

#if OPA_CPP14 == 1
template <typename T> struct opa_static_assert_impl_helper : std::false_type {};
#define OPA_STATIC_ASSERT_IMPL(T)                                              \
  static_assert(opa_static_assert_impl_helper<T>::value,                       \
                "this function has to be implemented for desired type");
#else
#define OPA_STATIC_ASSERT_IMPL(T)
#endif

#if OPA_CPP
#if !defined(SWIG)
namespace std {
template <class T1, class T2>
std::ostream &operator<<(std::ostream &os, const std::pair<T1, T2> &a) {
  os << "(";
  os << a.ST << ", " << a.ND;
  os << ")";
  return os;
}

#if OPA_CPP14 == 1
template <class T, class L>
std::ostream &operator<<(std::ostream &os, const std::unordered_set<T, L> &a) {
  os << "(";
  for (auto &x : a) os << x << ",";
  os << ")";
  return os;
}

template <class A, class B, class L>
std::ostream &operator<<(std::ostream &os,
                         const std::unordered_map<A, B, L> &a) {
  os << "(";
  for (auto &x : a) os << x << ",";
  os << ")";
  return os;
}

template <std::size_t N, typename... Args> class Opa_os_tuple_pipe {
public:
  static std::ostream &os_pipe_tuple(std::ostream &os,
                                     const std::tuple<Args...> &a) {
    Opa_os_tuple_pipe<N - 1, Args...>::os_pipe_tuple(os, a);
    os << "," << std::get<N>(a);
    return os;
  }
};

template <typename... Args> class Opa_os_tuple_pipe<0, Args...> {
public:
  static std::ostream &os_pipe_tuple(std::ostream &os,
                                     const std::tuple<Args...> &a) {
    os << std::get<0>(a);
    return os;
  }
};

template <typename... Args>
std::ostream &operator<<(std::ostream &os, const std::tuple<Args...> &a) {
  os << "<";
  Opa_os_tuple_pipe<sizeof...(Args) - 1, Args...>::os_pipe_tuple(os, a);
  os << ">";
  return os;
}

template <class T>
std::ostream &operator<<(std::ostream &os, const std::initializer_list<T> &a) {
  os << "(";
  for (auto &x : a) os << x << ",";
  os << ")";
  return os;
}

template <class T, size_t N>
std::ostream &operator<<(std::ostream &os, const std::array<T, N> &a) {
  os << "(";
  for (auto &x : a) os << x << ",";
  os << ")";
  return os;
}

#endif

template <class T, class L>
std::ostream &operator<<(std::ostream &os, const std::set<T, L> &a) {
  os << "(";
  for (auto &x : a) os << x << ",";
  os << ")";
  return os;
}

template <class A, class B, class L>
std::ostream &operator<<(std::ostream &os, const std::map<A, B, L> &a) {
  os << "(";
  for (auto &x : a) os << x << ",";
  os << ")";
  return os;
}

template <class T>
std::ostream &operator<<(std::ostream &os, const std::deque<T> &a) {
  os << "(";
  REP (i, a.size())
    os << a[i] << (i == a.size() - 1 ? "" : ",");
  os << ")";
  return os;
}

template <class T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &a) {
  os << "(";
  REP (i, a.size())
    os << a[i] << (i == a.size() - 1 ? "" : ",");
  os << ")";
  return os;
}

template <class A, class B>
std::pair<A, B> operator+(const std::pair<A, B> &a, const std::pair<A, B> &b) {
  return MP(a.ST + b.ST, a.ND + b.ND);
}

template <class A, class B>
std::pair<A, B> operator-(const std::pair<A, B> &a, const std::pair<A, B> &b) {
  return MP(a.ST - b.ST, a.ND - b.ND);
}
} // namespace std
#endif

template <class T> std::string to_c_code(const std::vector<T> &a) {
  std::ostringstream os;
  os << "{";
  REP (i, a.size()) {
    os << a[i] << ",";
    ;
  }
  os << "}";
  return os.str();
}

typedef std::vector<int> vi;
typedef std::pair<int, int> pii;

template <class T> static void checkmin(T &a, T b) {
  if (b < a) a = b;
}
template <class T> static void checkmax(T &a, T b) {
  if (b > a) a = b;
}
template <class T> static void out(T t[], int n) {
  REP (i, n)
    std::cout << std::hex << t[i] << " ";
  std::cout << std::endl;
}
template <class T> static void out(std::vector<T> t, int n = -1) {
  for (int i = 0; i < (n == -1 ? t.size() : n); ++i) std::cout << t[i] << " ";
  std::cout << std::endl;
}
static inline u64 count_bit_slow(u64 n) {
  return (n == 0) ? 0 : 1 + count_bit_slow(n & (n - 1));
}
static inline u64 low_bit(u64 n) { return (n ^ n - 1) & n; }
static inline u64 ctz(u64 n) { return (n == 0 ? -1 : ctz(n >> 1) + 1); }
static inline s8 log2_high_bit(u64 n) {
  s8 res = -1;
  while (n) ++res, n >>= 1;
  return res;
}

template <class T> std::string bitstring(const T &a, int n) {
  std::string res;
  res.resize(n);
  REP (i, n)
    res[i] = '0' + (a >> i & 1);
  return res;
}

static s64 toInt(std::string s) {
  s64 a;
  std::istringstream(s) >> a;
  return a;
}

template <class T> std::string toStr(const T &a) {
  std::ostringstream os;
  os << a;
  return os.str();
}

#if OPA_PIN == 0
template <class T> std::shared_ptr<const T> make_dummy_sptr(const T *a) {
  return std::shared_ptr<const T>(a, [](const T *) {});
}

template <class T> std::shared_ptr<T> make_dummy_sptr(T *a) {
  return std::shared_ptr<T>(a, [](T *) {});
}

template <class T> std::unique_ptr<T> make_dummy_uptr(T *a) {
  return std::unique_ptr<T>(a, [](T *) {});
}

template <class T> std::shared_ptr<T> make_sptr(T *a) {
  return std::shared_ptr<T>(a);
}
#endif

template <typename V, typename K, typename MP>
V FindWithDefault(const MP &mp, const K &key, const V &default_v = V()) {
  auto it = mp.find(key);
  if (it == mp.end()) return default_v;
  return it->second;
}

template <class T>
using SmallPQ = std::priority_queue<T, std::vector<T>, std::greater<T> >;

#endif

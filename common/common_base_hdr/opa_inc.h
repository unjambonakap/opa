#ifndef _H_OPA_INC
#define _H_OPA_INC


#if defined(__cplusplus)
#define __BEGIN_DECLS extern "C" {
#define __END_DECLS }
#define OPA_CPP 1
#define OPA_C 0

#if __cplusplus > 199711L
#define OPA_CPP11 1
#else
#define OPA_CPP11 0
#endif

#if __cplusplus > 201103L
#define OPA_CPP14 1
#else
#define OPA_CPP14 0
#endif

#else
#define __BEGIN_DECLS
#define __END_DECLS
#define OPA_CPP 0
#define OPA_C 1
#define OPA_CPP11 0

#define false 0
#define true 1
#endif

#ifdef __i386__
#define OPA_32 1
#else
#define OPA_32 0
#endif

#ifdef __PIN__
#define OPA_PIN 1
#else
#define OPA_PIN 0
#endif

#define OPA_64 (!OPA_32)

#ifdef __ASSEMBLY__
#define OPA_ASSEMBLY 1
#else
#define OPA_ASSEMBLY 0
#endif

#define LIST_COUNT(lst) (sizeof(lst) / sizeof((lst)[0]))
#define REP(i, n) for (int i = 0; i < int(n); ++i)
#define REP64(i, n) for (s64 i = 0; i < s64(n); ++i)
#define REPV(i, n) for (int i = (n)-1; (int)i >= 0; --i)
#define FOR(i, a, b) for (int i = (int)(a); i < (int)(b); ++i)
#define FORV(i, a, b) for (int i = (int)(b-1); i >= (int)(a); --i)
#define FOR64(i, a, b) for (s64 i = (s64)(a); i < (s64)(b); ++i)
#ifndef MIN
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#endif

#define FE(i, t)                                                               \
  for (__typeof((t).begin()) i = (t).begin(); i != (t).end(); ++i)
#define FEV(i, t)                                                              \
  for (__typeof((t).rbegin()) i = (t).rbegin(); i != (t).rend(); ++i)

#define _two(x) (1LL << (x))
#define ALL(a) (a).begin(), (a).end()

#define pb push_back
#define ST first
#define ND second
#define MP(x, y) std::make_pair(x, y)

#define INIT_FUNC_ATTR __attribute__((constructor))
#define OPA_XSTR(x) OPA_STR(x)
#define OPA_STR(x) #x

#define VA_NUM_ARGS_IMPL(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10, N, ...) N
#define VA_NUM_ARGS(...)                                                       \
  VA_NUM_ARGS_IMPL(__VA_ARGS__, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1)
#define OPA_DISPATCH(func, ...) OPA_DISPATCH_(func, VA_NUM_ARGS(__VA_ARGS__))
#define OPA_DISPATCH_(func, nargs) OPA_DISPATCH__(func, nargs)
#define OPA_DISPATCH__(func, nargs) func##nargs

#if OPA_ASSEMBLY == 0

#include <stdint.h>

#ifndef OPA_SKIP_LONG
typedef long long ll;
typedef unsigned long long ull;
typedef int64_t s64;
typedef uint64_t u64;
#endif

typedef unsigned int uint;
typedef int8_t s8;
typedef uint8_t u8;
typedef int16_t s16;
typedef uint16_t u16;
typedef int32_t s32;
typedef uint32_t u32;

typedef float f32;
typedef double f64;

#define OPA_RESOURCE_DECL(name, file) extern const char name[];

#else

#define END(f) .size f, .- f;

#define ENTRY(f)                                                               \
  .text;                                                                       \
  .globl f;                                                                    \
  .align 16;                                                                   \
  .type f, @function;                                                          \
  f:

#define OPA_RESOURCE_DECL(name, file)                                          \
  .section ".rodata";                                                          \
  .globl name;                                                                 \
  .type name, STT_OBJECT;                                                      \
  name:;                                                                       \
  .incbin file;                                                                \
  .byte 0;                                                                     \
  .size name, .- name;

#endif

#ifdef OPA_FONA_ARCH
#define LINE_TERMINAISON "\n"
#else
#define LINE_TERMINAISON "\n"
#endif

#define OPA_TRACE_ARGS(fmt, ...)
#define TRACE_STR(prefix, suffix) ((prefix ">" __FILE__ ":" OPA_XSTR(__LINE__) " >> " suffix LINE_TERMINAISON))
#define TRACE_STR1(msg) TRACE_STR("", msg)
#define TRACE_STR0() TRACE_STR("", "")

#define KTRACE(fmt, ...)                                                       \
  printk("%s:%s:%d, msg=" fmt "\n", __FILE__, __func__, __LINE__, __VA_ARGS__);
#define KTRACE0() printk("%s:%s:%d\n", __FILE__, __func__, __LINE__);
#define KTRACEX(fmt)                                                           \
  printk("%s:%s:%d, msg=" fmt "\n", __FILE__, __func__, __LINE__);
#include "opa_config.h"

/* Example of a dispatcher
#define OPA_ACTION(typ, func, data_t)                                          \
  case typ:                                                                    \
    func((data_t)msg->data);                                                   \
    break;

  switch (msg->type) {
    OPA_ACTION(opa_RprocControllerType_Boot, rproc_boot_proxy, const void *)
    OPA_ACTION(opa_RprocControllerType_Shutdown, rproc_shutdown_proxy,
               const void *)
  default:
    KTRACE("Invalid type %d", msg->type);
    break;
  }
#undef OPA_ACTION
*/
#endif

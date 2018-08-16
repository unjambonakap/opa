#pragma once
#include <opa_common.h>

OPA_NAMESPACE(opa, crypto, la)

template <typename T> void do_walsh(T *tb, int n) {
  s64 N = 1ull << n;
  REP64(i, n) {
    u64 step = 1ull << i;
    for (s64 j = 0; j < N; j += 2 * step) {
      REP (k, step) {
        T a = tb[j + k] + tb[j + k + step];
        T b = tb[j + k] - tb[j + k + step];
        tb[j + k] = a;
        tb[j + k + step] = b;
      }
    }
  }
}

template <typename T> void do_iwalsh(T *tb, int n) {
  do_walsh(tb, n);
  int N = 1 << n;
  REP (i, n) {
    //OPA_CHECK(tb[i] % N == 0, "bad ", i, tb[i], N);
    tb[i] >>= n;
  }
}

OPA_NAMESPACE_END(opa, crypto, la)

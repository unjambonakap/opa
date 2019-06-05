#pragma once

#include <cstdio>
#include <complex>

inline void test1(const std::complex<float> &tb) {
  printf("got %f %f\n", std::real(tb), std::imag(tb));

}

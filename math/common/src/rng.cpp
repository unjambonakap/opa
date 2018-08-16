#include <opa/math/common/rng.h>

OPA_NAMESPACE_DECL3(opa, math, common)

thread_local std::mt19937 rng(0);
thread_local std::mt19937_64 rng64(0);
void reseed_rng() {
  FILE *f = fopen("/dev/urandom", "rb");
  u64 data;
  OPA_CHECK0(fread(&data, 1, sizeof(data), f) == sizeof(data));
  fclose(f);
  rng.seed(data);
  rng64.seed(data);
}

Rng::~Rng() { fclose(m_urandom_file); }

Rng::Rng() { m_urandom_file = fopen("/dev/urandom", "rb"); }

std::string Rng::bytes(int n) const {
  std::string res(n, 0);
  OPA_CHECK0(fread(&res[0], 1, n, m_urandom_file) == n);
  return res;
}

std::string Rng2::bytes(int n) const {
  std::string res(n, 0);
  REP (i, n)
    res[i] = rng();
  return res;
}

OPA_NAMESPACE_DECL3_END

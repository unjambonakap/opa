#include <opa/math/common/fast_gf2.h>
#include <opa/math/common/rng.h>
#include <opa/math/common/Types.h>

OPA_NAMESPACE(opa, math, common)
GF_q<u32> gf128;
Gf128_t gcm_poly;
u8 cnt_bit[maxb];

void init_bitstuff() {
  cnt_bit[0] = 0;
  FOR(i, 1, maxb) { cnt_bit[i] = 1 ^ cnt_bit[i & i - 1]; }
}


void init_fast() {
  init_bitstuff();

  gcm_poly = PR_GF2.xpw(128) + PR_GF2.xpw(7) + PR_GF2.xpw(2) + PR_GF2.x() +
             PR_GF2.constant(1);
  gf128.init(&GF2, gcm_poly);
}

BitVec BitVec::rand(int n) {
  BitVec res(n);
  REP (i, n)
    res.set(i, rng() & 1);
  return res;
}

BitVec BitVec::FromBytes(const RangeCoverage &range, opa::stolen::StringRef s) {
  BitVec res(range.all().size());
  for (auto i : range.all()) res.set(i, s[i / 8] >> (i & 7) & 1);
  return res;
}

OPA_NAMESPACE_END(opa, math, common)

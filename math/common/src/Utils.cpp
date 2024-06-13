#include "common/Utils.h"

#include "common/Field.h"
#include "common/PolyRing.h"
#include <opa/math/common/Types.h>
#include <opa/math/common/fast_gf2.h>
#include <opa/math/common/float.h>
#include <opa/or/best_first_search.h>

using namespace std;
const int n_small_prime_test = 40;
const int n_rm_tries = 25;
const int factor_large_maxB = 2e9;
const int kSmallPrimeBound = 1e3;

DEFINE_int32(math_seed, 0, "");
DEFINE_bool(init_primes, true, "");
DEFINE_int32(primedb_maxv, 1e4, "");

OPA_NAMESPACE_DECL3(opa, math, common)

u32 *isc = nullptr;
inline bool get_isc(u32 x) { return isc[x >> 5] >> (x & 0x1f) & 1; }
inline void set_isc(u32 x) { isc[x >> 5] |= 1u << (x & 0x1f); }
int bitcount_tb[1 << BITCOUNT_BLK];

std::vector<u32> pl;

void init_primes(int bound) {
  int sz = bound / 32 + 1;
  isc = new u32[sz];
  memset(isc, 0, sizeof(isc[0]) * (sz));
  set_isc(0);
  set_isc(1);
  u32 ub = sqrt(bound) + 10;
  if (pl.size() > 0) return;
  for (u32 i = 2; i < bound; ++i)
    if (!get_isc(i)) {
      pl.push_back(i);
      if (i < ub)
        for (u32 j = i * i; j < bound; j += i) set_isc(j);
    }
}

void initMathCommon(int seed) {
  bitcount_tb[0] = 0;
  FOR (i, 1, 1 << BITCOUNT_BLK) bitcount_tb[i] = bitcount_tb[i >> 1] + (i & 1);

  if (seed == -1) seed = time(0);
  rng.seed(seed);
  rng64.seed(seed);
  bignum_init(seed);
  Float_init(seed);

  pl.clear();

  if (FLAGS_init_primes) init_primes(FLAGS_primedb_maxv);

  init_math_types();
  init_fast();
}

void init_math() { initMathCommon(FLAGS_math_seed); }

OPA_REGISTER_INIT(init_math, init_math);

bool isPrimeDB(u32 a) {
  assert(a < FLAGS_primedb_maxv);
  return !get_isc(a);
}

S64Factors factor_small(bignum a, bignum &rem) {
  assert(a > 0);

  S64Factors res;
  REP (i, pl.size()) {
    if (a == 1) break;
    int cnt = 0;
    while (a % pl[i] == 0) a.sdiv(pl[i]), ++cnt;

    if (cnt) res.push_back(MP((s64)pl[i], cnt));
  }
  rem = 1;
  if (a > 1 && a < bignum::fromu64(FLAGS_primedb_maxv) * FLAGS_primedb_maxv)
    res.push_back(MP(a.getu64(), 1));
  else
    rem = a;
  return res;
}

u32 nextPrimeSmall(u32 a) {
  for (; a < FLAGS_primedb_maxv; ++a)
    if (isPrimeDB(a)) return a;
  return 0;
}

bool testPrimeFermat(const bignum &n, int ntry) {
  bignum o = n - 1;
  REP (i, ntry) {
    bignum x = (n - 2).rand() + 2;
    x = x.faste(o, n);
    if (x != 1) return false;
  }
  return true;
}

bool testPrimeRM(const bignum &n, int ntry) {
  bignum o;
  int e;
  o = n - 1;
  while (o.get_bit(0) == 0) {
    ++e;
    o.srshift(1);
  }
  bignum minus1 = n - 1;

  REP (i, ntry) {
    bignum x = (n - 2).rand() + 2;
    x = x.powm(o, n);
    if (x == 1) return true;
    REP (i, e) {
      x = x * x % n;
      if (x == minus1) return true;
    }
  }
  return false;
}

bool testPrime(const bignum &n) {
  if (n < FLAGS_primedb_maxv) return isPrimeDB(n.getu32());
  init_primes(kSmallPrimeBound);
  int bound = std::min<int>(n_small_prime_test, pl.size());

  REP (i, bound)
    if (n % pl[i] == 0) return false;
  if (!testPrimeFermat(n, n_rm_tries)) return false;
  return testPrimeRM(n, n_rm_tries);
}

bool pollardP1(const bignum &n, int B, bignum &out) {
  bignum a = 2;
  FOR (i, 2, B + 1) {
    a = a.powm(i, n);
    out = n.gcd(a - 1);
    if (out > 1 && out < n) return true;
  }
  return false;
}

BGFactors factor_large(const bignum &a) {
  BGFactors res;
  map<bignum, int> factors;

  bignum rem;
  S64Factors small_factors = factor_small(a, rem);
  for (const auto &x : small_factors) factors[bignum::fromu64(x.ST)] += x.ND;

  queue<bignum> q;
  if (rem > 1) q.push(rem);

  while (q.size()) {
    bignum u = q.front();
    q.pop();

    if (testPrime(u)) {
      factors[u] += 1;
      continue;
    }

    bignum d;
    if (pollardP1(u, factor_large_maxB, d)) {
      q.push(d);
      q.push(u / d);
      continue;
    }
    cout << "FAIL ON >>" << u << endl;
    assert(0);
  }
  for (auto x : factors) res.push_back(x);
  return res;
}

u64 u32_faste(u32 a, u32 p, u32 mod) {
  u64 x = 1;
  for (; p; p >>= 1, a = (u64)a * a % mod)
    if (p & 1) x = (u64)x * a % mod;
  return x;
}

u32 u32_faste(u32 a, u32 p) {
  u32 x = 1;
  for (; p; p >>= 1, a = a * a)
    if (p & 1) x = x * a;
  return x;
}

std::vector<u32> genRand(u32 n, u32 sz) {

  std::vector<u32> tb(sz);
  for (u32 i = 0; i < sz; ++i) tb[i] = rng() % n;
  return tb;
}

u32 u32_gcd(u32 a, u32 b) {
  if (a < b) std::swap(a, b);
  return b ? u32_gcd(b, a % b) : a;
}

u32 u32_egcd(u32 a, u32 b, int &u, int &v) {
  if (a < b) return u32_egcd(b, a, v, u);
  return _u32_egcd(a, b, 1, 0, 0, 1, u, v);
}

u32 _u32_egcd(u32 a, u32 b, int ua, int va, int ub, int vb, int &u, int &v) {
  if (b == 0) {
    u = ua;
    v = va;
    return a;
  }
  u32 m = a / b;
  return _u32_egcd(b, a % b, ub, vb, ua - m * ub, va - m * vb, u, v);
}
u32 u32_lcm(u32 a, u32 b) { return a / u32_gcd(a, b) * b; }

bignum gen_prime(const bignum &n) {
  while (true) {
    bignum x = n.rand();
    if (testPrime(x)) return x;
  }
}

bignum next_prime(const bignum &n) {
  bignum v = n + 1;
  for (; !testPrime(v); v += 1)
    ;
  return v;
}

bignum crt_coprime(const std::vector<std::pair<bignum, bignum> > &lst, const bignum &mod) {

  bignum res = 0;
  bignum n = 1;
  const bignum *pmod = mod == -1 ? &n : &mod;
  for (auto &x : lst) n *= x.ND;

  for (auto &x : lst) {
    bignum nx = n / x.ND;
    res = (res + x.ST * nx * (nx % x.ND).inv(x.ND)) % *pmod;
  }
  return res;
}

void CRT::add(const bignum &q, int pw, const bignum &v) {
  if (m_mp.count(q)) {
    auto &x = m_mp[q];
    if (pw <= x.ST) {
      bignum qpw = q.pow(pw);
      OPA_CHECK0(x.ND % qpw == v);
      return;
    }
    bignum qpw = q.pow(x.ST);
    OPA_CHECK0(v % qpw == x.ND);
  }
  m_mp[q] = MP(pw, v);
}

bignum CRT::solve(const bignum &mod) const {
  std::vector<std::pair<bignum, bignum> > tb;
  for (auto &x : m_mp) tb.emplace_back(x.ND.ND, x.ST.pow(x.ND.ST));
  return crt_coprime(tb, mod);
}

bignum CRT::get_bound() const {
  bignum res = 1;
  for (auto &x : m_mp) res = res * x.ST.pow(x.ND.ST);
  return res;
}

OPA_NAMESPACE_DECL3_END

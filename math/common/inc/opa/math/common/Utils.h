#pragma once

#include <opa/math/common/base.h>
#include <opa/math/common/bignum.h>
#include <opa/math/common/rng.h>

OPA_NM_MATH_COMMON

typedef std::vector<std::pair<s64, int> > S64Factors;
typedef std::vector<std::pair<bignum, int> > BGFactors;

inline BGFactors to_bgfactors(const S64Factors &factors) {
  BGFactors res;
  for (auto &factor : factors) res.emplace_back(factor.first, factor.second);
  return res;
}

inline int get_permutation_signature(const std::vector<int> &tb) {
  int m = 0;
  std::vector<int> seen(tb.size(), 0);
  REP (i, tb.size()) {
    if (seen[i]) continue;
    int cnt = 0;
    for (int x = i; !seen[x]; seen[x] = 1, x = tb[x], ++cnt)
      ;
    m += cnt - 1;
  }
  return (m & 1) ? -1 : 1;
}

extern std::vector<u32> pl;
bignum gen_prime(const bignum &n);
bool testPrime(const bignum &n);
void initMathCommon(int seed = -1);
bool isPrimeDB(u32 a);
u32 nextPrimeSmall(u32 a);
inline static u64 u32_faste(u32 a, u32 p, u32 mod) {
  u64 x = 1;
  for (; p; p >>= 1, a = (u64)a * a % mod)
    if (p & 1) x = (u64)x * a % mod;
  return x;
}

inline static u32 u32_faste(u32 a, u32 p) {
  u32 x = 1;
  for (; p; p >>= 1, a = a * a)
    if (p & 1) x = x * a;
  return x;
}
bignum next_prime(const bignum &n);

// <remained, prime power>
bignum crt_coprime(const std::vector<std::pair<bignum, bignum> > &lst, const bignum &mod = -1);

class CRT {
public:
  void reset() { m_mp.clear(); }

  void add(const bignum &q, int pw, const bignum &v);
  bignum solve(const bignum &mod = -1) const;
  bignum get_bound() const;

  const std::map<bignum, std::pair<int, bignum> > &mp() { return m_mp; };

private:
  std::map<bignum, std::pair<int, bignum> > m_mp;
};

class CRT_Batched {
public:
  void init(const std::vector<bignum> &primes) {
    n = 1;
    for (auto &x : primes) n *= x;

    for (auto &x : primes) {
      bignum nx = n / x;
      coeffs.push_back(nx * (nx % x).inv(x));
    }
  }

  bignum solve(const std::vector<bignum> &vec) const {
    bignum res = 0;
    OPA_CHECK0(vec.size() == coeffs.size());
    REP (i, vec.size()) {
      res += vec[i] * coeffs[i];
    }
    return res % n;
  }

  bignum n;
  std::vector<bignum> coeffs;
};

template <class T> T t_faste_u32(T a, u32 p) {
  T x = T::identity();
  for (; p; p >>= 1, a = a * a)
    if (p & 1) x = x * a;

  return x;
}

template <class T> T t_faste(T a, bignum p) {
  T x = a.identity();
  for (; p != 0; p.srshift(1), a = a * a)
    if (p.get_bit(0)) x = x * a;
  return x;
}

template <class T> T t_faste(T a, bignum p, const bignum &mod) {
  T x = T::identity();
  for (; p != 0; p.srshift(1), a = a * a % mod)
    if (p.get_bit(0)) x = x * a % mod;
  return x;
}

// return pair <factor,power>
S64Factors factor_small(bignum a, bignum &rem);
BGFactors factor_large(const bignum &a);
inline S64Factors factor_small(u32 a) {
  bignum rem;
  S64Factors res;
  res = factor_small(a, rem);
  OPA_CHECK0(rem == 1);
  return res;
}

inline u32 u32_gcd(u32 a, u32 b) {
  if (a < b) std::swap(a, b);
  return b ? u32_gcd(b, a % b) : a;
}

inline u32 _u32_egcd(u32 a, u32 b, s32 ua, s32 va, s32 ub, s32 vb, s32 &u, s32 &v) {
  if (b == 0) {
    u = ua;
    v = va;
    return a;
  }
  u32 m = a / b;
  return _u32_egcd(b, a - m * b, ub, vb, ua - m * ub, va - m * vb, u, v);
}

inline u32 u32_egcd(u32 a, u32 b, s32 &u, s32 &v) {
  if (a < b) return u32_egcd(b, a, v, u);
  return _u32_egcd(a, b, 1, 0, 0, 1, u, v);
}

inline u32 u32inv_egcd(u32 a, u32 p) {
  s32 u, v;
  u32_egcd(a, p, u, v);
  if (u < 0) u += p;
  if (u >= p) u -= p;
  return u;
}

constexpr inline u64 _u64_egcd(u64 a, u64 b, s64 ua, s64 va, s64 ub, s64 vb, s64 &u, s64 &v) {
  if (b == 0) {
    u = ua;
    v = va;
    return a;
  }
  u64 m = a / b;
  return _u64_egcd(b, a - b * m, ub, vb, ua - m * ub, va - m * vb, u, v);
}

constexpr inline u64 u64_egcd(u64 a, u64 b, s64 &u, s64 &v) {
  if (a < b) return u64_egcd(b, a, v, u);
  return _u64_egcd(a, b, 1, 0, 0, 1, u, v);
}

constexpr inline u64 u64inv_egcd(u64 a, u64 p) {
  s64 u, v;
  u64_egcd(a, p, u, v);
  if (u < 0) u += p;
  if (u >= p) u -= p;
  return u;
}

inline u32 u32_lcm(u32 a, u32 b) { return a / u32_gcd(a, b) * b; }

std::vector<u32> genRand(u32 n, u32 sz);

constexpr int BITCOUNT_BLK = 16;
extern int bitcount_tb[1 << BITCOUNT_BLK];
inline int count_bit(u16 a) { return bitcount_tb[a]; }
inline int count_bit(u32 a) { return count_bit(u16(a)) + count_bit(u16(a >> 16)); }
inline int count_bit(u64 a) { return count_bit(u32(a)) + count_bit(u32(a >> 32)); }

OPA_NM_MATH_COMMON_END

#pragma once

#include <opa/crypto/lfsr.h>

OPA_NM_CRYPTO

class LFSR_GF2_small : public LFSR_Base<u32> {
public:
  LFSR_GF2_small() {}
  LFSR_GF2_small(u64 poly, int n, u64 state) { init(poly, n, state); }

  void init(u64 poly, int n, u64 state) {
    this->poly = poly;
    this->n = n;
    this->mask=  (1ull<<n) - 1;
    this->state = state;
  }

  static u64 poly_to_state(const OPA_MATH::Poly<u32> &p) {
    u64 res = 0;
    REP (i, p.size())
      res |= u64(p.get_safe(i)) << i;
    return res;
  }

  static LFSR_GF2_small from_lfsr(const LFSR<u32> &lfsr) {
    int n = lfsr.get_poly().deg();
    OPA_CHECK0(n < 64);
    u64 poly = poly = poly_to_state(lfsr.get_poly());
    u64 state = poly_to_state(lfsr.get_state());
    ;
    return LFSR_GF2_small(poly, n, state);
  }

  void set_rand() { set_state(opa::math::common::rng()); }

  static std::shared_ptr<LFSR_GF2_small> rand(int n) {
    OPA_MATH::Poly<u32> poly;
    poly = OPA_MATH::find_primitive_poly(OPA_MATH::GF2, n);
    OPA_MATH::PolyRing<u32> pr(&OPA_MATH::GF2);
    OPA_MATH::Poly<u32> init_state = pr.rand(n - 1);

    return std::make_shared<LFSR_GF2_small>(poly_to_state(poly), n,
                                            poly_to_state(init_state));
  }

  void set_state(u64 state) { this->state = state & ((1ull << this->n) - 1); }

  u32 get_next() {
    u64 r = state >> (n - 1) & 1;
    state <<= 1;
    state ^= (-r) & poly;
    return r;
  }

  u32 get_next_non_galois() {
    u64 r = state >> (n - 1) & 1;
    u32 nb = 1;
    u32 masked = state & poly;
    REP(i,n) nb ^= (masked >> i) & 1;
    state  = (state << 1 | nb) & mask;
    return r;
  }
  int size() const { return n; }

public:
  u64 mask;
  u64 state;
  u64 poly;
  u32 n;
};

OPA_NM_CRYPTO_END

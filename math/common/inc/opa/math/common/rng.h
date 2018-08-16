#pragma once

#include <opa_common.h>
#include <opa/utils/misc.h>

OPA_NAMESPACE_DECL3(opa, math, common)

extern thread_local std::mt19937 rng;
extern thread_local std::mt19937_64 rng64;
void reseed_rng();


class Rng {
 public:
   Rng();
   ~Rng();
  std::string bytes(int n) const;

  static Rng *get() { return opa::utils::Singleton2<Rng>::get(); }
 private:
  FILE *m_urandom_file;
};

class Rng2 {
 public:
  std::string bytes(int n) const;

  static Rng2 *get() { return opa::utils::Singleton2<Rng2>::get(); }
 private:
};


OPA_NAMESPACE_DECL3_END

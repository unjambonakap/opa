#include "misc.h"
#include <opa/math/common/Matrix.h>
#include <opa/math/common/Zn.h>

using namespace opa::math::common;

OPA_NAMESPACE_DECL2(opa, crypto)

class HillCipherImpl : public HillCipher {
public:
  HillCipherImpl(u32 n, u32 mod) : zn(mod), key(&zn, n, n) { this->n = n; }

  virtual bool setKey(const std::vector<u32> &m) {
    assert(n * n == m.size());

    key.initialize(&zn, n, n);
    ikey.initialize(&zn, n, n);

    REP(i, n) REP(j, n) key(i, j) = m[i * n + j];

    bool can = key.invert(&ikey);
    assert(can);
    return true;
  }
  virtual std::vector<u32> decrypt(const std::vector<u32> &cipher) {
    return ikey.eval(cipher);
  }
  virtual std::vector<u32> encrypt(const std::vector<u32> &plain) {
    return key.eval(plain);
  }

private:
  u32 n;
  Zn zn;
  Matrix<u32> key;
  Matrix<u32> ikey;
};

HillCipher *HillCipher::get(u32 n, u32 mod) {
  return new HillCipherImpl(n, mod);
}

OPA_NAMESPACE_DECL2_END

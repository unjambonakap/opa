#pragma once

#include <opa/crypto/cracker.h>

OPA_NM_CRYPTO_CRACKER

class ZipState : public CrackerState {
public:
  ZipState() {}
  virtual void reset(const u8 *src, int n);

  virtual void decrypt(const u8 *src, u8 *dest, int n);
  virtual void encrypt(const u8 *src, u8 *dest, int n);

private:
  u32 m_keys[3];
};

OPA_NM_CRYPTO_CRACKER_END

#pragma once

#include <opa/crypto/base.h>

typedef struct aes_key_st AES_KEY;

OPA_NM_CRYPTO

class Aes {
 public:
  Aes(opa::stolen::StringRef key, bool encrypt) { init(key, encrypt); }
  Aes() {}
  void init(opa::stolen::StringRef key, bool encrypt);

  std::string encrypt_raw(opa::stolen::StringRef data) const;
  std::string encrypt_ecb(opa::stolen::StringRef data) const;
  std::string encrypt_cbc(opa::stolen::StringRef data,
                          opa::stolen::StringRef iv) const;
  std::string encrypt_ctr(opa::stolen::StringRef data) const;

  std::string decrypt_raw(opa::stolen::StringRef data) const;
  std::string decrypt_ecb(opa::stolen::StringRef data) const;
  std::string decrypt_cbc(opa::stolen::StringRef data,
                          opa::stolen::StringRef iv) const;
  std::string decrypt_ctr(opa::stolen::StringRef data) const;

  const AES_KEY* get_key() const { return m_aes_key.get(); }
  AES_KEY* get_key() { return m_aes_key.get(); }

 private:
  int m_bs;
  std::string m_key;
  SPTR(AES_KEY) m_aes_key;
  bool m_is_encrypt;
};

OPA_NM_CRYPTO_END

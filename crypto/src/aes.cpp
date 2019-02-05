#include <opa/crypto/aes.h>
#include <openssl/aes.h>

OPA_NM_CRYPTO
namespace {
const bool kEncrypt = true;
}

void Aes::init(opa::stolen::StringRef key, bool encrypt) {
  m_aes_key.reset(new AES_KEY);
  m_key = key;
  m_bs = key.size();
  m_is_encrypt = encrypt;
  if (encrypt) {
    AES_set_encrypt_key((const u8 *)m_key.data(), m_key.size() * 8, get_key());
  } else {
    AES_set_decrypt_key((const u8 *)m_key.data(), m_key.size() * 8, get_key());
  }
}

std::string Aes::encrypt_raw(opa::stolen::StringRef data) const {
  OPA_CHECK0(data.size() == m_key.size());
  std::string res(m_key.size(), 0);
  AES_encrypt((const u8 *)data.data(), (u8 *)res.data(), get_key());
  return res;
}

std::string Aes::encrypt_ecb(opa::stolen::StringRef data) const {
  std::string res(data.size(), 0);
  OPA_CHECK0(data.size() % m_bs == 0);
  for (int i = 0; i < data.size(); i += m_bs) {
    AES_ecb_encrypt((const u8 *)data.data() + i, (u8 *)res.data() + i,
                    get_key(), kEncrypt);
  }
  return res;
}

std::string Aes::encrypt_cbc(opa::stolen::StringRef data,
                             opa::stolen::StringRef iv) const {
  std::string res(data.size(), 0);
  std::string niv = iv;
  AES_cbc_encrypt((const u8 *)data.data(), (u8 *)res.data(), data.size(),
                  get_key(), (u8 *)niv.data(), kEncrypt);
  return res;
}

std::string Aes::encrypt_ctr(opa::stolen::StringRef data) const {
  OPA_CHECK0(false);
  // std::string res(data.size());
  // AES_ctr128_encrypt((const u8 *)data.data(), (u8 *)res.data(), get_key(),
  //                   kEncrypt);
  // return res;
}

// Decryption
std::string Aes::decrypt_raw(opa::stolen::StringRef data) const {
  OPA_CHECK0(data.size() == m_key.size());
  std::string res(m_key.size(), 0);
  AES_decrypt((const u8 *)data.data(), (u8 *)res.data(), get_key());
  return res;
}
std::string Aes::decrypt_ecb(opa::stolen::StringRef data) const {
  std::string res(data.size(), 0);
  OPA_CHECK0(data.size() % m_bs == 0);
  for (int i = 0; i < data.size(); i += m_bs) {
    AES_ecb_encrypt((const u8 *)data.data() + i, (u8 *)res.data() + i,
                    get_key(), !kEncrypt);
  }
  return res;
}

std::string Aes::decrypt_cbc(opa::stolen::StringRef data,
                             opa::stolen::StringRef iv) const {
  std::string res(data.size(), 0);
  std::string niv = iv;
  AES_cbc_encrypt((const u8 *)data.data(), (u8 *)res.data(), data.size(),
                  get_key(), (u8 *)niv.data(), !kEncrypt);
  return res;
}

std::string Aes::decrypt_ctr(opa::stolen::StringRef data) const {
  OPA_CHECK0(false);
  // std::string res(data.size());
  // AES_ctr128_encrypt((const u8 *)data.data(), (u8 *)res.data(), get_key(),
  //                   !kEncrypt);
  // return res;
}

OPA_NM_CRYPTO_END

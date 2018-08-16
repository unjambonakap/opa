#include "hash.h"
#include <opa/utils/string.h>

#include <openssl/sha.h>
#include <openssl/md5.h>

using namespace opa::utils;

OPA_NM_CRYPTO

OPA_CLASS_STORE_REGISTER_BY_CL(Sha256);
OPA_CLASS_STORE_REGISTER_BY_CL(Sha512);
OPA_CLASS_STORE_REGISTER_BY_CL(Sha1);
OPA_CLASS_STORE_REGISTER_BY_CL(Md5);

const int SHA512_BLOCK_SIZE = SHA512_CBLOCK;
const int SHA512_OUTPUT_SIZE = SHA512_DIGEST_LENGTH;
const int SHA512_APPEND_LEN_SIZE = 16;

const int SHA256_BLOCK_SIZE = SHA256_CBLOCK;
const int SHA256_OUTPUT_SIZE = SHA256_DIGEST_LENGTH;
const int SHA256_APPEND_LEN_SIZE = 8;

const int SHA1_BLOCK_SIZE = SHA_CBLOCK;
const int SHA1_OUTPUT_SIZE = SHA_DIGEST_LENGTH;
const int SHA1_APPEND_LEN_SIZE = 8;

const int MD5_BLOCK_SIZE = MD5_CBLOCK;
const int MD5_OUTPUT_SIZE = MD5_DIGEST_LENGTH;
const int MD5_APPEND_LEN_SIZE = 8;

Hash::Hash(int block_size, int output_size, int append_len_size) {
  m_block_size = block_size;
  m_output_size = output_size;
  m_append_len_size = append_len_size;
}
void Hash::reset() const { do_reset(); }

const Hash &Hash::update_std(const std::string &s) const {
  return update((const u8 *)s.c_str(), s.size());
}

const Hash &Hash::update(const u8 *data, int n) const {
  do_update(data, n);
  return *this;
}

void Hash::update_raw_std(const std::string &s) {
  return update_raw((const u8 *)s.c_str(), s.size());
}

void Hash::update_raw(const u8 *data, int n) {
  OPA_CHECK(n % get_block_size() == 0, "Got ", n, get_block_size());
  for (int i = 0; i < n; i += get_block_size())
    do_update_raw(data + i);
}

void Hash::set_context(const u8 *data, int n, u64 len) {
  OPA_CHECK(n % get_output_size() == 0, "Got ", n, get_output_size());
  return do_set_context(data, len);
}
void Hash::set_context_std(const std::string &s, u64 len) {
  return set_context((const u8 *)s.c_str(), s.size(), len);
}

std::string Hash::get_context(u64 &res) const { return do_get_context(res); }

std::string Hash::get() const { return do_get(); }
std::string Hash::get_hex() const { return opa::utils::b2h(get()); }

// Sha256 funcs

Sha256::Sha256()
    : Hash(SHA256_BLOCK_SIZE, SHA256_OUTPUT_SIZE, SHA256_APPEND_LEN_SIZE) {
  m_ctx = new SHA256_CTX();
  reset();
}
Sha256::~Sha256() { delete m_ctx; }

void Sha256::do_reset() const { SHA256_Init(m_ctx); }

void Sha256::do_update(const u8 *data, int n) const {
  SHA256_Update(m_ctx, data, n);
}

void Sha256::do_update_raw(const u8 *data) { SHA256_Transform(m_ctx, data); }

std::string Sha256::do_get_context(u64 &res) const {
  const u32 *tb = (const u32 *)m_ctx->h;
  const int nblocks = get_output_size() / sizeof(tb[0]);
  u32 tsf[nblocks];

  REP (i, nblocks)
    tsf[i] = htobe32(tb[i]);
  res = ((u64)m_ctx->Nh) << 29 | m_ctx->Nl >> 3;
  return std::string((const char *)tsf, get_output_size());
}

void Sha256::do_set_context(const u8 *data, u64 len) {
  const u32 *tb = (const u32 *)data;
  const int nblocks = get_output_size() / sizeof(tb[0]);
  u32 *dest = m_ctx->h;

  REP (i, nblocks)
    dest[i] = be32toh(tb[i]);
  m_ctx->Nh = len >> 29;
  m_ctx->Nl = len << 3 & (u32)-1;
}

std::string Sha256::do_get() const {
  std::string res(get_output_size(), (char)0);
  SHA256_Final((u8 *)res.c_str(), m_ctx);
  return res;
}

// Sha512 funcs

Sha512::Sha512()
    : Hash(SHA512_BLOCK_SIZE, SHA512_OUTPUT_SIZE, SHA512_APPEND_LEN_SIZE) {
  m_ctx = new SHA512_CTX();
  reset();
}
Sha512::~Sha512() { delete m_ctx; }

void Sha512::do_reset() const { SHA512_Init(m_ctx); }

void Sha512::do_update(const u8 *data, int n) const {
  SHA512_Update(m_ctx, data, n);
}

void Sha512::do_update_raw(const u8 *data) { SHA512_Transform(m_ctx, data); }

std::string Sha512::do_get_context(u64 &res) const {
  const u64 *tb = (const u64 *)m_ctx->h;
  const int nblocks = get_output_size() / sizeof(tb[0]);
  u64 tsf[nblocks];

  REP (i, nblocks)
    tsf[i] = htobe64(tb[i]);
  res = ((u64)m_ctx->Nh) << 59 | m_ctx->Nl >> 3;
  return std::string((const char *)tsf, get_output_size());
}

void Sha512::do_set_context(const u8 *data, u64 len) {
  const u64 *tb = (const u64 *)data;
  const int nblocks = get_output_size() / sizeof(tb[0]);
  u64 *dest = (u64*)m_ctx->h;

  REP (i, nblocks)
    dest[i] = be64toh(tb[i]);
  m_ctx->Nh = len >> 59;
  m_ctx->Nl = len << 3;
}

std::string Sha512::do_get() const {
  std::string res(get_output_size(), (char)0);
  SHA512_Final((u8 *)res.c_str(), m_ctx);
  return res;
}

// Sha1 funcs

Sha1::Sha1() : Hash(SHA1_BLOCK_SIZE, SHA1_OUTPUT_SIZE, SHA1_APPEND_LEN_SIZE) {
  m_ctx = new SHA_CTX();
  reset();
}

Sha1::~Sha1() { delete m_ctx; }

void Sha1::do_reset() const { SHA1_Init(m_ctx); }

void Sha1::do_update(const u8 *data, int n) const {
  SHA1_Update(m_ctx, data, n);
}

void Sha1::do_update_raw(const u8 *data) { SHA1_Transform(m_ctx, data); }

std::string Sha1::do_get_context(u64 &res) const {
  const u32 *tb = (const u32 *)&m_ctx->h0;
  const int nblocks = get_output_size() / sizeof(tb[0]);
  u32 tsf[nblocks];

  REP (i, nblocks)
    tsf[i] = htobe32(tb[i]);
  res = ((u64)m_ctx->Nh) << 29 | m_ctx->Nl >> 3;
  return std::string((const char *)tsf, get_output_size());
}

void Sha1::do_set_context(const u8 *data, u64 len) {
  const u32 *tb = (const u32 *)data;
  const int nblocks = get_output_size() / sizeof(tb[0]);
  u32 *dest = &m_ctx->h0;

  REP (i, nblocks)
    dest[i] = be32toh(tb[i]);
  m_ctx->Nh = len >> 29;
  m_ctx->Nl = len << 3 & (u32)-1;
}

std::string Sha1::do_get() const {
  std::string res(get_output_size(), (char)0);
  SHA1_Final((u8 *)res.c_str(), m_ctx);
  return res;
}

// Md5 funcs

Md5::Md5() : Hash(MD5_CBLOCK, MD5_DIGEST_LENGTH, MD5_APPEND_LEN_SIZE) {
  m_ctx = new MD5_CTX();
  reset();
}

Md5::~Md5() { delete m_ctx; }

void Md5::do_reset() const { MD5_Init(m_ctx); }

void Md5::do_update(const u8 *data, int n) const { MD5_Update(m_ctx, data, n); }

void Md5::do_update_raw(const u8 *data) { MD5_Transform(m_ctx, data); }

std::string Md5::do_get_context(u64 &res) const {
  const u32 *tb = (const u32 *)&m_ctx->A;
  const int nblocks = get_output_size() / sizeof(tb[0]);
  u32 tsf[nblocks];

  REP (i, nblocks)
    tsf[i] = htobe32(tb[i]);
  res = ((u64)m_ctx->Nh) << 29 | m_ctx->Nl >> 3;
  return std::string((const char *)tsf, get_output_size());
}

void Md5::do_set_context(const u8 *data, u64 len) {
  const u32 *tb = (const u32 *)data;
  const int nblocks = get_output_size() / sizeof(tb[0]);
  u32 *dest = &m_ctx->A;

  REP (i, nblocks)
    dest[i] = be32toh(tb[i]);
  m_ctx->Nh = len >> 29;
  m_ctx->Nl = len << 3 & (u32)-1;
}

std::string Md5::do_get() const {
  std::string res(get_output_size(), (char)0);
  MD5_Final((u8 *)res.c_str(), m_ctx);
  return res;
}

OPA_NM_CRYPTO_END

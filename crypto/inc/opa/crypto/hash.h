#pragma once

#include <opa/crypto/base.h>
#include <opa_common.h>
#include <opa/utils/serialize.h>
#include <string>

typedef struct SHA256state_st SHA256_CTX;
typedef struct SHA512state_st SHA512_CTX;
typedef struct SHAstate_st SHA_CTX;
typedef struct MD5state_st MD5_CTX;

OPA_NM_CRYPTO

class Hash
    : public opa::utils::MapperFunc<opa::stolen::StringRef, std::string> {
 public:
  Hash(int block_size, int output_size, int append_len_size);
  virtual ~Hash() {}

  void reset() const;

  int get_block_size() const { return m_block_size; }
  int get_output_size() const { return m_output_size; }
  int get_append_len_size() const { return m_append_len_size; }

  const Hash &update(const u8 *data, int n) const;
  const Hash &update_std(const std::string &s) const;

  void update_raw(const u8 *data, int n);
  void update_raw_std(const std::string &s);

  void set_context(const u8 *data, int n, u64 len);
  void set_context_std(const std::string &s, u64 len);

  std::string get() const;
  std::string get_hex() const;
  virtual void operator()(const opa::stolen::StringRef &in,
                          std::string &out) const {
    reset();
    update((const u8 *)in.data(), in.size());
    out = get();
  }

  std::string get_context(u64 &res) const;

 protected:
  virtual void do_update(const u8 *data, int n) const = 0;
  virtual void do_update_raw(const u8 *data) = 0;
  virtual std::string do_get_context(u64 &res) const = 0;
  virtual void do_set_context(const u8 *data, u64 len) = 0;
  virtual std::string do_get() const = 0;
  virtual void do_reset() const = 0;

 private:
  int m_block_size;
  int m_output_size;
  int m_append_len_size;
};

class Sha256 : public Hash {
 public:
  Sha256();
  virtual ~Sha256();
  OPACS_GETTER_BY_CL(Sha256);

 protected:
  virtual void do_update(const u8 *data, int n) const;
  virtual void do_update_raw(const u8 *data);
  virtual std::string do_get_context(u64 &res) const;
  virtual void do_set_context(const u8 *data, u64 n);
  virtual std::string do_get() const;
  virtual void do_reset() const;

 private:
  mutable SHA256_CTX *m_ctx;
};

class Sha512 : public Hash {
 public:
  Sha512();
  virtual ~Sha512();
  OPACS_GETTER_BY_CL(Sha512);

 protected:
  virtual void do_update(const u8 *data, int n) const;
  virtual void do_update_raw(const u8 *data);
  virtual std::string do_get_context(u64 &res) const;
  virtual void do_set_context(const u8 *data, u64 n);
  virtual std::string do_get() const;
  virtual void do_reset() const;

 private:
  mutable SHA512_CTX *m_ctx;
};

class Sha1 : public Hash {
 public:
  Sha1();
  virtual ~Sha1();

  OPACS_GETTER_BY_CL(Sha1);

 protected:
  virtual void do_update(const u8 *data, int n) const;
  virtual void do_update_raw(const u8 *data);
  virtual std::string do_get_context(u64 &res) const;
  virtual void do_set_context(const u8 *data, u64 n);
  virtual std::string do_get() const;
  virtual void do_reset() const;

 private:
  mutable SHA_CTX *m_ctx;
};

class Md5 : public Hash {
 public:
  Md5();
  virtual ~Md5();
  OPACS_GETTER_BY_CL(Md5);

 protected:
  virtual void do_update(const u8 *data, int n) const;
  virtual void do_update_raw(const u8 *data);
  virtual std::string do_get_context(u64 &res) const;
  virtual void do_set_context(const u8 *data, u64 n);
  virtual std::string do_get() const;
  virtual void do_reset() const;

 private:
  mutable MD5_CTX *m_ctx;
};

OPA_NM_CRYPTO_END

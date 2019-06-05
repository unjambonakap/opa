#pragma once

#include <opa_common.h>
#include <opa/utils/serialize.h>
#include <opa/stolen/StringRef.h>

#define OPA_NM_CRYPTO                                                          \
  namespace opa {                                                              \
  namespace crypto {

#define OPA_NM_CRYPTO_END                                                      \
  }                                                                            \
  }

#define OPA_NM_CRYPTO_CRACKER                                                  \
  namespace opa {                                                              \
  namespace crypto {                                                           \
  namespace cracker {

#define OPA_NM_CRYPTO_CRACKER_END                                              \
  }                                                                            \
  }                                                                            \
  }

OPA_NM_CRYPTO_CRACKER

class CrackerChecker : public opa::utils::BaseStorable {
public:
  virtual ~CrackerChecker(){}
  virtual bool operator()(const std::string &str) const = 0;
  std::shared_ptr<CrackerChecker> non_owned_sptr() {
    return std::shared_ptr<CrackerChecker>(this, [](CrackerChecker *) {});
  }
};
OPA_DECL_SPTR(CrackerChecker, CrackerCheckerSptr);

struct Pattern : public opa::utils::ProtobufParams {
  virtual ~Pattern(){}
  std::string init;
  std::vector<int> mp;
  std::string charset;
  int shard_at=-1;
  std::vector<std::vector<int>> per_char_vals;
  OPA_TGEN_IMPL(init, mp, charset, shard_at, per_char_vals);
};

struct Res : public opa::utils::ProtobufParams {
  std::vector<std::string> tb;
  OPA_TGEN_IMPL(tb);
};

struct ShardParams : public opa::utils::ProtobufParams {
  Pattern pattern;
  std::shared_ptr<CrackerChecker> checker;
  std::string cur;
  int pos;
  int n_res;

  OPA_TGEN_IMPL(pattern, checker, cur, pos, n_res);
};

OPA_NM_CRYPTO_CRACKER_END

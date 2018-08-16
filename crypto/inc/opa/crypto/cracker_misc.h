#pragma once
#include <opa_common.h>
#include <opa/crypto/cracker_job.h>

OPA_NM_CRYPTO_CRACKER

struct MapperAndCheckerParams : public opa::utils::ProtobufParams {
  CrackerCheckerSptr checker;
  opa::utils::StrStrMapperFuncSptr mapper;
  OPA_TGEN_IMPL(checker, mapper);
};

class MapperAndChecker : public CrackerChecker {
public:
  MapperAndChecker() {}
  virtual bool operator()(const std::string &str) const override {
    std::string res;
    (*m_params.mapper)(str, res);
    return (*m_params.checker)(res);
  }

  OPA_SETTER(MapperAndCheckerParams, m_params, params);
  OPA_TGEN_IMPL(m_params);
  OPACS_GETTER_BY_CL(MapperAndChecker);

private:
  MapperAndCheckerParams m_params;
};

struct PatternCheckerParams : public opa::utils::ProtobufParams {
  std::vector<std::pair<int, char> > tb;

  void from_string(int pos, const opa::stolen::StringRef &str) {
    tb.clear();
    REP (i, str.size()) { tb.pb(MP(i + pos, str[i])); }
  }
  OPA_TGEN_IMPL(tb);
};

class PatternChecker : public CrackerChecker {
public:
  virtual bool operator()(const std::string &str) const override {
    for (const auto &x : m_params.tb) {
      if (str.size() <= x.ST)
        return false;
      if (str[x.ST] != x.ND)
        return false;
    }
    return true;
  }

  OPACS_GETTER_BY_CL(PatternChecker);
  OPA_SETTER(PatternCheckerParams, m_params, params);
  OPA_TGEN_IMPL(m_params);

private:
  PatternCheckerParams m_params;
};

OPA_NM_CRYPTO_CRACKER_END

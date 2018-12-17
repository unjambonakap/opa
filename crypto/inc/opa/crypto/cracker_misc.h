#pragma once
#include <opa/crypto/cracker_job.h>
#include <opa_common.h>

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

  bool check(const std::string &str) const {
    OPA_CHECK0(!tb.empty());
    if (tb[0].ST == -1) {
      int mv = tb.back().ST + 1;
      REP (i, mv + 1) {
        REP (j, tb.size()) {
          if (str[i + tb[j].ST + 1] != tb[j].ND) goto bad;
        }
        return true;
      bad:;
      }
      return false;
    } else {
      for (const auto &x : tb) {
        if (str.size() <= x.ST) return false;
        if (str[x.ST] != x.ND) return false;
      }
      return true;
    }
  }
  OPA_TGEN_IMPL(tb);
};

struct MultiplePatternCheckerParams : public opa::utils::ProtobufParams {
  std::vector<PatternCheckerParams> cx;
  OPA_SETTER(std::vector<PatternCheckerParams>, cx, cx);
  OPA_TGEN_IMPL(cx);
};

class PatternChecker : public CrackerChecker {
public:
  virtual bool operator()(const std::string &str) const override {
    return m_params.check(str);
  }

  OPACS_GETTER_BY_CL(PatternChecker);
  OPA_SETTER(PatternCheckerParams, m_params, params);
  OPA_TGEN_IMPL(m_params);

private:
  PatternCheckerParams m_params;
};

class MultiplePatternChecker : public CrackerChecker {
public:
  virtual bool operator()(const std::string &str) const override {
    for (auto &c : m_params.cx)
      if (c.check(str)) return true;
    return false;
  }

  OPACS_GETTER_BY_CL(MultiplePatternChecker);
  OPA_SETTER(MultiplePatternCheckerParams, m_params, params);
  OPA_TGEN_IMPL(m_params);

private:
  MultiplePatternCheckerParams m_params;
};

OPA_NM_CRYPTO_CRACKER_END

#pragma once
#include <opa/crypto/base.h>
#include <opa/utils/DataStruct.h>
#include <opa/threading/auto_job.h>
#include <opa/threading/runner.h>

OPA_NM_CRYPTO_CRACKER

extern std::string CHARSET_UPPER;
extern std::string CHARSET_LOWER;
extern std::string CHARSET_ALPHA;
extern std::string CHARSET_NUM;
extern std::string CHARSET_ALPHANUM;
extern std::string CHARSET_ALL;

struct CrackerParams : public opa::utils::ProtobufParams {
  std::vector<Pattern> patterns;
  std::shared_ptr<CrackerChecker> checker;
  bool single_res;
  OPA_TGEN_IMPL(checker);
};

class Cracker : public opa::threading::AutoJob<ShardParams, Res> {
public:
  void init(const CrackerParams &params) {
    m_params = params;
    m_unit.single_res = params.single_res;
    m_unit.checker = params.checker;
  }

  opa::threading::Job &get_job() { return *this; }

  bool go(int pos, bool gen) {
    if (gen && pos == std::min((int)pattern().mp.size() - 1, 1)) {
      bool more;
      m_unit.pos = pos;
      cb()(m_unit, more);
      return !more;
    }

    if (pos == pattern().mp.size()) {
      if (checker()(cur())) {
        m_res.tb.pb(cur());
        return m_unit.single_res;
      }
      return false;
    }

    REP (i, pattern().charset.size()) {
      cur()[pattern().mp[pos]] = pattern().charset[i];
      if (go(pos + 1, gen))
        return true;
    }
    return false;
  }

  virtual void auto_worker_do_work(const ShardParams &data,
                                   Res &out_res) override {
    m_unit = data;
    m_res.tb.clear();
    go(m_unit.pos, false);
    out_res = m_res;
  }

  virtual void auto_server_get_work() override {
    for (auto &x : m_params.patterns) {
      pattern() = x;
      cur() = pattern().init;
      if (go(0, true))
        return;
    }
  }

  virtual void
  auto_server_set_work_result(const Res &res,
                              opa::threading::DataId data_id) override {
    m_res.tb.insert(m_res.tb.begin(), ALL(res.tb));
  }

  virtual bool server_want_more_results() const override {
    return res().tb.size() == 0 || !m_params.single_res;
  }

  OPA_ACCESSOR(Pattern, m_unit.pattern, pattern);
  OPA_ACCESSOR(std::string, m_unit.cur, cur);
  OPA_ACCESSOR(Res, m_res, res);
  OPA_ACCESSOR(CrackerChecker, *m_unit.checker, checker);

  OPA_TGEN_IMPL(m_params);
  OPA_CLOUDY_JOB_DECL;

private:
  Res m_res;
  CrackerParams m_params;
  ShardParams m_unit;
};

OPA_NM_CRYPTO_CRACKER_END

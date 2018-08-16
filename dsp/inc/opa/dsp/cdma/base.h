#pragma once

#include <opa_common.h>
#include <opa/dsp/inc.h>
#include <opa/dsp/sym/base.h>
#include <opa/dsp/cdma/encoder.h>

OPA_NAMESPACE(opa, dsp)

class CDMAMultiplexer : public opa::utils::Initable {
public:
  virtual void init(double pw) {
    opa::utils::Initable::init();
    m_pw = pw;
  }

  void add_user(CDMAEncoderSPTR user) { m_users.pb(user); }

  bool has_next() const {
    for (auto &user : m_users)
      if (!user->has_next())
        return false;
    return true;
  }

  double get_next() {
    check_init();
    OPA_CHECK0(has_next());
    double res = 0;
    double pw_per_user = m_pw / m_users.size();
    for (auto &user : m_users) {
      DataType v = user->get_next();
      double constellation_point = user->get_mapped(v);
      res += pw_per_user * constellation_point;
    }
    return res;
  }

  OPA_ACCESSOR(std::vector<CDMAEncoderSPTR>, m_users, users);

private:
  std::vector<CDMAEncoderSPTR> m_users;
  double m_pw;
};

OPA_NAMESPACE_END(opa, dsp)

#pragma once

#include <opa/dsp/inc.h>

OPA_NAMESPACE(opa, dsp, gps)

constexpr int PRN_PERIOD = 1023;

class Satellite : public opa::utils::Initable {
public:
  struct Params {
    std::vector<u32> l1ca_seq;
  };
  Satellite(const Params &params) { init(params); }
  Satellite() {}

  virtual void init(const Params &params);

  int get_l1(int x) const;

  Params m_params;
};

class Satellites {
public:
  OPA_SINGLETON_GETTER(Satellites);

  int nsats() const { return m_sats.size(); }
  const Satellite &get(int i) const { return m_sats[i]; }

private:
  std::vector<Satellite> m_sats;
  Satellites();
};

OPA_NAMESPACE_END(opa, dsp, gps)

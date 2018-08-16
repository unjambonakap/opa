#pragma once

#include <opa_common.h>
#include <opa/crypto/lfsr.h>

OPA_NAMESPACE(opa, dsp)

typedef u32 DataType;
class DataProvider {
 public:
  virtual bool has_next() const = 0;
  virtual DataType get_next() {
    ++m_sym_count;
    return _get_next();
  }

  DataProvider *base() { return this; }
  OPA_ACCESSOR_R(int, m_sym_count, sym_count);

 protected:
  virtual DataType _get_next() = 0;

 private:
  int m_sym_count = 0;
};
OPA_DECL_SPTR(DataProvider, DataProviderSptr);

class RandomDataProvider : public DataProvider {
 public:
  virtual bool has_next() const override { return true; }
  virtual DataType _get_next() override { return opa::math::common::rng() % 2; }
};

class SimpleDataProvider : public DataProvider {
 public:
  virtual bool has_next() const override { return m_data.size() > 0; }

  virtual DataType _get_next() override {
    DataType res = m_data.front();
    m_data.pop_front();
    return res;
  }

  OPA_ACCESSOR(std::deque<DataType>, m_data, data);

 private:
  std::deque<DataType> m_data;
};

class CDMAEncoder : public opa::utils::Initable {

 public:
  struct DecoderResult {
    std::vector<std::pair<DataType, double> > prob_res;
  };
  struct Params {
    SPTR(opa::crypto::LFSR<DataType>) lfsr;
    u32 chip_size;
    DataProviderSptr data_provider;
  };

  virtual void init(const Params &params) {
    opa::utils::Initable::init();
    m_params = params;
    m_chip_pos = chip_size();
  }
  bool has_next() const { return data_provider()->has_next(); }

  DataType get_next() {
    OPA_CHECK0(has_next());
    if (m_chip_pos == chip_size())
      m_chip_pos = 0, m_cur_sym = data_provider()->get_next();
    ++m_chip_pos;
    DataType res = m_cur_sym;
    DataType lfsr_entry = get_next_key();
    res ^= lfsr_entry;
    return res;
  }

  DataType get_next_key() {
    DataType res = lfsr()->get_next();
    return res;
  }

  double get_mapped(DataType v) const { return -int(v) * 2 + 1; }

  DecoderResult decode(const std::vector<double> &data) {
    OPA_CHECK0(data.size() == chip_size());
    std::vector<DataType> keyed;
    double v = 0;
    REP(i, chip_size()) { v += data[i] * get_mapped(get_next_key()); }
    OPA_DISP0(v);
    DecoderResult res;
    return res;
  }
  OPA_ACCESSOR_R(SPTR(opa::crypto::LFSR<DataType>), m_params.lfsr, lfsr);
  OPA_ACCESSOR_R(u32, m_params.chip_size, chip_size);
  OPA_ACCESSOR(DataProviderSptr, m_params.data_provider, data_provider);

 private:
  Params m_params;
  DataType m_cur_sym;
  int m_chip_pos;
};
OPA_DECL_SPTR(CDMAEncoder, CDMAEncoderSPTR);

OPA_NAMESPACE_END(opa, dsp)

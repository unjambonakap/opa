#pragma once
#include <opa/crypto/la/base.h>
#include <opa/utils/misc.h>
#include <opa_common.h>

OPA_NAMESPACE(opa, crypto, la)
class SimpleRelDriver;

class RelationDriver : public opa::utils::Initable {
public:
  virtual ~RelationDriver() {}
  virtual void reset() = 0;
  virtual bool should_add_rel(const Relation::RelationCost &cost) = 0;
  virtual void add_rel(const Relation &rel) = 0;

  virtual void set_rels(Relations &rels) const =0;
  OPA_ACCESSOR_R(int, m_input_lim, input_lim);
  SimpleRelDriver *to_simple() { return (SimpleRelDriver *)this; }

protected:
  int m_input_lim;
  mutable Relations rels;
};

class SimpleRelDriver : public RelationDriver {
public:
  SimpleRelDriver(double thresh, int count) { init(thresh, count); }

  virtual void init(double thresh, int count) {
    RelationDriver::init();
    m_data.init(count);
    m_thresh = thresh;
    m_count = count;
  }

  ~SimpleRelDriver() { fini(); }
  virtual void fini() override { m_data.fini(); }

  virtual bool should_add_rel(const Relation::RelationCost &cost) override {
    if (cost.bias() < m_thresh)
      return false;
    return m_data.should_add(cost);
  }
  virtual void add_rel(const Relation &rel) override {
    m_data.add(rel.cost, rel);
  }

  virtual void reset() override { m_data.clear(); }
  OPA_ACCESSOR(double, m_thresh, thresh);
  virtual void set_rels(Relations &rels) const override {
    m_data.get_items(rels.tb);
  }

private:

  double m_thresh;
  int m_count;
  opa::utils::KBestContainer<Relation::RelationCost, Relation> m_data;
};

OPA_NAMESPACE_END(opa, crypto, la)

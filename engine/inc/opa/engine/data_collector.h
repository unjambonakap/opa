#pragma once

#include <opa/engine/conf.h>
#include <opa/engine/common.h>

OPA_NAMESPACE_DECL2(opa, engine)
class Game;
class DataCollectorInternal;

class DataCollector : public opa::utils::Initable {
public:
  DataCollector();
  virtual void init(const Game *game);
  virtual void fini();
  ~DataCollector();

  void update();
  OPA_ACCESSOR_R(double, m_fps, fps);

private:
  const Game *m_game;
  double m_fps;
  std::unique_ptr<DataCollectorInternal> m_internal;
};
OPA_DECL_SPTR(DataCollector, DataCollectorPtr)

OPA_NAMESPACE_DECL2_END

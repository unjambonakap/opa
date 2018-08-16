#include "data_collector.h"
#include "Renderer.h"
#include <experimental/optional>

using namespace opa::math::game;
OPA_NAMESPACE_DECL2(opa, engine)

class DataCollectorInternal {
public:
  std::experimental::optional<u64> t;
};

void DataCollector::init(const Game *game) {
  opa::utils::Initable::init();
  m_internal.reset(new DataCollectorInternal);
  m_game = game;
  m_fps = 0;
}

DataCollector::DataCollector() {}
DataCollector::~DataCollector() {}

void DataCollector::fini() { opa::utils::Initable::fini(); }

void DataCollector::update() {
  u64 tmp = m_game->time().us();

  if (m_internal->t) {
    m_fps = m_fps * 0.9 + (1e6 / (tmp - *m_internal->t)) * 0.1;
  }
  m_internal->t = tmp;
}

OPA_NAMESPACE_DECL2_END

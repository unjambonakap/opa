#pragma once

#include <opa/engine/conf.h>
#include <opa/engine/common.h>

OPA_NAMESPACE_DECL2(opa, engine)
class Renderer;

template <class T> class DelayedContainer {
public:
  void init(int n) {
    this->n = n;
    objs.resize(n);
    pos = -n;
  }

  T &next_write() {
    int cpos = pos;
    pos = (pos + 1) % n;
    return objs[(cpos + n) % n];
  }
  bool can_read() { return pos >= 0; }

  T &next_read() { return objs[pos]; }

  int n;
  int pos;
  std::vector<T> objs;
};

class RenderStreamer : public opa::utils::Initable {
public:
  RenderStreamer(){}
  virtual void init(const Renderer *renderer);
  virtual void fini();

  void update();

private:
  DelayedContainer<u32> m_container;
  const Renderer *renderer;
  std::ofstream ofs;
  Image img;
};
OPA_DECL_SPTR(RenderStreamer, RenderStreamerPtr)

OPA_NAMESPACE_DECL2_END

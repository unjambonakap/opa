#pragma once

#include <opa/engine/Game.h>
#include <opa/engine/Strip.h>
#include <opa/engine/conf.h>

OPA_NAMESPACE_DECL2(opa, engine)
class SceneObj;

class DrawData {
public:
  DrawData() {}
  DrawData(const glm::mat4 &proj, const SceneObj *base)
      : m_proj(proj), m_base(base) {}

  OPA_ACCESSOR(glm::mat4, m_proj, proj)
  const SceneObj *base() const { return m_base; }

private:
  glm::mat4 m_proj;
  const SceneObj *m_base;
};

class Drawer : public opa::utils::Initable {
public:
  virtual void init() override { opa::utils::Initable::init(); }
  virtual void draw(const DrawData &data) = 0;
};
OPA_DECL_SPTR(Drawer, DrawerPtr);

class FuncDrawer : public Drawer {
public:
  typedef std::function<void(const DrawData &)> DrawFunc;
  virtual void init(DrawFunc func) {
    Drawer::init();
    m_func = func;
  }
  virtual void draw(const DrawData &data) override { m_func(data); }

private:
  DrawFunc m_func;
};

class DrawerHandler : public opa::utils::Initable {
public:
  virtual void init(SceneObj *obj);
  void draw(const DrawData &data, bool from_cam = false);
  void register_drawer(DrawerPtr drawer);
  OPA_ACCESSOR_PTR(SceneObj, m_obj, obj);

private:
  std::set<DrawerPtr> m_drawers;
  SceneObj *m_obj;
};

class TriangleContainerDrawer : public Drawer {
public:
  struct Params {
    TriangleContainer *tr;
    Params(TriangleContainer *tr) { this->tr = tr; }
    Params() {}
  };
  virtual void init(const Params &params);
  virtual void draw(const DrawData &data) override;

private:
  Params m_params;
};

OPA_NAMESPACE_DECL2_END

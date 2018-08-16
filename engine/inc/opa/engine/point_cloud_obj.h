#pragma once

#include <opa/engine/Drawer.h>
#include <opa/engine/Game.h>
#include <opa/engine/Intersection.h>
#include <opa/engine/Points.h>
#include <opa/engine/Scene.h>
#include <opa/engine/conf.h>

OPA_NAMESPACE_DECL2(opa, engine)
class SceneObj;

struct PointCloudDrawerParams {
  PointCloud *obj;
  PointCloudDrawerParams(PointCloud *obj) { this->obj = obj; }
  PointCloudDrawerParams() {}
};

class PointCloudDrawer : public Drawer {
public:
  virtual void init(const PointCloudDrawerParams &params);
  virtual void draw(const DrawData &data) override;

private:
  PointCloudDrawerParams m_params;
};

struct PointCloudIsecParams {
  PointCloud *obj;
  PointCloudIsecParams() {}
  PointCloudIsecParams(PointCloud *obj) { this->obj = obj; }
};

class PointCloudIsec : public Isec {
public:
  virtual bool find_intersection(const Ray &ray,
                                 IntersectionResult &res) override;
  virtual void init(const PointCloudIsecParams &params) {
    Isec::init();
    this->m_params = params;
  }

private:
  PointCloudIsecParams m_params;
};
OPA_DECL_SPTR(PointCloudIsec, PointCloudIsecPtr);

class PointCloudSceneObj : public SceneObj {
public:
  OPA_ACCESSOR(PointCloud, m_pc, pc);
  void init(Game *game);

protected:
  PointCloud m_pc;
};
OPA_DECL_SPTR(PointCloudSceneObj, PointCloudSceneObjPtr);

OPA_NAMESPACE_DECL2_END

#include "point_cloud_obj.h"

OPA_NAMESPACE_DECL2(opa, engine)

using namespace opa::math::game;
using namespace opa::utils;

void PointCloudSceneObj::init(Game *game) {
  SceneObj::init(game);
  auto x = new PointCloudIsec;
  x->init(PointCloudIsecParams(&this->pc()));
  auto y = new PointCloudDrawer;
  y->init(PointCloudDrawerParams(&this->pc()));

  isec().register_isec(IsecPtr(x));
  drawer().register_drawer(DrawerPtr(y));
}

void PointCloudDrawer::draw(const DrawData &data) { m_params.obj->draw(data); }

void PointCloudDrawer::init(const PointCloudDrawerParams &params) {

  Drawer::init();
  m_params = params;
}

bool PointCloudIsec::find_intersection(const Ray &ray,
                                       IntersectionResult &res) {
  res.found = false;
  return false;
}

OPA_NAMESPACE_DECL2_END

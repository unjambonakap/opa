#include "Drawer.h"
#include "Scene.h"
#include "SceneTreeHelper.h"
#include "Strip.h"

OPA_NAMESPACE_DECL2(opa, engine)

void DrawerHandler::init(SceneObj *obj) { m_obj = obj; }

void DrawerHandler::draw(const DrawData &data, bool from_cam) {
  obj()->check_register();
  if (obj()->status().hide) return;

  DrawData local_data = data;
  if (!from_cam) {
    local_data.proj() = data.proj() * obj()->pos_obj().mat_to_world();
  }

  for (auto x : obj()->tree_struct()->children)
    x->obj->drawer().draw(local_data);
  for (auto &x : m_drawers)
    x->draw(local_data);
}

void DrawerHandler::register_drawer(DrawerPtr drawer) {
  m_drawers.insert(drawer);
}

void TriangleContainerDrawer::init(const Params &params) {
  Drawer::init();
  m_params = params;
}

void TriangleContainerDrawer::draw(const DrawData &data) {

  for (auto x : m_params.tr->strips().used())
    m_params.tr->strips().get(x).draw(data);

  for (auto x : m_params.tr->triangles().used())
    m_params.tr->triangles().get(x).draw(data);
}

OPA_NAMESPACE_DECL2_END

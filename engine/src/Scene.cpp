#include "Scene.h"
#include "Camera.h"
#include "PosAccessor.h"
#include "SceneTreeHelper.h"

using namespace std;

OPA_NAMESPACE_DECL2(opa, engine)

SceneObj::SceneObj() {
  m_game = 0;
  m_tree_struct = 0;
}

Scene::Scene() {}
void Scene::init(Game *game) {
  status().drawable = 0;
  SceneObj::init(game);
}

void Scene::draw() {
  DrawData data = camera()->get_draw_data();
  data.proj() = data.proj() * camera()->obj()->pos()->get_mat_from(glm::mat4(1.));
  drawer().draw(data, true);
}

void SceneObj::set_tree_struct(TreeStruct *ts) { m_tree_struct = ts; }
void SceneObj::init(Game *game) {
  opa::utils::Initable::init();
  m_game = game;

  m_isec.init(this);
  m_drawer.init(this);

  m_pos.reset(new PosAccessor());
  m_pos->init(this);
}

void SceneObj::check_register() { game()->check_register(this); }

OPA_NAMESPACE_DECL2_END

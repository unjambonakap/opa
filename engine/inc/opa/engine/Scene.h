#pragma once

#include <opa/engine/Drawer.h>
#include <opa/engine/Game.h>
#include <opa/engine/Intersection.h>
#include <opa/engine/PosAccessor.h>
#include <opa/engine/PosObject.h>
#include <opa/engine/SceneTreeHelper.h>
#include <opa/engine/Strip.h>
#include <opa/engine/conf.h>

OPA_NAMESPACE_DECL2(opa, engine)

class SceneObj;
OPA_DECL_SPTR(SceneObj, SceneObjPtr);

struct SceneObjStatus {
  int ghost : 1;
  int drawable : 1;
  int highlight : 1;
  int weak : 1;
  int hide : 1;

  SceneObjStatus() {
    ghost = 0;
    drawable = 1;
    highlight = 0;
    weak = 0;
    hide = 0;
  }
};

class SceneObj : public opa::utils::Initable {
public:
  SceneObj();

  virtual void init(Game *game);

  OPA_ACCESSOR(PosObject, m_pos_obj, pos_obj);
  OPA_ACCESSOR(DrawerHandler, m_drawer, drawer);
  OPA_ACCESSOR(IntersectionHandler, m_isec, isec);

  OPA_ACCESSOR(SceneObjStatus, m_status, status);
  OPA_ACCESSOR_PTR(Game, m_game, game);
  OPA_ACCESSOR_PTR(TreeStruct, m_tree_struct, tree_struct);
  OPA_ACCESSOR_PTR(PosAccessor, m_pos.get(), pos);
  OPA_ACCESSOR(ObjId, m_id, id);
  OPA_ACCESSOR_PTR(void, m_private, priv);
  OPA_ACCESSOR(std::string, m_debug, debug);
  OPA_ACCESSOR(bool, m_fixed, fixed);

  void set_tree_struct(TreeStruct *ts);

  //void draw(const DrawData &data);
  void check_register();
  virtual opa::math::game::Mesh *get_mesh(double precision = 0.) const {
    OPA_CHECK0(false);
  };

protected:
  PosAccessorPtr m_pos;
  PosObject m_pos_obj;

  TreeStruct *m_tree_struct;
  std::string m_debug;

  bool m_fixed = false;
  DrawerHandler m_drawer;
  IntersectionHandler m_isec;
  SceneObjStatus m_status;
  Game *m_game;
  ObjId m_id;
  void *m_private;

private:
  virtual void init() override {}
};
OPA_DECL_SPTR(SceneObj, SceneObjPtr)

class Scene : public SceneObj {
public:
  Scene();
  virtual void init(Game *game) override;

  void draw();
  OPA_ACCESSOR_PTR(Camera, m_game->camera(), camera);

private:
};
OPA_DECL_SPTR(Scene, ScenePtr)

OPA_NAMESPACE_DECL2_END

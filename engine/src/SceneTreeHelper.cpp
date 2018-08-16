#include "SceneTreeHelper.h"
#include <opa/engine/Scene.h>

using namespace opa::math::game;

OPA_NAMESPACE_DECL2(opa, engine)

TreeStruct::TreeStruct() {
  obj = 0;
  par = 0;
  depth = -1;
}

SceneTreeHelper::SceneTreeHelper() {}

void SceneTreeHelper::init() {}

void SceneTreeHelper::add(SceneObj *obj, SceneObj *par) {
  TreeStruct *st = create(obj);
  TreeStruct *par_st = nullptr;
  if (par != nullptr) {
    par_st = get(par);
    OPA_ASSERT0(par_st != nullptr);
    par_st->children.push_back(st);
  }

  st->par = par_st;
  st->depth = (par_st ? par_st->depth + 1 : 0);
  obj->set_tree_struct(st);
}

void SceneTreeHelper::remove_from_map_only(const std::vector<SceneObj *> &objs) {
  for (auto &obj : objs) m_map.erase(obj);
}

void SceneTreeHelper::remove(SceneObj *obj) {
  TreeStruct *st = get(obj);
  OPA_ASSERT0(st != nullptr);
  TreeStruct *par_st = get(parent(obj));
  par_st->children.resize(std::remove(ALL(par_st->children), st) -
                          par_st->children.begin());
  m_map.erase(obj);
}

SceneTreeHelper::TreePath SceneTreeHelper::get_path(const SceneObj *xa,
                                                    const SceneObj *xb) const {
  TreeStruct *a = get(xa);
  TreeStruct *b = get(xb);
  SceneObj *cur;
  TreePath p1, p2;
  while (a->depth < b->depth)
    p2.push_back(PathNode(b, TreeDir::Down)), b = b->par;
  while (b->depth < a->depth)
    p1.push_back(PathNode(a, TreeDir::Up)), a = a->par;

  while (a != b) {
    p1.push_back(PathNode(a, TreeDir::Up)), a = a->par;
    p2.push_back(PathNode(b, TreeDir::Down)), b = b->par;
  }

  reverse(ALL(p2));
  p1.insert(p1.end(), ALL(p2));
  return p1;
}

SceneObj *SceneTreeHelper::parent(const SceneObj *a) const {
  TreeStruct *x = get(a);
  if (x->par == nullptr) return nullptr;
  return x->par->obj;
}

TreeStruct *SceneTreeHelper::get(const SceneObj *obj) const {
  if (!m_map.count(obj)) return nullptr;
  return m_map.find(obj)->second;
}

TreeStruct *SceneTreeHelper::create(SceneObj *obj) {
  OPA_CHECK0(m_map.count(obj) == 0);
  TreeStruct *st = new TreeStruct;
  st->par = nullptr;
  st->depth = -1;
  st->obj = obj;
  m_map[obj] = st;
  return st;
}

void SceneTreeHelper::list_children_internal(
  SceneObj *obj, std::vector<SceneObj *> &lst) const {
  for (auto &e : obj->tree_struct()->children) {
    lst.push_back(e->obj);
    list_children_internal(e->obj, lst);
  }
}

std::vector<SceneObj *> SceneTreeHelper::list_children(SceneObj *obj) const {
  std::vector<SceneObj *> res;
  list_children_internal(obj, res);
  return res;
}

SceneObj *SceneTreeHelper::get_non_weak(SceneObj *a) const {
  while (a->status().weak) a = parent(a);
  return a;
}

OPA_NAMESPACE_DECL2_END

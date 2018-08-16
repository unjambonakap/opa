#pragma once

#include <opa/engine/conf.h>

OPA_NAMESPACE_DECL2(opa, engine)
class SceneObj;

struct TreeStruct {
  SceneObj *obj;
  TreeStruct *par;
  s32 depth;
  std::vector<TreeStruct *> children;
  TreeStruct();
};

class SceneTreeHelper : public opa::utils::Initable {
public:
  enum TreeDir { Down, Up };

  struct PathNode {
    TreeStruct *x;
    TreeDir d;
    PathNode(TreeStruct *x, TreeDir d) : x(x), d(d) {}
  };

  typedef std::vector<PathNode> TreePath;

  SceneTreeHelper();
  virtual void init();
  void add(SceneObj *obj, SceneObj *par);
  void remove(SceneObj *obj);
  void remove_from_map_only(const std::vector<SceneObj*> &objs);
  std::vector<SceneObj *> list_children(SceneObj *obj) const;
  SceneObj *get_non_weak(SceneObj *a) const;

  TreePath get_path(const SceneObj *a, const SceneObj *b) const;
  SceneObj *parent(const SceneObj *a) const;
  TreeStruct *get(const SceneObj *obj) const;

private:
  void list_children_internal(SceneObj *obj,
                              std::vector<SceneObj *> &lst) const;
  TreeStruct *create(SceneObj *obj);

  std::unordered_map<const SceneObj *, TreeStruct *> m_map;
};
OPA_DECL_SPTR(SceneTreeHelper, SceneTreeHelperPtr)

OPA_NAMESPACE_DECL2_END

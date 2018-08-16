#pragma once
#include <opa/engine/conf.h>

OPA_NAMESPACE_DECL2(opa, engine)
class SceneObj;

class PosAccessor : public opa::utils::Initable {
public:
  PosAccessor();

  virtual void init(SceneObj *obj);
  void do_mov(const glm::vec3 &dp, const SceneObj *base = nullptr);
  void do_rot(const Rot &v, const SceneObj *base = nullptr);

  Rot get_rot_from(const Rot &v, const SceneObj *base = nullptr);
  Rot get_rot_to(const Rot &v, const SceneObj *base = nullptr);

  Dir get_front(const SceneObj *base = nullptr) const;
  Dir get_left(const SceneObj *base = nullptr) const;
  Dir get_up(const SceneObj *base = nullptr) const;

  Dir get_dir_from(const Dir &dir, const SceneObj *base = nullptr) const;
  Dir get_dir_to(const Dir &dir, const SceneObj *base = nullptr) const;
  Pos get_pos_from(const Pos &pos, const SceneObj *base = nullptr) const;
  Pos get_pos_to(const Pos &pos, const SceneObj *base = nullptr) const;

  Pos get_self_pos(const SceneObj *base = nullptr) const;
  OPA_ACCESSOR_PTR(SceneObj, m_obj, obj);

  glm::mat4 get_mat_to(const glm::mat4 &v,
                       const SceneObj *base = nullptr) const;
  glm::mat4 get_mat_from(const glm::mat4 &v,
                         const SceneObj *base = nullptr) const;

private:
  Rot rot_apply_path(const SceneObj *from, const SceneObj *to,
                     const Rot &v) const;
  Pos pos_apply_path(const SceneObj *from, const SceneObj *to,
                     const Pos &v) const;
  Dir dir_apply_path(const SceneObj *from, const SceneObj *to,
                     const Dir &v) const;
  glm::mat4 mat_apply_path(const SceneObj *from, const SceneObj *to,
                           const glm::mat4 &v) const;

  const SceneObj *normalize(const SceneObj *base) const;
  virtual void init() {}
  SceneObj *m_obj;
};
OPA_DECL_SPTR(PosAccessor, PosAccessorPtr);

OPA_NAMESPACE_DECL2_END

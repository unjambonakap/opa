#pragma once
#include <opa/engine/conf.h>

OPA_NAMESPACE_DECL2(opa, engine)

class PosObject {
public:
  PosObject();

  const glm::mat4 &mat() const;
  glm::mat4 mat_to_world() const;
  glm::mat4 mat_to_local() const;

  OPA_ACCESSOR(glm::vec3, m_pos, pos);
  OPA_ACCESSOR(glm::quat, m_rot, rot);
  OPA_ACCESSOR(double, m_scale, scale);

  Rot rot_to_world() const;
  Rot rot_to_local() const;

  void mov_rel(const glm::vec3 &v);
  void mov_abs(const glm::vec3 &v);
  void rot_rel(const glm::quat &r);
  void rot_abs(const glm::quat &r);

  glm::vec3 get_front() const;
  glm::vec3 get_left() const;
  glm::vec3 get_up() const;

  glm::vec3 get_world_dir(const Dir &dir) const;
  glm::vec3 get_local_dir(const Dir &dir) const;
  glm::vec3 get_local_pos(const Pos &pos) const;
  glm::vec3 get_world_pos(const Pos &pos) const;

private:
  double m_scale = 1;
  glm::vec3 m_pos;
  glm::quat m_rot;
  mutable glm::mat4 m_mat;
};
OPA_DECL_SPTR(PosObject, PosObjectPtr);

OPA_NAMESPACE_DECL2_END

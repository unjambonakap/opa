#include "PosObject.h"

using namespace opa::math::game;

OPA_NAMESPACE_DECL2(opa, engine)

const glm::mat4 &PosObject::mat() const {
  m_mat = glm::translate(glm::mat4(1.), m_pos) * glm::mat4_cast(m_rot) *
          glm::scale(Pos(m_scale, m_scale, m_scale));
  return m_mat;
}

PosObject::PosObject() {
  m_pos = glm::vec3(0);
  m_rot = glm::quat(1., 0, 0, 0);
}

void PosObject::mov_rel(const glm::vec3 &v) { m_pos += m_rot * v; }
void PosObject::mov_abs(const glm::vec3 &v) { m_pos += v; }

void PosObject::rot_rel(const glm::quat &r) { m_rot = m_rot * r; }
void PosObject::rot_abs(const glm::quat &r) { m_rot = r * m_rot; }

glm::vec3 PosObject::get_front() const { return m_rot * vec_x; }
glm::vec3 PosObject::get_left() const { return m_rot * vec_y; }
glm::vec3 PosObject::get_up() const { return m_rot * vec_z; }
glm::vec3 PosObject::get_world_dir(const Dir &dir) const { return m_rot * dir; }
glm::vec3 PosObject::get_local_dir(const Dir &dir) const {
  return glm::inverse(m_rot) * dir;
}
glm::vec3 PosObject::get_local_pos(const Pos &pos) const {
  return get_local_dir((pos - m_pos) / m_scale);
}
glm::vec3 PosObject::get_world_pos(const Pos &pos) const {
  return get_world_dir(pos * m_scale) + m_pos;
}

Rot PosObject::rot_to_world() const { return m_rot; }
Rot PosObject::rot_to_local() const { return glm::inverse(m_rot); }

glm::mat4 PosObject::mat_to_world() const { return mat(); }
glm::mat4 PosObject::mat_to_local() const { return glm::inverse(mat()); }

OPA_NAMESPACE_DECL2_END

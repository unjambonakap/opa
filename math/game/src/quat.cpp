#include "opa/math/game/quat.h"

#include <opa/math/common/Utils.h>
#include <opa/utils/string.h>

using namespace std;
static const float eps = 1e-6;

OPA_NAMESPACE_DECL3(opa, math, game)

glm::quat quat_from_vec_rot_safe(const glm::vec3 &vec, float angle) {
  return quat_from_vec_rot(glm::normalize(vec), angle);
}

glm::quat quat_from_vec_rot(const glm::vec3 &vec, float angle) {
  double cx = cos(angle / 2);
  double sx = sin(angle / 2);
  glm::quat res = glm::quat(cx, sx * vec);
  return res;
}

glm::quat quat_tsf_vec_safe(const glm::vec3 &from, const glm::vec3 &to) {
  return quat_tsf_vec(glm::normalize(from), glm::normalize(to));
}

glm::quat quat_tsf_vec(const glm::vec3 &from, const glm::vec3 &to) {
  float dot = glm::dot(from, to);
  if (OPA_FLOAT_EQ(fabs(dot), 1, eps)) {
    if (dot > 0)
      return glm::quat(1., 0, 0, 0);

    glm::vec3 pivot = vec_make_ortho(vec_rand_uni(), from);
    return quat_from_vec_rot(pivot, OPA_PI);

  } else {
    glm::vec3 dir = glm::normalize(glm::cross(from, to));
    float angle = vec_get_angle(from, to);
    glm::quat res = quat_from_vec_rot(dir, angle);
    return res;
  }
}

glm::quat quat_look_at_safe(const glm::vec3 &front, const glm::vec3 &up) {
  glm::vec3 t1, t2;
  t1 = glm::normalize(front);
  t2 = glm::normalize(up);

  return quat_look_at(t1, vec_make_ortho(t2, t1));
}

glm::quat quat_xy(const Dir &x, const Dir &y) {
  return quat_look_at(x, glm::cross(x, y));
}

glm::quat quat_look_at(const glm::vec3 &front, const glm::vec3 &up) {
  glm::quat res = quat_tsf_vec(vec_x, front);
  glm::vec3 cur_up = res * vec_z;
  float angle = vec_get_angle(cur_up, up);

  glm::vec3 dir = glm::cross(cur_up, up);
  if (glm::dot(dir, front) < 0)
    angle = -angle;
  auto res2 = quat_from_vec_rot(front, angle);

  return res2 * res;
}

OPA_NAMESPACE_DECL3_END

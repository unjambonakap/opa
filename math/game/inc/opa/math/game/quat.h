#pragma once

#include <opa/math/game/base.h>
#include <opa/math/game/conf.h>

OPA_NAMESPACE_DECL3(opa, math, game)

glm::quat quat_from_vec_rot_safe(const glm::vec3 &vec, float angle);
glm::quat quat_from_vec_rot(const glm::vec3 &vec, float angle);
glm::quat quat_tsf_vec_safe(const glm::vec3 &from, const glm::vec3 &to);
glm::quat quat_tsf_vec(const glm::vec3 &from, const glm::vec3 &to);
glm::quat quat_look_at_safe(const glm::vec3 &front, const glm::vec3 &up);
glm::quat quat_xy(const Dir &x, const Dir &y);
glm::quat quat_look_at(const glm::vec3 &front,
                       const glm::vec3 &up); // front and up ortho
inline glm::quat quat_from_vec4(const glm::vec4 &v) {
  return glm::quat(v.w, v.x, v.y, v.z);
}

inline glm::quat rot_to_box_space(const BoxSpec &box) {
  return glm::inverse(
    quat_look_at(glm::normalize(box.vec(0)), glm::normalize(box.vec(2))));
}

static inline glm::quat quat_from_euler(const glm::vec3 yaw_pitch_roll) {
  return glm::toQuat(
    glm::yawPitchRoll(yaw_pitch_roll.x, yaw_pitch_roll.y, yaw_pitch_roll.z));
}

static inline glm::vec3 quat_to_euler(const glm::quat rot) {
  return glm::eulerAngles(rot);
}

static inline double quat_dist(const Rot &a, const Rot &b) {
  return glm::length(quat_to_euler(glm::inverse(b) * a));
}

const opa::OR::C1Manifold CircleManifold = opa::OR::C1Manifold(2 * OPA_PI);

class SphereManifold : public opa::OR::NManifold<Pos2> {
public:
  SphereManifold()
      : opa::OR::NManifold<Pos2>({ &CircleManifold, &CircleManifold }) {}

  virtual double compute_dist(const PointType &a,
                              const PointType &b) const override {
    return glm::length(a - b);
  }
};

class RotManifold : public opa::OR::NManifold<Pos> {
public:
  RotManifold()
      : opa::OR::NManifold<Pos>(
          { &CircleManifold, &CircleManifold, &CircleManifold }) {}

  virtual double compute_dist(const PointType &a,
                              const PointType &b) const override {
    return glm::length(a - b);
  }

  static Rot to_rot(const PointType &pt) { return quat_from_euler(pt); }
  static Pos from_rot(const Rot &rot) { return quat_to_euler(rot); }
};

OPA_NAMESPACE_DECL3_END

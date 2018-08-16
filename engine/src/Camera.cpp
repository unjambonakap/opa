#include "Camera.h"

using namespace opa::math::game;

OPA_NAMESPACE_DECL2(opa, engine)

Camera::Camera() { m_params.zoom = 0; }

void Camera::update_viewport(int w, int h) {
  m_params.w = w;
  m_params.h = h;
}

float Camera::get_y_aperture() const {
  return std::min(PI / 2, PI / 4 * pow(0.8, -m_params.zoom));
}

float Camera::get_x_aperture() const {
  return get_y_aperture() * get_aspect_ratio();
}
float Camera::get_aspect_ratio() const { return 1.f * m_params.w / m_params.h; }

void Camera::setup_perspective() const {
  m_perspective =
    glm::perspective<float>(get_y_aperture(), get_aspect_ratio(),
                            m_params.near_field, m_params.far_field) *
    glm::lookAt(vec_0, vec_x, vec_z);
}

void Camera::setup() const {
  setup_perspective();
  m_iperspective = glm::inverse(m_perspective);
}

DrawData Camera::get_draw_data() {
  setup_perspective();
  return DrawData(m_perspective, obj().get());
}

Pos Camera::to_screen_pos(const Pos &pos) const {
  return ApplyMat(m_perspective, pos);
}

Pos Camera::from_screen_pos(const Pos &pos, const PosObject &campos) const {
  return ApplyMat(campos.mat_to_world() * m_iperspective, pos);
}
Pos Camera::to_screen_pos(const Pos &pos, const PosObject &campos) const {
  return ApplyMat(m_perspective * campos.mat_to_local(), pos);
}

// todo: return ray
Ray Camera::get_click_dir(const glm::vec2 &click) const {
  Pos screen_pos = Pos((click / Pos2(viewport()) - Pos2(0.5, 0.5)) * 2, 0);
  screen_pos.y *= -1;

  Ray res;
  res.pos = obj()->pos()->get_self_pos();
  res.dir = glm::normalize(this->from_screen_pos(screen_pos, obj()->pos_obj()) -
                           res.pos);
  return res;
}

PointVec Camera::get_view_frustrum_points_world(const PosObject &campos) const {
  PointVec res;
  REP (i, 8) {
    Pos pt(OPA_BITSIGN(i >> 2), OPA_BITSIGN(i >> 1), OPA_BITSIGN(i));
    res.push_back(from_screen_pos(pt, campos));
  }

  return res;
}

Ray Camera::get_click_dir_with_pos(const glm::vec2 &click,
                                   const PosObject &campos) const {
  Pos screen_pos =
    Pos((-click / Pos2(viewport()) + glm::vec2(0.5, 0.5)) * 2, 0);

  Ray res;
  res.pos = obj()->pos()->get_self_pos();
  res.dir = glm::normalize(this->from_screen_pos(screen_pos, campos) - res.pos);
  return res;
}

std::vector<math::game::HyperPlaneSpec>
Camera::hps_cam_space(const PosObject &campos) const {
  std::vector<math::game::HyperPlaneSpec> res;

  if (0) {
    Mat4 mat_norm(glm::transpose(campos.mat_to_world() * this->iperspective()));
    REP (i, 4) {
      Dir normal = Dir(i & 1, (1 ^ i) & 1, 0) * OPA_BITSIGN(i >> 1);
      Dir normal_other = ApplyMat(mat_norm, normal);
      double nv = 1 / glm::length(normal_other);
      OPA_DISP0(mat_norm, normal, normal_other, nv);
      res.push_back(HyperPlaneSpec{ PlaneSpec(normal_other * nv, nv) });
    }
  } else {
    std::vector<Pos> pos_near;
    std::vector<Pos> pos_far;
    Mat4 to_world = campos.mat_to_world() * this->iperspective();
    REP (i, 4) {
      pos_near.push_back(ApplyMat(
        to_world, Pos(OPA_BITSIGN(i >> 1), OPA_BITSIGN((i + 1) >> 1), -1)));
      pos_far.push_back(ApplyMat(
        to_world, Pos(OPA_BITSIGN(i >> 1), OPA_BITSIGN((i + 1) >> 1), 1)));
    }
    REP (i, 4) {
      Pos normal = math::game::get_plane_normal(
        pos_near[i], pos_near[(i + 1) % 4], pos_far[i]);
      if (glm::dot(pos_near[i] - pos_near[(i + 2) % 4], normal) < 0) normal = -normal;
      res.push_back(
        HyperPlaneSpec{ PlaneSpec(normal, glm::dot(normal, pos_near[i])) });
    }
  }

  return res;
}

OPA_NAMESPACE_DECL2_END

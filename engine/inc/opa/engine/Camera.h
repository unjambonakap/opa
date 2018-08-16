#pragma once

#include <opa/engine/Drawer.h>
#include <opa/engine/Scene.h>
#include <opa/engine/conf.h>

OPA_NAMESPACE_DECL2(opa, engine)
struct CameraParams {
  double zoom;
  int w, h;
  double near_field = 0.1;
  double far_field = 10000;
};
class Camera {
public:
  Camera();

  void update_viewport(int w, int h);
  DrawData get_draw_data();
  OPA_ACCESSOR(SceneObjPtr, m_obj, obj);
  OPA_ACCESSOR(CameraParams, m_params, params);

  OPA_ACCESSOR(double, m_params.zoom, zoom);
  OPA_ACCESSOR(double, m_params.near_field, near_field);
  OPA_ACCESSOR(double, m_params.far_field, far_field);
  OPA_ACCESSOR(Mat4, m_perspective, perspective);
  OPA_ACCESSOR(Mat4, m_iperspective, iperspective);

  Ray get_click_dir(const glm::vec2 &click) const;
  Ray get_click_dir_with_pos(const glm::vec2 &click,
                             const PosObject &campos) const;
  float get_y_aperture() const;
  float get_x_aperture() const;
  float get_aspect_ratio() const;

  Pos delta_pix_w() const { return Pos(1, 0, 0) / m_params.w; }
  Pos delta_pix_h() const { return Pos(0, 1, 0) / m_params.h; }

  void setup_perspective() const;
  void setup() const;
  Pos to_screen_pos(const Pos &pos) const;
  Pos to_screen_pos(const Pos &pos, const PosObject &campos) const;
  Pos from_screen_pos(const Pos &pos, const PosObject &campos) const;
  IPos2 viewport() const { return IPos2(m_params.w, m_params.h); }

  std::vector<math::game::HyperPlaneSpec>
  hps_cam_space(const PosObject &campos) const;

  PointVec get_view_frustrum_points_world(const PosObject &campos) const;

private:
  CameraParams m_params;
  SceneObjPtr m_obj;
  mutable glm::mat4 m_perspective, m_iperspective;
};
OPA_DECL_SPTR(Camera, CameraPtr)

OPA_NAMESPACE_DECL2_END

#include "RotatorService.h"
#include "Camera.h"
#include "Intersection.h"
#include "Scene.h"

using namespace opa::math::game;
using namespace std;

OPA_NAMESPACE_DECL2(opa, engine)
const ServiceKey RotatorService::KEY = "RotatorService";

void RotatorService::start() {
  m_rotate_listener =
    std::make_shared<InputCaller>([this](InputHandler *) { return false; });
  m_target_obj = 0;

  m_game->updater().add(UpdateMark(UpdateMark::InputUpdate + 1),
                        UpdaterPtr(new Updater([this]() { this->update(); })));
}

void RotatorService::stop() {}

void RotatorService::update() {
  const PosObject &pos = camera()->obj()->pos_obj();

  if (input_handler()->tracker().just_pressed(InputEnum::Key_R)) {
    if (m_target_obj)
      m_target_obj = 0;
    else {
      IntersectionResult res;
      Ray ray = camera()->get_click_dir(input_handler()->mouse().cur());

      if (scene()->isec().find_intersection(ray, res)) {
        m_target_obj = res.obj;
        OPA_CHECK0(m_target_obj != nullptr);
        while (m_target_obj->fixed()) {
          m_target_obj = game()->tree()->parent(m_target_obj);
          OPA_CHECK0(m_target_obj != nullptr);
        }
      }
    }
  }

  if (m_target_obj) {
    float speed = 0.05;
    auto df = input_handler()->mouse().diff();
    df.x = OPA_DRANGE_REMAP_CENTERED(0.5, df.x * speed);
    df.y = OPA_DRANGE_REMAP_CENTERED(0.5, df.y * speed);
    df = -df;

    auto v1 = camera()->obj()->pos()->get_dir_to(-vec_x, m_target_obj);
    auto v2 = camera()->obj()->pos()->get_dir_to(
      -vec_x + vec_y * df.x + vec_z * df.y, m_target_obj);
    auto quat2 = quat_tsf_vec_safe(v1, v2);

    m_target_obj->pos_obj().rot_rel(quat2);
  }
}

OPA_NAMESPACE_DECL2_END

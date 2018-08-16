#include "CameraObserver.h"
#include "Camera.h"
#include "InputHandler.h"
#include "Scene.h"
#include "common.h"

using namespace opa::math::game;
using namespace std;

OPA_NAMESPACE_DECL2(opa, engine)

const ServiceKey CameraObserver::KEY = "CameraObserver";

ModeManager::ModeManager() { activate_mode(CameraMode::Normal); }

bool ModeManager::activate_mode(CameraMode mode) {
  if (mode == CameraMode::Rotate) {
    if (in_mode(CameraMode::Pivot))
      return false;
    m_modes.erase(CameraMode::Normal);
  } else if (mode == CameraMode::Pivot) {
    if (in_mode(CameraMode::Rotate))
      return false;
    m_modes.erase(CameraMode::Normal);
  }

  m_modes.insert(mode);
  return true;
}

void ModeManager::deactivate_mode(CameraMode mode) {
  m_modes.erase(mode);
  if (m_modes.size() == 0)
    activate_mode(CameraMode::Normal);
}

bool ModeManager::in_mode(CameraMode mode) { return m_modes.count(mode); }

CameraObserver::CameraObserver() { m_game = 0; }

void CameraObserver::start() {
  input_handler()->check_init();

  m_key_cb = std::make_shared<InputCaller>(
    [this](InputHandler *) { return this->handle_key_event(); });

  m_mouse_key_cb = std::make_shared<InputCaller>(
    [this](InputHandler *) { return this->handle_mouse_key_event(); });

  FE_LIST(InputEnum, x, input_handler()->input_cbs()[x].add(m_key_cb),
          {
            InputEnum::Key_UP, InputEnum::Key_DOWN, InputEnum::Key_LEFT,
            InputEnum::Key_RIGHT,
          });
  input_handler()->input_cbs()[InputEnum::Key_C].add(std::make_shared<InputCaller>(
    [this](InputHandler *) { return this->handle_center_req(); }));

  FE_LIST(InputEnum, x, input_handler()->input_cbs()[x].add(m_mouse_key_cb),
          { InputEnum::Button_LEFT, InputEnum::Button_RIGHT });

  m_game->updater().add(
    UpdateMark(UpdateMark::InputUpdate + 1),
    UpdaterPtr(new Updater([this]() { this->update(); })));
}

void CameraObserver::stop() {}

bool CameraObserver::handle_key_event() {
  if (!m_mode.in_mode(CameraMode::Normal))
    return false;

  if (input_handler()->action().pressed()) {
    double speed = 0.1;
    Pos mov_vec(0.);
    switch (input_handler()->input()->type()) {
      OPA_CASE(InputEnum::Key_UP, mov_vec += vec_x;)
      OPA_CASE(InputEnum::Key_DOWN, mov_vec -= vec_x;)
      OPA_CASE(InputEnum::Key_LEFT, mov_vec += vec_y;)
      OPA_CASE(InputEnum::Key_RIGHT, mov_vec -= vec_y;)
    default:
      OPA_CHECK(false, "Unexpected key >> %s\n",
                input_handler()->input()->name().c_str());
      break;
    }
    camera()->obj()->pos_obj().mov_rel(mov_vec * speed);
    OPA_DISP0(mov_vec, speed, camera()->obj()->pos_obj().pos(), vec_x, vec_y);
    return true;
  }
  return false;
}

bool CameraObserver::handle_center_req() {
  if (!input_handler()->action().pressed())
    return false;
  SceneObj *main_obj = scene();
  if (input_handler()->mods().shift) {
    IntersectionResult res;
    Ray ray = camera()->get_click_dir(input_handler()->mouse().cur());

    if (scene()->isec().find_intersection(ray, res, true)) {
      main_obj = game()->tree()->get_non_weak(res.obj);
    } else
      return true;
  }
  std::vector<SceneObj *> objs = game()->tree()->list_children(main_obj);
  objs.push_back(main_obj);

  OPA_DISP("Got %d objs", objs.size());
  std::vector<opa::math::game::MeshData> meshes;
  for (auto &e : objs) {
    if (e->status().drawable) {
      meshes.push_back(
        { e->get_mesh(),
          e->pos()->get_mat_to(glm::mat4(1.), camera()->obj().get()) });
    }
  }

  Pos new_cam_pos(0.);
  Rot new_cam_rot;
  OPA_TRACEM("Computing new  center");
  compute_new_center(camera()->get_draw_data().proj(), meshes, &new_cam_pos,
                     &new_cam_rot);
  OPA_DISPL("Got new center ", new_cam_pos, new_cam_rot);
  camera()->obj()->pos_obj().mov_rel(new_cam_pos);
  camera()->obj()->pos_obj().rot_rel(new_cam_rot);
  return true;
}

bool CameraObserver::handle_mouse_key_event() {
  if (input_handler()->action().pressed()) {
    if (m_mode.activate_mode(CameraMode::Pivot)) {
      return true;
    }
  } else if (input_handler()->action().released()) {
    m_mode.deactivate_mode(CameraMode::Pivot);
    m_pivot.deactivate();
    return true;
  }
  return false;
}

void CameraObserver::update() {

  if (input_handler()->mods().ctrl == 1)
    m_mode.activate_mode(CameraMode::Rotate);
  else
    m_mode.deactivate_mode(CameraMode::Rotate);

  input_handler()->enable_cursor(!m_mode.in_mode(CameraMode::Rotate));
  PosObject &pos = camera()->obj()->pos_obj();

  if (m_mode.in_mode(CameraMode::Pivot)) {
    if (!m_pivot.is_activated()) {
      IntersectionResult res;
      Ray ray = camera()->get_click_dir(input_handler()->mouse().cur());

      if (scene()->isec().find_intersection(ray, res, true)) {
        Pos inter_pos = res.get_world_pos();
        cout << "find intersect at " << inter_pos << endl;
        cout << "loc space >> " << res.pos << endl;
        auto quat = quat_tsf_vec_safe(pos.get_front(), inter_pos - ray.pos);
        pos.rot_abs(quat);

        input_handler()->center_cursor();
        m_pivot.activate(res.obj, inter_pos);
      } else {
        puts("no reuslt");

        Pos inter_pos = Pos(0);
        auto quat = quat_tsf_vec_safe(pos.get_front(), inter_pos - ray.pos);
        pos.rot_abs(quat);

        input_handler()->center_cursor();
        m_pivot.activate(scene(), inter_pos);
      }

    } else {

      double speed = 0.1;
      auto df = input_handler()->mouse().diff();
      df.x = OPA_DRANGE_REMAP_CENTERED(0.5, df.x * speed);
      df.y = OPA_DRANGE_REMAP_CENTERED(0.5, df.y * speed);
      if (glm::length(df) > 1e-6) {
        auto quat =
          quat_tsf_vec_safe(vec_x, vec_x + vec_y * df.x + vec_z * df.y);
        PosObject &pos = camera()->obj()->pos_obj();
        pos.rot_rel(quat);
        double dist = glm::length(pos.pos() - m_pivot.center());
        pos.pos() = m_pivot.center() - dist * pos.get_front();
      }
    }
  }

  if (m_mode.in_mode(CameraMode::Rotate)) {

    if (1) {
      double speed = 0.01;
      auto df = input_handler()->mouse().diff();
      df.x = OPA_DRANGE_REMAP_CENTERED(0.5, -df.x * speed);
      df.y = OPA_DRANGE_REMAP_CENTERED(0.5, -df.y * speed);
      auto quat = quat_tsf_vec_safe(vec_x, vec_x + vec_y * df.x + vec_z * df.y);
      OPA_DISP0(df, speed, quat); 
      pos.rot_rel(quat);
    } else {
      double rot_speed = 0.05;

      auto df = input_handler()->mouse().diff();
      df.x = OPA_DRANGE_REMAP_CENTERED(PI / 8, -df.x * rot_speed);
      df.y = OPA_DRANGE_REMAP_CENTERED(PI / 8, df.y * rot_speed);

      glm::quat rot;
      rot = quat_from_vec_rot(vec_z, df.x);
      cout << rot * vec_x << " for angle " << df.x << endl;
      rot *= quat_from_vec_rot(vec_y, df.y);
      pos.rot_rel(rot);
    }
  }

  {
    double zoom_rate = -0.1;

    if (input_handler()->scroll().y < -0.5)
      zoom_rate *= -1;
    else if (input_handler()->scroll().y < 0.5)
      zoom_rate = 0;
    camera()->zoom() += zoom_rate;
    camera()->zoom() = std::max<double>(camera()->zoom(), -5);
    camera()->zoom() = std::min<double>(camera()->zoom(), 5);
  }
}

OPA_NAMESPACE_DECL2_END

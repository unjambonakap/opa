#pragma once

#include <opa/engine/conf.h>
#include <opa/engine/Game.h>
#include <opa/engine/Service.h>
#include <opa/engine/InputHandler.h>

OPA_NAMESPACE_DECL2(opa, engine)

enum CameraMode : int {
    Normal,
    Rotate,
    Pivot,
};

class PivotHelper {
  public:
    PivotHelper() { m_obj = 0; }

    bool is_activated() { return m_obj != 0; }
    void activate(SceneObj *obj, const Pos &center) {
        m_obj = obj;
        m_center = center;
    }
    void deactivate() { m_obj = 0; }

    OPA_ACCESSOR_R(Pos, m_center, center);

  private:
    SceneObj *m_obj;
    Pos m_center;
};

class ModeManager {
  public:
    ModeManager();
    bool activate_mode(CameraMode mode);
    void deactivate_mode(CameraMode mode);
    bool in_mode(CameraMode mode);

  private:
    std::set<CameraMode> m_modes;
};

class CameraObserver : public Service {
  public:
    static const ServiceKey KEY;

    virtual void start() override;
    virtual void stop() override;
    CameraObserver();

    bool handle_key_event();
    bool handle_center_req();
    bool handle_mouse_key_event();
    void update();

  private:
    virtual void init() override {}
    ModeManager m_mode;

    PivotHelper m_pivot;

    InputCallerPtr m_key_cb;
    InputCallerPtr m_mouse_key_cb;
};

OPA_NAMESPACE_DECL2_END

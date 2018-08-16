#pragma once

#include <opa/engine/conf.h>
#include <opa/engine/Game.h>
#include <opa/engine/Service.h>
#include <opa/engine/InputHandler.h>


OPA_NAMESPACE_DECL2(opa, engine)

class RotatorService : public Service {
  public:
    static const ServiceKey KEY;

    virtual void start() override;
    virtual void stop() override;
    void update();

  private:
    InputCallerPtr m_rotate_listener;

    SceneObj *m_target_obj;
};
OPA_NAMESPACE_DECL2_END

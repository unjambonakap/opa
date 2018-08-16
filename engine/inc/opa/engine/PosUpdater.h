#pragma once

#include <opa/engine/conf.h>
#include <opa/engine/Game.h>
#include <opa/engine/Service.h>
#include <opa/engine/Camera.h>

OPA_NAMESPACE_DECL2(opa, engine)

class PosUpdaterSub {
  public:
    virtual void notify(Camera *data) = 0;
};
OPA_DECL_SPTR(PosUpdaterSub, PosUpdaterSubPtr);

typedef PubSub<Camera *, PosUpdaterSub> PosUpdater;
OPA_DECL_SPTR(PosUpdater, PosUpdaterPtr);

OPA_NAMESPACE_DECL2_END

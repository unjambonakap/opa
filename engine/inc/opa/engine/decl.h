#pragma once

#include <opa/engine/conf.h>

OPA_NAMESPACE_DECL2(opa, engine)

class Camera;
class Renderer;
class InputHandler;
class Scene;
class SceneObj;
class ServiceManager;
class SceneTreeHelper;
class RenderStreamer;
class DataCollector;

OPA_DECL_SPTR(Camera, CameraPtr)
OPA_DECL_SPTR(Renderer, RendererPtr)
OPA_DECL_SPTR(InputHandler, InputHandlerPtr)
OPA_DECL_SPTR(Scene, ScenePtr)
OPA_DECL_SPTR(ServiceManager, ServiceManagerPtr)
OPA_DECL_SPTR(SceneTreeHelper, SceneTreeHelperPtr)
OPA_DECL_SPTR(SceneObj, SceneObjPtr)
OPA_DECL_SPTR(RenderStreamer, RenderStreamerPtr)
OPA_DECL_SPTR(DataCollector, DataCollectorPtr)


OPA_NAMESPACE_DECL2_END

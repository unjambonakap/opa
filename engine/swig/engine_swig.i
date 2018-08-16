%module(directors="1") opa_engine_swig

%include "std_string.i"
%include "typemaps.i"
%include "opa.i"
%include "glm_base.i"
%include "math_game_swig.i"


%include "opa/engine/Inputs_vals_swig.h"


%{
#include "opa/utils/debug.h"
#include "opa/engine/conf.h"
#include "opa/engine/decl.h"
#include "opa/engine/common.h"
#include "opa/engine/PosObject.h"
#include "opa/engine/Camera.h"
#include "opa/engine/data_collector.h"

#include "opa/engine/Game.h"
#include "opa/engine/Intersection.h"
#include "opa/engine/Inputs.h"
#include "opa/engine/InputHandler.h"
#include "opa/engine/primitives/Base.h"
#include "opa/engine/earth_render.h"
#include "opa/engine/Points.h"
#include "opa/engine/point_cloud_obj.h"

using namespace opa;
using namespace opa::engine;
using namespace opa::math::game;
%}
namespace std {
  %template(vec_autokd_node) std::vector<const opa::engine::AutoKDNode*>;
}

%include "opa/engine/conf.h"
%include "opa/engine/decl.h"



%feature("director") opa::engine::Updater_Swig;
%feature("nodirector") opa::engine::Updater_Swig::swig_call;
%feature("nodirector") opa::engine::Updater_Swig::to_updater;
%feature("director") opa::engine::InputCaller_Swig;
%feature("nodirector") opa::engine::InputCaller_Swig::swig_call;
%feature("nodirector") opa::engine::InputCaller_Swig::to_updater;

%shared_ptr(opa::engine::Scene);
%shared_ptr(opa::engine::SceneObj);
%shared_ptr(opa::engine::PrimitiveSceneObj);
%shared_ptr(opa::engine::Buffer);
%shared_ptr(opa::engine::Attribute);
%shared_ptr(opa::engine::Square);
%shared_ptr(opa::engine::Cube);
%shared_ptr(opa::engine::Icosahedron);
%shared_ptr(opa::engine::VBO);
%shared_ptr(opa::engine::TrVBO);
%shared_ptr(opa::engine::FanVBO);
%shared_ptr(opa::engine::Image);
%shared_ptr(opa::engine::Texture);
%shared_ptr(opa::engine::MultiColorTexture);
%shared_ptr(opa::engine::TrianglesSceneObj);

%shared_ptr(opa::engine::DataCollector);
%shared_ptr(opa::engine::Game);
%shared_ptr(opa::engine::Updater);
%shared_ptr(opa::engine::InputCaller);
%shared_ptr(opa::engine::InputHandler);
%shared_ptr(opa::engine::IntersectionHandler);
%shared_ptr(opa::engine::Isec);
%shared_ptr(opa::engine::SphereIsec);
%shared_ptr(opa::engine::TriangleContainerIsec);
%shared_ptr(opa::engine::IntersectionResult);

%shared_ptr(opa::engine::PointCloudSceneObj);
%shared_ptr(opa::engine::PointCloudIsec);
%shared_ptr(opa::engine::PointCloudDrawer);
%shared_ptr(opa::engine::PointCloud);

%typemap(in, numinputs=0) (std::vector<const opa::engine::AutoKDNode *> &nodes_out) (std::vector<const opa::engine::AutoKDNode*> data) { 
  $1=&data;
}
%apply std::vector<const opa::engine::AutoKDNode*>& nodes_out { std::vector<const opa::engine::AutoKDNode*>& };

//%typemap(argout) std::vector<const opa::engine::AutoKDNode*>& "";
//%apply std::vector<const opa::engine::AutoKDNode*>& { std::vector<const opa::engine::AutoKDNode*>& };

%typemap(argout) (std::vector<const opa::engine::AutoKDNode*>& nodes_out) {
  %append_output(swig::from(*$1));
}

%typemap(in, numinputs=0) (opa::engine::IntersectionResult& res) (opa::engine::IntersectionResult data) { 
  $1=&data;
}
%apply opa::engine::IntersectionResult& res { opa::engine::IntersectionResult& };

%typemap(argout) (opa::engine::IntersectionResult& res) {
  %append_output(SWIG_NewPointerObj(new opa::engine::IntersectionResultPtr(new opa::engine::IntersectionResult(*$1)), $descriptor(std::shared_ptr<opa::engine::IntersectionResult>*), SWIG_POINTER_OWN));
}

%ignore opa::engine::operator<<;

%include "opa/utils/debug.h"
%include "opa/engine/common.h"
%include "opa/engine/PosObject.h"
%include "opa/engine/data_collector.h"

%include "opa/engine/Game.h"
%include "opa/engine/Intersection.h"
%include "opa/engine/Scene.h"
%include "opa/engine/Camera.h"
%include "opa/engine/Inputs.h"
%include "opa/engine/InputHandler.h"
%include "opa/engine/primitives/Base.h"
%include "opa/engine/earth_render.h"
%include "opa/engine/Points.h"
%include "opa/engine/point_cloud_obj.h"



namespace opa{
namespace engine{

%template(RankCaller_Updater) opa::engine::RankCaller<UpdaterPtr>;
%template(RankCaller_InputCaller) opa::engine::RankCaller<InputCallerPtr>;
}
}

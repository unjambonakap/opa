%module opa_math_game_swig

%include "opa.i"
%include "opa_common_swig.i"
%include "or_swig.i"
%include "glm_base.i"
%include "math_common_swig_base.i"

%{
#include "opa/algo/graph.h"
#include "opa/math/game/base.h"
#include "opa/math/game/graph.h"
#include "opa/math/game/mesh.h"
#include "opa/math/game/mesh_util.h"
#include "opa/math/game/point_cloud.h"
#include "opa/math/game/geo_2d.h"
#include "opa/math/game/intersection.h"
#include "opa/math/game/mesh_primitives.h"
#include "opa/math/game/mesh_atomize.h"
#include "opa/math/game/base_impl.h"
#include "opa/math/game/swig_decl.h"


using namespace opa::math::game;
using namespace opa;
%}

%shared_ptr(opa::algo::FastGraph);

%shared_ptr(opa::math::game::Mesh);
%shared_ptr(opa::math::game::FaceCollection);
%shared_ptr(opa::math::game::MeshBuilder);
%shared_ptr(opa::math::game::GraphWithPosData);
%shared_ptr(opa::math::game::FaceGraph);
%shared_ptr(opa::math::game::FaceGraph2);

%rename(from_) opa::algo::EdgeData::from;

%include "opa/algo/graph.h"
%include "opa/math/game/base.h"

namespace opa{
namespace math{
namespace game{
%template(BoxAASpec_2DPython) opa::math::game::BoxAASpec_Gen<2, opa::Pos2>;
%template(BoxAASpec_3DPython) opa::math::game::BoxAASpec_Gen<3, opa::Pos>;
%template(BoxSpec_2DPython) opa::math::game::BoxSpec_Gen<2, opa::Pos2, opa::Mat2>;
%template(BoxSpec_3DPython) opa::math::game::BoxSpec_Gen<3, opa::Pos, opa::Mat3>;
}
}
}
%template(DataWhitenerPos) opa::math::common::DataWhitener<opa::Pos, 3>;

%include "opa/math/game/swig_decl.h"
%include "opa/math/game/base_impl.h"
%include "opa/math/game/graph.h"
%include "opa/math/game/mesh.h"
%include "opa/math/game/mesh_util.h"
%include "opa/math/game/point_cloud.h"
%include "opa/math/game/geo_2d.h"
%include "opa/math/game/mesh_primitives.h"
%include "opa/math/game/mesh_atomize.h"
%include "opa/math/game/intersection.h"
%ignore opa::math::game::C1Manifold;
%ignore opa::math::game::RotManifold;



namespace std {
  %template(vec_graph) std::vector<std::shared_ptr<opa::algo::FastGraph>>;
  %template(vec_meshptr) std::vector<std::shared_ptr<opa::math::game::Mesh>>;
  %template(vec_box) std::vector<opa::math::game::BoxSpec>;
  %template(vec_round_debug) std::vector<opa::math::game::atom::InitRoundDebug>;
  %template(vec_cluster_data) std::vector<opa::math::game::atom::ClusterData>;
  %template(vec_edge_data) std::vector<opa::algo::EdgeData>;
  %template(vec_plane) std::vector<opa::math::game::PlaneSpec>;
  %template(vec_hplane) std::vector<opa::math::game::HyperPlaneSpec>;
}


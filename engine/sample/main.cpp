#include <gdal.h>
#include <gdal_priv.h>
#include <opa/engine/CameraObserver.h>
#include <opa/engine/InputHandler.h>
#include <opa/engine/Renderer.h>
#include <opa/math/game/quat.h>
#include <opa_common.h>

#include <opa/engine/Camera.h>
#include <opa/engine/CameraObserver.h>
#include <opa/engine/InputHandler.h>
#include <opa/engine/Renderer.h>
#include <opa/engine/RotatorService.h>
#include <opa/engine/Scene.h>
#include <opa/engine/SceneTreeHelper.h>
#include <opa/engine/Service.h>
#include <opa/engine/Sphere.h>
#include <opa/engine/gdal_obj.h>
#include <opa/engine/graph/graph.h>
#include <opa/engine/primitives/Base.h>
#include <opa/engine/procgen/Texture.h>

#include <opa/engine/mesh/MeshDrawer.h>
#include <opa/engine/mesh/MeshIsec.h>
#include <opa/engine/mesh/mesh.h>
#include <opa/math/common/Utils.h>
#include <opa/math/game/intersection.h>
#include <opa/math/game/mesh_export.h>
#include <opa/math/game/mesh_loader.h>
#include <opa/math/game/mesh_util.h>
#include <opa/math/game/point_cloud.h>

#include <opa/engine/graph/graph.h>

#include <Poco/Exception.h>
#include <Poco/File.h>
#include <Poco/Foundation.h>
#include <Poco/Glob.h>
#include <Poco/Path.h>

DEFINE_string(action, "", "");
DEFINE_string(gdal_glob, "*.tif", "");
DEFINE_string(stl_file, "", "");
DEFINE_string(outfile, "", "");

using namespace opa::engine;
using namespace opa::math::game;
using namespace std;
using namespace opa;

typedef void (*action_func)();
std::unordered_map<string, action_func> actions;

FaceCollection load_stl() {

  OPA_CHECK0(!FLAGS_stl_file.empty());
  return FaceCollection()
    .load_stl(FLAGS_stl_file)
    .filter([](const PointVec &tr) {
      Pos target(15.0793, 5.95295, -10.8796);
      double v = 1e100;
      for (auto p : { tr[0], tr[1], tr[2] }) {
        v = std::min<double>(v, glm::length(p - target));
      }
      return (v > 1e-3);
    });
}

void game_test() {
  Game game;
  game.init();

  auto polyhedra = std::make_shared<DisplayMesh>();
  load_stl().to_mesh(polyhedra.get());
  {
    polyhedra->init();
    auto tex = new MultiColorTexture();
    int nc = 1000;
    tex->init(nc);
    polyhedra->set_texture(TexturePtr(tex));

    auto mesh_obj = std::make_shared<MeshObj>();
    mesh_obj->init(&game, polyhedra);
    mesh_obj->pos_obj().pos() = glm::vec3(0, 0, 0);
    mesh_obj->pos_obj().scale() = 0.01;
    game.register_obj(mesh_obj, game.scene());
  }

  auto sphere = std::make_shared<Sphere>();
  sphere->init(&game, 1);
  sphere->build(1000);
  SphereSpec res = get_min_enclosing_sphere(polyhedra->point_cloud());
  sphere->pos_obj().pos() = Pos(1, 1, 1);
  game.register_obj(sphere, game.scene());

  game.run();
}

void gdal_test() {
  Game game;
  game.init();

  if (0) {
    OPA_CHECK0(!FLAGS_gdal_glob.empty());
    std::set<std::string> gdal_files;

    try {
      auto p1 = Poco::Path(FLAGS_gdal_glob);
      Poco::Glob::glob(p1, gdal_files);
    } catch (const Poco::AssertionViolationException &e) {
      cout << e.message() << endl;
      throw;
    }
    GdalObj *obj = new GdalObj;
    game.register_obj(SceneObjPtr(obj), game.scene());
    obj->init(&game, vector<string>(ALL(gdal_files)));

    obj->pos_obj().rot() = quat_look_at_safe(vec_x + vec_y + vec_z, vec_z);
    obj->pos_obj().pos() = Pos(8, 8, 8);
  }
  game.run();
}

void graph_test() {
  Game game;
  game.init();

  graph::DisplayGraph<VertexBase, EdgeBase> graph_obj;
  graph_obj.init();
  SceneObj *obj = graph_obj.create_scene_obj();
  obj->init(&game);

  auto v1 = graph_obj.create_vertex().id;
  auto v2 = graph_obj.create_vertex().id;
  auto v3 = graph_obj.create_vertex().id;
  auto v4 = graph_obj.create_vertex().id;
  auto v5 = graph_obj.create_vertex().id;

  auto e1 = graph_obj.add_edge(v1, v2).id;
  auto e2 = graph_obj.add_edge(v2, v3).id;
  auto e3 = graph_obj.add_edge(v2, v4).id;
  auto e4 = graph_obj.add_edge(v2, v5).id;

  graph_obj.vertex_attr(v1).pos = Pos4(4, 4, 3, 0.2);
  graph_obj.vertex_attr(v1).col = Col(0, 0, 1);
  graph_obj.vertex_attr(v2).pos = Pos4(4, 4, 2, 0.3);
  graph_obj.vertex_attr(v2).col = Col(0, 1, 0);
  graph_obj.vertex_attr(v3).pos = Pos4(2 ,4, 1, 0.1);
  graph_obj.vertex_attr(v3).col = Col(1, 0, 0);
  graph_obj.vertex_attr(v4).pos = Pos4(5, 3, 1, 0.1);
  graph_obj.vertex_attr(v4).col = Col(0, 1, 1);
  graph_obj.vertex_attr(v5).pos = Pos4(5, 5, 1, 0.1);
  graph_obj.vertex_attr(v5).col = Col(1, 1, 1);

  graph_obj.update_edge(e1);
  graph_obj.update_edge(e2);
  graph_obj.update_edge(e3);
  graph_obj.update_edge(e4);
  graph_obj.edge_attr(e1).col = Col(0, 0, 1);
  graph_obj.edge_attr(e2).col = Col(0, 1, 1);
  graph_obj.edge_attr(e3).col = Col(1, 0, 1);
  graph_obj.edge_attr(e4).col = Col(1, 1, 1);

  game.register_obj(SceneObjPtr(obj), game.scene());
  game.run();
}

void test() {
  Game game;
  game.init();

  game.run();
}

int main(int argc, char *argv[]) {
  opa::init::opa_init(argc, argv);
  actions["gdal"] = gdal_test;
  actions["game"] = game_test;
  actions["graph"] = graph_test;
  actions["test"] = test;

  string action = FLAGS_action;
  OPA_CHECK0(!action.empty());
  OPA_CHECK0(actions.count(action));
  actions[action]();
  return 0;
}

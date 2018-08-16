#include <opa/math/game/point_cloud.h>
#include <opa/math/game/quat.h>
#include <opa_common.h>

#include <opa/math/common/Utils.h>
#include <opa/math/game/conf.h>
#include <opa/math/game/intersection.h>
#include <opa/math/game/mesh.h>
#include <opa/math/game/mesh_export.h>
#include <opa/math/game/mesh_loader.h>
#include <opa/math/game/mesh_util.h>
#include <opa/utils/files.h>

#include <Poco/Exception.h>
#include <Poco/File.h>
#include <Poco/Foundation.h>
#include <Poco/Glob.h>
#include <Poco/Path.h>

DEFINE_string(action, "", "");
DEFINE_string(stl_file, "", "");
DEFINE_string(outfile, "", "");

using namespace opa;
using namespace opa::utils;
using namespace opa::math::game;
using namespace std;

typedef void (*action_func)();
std::unordered_map<string, action_func> actions;

FaceCollection load_stl() {

  OPA_CHECK0(!FLAGS_stl_file.empty());
  return FaceCollection(true)
    .load_stl(FLAGS_stl_file)
    .filter([](const PointVec &tr) {
      Pos target(15.0793, 5.95295, -10.8796);
      double v = 1e100;
      for (auto p : tr) {
        v = std::min<double>(v, glm::length(p - target));
      }
      return (v > 1e-3);
    });
}

void export_test() {
  auto m = load_stl().to_mesh();
  proto::MeshFaceGraphList face_graph;

  m->get_face_graph(&face_graph);
  std::ofstream(FLAGS_outfile, std::ofstream::binary)
    << face_graph.SerializeAsString();
}

void clean_stl() {
  auto m = load_stl();
  m.write_stl(FLAGS_outfile);
}

void connected_test() {
  auto m = load_stl().to_mesh();

  auto res = m->split_connected();
  OPA_DISP0(res.size());
  proto::MeshFaceGraphList face_graph;
  dump_vector_to_files<SPTR(Mesh)>(
    FLAGS_outfile, res, [&](BufferWriter &writer, const SPTR(Mesh) & mesh) {
      writer.put(mesh->get_face_graph(&face_graph).DebugString());
    });
}

void sphere_test() {
  auto m = load_stl().to_mesh();
  auto pts = m->point_cloud();
  for (auto &p : pts) {
    OPA_DISP0(p);
  }

  SphereSpec res = get_min_enclosing_sphere(pts);
  OPA_DISP0(get_sphere_radius(pts, res.center));
  OPA_DISP0(res.center, res.radius);
}

void graph_test() {
  auto m = load_stl().to_mesh();
  SPTR(opa::algo::FastGraph) graph = m->compute_graph();
  OPA_DISP0(graph->n);
  for (auto &cc : graph->split_cc()) {
    auto walk = cc->get_cover_walk();
    OPA_DISP0(walk, cc->n, cc->n_edges(), walk.size());
  }
}

void cut_plan_test() {
  /*
  auto m = load_stl();
  Pos plane{ 1, 0, 0 };
  double pos = 0;
  auto mesh = m.filter([&](const Triangle &tr) {
    int a = std::signbit(point_signed_dist_to_plane(std::get<0>(tr), plane) -
  pos);
    int b = std::signbit(point_signed_dist_to_plane(std::get<1>(tr), plane) -
  pos);
    int c = std::signbit(point_signed_dist_to_plane(std::get<2>(tr), plane) -
  pos);
    return a == 0 && a == b && b == c;
  }).to_mesh();

  OPA_TRACEM("Splitting connected");
  auto con_list = mesh->split_connected();
  OPA_TRACEM("Done split");
  dump_vector_to_files<UPTR(Mesh)>(FLAGS_outfile, con_list,
                       [](BufferWriter &writer, const UPTR(Mesh)& mesh) {
                         mesh->to_triangles()->write_stl(writer);
                       });
                       */

  PointVec tb = {
    { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 },
  };

  std::vector<Mesh> meshes;
  if (0) {
    meshes.emplace_back();
    MeshBuilder::PolyhedraFromCloud(&meshes.back(), vec_point_to_info_data(tb));
  } else {
    auto m = load_stl().to_mesh();
    auto res = m->split_connected();
    for (auto &e : res) meshes.push_back(*e);
  }
  Pos plane_center = Pos(0.5, 0.5, -0);

  HyperPlaneSpec hp =
    HyperPlaneSpec{ PlaneSpec::FromPointAndVec(plane_center, Pos(1, 3, -5)) };

  if (0) {
    hp =
      HyperPlaneSpec{ PlaneSpec::FromPointAndVec(plane_center, Pos(0, 0, -1)) };
  }

  utils::FilenameSharder sharder;
  sharder.set_pattern("mesh_{}").set_number_elem(meshes.size()).build();
  proto::MeshList meshes_proto;

  {
    PointVec pc;
    for (auto &mesh : meshes) mesh.add_to_point_cloud(pc);
    Dir front_dir = hp.plane.get_rand_plane_vec();
    Point2Vec pc2 = ApplyOp<Pos, Pos2>(
      pc, [&](const Pos &pos) { return hp.plane.proj(pos, front_dir); });

    // auto debug_trs = TriangleCollection();
    //// mesh.to_triangles();
    // debug_trs.add_plane(hp.plane, plane_center, 10);
    // debug_trs.write_stl("/tmp/res.out");
    {
      auto tx = compute_aabb_box2d(pc2);
      OPA_DISP0(tx.low, tx.high);
    }
    Exporter()
      .to_mesh(
        FaceCollection().add_box(compute_aabb_box2d(pc2).to_box().expand(1.1),
                                 hp.plane, front_dir),
        meshes_proto.add_mesh())
      ->set_name("cut_plan");
  }

  REP (i, meshes.size()) {
    std::string name = sharder.get(i);
    Exporter()
      .to_mesh(meshes[i], meshes_proto.add_mesh())
      ->set_name(name + "_orig");
  }
  {
    std::ofstream tmpfile("/tmp/meshes.out", std::ios::binary);
    meshes_proto.SerializeToOstream(&tmpfile);
  }

  REP (i, meshes.size()) {
    if (i != 3) continue;
    std::string name = sharder.get(i);
    auto res = mesh_cut_plan(meshes[i], hp);
    Exporter().to_mesh_list(res, &meshes_proto, name + "_proc");
    OPA_TRACES(res.size());
  }

  {
    std::ofstream tmpfile("/tmp/meshes.out", std::ios::binary);
    meshes_proto.SerializeToOstream(&tmpfile);
  }
}

int main(int argc, char *argv[]) {
  opa::init::opa_init(argc, argv);
  actions["export"] = export_test;
  actions["connected"] = connected_test;
  actions["sphere"] = sphere_test;
  actions["cut_plan"] = cut_plan_test;
  actions["clean_stl"] = clean_stl;
  actions["graph"] = graph_test;

  string action = FLAGS_action;
  OPA_CHECK0(!action.empty());
  OPA_CHECK0(actions.count(action));
  actions[action]();
  return 0;
}

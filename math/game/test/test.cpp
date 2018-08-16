#include <gtest/gtest.h>
#include <opa/algo/base.h>
#include <opa/algo/graph.h>
#include <opa/math/game/geo_2d.h>
#include <opa/math/game/intersection.h>
#include <opa/math/game/mesh_atomize.h>
#include <opa/math/game/mesh_primitives.h>
#include <opa/math/game/mesh_util.h>
#include <opa/math/game/point_cloud.h>
#include <opa/math/game/quat.h>
#include <opa/or/adaptative_search.h>
#include <opa_common.h>

#include <lemon/concepts/graph.h>
#include <lemon/dijkstra.h>
#include <lemon/list_graph.h>
#include <lemon/matching.h>
#include <lemon/smart_graph.h>

using namespace opa::math::game;
using namespace opa::math::game::atom;
using namespace std;
using namespace opa;
using namespace opa::OR;

const float eps = 1e-6;
#define COMP_FLOAT(a, b) ASSERT_TRUE(OPA_FLOAT_EQ(a, b, eps))
#define COMP_VEC(a, b) ASSERT_TRUE(OPA_VEC_EQ(a, b, eps))

TEST(TsfVec, Test1) {
  glm::vec3 a, b;
  a = vec_rand_uni();
  b = vec_rand_uni();
  auto res = quat_tsf_vec(a, b);

  COMP_VEC(b, res * a);
}

TEST(MakeOrtho, Test1) {
  glm::vec3 a, b;
  a = vec_rand_uni();
  b = vec_rand_uni();
  auto res = vec_make_ortho(a, b);
  COMP_FLOAT(glm::dot(res, b), 0);
  COMP_FLOAT(glm::length(res), 1);
}

TEST(VecRot, Test1) {
  glm::vec3 a, b, rot_vec;
  a = vec_rand_uni();
  b = vec_rand_uni();
  rot_vec = glm::normalize(glm::cross(a, b));

  float angle = vec_get_angle(a, b);
  auto res = quat_from_vec_rot(rot_vec, angle);
  COMP_VEC(b, res * a);
}

TEST(VecRot, Test2) {
  glm::vec3 a, b, rot_vec;
  a = vec_rand_uni();
  b = vec_rand_uni();
  rot_vec = glm::normalize(glm::cross(a, b));

  float angle = vec_get_angle(a, b);
  auto res = quat_from_vec_rot(rot_vec, -angle);
  COMP_VEC(a, res * b);
}

TEST(LookAt, Test1) {
  glm::vec3 a, b;
  a = vec_rand_uni();
  b = vec_rand_uni();
  b = vec_make_ortho(b, a);
  auto res = quat_look_at(a, b);

  glm::vec3 x = vec_x;
  glm::vec3 z = vec_z;
  // cout << res *x << ' ' << res *z << endl;
  // cout << a << ' ' << b << endl;
  COMP_VEC(a, res * x);
  COMP_VEC(b, res * z);
}

TEST(TestQuat, Test1) {
  {
    glm::vec3 a, b;
    a = vec_rand_uni();
    b = vec_rand_uni();
    b = vec_make_ortho(b, a);
    auto res = quat_look_at(a, b);
    glm::vec3 res_xyx = glm::vec3(res.x, res.y, res.z);
    res_xyx = glm::normalize(res_xyx);

    OPA_DISP_VARS3(a, b, res_xyx);
  }

  {
    auto tmp = quat_from_vec_rot(vec_y, PI / 4);
    cout << tmp * vec_x << endl;

    cout << quat_from_vec_rot(vec_z, PI / 4) * vec_x << endl;
  }
}

TEST(Intersection, LineSphere1) {
  glm::vec3 p, d, res;
  float r;

  r = 2;
  p = glm::vec3(3, 0, 0);
  d = glm::normalize(glm::vec3(-1, 0.5, 0));
  ASSERT_TRUE(line_sphere_intersection(p, d, r, res));
  cout << res << endl;
}

TEST(Intersection, PointInTr) {
  glm::vec3 t2 = 2 * vec_x;
  glm::vec3 t3 = 2 * vec_y;
  bool is_tr = true;
  ASSERT_TRUE(point_in_parallelogram(glm::vec3(0, 0, 0), t2, t3, is_tr));
  ASSERT_TRUE(point_in_parallelogram(glm::vec3(0, 2, 0), t2, t3, is_tr));
  ASSERT_TRUE(point_in_parallelogram(glm::vec3(2, 0, 0), t2, t3, is_tr));
  ASSERT_TRUE(point_in_parallelogram(glm::vec3(1, 1, 0), t2, t3, is_tr));
  ASSERT_TRUE(point_in_parallelogram(glm::vec3(0.5, 0.5, 0), t2, t3, is_tr));

  ASSERT_FALSE(point_in_parallelogram(glm::vec3(2, 1, 0), t2, t3, is_tr));
  ASSERT_FALSE(point_in_parallelogram(glm::vec3(1, 1.01, 0), t2, t3, is_tr));
}

TEST(Intersection, PointInParallelogram) {
  glm::vec3 t2 = 2 * vec_x;
  glm::vec3 t3 = 2 * vec_y;
  bool is_tr = false;
  ASSERT_TRUE(point_in_parallelogram(glm::vec3(0, 0, 0), t2, t3, is_tr));
  ASSERT_TRUE(point_in_parallelogram(glm::vec3(0, 2, 0), t2, t3, is_tr));
  ASSERT_TRUE(point_in_parallelogram(glm::vec3(2, 0, 0), t2, t3, is_tr));
  ASSERT_TRUE(point_in_parallelogram(glm::vec3(1, 1, 0), t2, t3, is_tr));
  ASSERT_TRUE(point_in_parallelogram(glm::vec3(2, 2, 0), t2, t3, is_tr));
  ASSERT_TRUE(point_in_parallelogram(glm::vec3(1, 2, 0), t2, t3, is_tr));
  ASSERT_TRUE(point_in_parallelogram(glm::vec3(0.5, 0.5, 0), t2, t3, is_tr));

  ASSERT_FALSE(point_in_parallelogram(glm::vec3(2.1, 2, 0), t2, t3, is_tr));
  ASSERT_FALSE(point_in_parallelogram(glm::vec3(2.1, 0, 0), t2, t3, is_tr));
}

TEST(Intersection, LineTr1) {
  glm::vec3 t2 = 2 * vec_x;
  glm::vec3 t3 = 2 * vec_y;
  glm::vec3 t1 = glm::vec3();
  glm::vec3 res;
  bool is_tr = true;

  glm::vec3 p = glm::vec3(-5, -5, -5);
  ASSERT_TRUE(line_parallelogram_intersection(p, glm::vec3(0, 0, 0) - p, t1, t2,
                                              t3, res, is_tr));
  ASSERT_TRUE(line_parallelogram_intersection(p, glm::vec3(2, 0, 0) - p, t1, t2,
                                              t3, res, is_tr));
  ASSERT_TRUE(line_parallelogram_intersection(p, glm::vec3(0, 2, 0) - p, t1, t2,
                                              t3, res, is_tr));
  ASSERT_TRUE(line_parallelogram_intersection(p, glm::vec3(1, 1, 0) - p, t1, t2,
                                              t3, res, is_tr));
  ASSERT_FALSE(line_parallelogram_intersection(p, glm::vec3(1.1, 1, 0) - p, t1,
                                               t2, t3, res, is_tr));
}

TEST(Intersection, LinePr1) {
  glm::vec3 t2 = 2 * vec_x;
  glm::vec3 t3 = 2 * vec_y;
  glm::vec3 t1 = glm::vec3();
  glm::vec3 res;
  bool is_tr = false;

  glm::vec3 p = glm::vec3(-5, -5, -5);
  ASSERT_TRUE(line_parallelogram_intersection(p, glm::vec3(0, 0, 0) - p, t1, t2,
                                              t3, res, is_tr));
  ASSERT_TRUE(line_parallelogram_intersection(p, glm::vec3(2, 0, 0) - p, t1, t2,
                                              t3, res, is_tr));
  ASSERT_TRUE(line_parallelogram_intersection(p, glm::vec3(0, 2, 0) - p, t1, t2,
                                              t3, res, is_tr));
  ASSERT_TRUE(line_parallelogram_intersection(p, glm::vec3(1, 1, 0) - p, t1, t2,
                                              t3, res, is_tr));
  ASSERT_TRUE(line_parallelogram_intersection(p, glm::vec3(1.5, 1.9, 0) - p, t1,
                                              t2, t3, res, is_tr));
  ASSERT_TRUE(line_parallelogram_intersection(p, glm::vec3(2, 2, 0) - p, t1, t2,
                                              t3, res, is_tr));
  ASSERT_FALSE(line_parallelogram_intersection(p, glm::vec3(2.1, 2, 0) - p, t1,
                                               t2, t3, res, is_tr));
  ASSERT_FALSE(line_parallelogram_intersection(p, glm::vec3(-0.1, 0, 0) - p, t1,
                                               t2, t3, res, is_tr));
  ASSERT_FALSE(line_parallelogram_intersection(p, glm::vec3(0, -0.1, 0) - p, t1,
                                               t2, t3, res, is_tr));
}

TEST(Projection, Test1) {
  auto m = glm::perspective((float)PI / 2.f, 1.f, 0.01f, 100.0f) *
           glm::lookAt(vec_0, vec_x, -vec_z);
  auto m2 = glm::perspective((float)PI / 2.f, 1.f, 0.01f, 100.0f);
  glm::vec3 pt;
  pt = { 0, 0, -1 };
  OPA_DISP0(ApplyMat(m, pt), ApplyMat(m2, pt));
  pt = { 0, 1, -1 };
  OPA_DISP0(ApplyMat(m, pt), ApplyMat(m2, pt));
  pt = { 1, 0, -1 };
  OPA_DISP0(ApplyMat(m, pt), ApplyMat(m2, pt));
  pt = { -1, 0, -1 };
  OPA_DISP0(ApplyMat(m, pt), ApplyMat(m2, pt));
  pt = { -1, -1, -1 };
  OPA_DISP0(ApplyMat(m, pt), ApplyMat(m2, pt));

  pt = { 1, 0, 0 };
  glm::mat4 t1 = glm::mat4_cast(quat_look_at_safe(vec_x, vec_z));
  glm::mat4 t2 = glm::lookAt(vec_0, vec_x, vec_z);
  puts("FUU");
  OPA_DISP0(ApplyMat(t1, pt), ApplyMat(t2, pt));
}

TEST(Center, Test1) {
  PointVec tb = {
    { -1, 1, 0 },
    { 1, 1, 0 },
    { 3, -1, 0 },
    { -1, -1, 0 },
  };
  Mat4 proj_mat = glm::perspective((float)PI / 4.f, 2.f, 0.01f, 100.0f);
  OPA_DISP0(compute_perspective_fovx(glm::inverse(proj_mat)),
            compute_perspective_fovy(glm::inverse(proj_mat)));

  Pos new_center;
  Rot new_rot;
  compute_new_center(proj_mat, tb, &new_center, &new_rot, 0.9);

  auto trans_mat = glm::translate(glm::mat4(), -new_center);
  auto xz_mat = glm::inverse(glm::mat4_cast(new_rot)) * trans_mat;
  auto view_mat = glm::lookAt(vec_0, vec_z, vec_y) * xz_mat;
  OPA_DISP0(new_center, new_rot);

  for (auto &e : tb) {
    auto p1 = ApplyMat(view_mat, e);
    auto p2 = ApplyMat(proj_mat, p1);
    auto p3 = ApplyMat(xz_mat, e);
    OPA_DISP0(e, p1, p2, p3);
  }
}

TEST(EulerAngles, Test1) {
  auto r1 = quat_look_at_safe(vec_y, vec_z);
  Rot r = r1 * glm::inverse(quat_look_at_safe(vec_x, vec_z));

  auto e = quat_to_euler(r);
  OPA_DISP0(e);
  OPA_DISP0(glm::toQuat(glm::eulerAngleYXZ(e.y, e.x, e.z)) * vec_x);
  ;
  OPA_DISP0(r * vec_x);
  ;
}

TEST(ConvexHull, Test1) {
  Point2Vec tb = {
    { -1, 1 }, { 1, 1 }, { 3, 0 }, { -1, -1 }, { 1, -1 }, { 0, 0 },
  };
  std::vector<int> hull_ids = compute_convex_hull(tb);
  Point2Vec lst;
  for (auto x : hull_ids) lst.push_back(tb[x]);
  OPA_DISP0(lst);

  OPA_DISP0(inside_polygon_unsafe(lst, { 0, 0 }));
  OPA_DISP0(inside_polygon_unsafe(lst, { 4, 0 }));
  OPA_DISP0(inside_polygon_unsafe(lst, { 0, 0.9 }));
  OPA_DISP0(inside_polygon_unsafe(lst, { 0, 1.1 }));
}

TEST(BoundingBox2D, Test1) {
  Point2Vec tb = {
    { -1, 1 }, { 1, 1 }, { 3, 0 }, { -1, -1 }, { 1, -1 },
  };
  Rot2 r = d2_rot(PI / 4);
  // r = d2_rot(0);
  for (auto &pt : tb) pt = r * pt;
  OPA_DISP0(r);
  Box2DSpec res = compute_best_box2d(tb);
  OPA_DISPL0(res.area(), res.corner, res.v[0], res.v[1]);
  for (auto &pt : tb) {
    OPA_CHECK(res.in(pt), pt);
  }
}

TEST(BoundingBox2D, TestRand) {
  Point2Vec tb;
  int npoints = 100;
  REP (i, npoints)
    tb.push_back(vec2_rand_uni());
  Box2DSpec res = compute_best_box2d(tb);
  Box2DSpec res2 = compute_best_box2d_dumb(tb);
  OPA_DISPL0(res.area(), res.corner);
  OPA_DISPL0(res2.area(), res2.corner);
}

TEST(BoundingBox, TestRand) {
  PointVec tb;
  int npoints = 30;
  REP (i, npoints)
    tb.push_back(vec_rand_uni());
  REP (i, 10) {
    BoxSpec res2 = compute_best_box(tb);
    OPA_DISPL0(res2.area(), res2.corner);

    BoxSpec res = compute_best_box_vnaze(tb);
    OPA_DISPL0(res.area(), res.corner);

    BoxSpec res3 = compute_best_box_dumb(tb);
    OPA_DISPL0(res3.area(), res3.corner);
  }
}

TEST(BoundingBox, Test1) {
  PointVec tb = {
    { -1, 1, 0 }, { 1, 1, 0 }, { 3, -1, 0 }, { -1, -1, 0 }, { 1, 1, 1 },
  };
  Rot r = quat_from_euler(Pos(0.2, 0.4, 0.6));
  // Rot r;
  for (auto &pt : tb) pt = r * pt;
  OPA_DISP0(r);
  BoxSpec res = compute_best_box(tb);
  OPA_DISPL0(res.area(), res.corner, res.v[0], res.v[1], res.v[2]);
  for (auto &pt : tb) {
    OPA_CHECK(res.in(pt), pt);
  }
  BoxAASpec aabb = compute_aabb_box(tb, Rot());
  OPA_DISP0(aabb.area());
}

TEST(AdaptativeSearch, Test1) {
  OR::C1Manifold manifold;
  typedef OR::ConstRManagedDesc<double, OR::C1Manifold> CurManagedDesc;
  CurManagedDesc managed_desc(manifold, 10, 5);
  OR::RoundStoppingCriteria stop_criteria(30);
  typedef OR::ManagedSimpleAdaptativeStrategy<double, OR::C1Manifold,
                                              CurManagedDesc>
    Strategy;
  auto strategy = Strategy(managed_desc, stop_criteria, manifold);

  auto search = OR::AdaptativeSearcher<double, Strategy>(
    [](const double &a) { return std::abs(a - PI / 4); }, strategy);
  double res = search.find_minimum().point;
  OPA_DISP("found result ", res);
}

TEST(AdaptativeSearch, Test2D) {
  SphereManifold manifold;
  typedef OR::NManifold<Pos2> BaseManifold;
  typedef OR::ConstNDimManagedDesc<Pos2, BaseManifold> CurManagedDesc;
  CurManagedDesc managed_desc(manifold, 10, 10);
  OR::RoundStoppingCriteria stop_criteria(30);

  typedef OR::ManagedSimpleAdaptativeStrategy<Pos2, BaseManifold,
                                              CurManagedDesc>
    Strategy;
  auto strategy = Strategy(managed_desc, stop_criteria, manifold);

  auto search = OR::AdaptativeSearcher<Pos2, Strategy>(
    [](const Pos2 &a) {
      Pos2 tmp = Pos2(cos(a.x) * sin(a.y), cos(a.y));
      OPA_DISP0(tmp, a, glm::length2(tmp - Pos2(1, 1)));
      return glm::length2(tmp - Pos2(1, 1));
    },
    strategy);
  Pos2 res = search.find_minimum().point;
  OPA_DISP("found result ", res);
}

TEST(AdaptativeSearch, Test2D_2) {
  SphereManifold manifold;
  typedef OR::NManifold<Pos2> BaseManifold;
  typedef OR::ConstNDimManagedDesc<Pos2, BaseManifold> CurManagedDesc;
  CurManagedDesc managed_desc(manifold, 10, 10, Pos2(), { 0 });
  OR::RoundStoppingCriteria stop_criteria(30);

  typedef OR::ManagedSimpleAdaptativeStrategy<Pos2, BaseManifold,
                                              CurManagedDesc>
    Strategy;
  auto strategy = Strategy(managed_desc, stop_criteria, manifold);

  auto search = OR::AdaptativeSearcher<Pos2, Strategy>(
    [](const Pos2 &a) {
      Pos2 tmp = Pos2(cos(a.x) * cos(a.y), sin(a.y));
      OPA_DISP0(tmp, a, glm::length2(tmp - Pos2(1, 1)));
      return glm::length2(tmp - Pos2(1, 1));
    },
    strategy);
  Pos2 res = search.find_minimum().point;
  OPA_DISP("found result ", res);
}

TEST(PointMatcher, Test1) {
  double eps = 1e-7;
  PointVec to_match = {
    vec_x,
    vec_y,
    vec_z,
    vec_x * (1 + eps),
    vec_y,
    vec_x + vec_z,
    vec_z * (1 - eps),
  };

  PointMatcher matcher;
  OPA_CHECK0(matcher.load_batch(to_match));
  std::vector<int> rmp;
  for (auto &pt : to_match) {
    rmp.emplace_back(matcher.query(pt));
  }

  algo::UnionJoin uj(to_match.size());
  uj.merge(0, 3);
  OPA_TRACES(uj.get_clusters());
  uj.merge(1, 4);
  OPA_TRACES(uj.get_clusters());
  uj.merge(2, 6);

  OPA_TRACES(to_match);
  OPA_TRACES(uj.get_clusters());
  OPA_TRACES(rmp);
  OPA_CHECK0(uj == algo::UnionJoin::FromVec(rmp));
}

TEST(FaceOrientation, Test1) {
  PointVec tb = { { 0.000000, -89.549263, -170.910568 },
                  { 0.000000, -79.679008, -164.627838 },
                  { 0.000000, -66.848038, -166.458832 },
                  { -0.000000, -56.284607, -171.097809 },
                  { 0.000000, -49.359703, -175.914886 },
                  { 0.000000, -47.715630, -179.802826 },
                  { 0.000000, -44.085541, -197.595215 },
                  { 0.000000, -46.719570, -205.607208 },
                  { 0.000000, -53.499737, -215.031876 },
                  { 0.000000, -69.314308, -223.780746 },
                  { 0.000000, -71.586792, -224.213013 },
                  { 0.000000, -78.406769, -222.738464 },
                  { 0.000000, -90.256172, -219.353745 },
                  { 0.000000, -92.253845, -217.648621 },
                  { 0.000000, -101.535362, -200.770020 },
                  { 0.000000, -101.534790, -187.691223 },
                  { 0.000000, -97.825287, -177.743454 } };

  tb = { { 0.000000, 3.039902, 2.596634 },
         { 0.000000, 1.932470, 6.879122 },
         { 0.000000, 0.800208, 8.872557 } };
  tb = { { 0.000000, 12.164412, 5.016092 },
         { 0.000000, 9.401305, 16.223434 },
         { 0.000000, 0.800208, 8.872557 },
         { 0.000000, 1.932470, 6.879122 },
         { 0.000000, 3.039902, 2.596634 } };
  //  ../opa/math/game/src/mesh_util.cpp:mesh_cut_plan:980,
  //  msg=tr=((0.000000,9.401305,16.223434),(0.000000,0.800208,8.872557),(0.000000,3.039902,2.596634)),face_normal(tr)=(1.000000,0.000000,-0.000000),
  //  ../opa/math/game/src/mesh_util.cpp:mesh_cut_plan:980,
  //  msg=tr=((0.000000,12.164412,5.016092),(0.000000,9.401305,16.223434),(0.000000,3.039902,2.596634)),face_normal(tr)=(1.000000,0.000000,0.000000),

  // PointVec res = tr_to_vec(realign_tr(vec_to_tr(tb), -vec_x));
  PointVec res2 = realign_face(tb, -vec_x);
  OPA_DISP0(res2, face_normal(res2), tb, face_normal(tb));
}

TEST(FaceGraph, Test1) {

  BoxAASpec box;
  box.update(Pos(0));
  box.update(Pos(1));
  auto fc = FaceCollection();
  add_icosahedron(&fc);
  auto mesh = fc.to_mesh();
  OPA_DISP0(mesh->str());

  auto face_graph = mesh->compute_face_graph();
  lemon::Dijkstra<lemon::SmartGraph, lemon::SmartGraph::ArcMap<double> > djk(
    *face_graph->graph, *face_graph->edge_map);
  lemon::SmartGraph::NodeMap<double> nm(*face_graph->graph);
  djk.distMap(nm);
  djk.init();

  djk.addSource(face_graph->graph->nodeFromId(0));
  djk.start();
  for (lemon::SmartGraph::NodeIt n(*face_graph->graph); n != lemon::INVALID;
       ++n) {
    OPA_DISP0(face_graph->graph->id(n), nm[n]);
  }
  djk.init();
  djk.addSource(face_graph->graph->nodeFromId(10));
  djk.start();
  for (lemon::SmartGraph::NodeIt n(*face_graph->graph); n != lemon::INVALID;
       ++n) {
    OPA_DISP0(face_graph->graph->id(n), nm[n]);
  }
}

TEST(FaceGraph, Test2) {
  BoxAASpec box;
  box.update(Pos(0));
  box.update(Pos(1));
  auto fc = FaceCollection();
  add_icosahedron(&fc);
  auto mesh = fc.to_mesh();
  OPA_DISP0(mesh->str());

  auto face_graph = mesh->compute_face_graph2();
  algo::FastGraph::DijkstraCtx ctx;
  ctx.edge_costs = face_graph->edge_costs;

  ctx.sources = { 0 };
  face_graph->graph->dijkstra(ctx);
  OPA_DISP0(ctx.dists);

  ctx.sources = { 10 };
  face_graph->graph->dijkstra(ctx);
  OPA_DISP0(ctx.dists);
}

TEST(BoxRot, Test1) {
  REP (i, 10) {
    PointVec pv;
    REP (k, 100)
      pv.push_back(vec_rand_uni());
    BoxSpec box = compute_best_box(pv);
    Rot rot = rot_to_box_space(box);
    BoxAASpec aabox = BoxAASpec::FromBoxSpace(box);
    for (auto &x : pv)
      OPA_CHECK(aabox.in(rot * x), rot * x, aabox.low, aabox.high);
    OPA_DISP0(aabox.area(), box.area());
  }
}

TEST(Atomize, Test1) {
  auto mesh =
    SPTR(Mesh)(FaceCollection(true)
                 .load_stl("/home/benoit/programmation/robotics/test.stl")
                 .to_mesh());
  auto res = atom::atomize_mesh_K(*mesh, 2);

  for (auto &x : res) {
    OPA_DISP0(x.str());
  }
}

TEST(Atomize, Test2) {
  auto mesh = SPTR(Mesh)(
    FaceCollection(true)
      .load_stl("/home/benoit/programmation/robotics/obj1.clean2.stl")
      .to_mesh());
  mesh->whiten();
  auto cmesh = mesh->split_connected();
  int id = 0;
  REP (i, cmesh.size()) {
    if (cmesh[id]->vertices().size() < cmesh[i]->vertices().size()) id = i;
  }
  auto curmesh = cmesh[id];
  KAtomizer k(*curmesh, 5);
  auto asst = k.initial_assignment();

  auto search_params = KAtomizer_SearchParams();
  auto ts = TimeStoppingCriteria(int(20e6));
  search_params.stop_criteria = &ts;
  search_params.T0 = 0.5;
  search_params.state = 0;
  k.improve_assignment(&asst, search_params);

  for (auto &x : asst.cluster_data) {
    OPA_DISP0(x.box.str());
  }
}

TEST(QuatLookAt, Test1) {
  Pos p1 = glm::normalize(Pos{ -0.229824, 0.688104, 0.688254 });
  Pos p2 = vec_make_ortho_safe({ 0.688472, 0.614790, -0.384758 }, p1);
  Rot rot = quat_look_at_safe(p1, p2);
  Rot irot = glm::inverse(rot);
  OPA_DISP0(p1, p2, rot);
  OPA_TRACES(irot * p1);
  OPA_TRACES(irot * p2);
  OPA_TRACES(p1 * irot);
  OPA_TRACES(p2 * irot);

  OPA_TRACES(rot * p1);
  OPA_TRACES(rot * p2);
  OPA_TRACES(rot * vec_x);
  OPA_TRACES(rot * vec_z);
}

TEST(LineLine, Test1) {
  opa::Pos2 p1, p2;

  OPA_CHECK0(line_line_intersection_2d(
               LineSpec2D{ opa::Pos2(0.8, -0.55), opa::Pos2(1.2, 0.4) },
               LineSpec2D{ opa::Pos2(-1.35, -0.23), opa::Pos2(1.42, 0.22) }, p1,
               p2) == 1);
  OPA_DISP0(p1, p2);

  OPA_CHECK0(line_line_intersection_2d(
               LineSpec2D{ opa::Pos2(1, 2), opa::Pos2(2, 3) },
               LineSpec2D{ opa::Pos2(1.5, 2.5), opa::Pos2(2.5, 3.5) }, p1,
               p2) == 2);
  OPA_DISP0(p1, p2);

  OPA_CHECK0(line_line_intersection_2d(
               LineSpec2D{ opa::Pos2(1, 2), opa::Pos2(2, 3) },
               LineSpec2D{ opa::Pos2(1, 2), opa::Pos2(0.5, 1.5) }, p1,
               p2) == 1);
  OPA_DISP0(p1, p2);

  OPA_CHECK0(
    line_line_intersection_2d(
      LineSpec2D{ Pos2{ 0.809628, -0.553680 }, Pos2{ 1.201028, 0.411268 } },
      LineSpec2D{ Pos2{ 0.265543, -1.394347 }, Pos2{ 1.429934, 0.232532 } }, p1,
      p2) == 0);
  OPA_DISP0(p1, p2);

  OPA_CHECK0(
    line_line_intersection_2d(
      LineSpec2D{ Pos2{ 0.809628, -0.553680 }, Pos2{ 1.201028, 0.411268 } },
      LineSpec2D{ Pos2{ -1.362726, -0.231046 }, Pos2{ 1.429934, 0.232532 } },
      p1, p2) == 1);
  OPA_DISP0(p1, p2);

  OPA_CHECK0(
    line_line_intersection_2d(
      LineSpec2D{ Pos2{ 0.809628, -0.553680 }, Pos2{ 1.201028, 0.411268 } },
      LineSpec2D{ Pos2{ 1.429934, 0.232532 }, Pos2{ -1.362726, -0.231046 } },
      p1, p2) == 1);
  OPA_DISP0(p1, p2);
}

TEST(LineLine, Test2) {
  opa::Pos2 p1, p2;

  OPA_CHECK0(line_line_intersection_2d_v2(
               LineSpec2D{ opa::Pos2(0.8, -0.55), opa::Pos2(1.2, 0.4) },
               LineSpec2D{ opa::Pos2(-1.35, -0.23), opa::Pos2(1.42, 0.22) }, p1,
               p2) == 1);
  OPA_DISP0(p1, p2);

  OPA_CHECK0(line_line_intersection_2d_v2(
               LineSpec2D{ opa::Pos2(1, 2), opa::Pos2(2, 3) },
               LineSpec2D{ opa::Pos2(1.5, 2.5), opa::Pos2(2.5, 3.5) }, p1,
               p2) == 2);
  OPA_DISP0(p1, p2);

  OPA_CHECK0(line_line_intersection_2d_v2(
               LineSpec2D{ opa::Pos2(1, 2), opa::Pos2(2, 3) },
               LineSpec2D{ opa::Pos2(1, 2), opa::Pos2(0.5, 1.5) }, p1,
               p2) == 1);
  OPA_DISP0(p1, p2);

  OPA_CHECK0(
    line_line_intersection_2d_v2(
      LineSpec2D{ Pos2{ 0.809628, -0.553680 }, Pos2{ 1.201028, 0.411268 } },
      LineSpec2D{ Pos2{ 0.265543, -1.394347 }, Pos2{ 1.429934, 0.232532 } }, p1,
      p2) == 0);
  OPA_DISP0(p1, p2);

  OPA_CHECK0(
    line_line_intersection_2d_v2(
      LineSpec2D{ Pos2{ 0.809628, -0.553680 }, Pos2{ 1.201028, 0.411268 } },
      LineSpec2D{ Pos2{ -1.362726, -0.231046 }, Pos2{ 1.429934, 0.232532 } },
      p1, p2) == 1);
  OPA_DISP0(p1, p2);

  OPA_CHECK0(
    line_line_intersection_2d_v2(
      LineSpec2D{ Pos2{ 0.809628, -0.553680 }, Pos2{ 1.201028, 0.411268 } },
      LineSpec2D{ Pos2{ 1.429934, 0.232532 }, Pos2{ -1.362726, -0.231046 } },
      p1, p2) == 1);
  OPA_DISP0(p1, p2);
}

TEST(TrTrIntersection, Test1) {
  Tr2D tr0{
    Pos2{ 0.033604, 0.000743 },
    Pos2{ 0.809628, -0.553680 },
    Pos2{ 1.201028, 0.411268 },
  };
  Tr2D tr1{
    Pos2{ -1.362726, -0.231046 },
    Pos2{ 0.265543, -1.394347 },
    Pos2{ 1.429934, 0.232532 },
  };
  auto res = tr_tr_intersection(tr0, tr1);
  OPA_DISP0(res);
}

GTEST_API_ int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  opa::init::opa_init(argc, argv);
  return RUN_ALL_TESTS();
}

#include <gtest/gtest.h>

#include <fplll/fplll.h>
#include <glib/strings/strutil.h>
#include <opa/math/co/barvinok.h>
#include <opa/math/co/base.h>
#include <opa/math/common/matrix_utils.h>
#include <opa/math/common/utils_num.h>
#include <opa/math/game/proto/common.pb.h>
#include <opa/or/grid_search.h>
#include <opa/utils/misc.h>

using namespace opa;
DEFINE_int32(maxd, 2, "");
DEFINE_string(dim_list, "2", "");
DEFINE_string(dump_file, "/tmp/debug.out", "");
DEFINE_int32(nrepeat, 1, "");
DEFINE_int32(max_depth, 1, "");
DEFINE_bool(debug_check, false, "");
DEFINE_bool(online_mode, true, "");
DEFINE_bool(debug_check_more, false, "");
DEFINE_bool(check_answers, true, "");
DEFINE_bool(simple, false, "");
DEFINE_bool(check_bruteforce, true, "");
std::vector<u32> dim_list;

using namespace opa::math::common;
using namespace opa::math::co;
using namespace opa::math::game;
using namespace std;

void init_test() {
  auto res = glib::Split(FLAGS_dim_list, ",");
  dim_list = utils::transform_container(res, &utils::Conv::from_str2<u32>);
}
Barvinok::Configuration build_conf() {
  Barvinok::Configuration conf;
  conf.enumerate_ub = int(10);
  conf.max_depth = FLAGS_max_depth;
  conf.online_mode = FLAGS_online_mode;
  return conf;
}

std::vector<ZModule_t> gen_crosspoly_vertices(int n) {
  std::vector<ZModule_t> vertices;
  REP (i, n) {
    ZModule_t coords(n);
    coords[i] = 1;
    vertices.push_back(coords);
    coords[i] = -1;
    vertices.push_back(coords);
  }
  return vertices;
}

std::vector<ZModule_t> gen_ncube_vertices(int n) {
  std::vector<ZModule_t> vertices;
  REP (i, 1 << n) {
    std::vector<bignum> coords;
    REP (j, n) { coords.push_back(i >> j & 1); }
    vertices.emplace_back(coords);
  }
  return vertices;
}

Polyhedra_VertexRepr
gen_poly_from_zvertices(const std::vector<ZModule_t> &vertices) {
  std::vector<QModule_t> qvertices = lift_zmodules(vertices);
  Polyhedra_VertexRepr polyhedra;
  polyhedra.init(qvertices);
  return polyhedra;
}

Polyhedra_VertexRepr
gen_poly_from_vertices(const std::vector<QModule_t> &vertices) {
  Polyhedra_VertexRepr polyhedra;
  polyhedra.init(vertices);
  return polyhedra;
}

Z get_bruteforce(const Polyhedra_VertexRepr &polyhedra,
                 const std::vector<Hyperplane> *on_planes = nullptr) {
  int d = polyhedra.d;
  std::vector<utils::MinFinder<Z> > minv(d);
  std::vector<utils::MaxFinder<Z> > maxv(d);

  for (auto &v : polyhedra.vertices) {
    REP (i, d)
      minv[i].update(QF.floor(v[i])), maxv[i].update(QF.ceil(v[i]));
  }
  std::vector<std::vector<s32> > ranges;
  REP (i, d) {
    OPA_DISP("Bruteforce ", i, minv[i].get(), maxv[i].get());
    ranges.push_back(utils::Range<s32>::StepRange(minv[i].get().gets32(),
                                                  maxv[i].get().gets32() + 1, 1)
                       .tb());
  }

  int count = 0;
  OR::CrossProdGen<s32>(ranges, [&](const std::vector<s32> &v) {
    ZModule_t p(utils::transform_container(v, &bignum::froms32));
    if (on_planes != nullptr &&
        !std::all_of(ALL(*on_planes),
                     std::bind(&Hyperplane::is_on, std::placeholders::_1,
                               lift_zmodule(p))))
      return;
    if (polyhedra.contains(lift_zmodule(p))){
    OPA_DISP("Bruteforce point ", p);
    ++count;
    }
  });
  return count;
}

void set_mesh_face_from_points(proto::MeshFace *mf,
                               const std::vector<ZModule_t> &points) {
  for (auto &pt : points) {
    auto pos = mf->add_vertex();
    pos->set_x(pt[0].getfloat().to_double());
    pos->set_y(pt[1].getfloat().to_double());
    pos->set_z(pt[2].getfloat().to_double());
  }
}

void dump_descent(proto::MeshList &ml, const Barvinok::DescentData &data,
                  const string &desc) {
  int dim = 3;
  if (data.children.size() == 0) {
    OPA_DISP0(data.rd.cone.rays, data.rd.sgn, data.rd.exclude_face);
    REP (i, dim) {
      proto::Mesh *single_mesh = ml.add_mesh();
      single_mesh->set_name(
        utils::stdsprintf("%s-%d", desc.c_str(), data.rd.exclude_face[i]));
      set_mesh_face_from_points(single_mesh->add_face(),
                                { ZModule_t(dim),
                                  data.rd.cone.rays[(i + 1) % dim],
                                  data.rd.cone.rays[(i + 2) % dim] });
    }
  } else {
    proto::Mesh *single_mesh = ml.add_mesh();
    single_mesh->set_name("split_vec");
    OPA_DISP0(data.rd.cone.rays, data.rd.sgn, data.rd.exclude_face, data.w,
              data.hp_status);
    set_mesh_face_from_points(single_mesh->add_face(),
                              { ZModule_t(dim), data.w });
    for (auto &child : data.children) dump_descent(ml, child, desc);
  }
}

void check_steps(int dim, Barvinok &left, Barvinok &right,
                 const Barvinok::DescentData &data, int max_depth_l) {
  bool is_leaf = data.children.size() == 0;
  for (auto &child : data.children) {
    check_steps(dim, left, right, child, max_depth_l - 1);
  }

  left.reset_residue_helper(dim);
  right.reset_residue_helper(dim);
  left.conf.max_try = max_depth_l;

  Barvinok::DescentData dl, dr;
  left.rec(data.rd, dl);
  right.rec(data.rd, dr);
  OPA_CHECK0(dl.ok);
  OPA_CHECK0(dr.ok);

  OPA_DISP("CHECK STEP >> ", left.residue_helper.cone_ctx.resv,
           right.residue_helper.cone_ctx.resv, is_leaf);
  OPA_DISP("CHECK STEP2 >> ", left.residue_helper.res,
           right.residue_helper.res);

  if (left.residue_helper.cone_ctx.resv != right.residue_helper.cone_ctx.resv) {

    if (dim == 3) {
      OPA_DISP0(data.hp_status, data.w);
      proto::MeshList ml;
      dump_descent(ml, dr, "orig");
      dump_descent(ml, dl, "split");
      std::ofstream tmpfile(FLAGS_dump_file, std::ios::binary);
      ml.SerializeToOstream(&tmpfile);
    }

    OPA_DISP("failing on reduction of ", data.rd,
             left.residue_helper.cone_ctx.resv,
             right.residue_helper.cone_ctx.resv);
    OPA_CHECK0(!is_leaf);

    // Should be : bitstring(bit_exclude, 9)=
    // 011 000 100,bitstring(bit_sgn,
    // 3)=011,
    // L0: 101
    // L1: 001 001 010
    if (FLAGS_debug_check_more) {

      REP (bit_exclude, 1 << (dim * dim)) {
        REP (bit_sgn, 1 << dim) {
          left.reset_residue_helper(dim);
          REP (j, dim) {
            Barvinok::RecData nrec = dl.children[j].rd;
            REP (k, dim) {
              nrec.exclude_face[k] = (bit_exclude >> (j * dim + k) & 1);
            }
            nrec.sgn = (bit_sgn >> j & 1) ? 1 : -1;
            Barvinok::DescentData tmp;
            left.rec(nrec, tmp);
          }
          if (left.residue_helper.cone_ctx.resv ==
              right.residue_helper.cone_ctx.resv) {
            OPA_DISP("Should be ", bitstring(bit_exclude, 9),
                     bitstring(bit_sgn, 3));
          }
        }
      }
    }
    OPA_CHECK0(false);
  }
}

void check_poly(const Polytope_ConeRepr &cone) {
  Barvinok::Configuration conf_reduce;
  conf_reduce.enumerate_ub = int(10);
  conf_reduce.max_try = 1;
  conf_reduce.online_mode = FLAGS_online_mode;
  conf_reduce.max_depth = FLAGS_max_depth;
  conf_reduce.fixed_lambda = true;
  conf_reduce.debug = true;
  conf_reduce.lambda = generate_lambda_rand(cone.cones[0].vertex.size());

  Barvinok::Configuration conf_single = conf_reduce;
  conf_single.max_depth = 0;

  for (auto &x : cone.cones) {
    int dim = x.vertex.size();

    Barvinok barvinok_reduce;
    barvinok_reduce.init(conf_reduce);
    auto steps = barvinok_reduce.prepare_cone(x);

    Barvinok barvinok_single;
    barvinok_single.init(conf_single);

    barvinok_single.reset_residue_helper(dim);
    barvinok_reduce.reset_residue_helper(dim);
    barvinok_reduce.add_cone(x);
    barvinok_single.add_cone(x);
    if (barvinok_single.residue_helper.cone_ctx.resv !=
        barvinok_reduce.residue_helper.cone_ctx.resv) {
      Barvinok::DescentData dleft = barvinok_reduce.descent_data;
      Barvinok::DescentData dright = barvinok_single.descent_data;

      check_steps(x.vertex.size(), barvinok_reduce, barvinok_single, dleft,
                  conf_reduce.max_depth);
    }
  }
}

TEST(SimplicialCone, Parallelepiped) {
  int n = 3, m = 2;
  Matrix<Z> m1(&Ring_Z, n, m);
  m1.setCol(0, { 8, 4, 8 });
  m1.setCol(1, { 4, 8, 4 });

  SimplicialCone cone;
  cone.init(m1);
  OPA_DISP0(cone.list_parallelepiped());
}

TEST(SimplicialCone, Parallelepiped2) {
  int n = 2, m = 2;
  Matrix<Z> m1(&Ring_Z, n, m);
  m1.setCol(0, { 1, 5 });
  m1.setCol(1, { 1, 0 });
  QModule_t pt = QModule_t::zero(n);

  SimplicialCone cone;
  cone.init(m1);
  OPA_DISP0(cone.list_parallelepiped(pt, { false, false }));
  OPA_DISP0(cone.list_parallelepiped(pt, { true, false }));
  OPA_DISP0(cone.list_parallelepiped(pt, { true, true }));
  OPA_DISP0(cone.list_parallelepiped(pt, { false, true }));
}

TEST(SimplicialCone, Split) {
  int n = 2, m = 2;
  Matrix<Z> m1(&Ring_Z, n, m);
  m1.setCol(0, { 1, 5 });
  m1.setCol(1, { 0, 1 });
  SimplicialCone cone;
  cone.init(m1);
  auto res = cone.find_w();
  OPA_DISP0(res);
}

TEST(BasicMat, SNF) {
  int n = 3, m = 2;
  Matrix<bignum> m1(&Ring_Z, n, m);
  m1.setCol(0, { 8, 4, 8 });
  m1.setCol(1, { 4, 8, 4 });
  Matrix<bignum> p, d, q, ip, iq;
  snf(m1, &d, &p, &q);
  auto res = q * m1 * p;
  OPA_DISP0(p);
  OPA_DISP0(d);
  OPA_DISP0(q);
  OPA_DISP0(res);

  Matrix<bignum> lat_hnf;
  hnf(m1, &lat_hnf);

  OPA_CHECK(p.invert(&ip), p);
  OPA_CHECK(q.invert(&iq), q);
  OPA_DISP0(ip, iq, lat_hnf);
  std::vector<std::vector<bignum> > ranges;
  REP (i, n) {
    if (i < m) {
      ranges.push_back(
        utils::Range<bignum>::StepRange(0, d.get(i, i).abs(), 1).tb());
    } else
      ranges.push_back({ 0 });
  }

  auto v1 = m1.get_col(0).tovec();
  auto v2 = m1.get_col(1).tovec();

  OPA_DISP0(Ring_Z.dot(v1, v1), Ring_Z.dot(v2, v2));
  OR::CrossProdGen<bignum>(ranges, [&](const std::vector<bignum> &v) {
    auto tmp = iq.eval(v);
    std::vector<bignum> tmp_reduced = tmp;
    reduce_hnf_vec(lat_hnf, tmp_reduced, nullptr, true);
    OPA_DISP0(v, tmp, tmp_reduced, Ring_Z.dot(v1, tmp_reduced),
              Ring_Z.dot(v2, tmp_reduced));
  });
}

TEST(BarvinokRed, LLL1) {
  int n = 4;
  std::vector<ZModule_t> rays;
  rays.emplace_back(std::vector<bignum>{ -0x2f, -0x5e, -0x5e, -0x5e });
  rays.emplace_back(std::vector<bignum>{ -0x37, 0, 0, 0 });
  rays.emplace_back(std::vector<bignum>{ 0, -0x81, 0, 0 });
  rays.emplace_back(std::vector<bignum>{ 0, 0, 0x4d, -0x35 });
  SimplicialCone cone;
  cone.init(rays);

  auto redmat = cone.find_nw_matrix();
  OPA_DISP0(redmat);
  OPA_DISP0(cone.parallelepiped_count);

  REP (i, n) {
    ZModule_t repl_l(redmat.get_col(i));
    ZModule_t repl(cone.mat.eval(repl_l));
    gcd_simplify(Ring_Z, repl.elems);
    if (cone.is_on_cone_boundary(repl)) continue;
    OPA_DISP("Replacing ", i, repl);

    REP (j, n) {
      Matrix<Z> nmat = cone.mat.clone();
      nmat.set_col(j, repl);
      OPA_DISP0(j, hnf_det(nmat));
    }
  }
}

TEST(BarvinokRed, LLL2) {
  int n = 4;
  std::vector<std::vector<bignum> > data;
  data.push_back(std::vector<bignum>{ -0x2f, -0x5e, -0x5e, -0x5e });
  data.push_back(std::vector<bignum>{ -0x37, 0, 0, 0 });
  data.push_back(std::vector<bignum>{ 0, -0x81, 0, 0 });
  data.push_back(std::vector<bignum>{ 0, 0, 0x4d, -0x35 });
  Matrix<Z> mat(&Ring_Z, n, n);
  mat.set_rows(data);
  SimplicialCone cone;
  cone.init(mat);

  auto redmat = cone.find_nw_matrix();
  OPA_DISP0(redmat);
  OPA_DISP0(cone.parallelepiped_count);

  REP (i, n) {
    ZModule_t repl_l(redmat.get_col(i));
    ZModule_t repl(cone.mat.eval(repl_l));
    OPA_DISP("Replacing ", i, repl);
    gcd_simplify(Ring_Z, repl.elems);
    // if (cone.is_on_cone_boundary(repl)) continue;
    OPA_DISP("Replacing ", i, repl);

    REP (j, n) {
      Matrix<Z> nmat = cone.mat.clone();
      nmat.set_col(j, repl);
      OPA_DISP0(j, hnf_det(nmat));
    }
  }
}

TEST(Polyhedra, NCrossPoly) {
  int maxn = FLAGS_maxd;
  FOR (n, 2, maxn + 1) {
    Polyhedra_VertexRepr polyhedra =
      gen_poly_from_zvertices(gen_crosspoly_vertices(n));
    OPA_DISP0(polyhedra.faces);
    OPA_DISP0(polyhedra.cone_repr);
  }
}

TEST(Polyhedra, NCube) {
  int maxn = FLAGS_maxd;
  FOR (n, 2, maxn + 1) {
    Polyhedra_VertexRepr polyhedra =
      gen_poly_from_zvertices(gen_ncube_vertices(n));
    OPA_DISP0(polyhedra.faces);
    OPA_DISP0(polyhedra.cone_repr);
  }
}

TEST(Barvinok, Cube) {
  int maxn = FLAGS_maxd;
  Barvinok::Configuration conf;
  conf.enumerate_ub = int(1e5);
  conf.max_try = 1;
  conf.online_mode = true;

  FOR (n, 2, maxn + 1) {
    Polyhedra_VertexRepr polyhedra =
      gen_poly_from_zvertices(gen_ncube_vertices(n));

    Barvinok barvinok;
    barvinok.init(conf);
    Z res = barvinok.solve(polyhedra.cone_repr);
    OPA_DISP0(res);
    if (FLAGS_check_bruteforce) {
      OPA_CHECK_EQ0(res, get_bruteforce(polyhedra));
    }
  }
}

Z handle_poly(const Polyhedra_VertexRepr &polyhedra,
              const Barvinok::Configuration &conf) {
  Barvinok barvinok;
  barvinok.init(conf);
  auto res = barvinok.try_solve(polyhedra.cone_repr);
  if (FLAGS_debug_check) {
    check_poly(polyhedra.cone_repr);
  }
  if (!res.ok) {
    OPA_DISP("RESULT KAPPA ", res);
    OPA_CHECK0(false);

  } else {
    OPA_DISP0(res.val);
    if (FLAGS_check_bruteforce) {
      OPA_CHECK_EQ0(res.val, get_bruteforce(polyhedra));
    }
    OPA_DISP0(barvinok.descent_data.transitions);
  }
  return res.val;
}

Z handle_poly_with_planes(const Polyhedra_VertexRepr &polyhedra,
                          const std::vector<Hyperplane> &planes,
                          const Barvinok::Configuration &conf) {

  Polyhedra_VertexRepr target = compute_poly_proj(polyhedra, planes);
  OPA_DISP0(target.vertices);
  Z res = handle_poly(target, conf);
  Z check_res = get_bruteforce(polyhedra, &planes);
  OPA_CHECK_EQ0(res, check_res);
  OPA_DISP0(res);
  return res;
}

TEST(Barvinok, CrossPoly) {
  Barvinok::Configuration conf;
  conf.enumerate_ub = int(10);
  conf.max_try = 1;
  conf.online_mode = FLAGS_online_mode;
  conf.max_depth = FLAGS_max_depth;
  // conf.fixed_lambda = true;
  // conf.lambda =
  //  ZModule_t(utils::Range<bignum>::StepRange(1, n + 1, 1).tb());

  for (auto d : dim_list) {
    REP (checkstep, FLAGS_nrepeat) {
      std::vector<QModule_t> crosspoly_vertices =
        lift_zmodules(gen_crosspoly_vertices(d));

      if (!FLAGS_simple) {
        for (auto &x : crosspoly_vertices) x = x * QF.import(rng() % 8 + 1, 10);
      }
      Polyhedra_VertexRepr polyhedra =
        gen_poly_from_vertices(crosspoly_vertices);

      handle_poly(polyhedra, conf);
    }
  }
}

TEST(Barvinok, CrossPolyShift) {
  Barvinok::Configuration conf;
  conf.enumerate_ub = int(10);
  conf.max_depth = FLAGS_max_depth;
  conf.online_mode = FLAGS_online_mode;
  // conf.max_try = 1;
  // res=(-4ee / 91b),

  for (auto d : dim_list) {
    REP (checkstep, FLAGS_nrepeat) {
      std::vector<QModule_t> crosspoly_vertices =
        lift_zmodules(gen_crosspoly_vertices(d));

      if (!FLAGS_simple) {
        for (auto &x : crosspoly_vertices)
          x = x * QF.import(rng() % 10 + 5, 10);
      }
      Polyhedra_VertexRepr polyhedra =
        gen_poly_from_vertices(crosspoly_vertices);
      handle_poly(polyhedra, conf);
    }
  }
}

TEST(Barvinok, Test1) {
  Barvinok::Configuration conf;
  conf.enumerate_ub = int(1e5);
  int n = 2;

  std::vector<ZModule_t> vertices;
  vertices.emplace_back(std::vector<bignum>{ 1, 2 });
  vertices.emplace_back(std::vector<bignum>{ 4, 7 });
  vertices.emplace_back(std::vector<bignum>{ 7, 7 });
  vertices.emplace_back(std::vector<bignum>{ 6, 4 });
  std::vector<QModule_t> qvertices = lift_zmodules(vertices);
  Polyhedra_VertexRepr polyhedra;
  polyhedra.init(qvertices);
  OPA_DISP0(polyhedra.faces);
}

TEST(Barvinok, LowDim0) {
  auto conf = build_conf();

  std::vector<Hyperplane> proj_planes;
  proj_planes.emplace_back(lift_zmodule(std::vector<bignum>{ 1, 0, 1 }),
                           QF.import(2));

  std::vector<QModule_t> crosspoly_vertices =
    lift_zmodules(gen_crosspoly_vertices(3));

  for (auto &x : crosspoly_vertices) x = x * QF.import(3);
  Polyhedra_VertexRepr polyhedra = gen_poly_from_vertices(crosspoly_vertices);
  handle_poly_with_planes(polyhedra, proj_planes, conf);
}

TEST(Barvinok, LowDim1) {
  auto conf = build_conf();

  std::vector<Hyperplane> proj_planes;
  proj_planes.emplace_back(
    lift_zmodule(std::vector<bignum>{ -0x34, 0x3d, -0x1 }), QF.getE());

  std::vector<QModule_t> crosspoly_vertices =
    lift_zmodules(gen_crosspoly_vertices(3));

  for (auto &x : crosspoly_vertices) x = x * QF.import(40);
  Polyhedra_VertexRepr polyhedra = gen_poly_from_vertices(crosspoly_vertices);
  handle_poly_with_planes(polyhedra, proj_planes, conf);
}

TEST(Barvinok, LowDim2) {
  auto conf = build_conf();

  int dim = 5;
  std::vector<Hyperplane> proj_planes;
  proj_planes.emplace_back(
    std::vector<Q>{ QF.import(0x1, 0x2), QF.import(0x13, 0x4),
                    QF.import(-0xd, 0x2), QF.import(0x1), QF.import(0x0) },
    QF.getZ());
  proj_planes.emplace_back(std::vector<Q>{ QF.import(0x0), QF.import(0x3, 0x2),
                                           QF.import(0x0), QF.import(0x0),
                                           QF.import(0x1) },
                           QF.getZ());
  // proj_planes.emplace_back(std::vector<Q>{QF.import(0x1 , 0x2),
  // QF.import(-0x9 , 0x4), QF.import(0x1 , 0x2), QF.import(0x0),
  // QF.import(0x0), QF.import(0x1)}, QF.getZ());
  // proj_planes.emplace_back(std::vector<Q>{QF.import(-0x1 , 0x2),
  // QF.import(-0x7 , 0x4), QF.import(0xb , 0x2), QF.import(0x0),
  // QF.import(0x0), QF.import(0x0), QF.import(0x1)}, QF.getZ());

  std::vector<QModule_t> crosspoly_vertices =
    lift_zmodules(gen_crosspoly_vertices(dim));

  for (auto &x : crosspoly_vertices) x = x * QF.import(40);
  Polyhedra_VertexRepr polyhedra = gen_poly_from_vertices(crosspoly_vertices);
  handle_poly_with_planes(polyhedra, proj_planes, conf);
}

TEST(Lattice, Gen2) {
  int n = 5;
  auto mat = generate_lattice(n, n, 100, 100);
  std::vector<Z> tb = mat.get_col(0).tovec();

  auto m2 = expand_unimodular(&Ring_Z, tb, 0);
  OPA_DISP0(m2);
  OPA_DISP0(hnf_det(m2));
  OPA_DISP0(m2.evalT(tb));
  OPA_DISP0(m2.eval(tb));
  auto im2 = m2.inverse();
  OPA_DISP0(im2.evalT(tb));
  OPA_DISP0(im2.eval(tb));
  OPA_DISP0(m2);
  OPA_DISP0(im2);
  OPA_DISP0(tb);
}

/*
Barvinok::BarvinokPolytope polytope;

{
  Barvinok::BarvinokCone cone;
  cone.vertex = lift_zmodule(ZModule_t(std::vector<bignum>{ 1, 2 }));
  cone.rays.push_back(lift_zmodule(ZModule_t(std::vector<bignum>{ 5, 2 }));
  cone.rays.push_back(lift_zmodule(ZModule_t(std::vector<bignum>{ 3, 5 }));
  polytope.cones.push_back(cone);
}
{
  Barvinok::BarvinokCone cone;
  cone.vertex = lift_zmodule(ZModule_t(std::vector<bignum>{ 4, 7 }));
  cone.rays.push_back(lift_zmodule(ZModule_t(std::vector<bignum>{ -3, -5}));
  cone.rays.push_back(lift_zmodule(ZModule_t(std::vector<bignum>{ 3, 0 }));
  polytope.cones.push_back(cone);
}
{
  Barvinok::BarvinokCone cone;
  cone.vertex = lift_zmodule(ZModule_t(std::vector<bignum>{ 7, 7 }));
  cone.rays.push_back(lift_zmodule(ZModule_t(std::vector<bignum>{ -3, 0}));
  cone.rays.push_back(lift_zmodule(ZModule_t(std::vector<bignum>{ -1, -3}));
  polytope.cones.push_back(cone);
}
{
  Barvinok::BarvinokCone cone;
  cone.vertex = lift_zmodule(ZModule_t(std::vector<bignum>{ 6, 4 }));
  cone.rays.push_back(lift_zmodule(ZModule_t(std::vector<bignum>{ 1, 3}));
  cone.rays.push_back(lift_zmodule(ZModule_t(std::vector<bignum>{ -5, -2}));
  polytope.cones.push_back(cone);
}

Barvinok barvinok;
barvinok.init(conf);
Z res = barvinok.solve(n, polytope);
OPA_DISP0(barvinok.residue_helper.data);
OPA_DISP0(res);
}

TEST(Barvinok, Test2) {
Barvinok::Configuration conf;
conf.enumerate_ub = int(1e5);
int n = 2;

Barvinok::BarvinokPolytope polytope;

{
  Barvinok::BarvinokCone cone;
  cone.vertex = lift_zmodule(ZModule_t(std::vector<bignum>{ 0, 0 }));
  cone.rays.push_back(lift_zmodule(ZModule_t(std::vector<bignum>{ 1, 0}));
  cone.rays.push_back(lift_zmodule(ZModule_t(std::vector<bignum>{ 0, 1}));
  polytope.cones.push_back(cone);
}
{
  Barvinok::BarvinokCone cone;
  cone.vertex = lift_zmodule(ZModule_t(std::vector<bignum>{ 2, 0 }));
  cone.rays.push_back(lift_zmodule(ZModule_t(std::vector<bignum>{ -1, 0}));
  cone.rays.push_back(lift_zmodule(ZModule_t(std::vector<bignum>{ 0, 1}));
  polytope.cones.push_back(cone);
}
{
  Barvinok::BarvinokCone cone;
  cone.vertex = lift_zmodule(ZModule_t(std::vector<bignum>{ 2, 2 }));
  cone.rays.push_back(lift_zmodule(ZModule_t(std::vector<bignum>{ -1, 0}));
  cone.rays.push_back(lift_zmodule(ZModule_t(std::vector<bignum>{ 0, -1}));
  polytope.cones.push_back(cone);
}
{
  Barvinok::BarvinokCone cone;
  cone.vertex = lift_zmodule(ZModule_t(std::vector<bignum>{ 0, 2 }));
  cone.rays.push_back(lift_zmodule(ZModule_t(std::vector<bignum>{ 1, 0}));
  cone.rays.push_back(lift_zmodule(ZModule_t(std::vector<bignum>{ 0, -1}));
  polytope.cones.push_back(cone);
}

Barvinok barvinok;
barvinok.init(conf);
Z res = barvinok.solve(n, polytope);
OPA_DISP0(barvinok.residue_helper.data);
OPA_DISP0(res);
}

TEST(Barvinok, Test3) {
Barvinok::Configuration conf;
conf.enumerate_ub = int(1e5);
int n = 2;

Barvinok::BarvinokPolytope polytope;
std::vector<ZModule_t> vertices;
vertices.emplace_back({0, 0});
vertices.emplace_back({2, 0});
vertices.emplace_back({4, 2});
vertices.emplace_back({0, 2});
std::vector<QModule_t> qvertices = lift_zmodules(vertices);

Polyhedra_VertexRepr polyhedra;
polyhedra.init(qvertices);
OPA_DISP0(polyhedra.faces);

{
  Barvinok::BarvinokCone cone;
  cone.vertex = lift_zmodule(ZModule_t(std::vector<bignum>{ 0, 2 }));
  cone.rays.push_back(lift_zmodule(ZModule_t(std::vector<bignum>{ 1, 0}));
  cone.rays.push_back(lift_zmodule(ZModule_t(std::vector<bignum>{ 0, -1}));
  cone.vertex = Module_t(std::vector<bignum>{ 0, 0 });
  cone.rays.emplace_back(std::vector<bignum>{ 2, 0 });
  cone.rays.emplace_back(std::vector<bignum>{ 0, 2 });
  polytope.cones.push_back(cone);
}
{
  Barvinok::BarvinokCone cone;
  cone.vertex = Module_t(std::vector<bignum>{ 2, 0 });
  cone.rays.emplace_back(std::vector<bignum>{ -2, 0 });
  cone.rays.emplace_back(std::vector<bignum>{ 2, 2 });
  polytope.cones.push_back(cone);
}
{
  Barvinok::BarvinokCone cone;
  cone.vertex = Module_t(std::vector<bignum>{ 4, 2 });
  cone.rays.emplace_back(std::vector<bignum>{ -2, -2 });
  cone.rays.emplace_back(std::vector<bignum>{ -4, -0 });
  polytope.cones.push_back(cone);
}
{
  Barvinok::BarvinokCone cone;
  cone.vertex = Module_t(std::vector<bignum>{ 0, 2 });
  cone.rays.emplace_back(std::vector<bignum>{ 4, 0 });
  cone.rays.emplace_back(std::vector<bignum>{ 0, -2 });
  polytope.cones.push_back(cone);
}

Barvinok barvinok;
barvinok.init(conf);
Z res = barvinok.solve(n, polytope);
OPA_DISP0(barvinok.residue_helper.data);
OPA_DISP0(res);
}

TEST(Barvinok, TestSplit) {
Barvinok::Configuration conf;
conf.enumerate_ub = 1;
conf.max_depth = 1;
int n = 2;

std::vector<Module_t> vertices;
vertices.emplace_back(std::vector<bignum>{ 0, 0 });
vertices.emplace_back(std::vector<bignum>{ 1, 0 });
vertices.emplace_back(std::vector<bignum>{ 1, 5 });
Simplex simplex(vertices);

Barvinok barvinok;
barvinok.init(conf);
Z res = barvinok.solve_simplex(simplex);
OPA_DISP0(barvinok.residue_helper.data);
OPA_DISP0(res);
}

Z get_nosplit(const Simplex &simplex) {
Barvinok::Configuration conf;
conf.enumerate_ub = 1;
conf.max_depth = 0;

Barvinok barvinok;
barvinok.init(conf);
Z res = barvinok.solve_simplex(simplex);
return res;
}


TEST(Barvinok, Test3D_1) {
Barvinok::Configuration conf;
conf.enumerate_ub = 1;
conf.max_depth = 1;
int n = 3;

std::vector<Module_t> vertices;
vertices.emplace_back(std::vector<bignum>{ 0, 1, 1 });
vertices.emplace_back(std::vector<bignum>{ 1, 1, 1 });
vertices.emplace_back(std::vector<bignum>{ 1, 3, 5 });
vertices.emplace_back(std::vector<bignum>{ -3, 7, 0 });
Simplex simplex(vertices);

Barvinok barvinok;
barvinok.init(conf);
Z res = barvinok.solve_simplex(simplex);
OPA_DISP0(barvinok.residue_helper.data);
Z ans1 = get_nosplit(simplex);
Z ans2 = get_bruteforce(simplex);
OPA_DISP0(res, ans1, ans2);
}

TEST(Barvinok, GenLattice) {
int n = 3;
int bound = 100;
auto mat = generate_lattice(n, n, bound, 100);

std::vector<Module_t> vertices;
vertices.push_back(Module_t::zero(n));
REP (i, n)
  vertices.emplace_back(mat.get_col(i).tovec());
Simplex simplex(vertices);

Barvinok::Configuration conf;
conf.enumerate_ub = 1;
conf.max_depth = 1;
Barvinok barvinok;
barvinok.init(conf);
Z res = barvinok.solve_simplex(simplex);
OPA_DISP0(barvinok.residue_helper.data);
Z ans1 = get_nosplit(simplex);
Z ans2 = get_bruteforce(simplex);
OPA_DISP0(res, ans1, ans2);
}

template <class T>
Matrix<T> expand_unimodular_mdim(const Ring<T> *ring,
                               const std::vector<std::vector<T> > &tb) {
// not working
// TODO: find algorithm
//OPA_CHECK0(tb.size() > 0);
//int n = tb[0].size();
//Matrix<T> res = Matrix<T>::identity(ring, n);

//REP (i, tb.size()) {
//  std::vector<T> remapped = res.evalT(tb[i]);
//  remapped.resize(n - i);
//  int curn = n - i;
//  auto nxt = expand_unimodular(ring, remapped, curn - 1);
//  res.get_mutable(0, 0, curn, curn).self_mul(nxt);
//  OPA_DISP0(res, nxt);
//}
//return res;
}

TEST(Barvinok, LowDim1) {
int dim = 3;
auto mat = generate_lattice(3, 3, 10, 100);
auto proj_mat = expand_unimodular(&Ring_Z, mat.get_row(0).tovec(), 0);

Barvinok::BarvinokPolytope polytope;

{
  Barvinok::BarvinokCone cone;
  cone.vertex = Module_t(std::vector<bignum>{ 0, 0 });
  cone.rays.emplace_back(std::vector<bignum>{ 2, 0 });
  cone.rays.emplace_back(std::vector<bignum>{ 0, 2 });
  polytope.cones.push_back(cone);
}
{
  Barvinok::BarvinokCone cone;
  cone.vertex = Module_t(std::vector<bignum>{ 2, 0 });
  cone.rays.emplace_back(std::vector<bignum>{ -2, 0 });
  cone.rays.emplace_back(std::vector<bignum>{ 2, 2 });
  polytope.cones.push_back(cone);
}
{
  Barvinok::BarvinokCone cone;
  cone.vertex = Module_t(std::vector<bignum>{ 4, 2 });
  cone.rays.emplace_back(std::vector<bignum>{ -2, -2 });
  cone.rays.emplace_back(std::vector<bignum>{ -4, -0 });
  polytope.cones.push_back(cone);
}
{
  Barvinok::BarvinokCone cone;
  cone.vertex = Module_t(std::vector<bignum>{ 0, 2 });
  cone.rays.emplace_back(std::vector<bignum>{ 4, 0 });
  cone.rays.emplace_back(std::vector<bignum>{ 0, -2 });
  polytope.cones.push_back(cone);
}

Barvinok barvinok;
barvinok.init(conf);
Z res = barvinok.solve(n, polytope);
OPA_DISP0(barvinok.residue_helper.data);
OPA_DISP0(res);
}

TEST(Lattice, Gen1) {
auto mat = generate_lattice(5, 5, 100, 100);
OPA_DISP0(mat);
OPA_DISP0(hnf_det(mat));
}


TEST(Lattice, GenMul) {
int n = 5;
auto mat = generate_lattice(n, n, 100, 100);
auto vspace = mat.get_submatrix(0, 0, 2).clone();
gram_schmidt(vspace);
OPA_DISP0(vspace);

std::vector<std::vector<Z> > tb = vspace.as_rows();

auto m2 = expand_unimodular_mdim(&Ring_Z, tb);
OPA_DISP0(m2);
OPA_DISP0(hnf_det(m2));
auto im2 = m2.inverse();
OPA_DISP0(im2.evalT(tb[0]));
OPA_DISP0(im2.evalT(tb[1]));
OPA_DISP0(m2);
OPA_DISP0(im2);
OPA_DISP0(tb);
}

TEST(Taylor, Exp1) {
TaylorExp exp1(1);
OPA_DISP0(exp1.get_as_poly(10));
}

TEST(Taylor, X1Exp1) {
Taylor_x_over_1_eminusx t(10);
OPA_DISP0(t.get_as_poly(10));
OPA_DISP0(t.clist);
}
*/

GTEST_API_ int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  opa::init::opa_init(argc, argv);
  init_test();
  return RUN_ALL_TESTS();
}

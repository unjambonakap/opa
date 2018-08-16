#include "earth_render.h"

#include <opa/engine/Camera.h>
#include <opa/engine/PosObject.h>
#include <opa/engine/common.h>
#include <opa/engine/conf.h>
#include <opa/engine/primitives/Base.h>
#include <opa/math/game/conf.h>
#include <opa/math/game/geo_2d.h>
#include <opa/math/game/intersection.h>
#include <opa/math/game/mesh_util.h>
#include <opa/or/best_first_search.h>
#include <opa/utils/csv.h>
#include <opa_common.h>

double kEarthRadiusM = 6.371;
using namespace opa::math::game;

namespace opa {
namespace engine {
const FaceIndexList square_tr_id_from_coord = { { 0, 1, 3 }, { 3, 2, 0 } };
const FaceIndex square_coord_to_poly = { 0, 1, 3, 2 };
const int kMaxDepth = 19;
const int kMinDepth = 4;

namespace earth {
Pos2 uv_to_latlng(const Pos2 &uv) {
  double lng = (uv.x - 0.5) * 2 * PI;
  double lat = std::atan(std::exp(2 * PI * (0.5 - uv.y))) * 2 - PI / 2;
  return Pos2{ lat, lng };
}

Pos latlng_to_pos(const Pos2 &latlng) {
  double a = cos(latlng.x);
  return Pos(cos(latlng.y) * a, sin(latlng.y) * a, sin(latlng.x));
}

Pos get_normal(const Pos2 &uv) {
  Pos2 latlng = uv_to_latlng(uv);
  return latlng_to_pos(latlng);
}

Pos get_pos(const Pos2 &uv) { return get_normal(uv); }
Pos2 get_corner_latlng(const AutoKDNode *node, int id) {
  return earth::uv_to_latlng(
    Pos2(node->low + (node->high - node->low) * Pos2(id & 1, id >> 1)));
}

Pos uv_to_pos(const AutoKDNode *node, const Pos2 &uv) {
  return earth::get_pos(Pos2(node->low + (node->high - node->low) * uv));
}

Pos get_corner_pos(const AutoKDNode *node, int id) {
  return uv_to_pos(node, Pos2(id & 1, id >> 1));
}

} // namespace earth
bool file_exists(const std::string &fname) {
  return std::ifstream(fname).good();
}

Pos2 screen_space_to_pix_units(const EnvDesc &env, const Pos2 &pos) {
  return pos / 2 * Pos2(env.camera.viewport());
}

bool tr_in_screen(const Tr3D &tr) {

  Pos normal = -get_plane_normal(
    tr[0], tr[1], tr[2]); // flipping normal because of screen coordinates for
                          // texture (00 top left)
  double v = glm::dot(normal, tr[0]);
  double m = glm::dot(glm::abs(normal), Pos(1));
  return (v <= m);
}

TexturePtr TileResourceManager::create_tex(const AutoKDNode &node) const {
  Image img;
  if (m_mosaic) {
    return m_cp.get_rand();
  } else {
    TexturePtr tex(new Texture);
    img.init(256, 256, Image::ImageFormat::RGB);
    std::string fname = get_tile_path(node);
    if (!file_exists(fname)) {
      memset(img.buf(), 0x80, img.bytesize());
    } else {
      OPA_CHECK(img.load_png(fname), fname);
    }
    tex->init(img);
    return tex;
  }
}
AutoKDTree::AutoKDTree() {

  m_pool.emplace_back();
  root = &m_pool.back();
  root->depth = 0;
  root->low = Pos2(0, 0);
  root->high = Pos2(1, 1);
  root->is_leaf = false;
  root->coord = IPos2(0, 0);
  screen2d_box = opa::math::game::BoxAASpec_Gen2D().from_center_and_dim(
    Pos2(0, 0), Pos2(2, 2));
}

double AutoKDTree::required_precision(bool is_actif) const {
  return params.precision * (is_actif ? 1.1 : 0.9);
}

void AutoKDTree::prepare_node(AutoKDNode *node) {
  if (node->is_leaf || node->next[0]) return; // is already prepared

  if ((node->is_leaf = (node->depth == kMaxDepth))) return;
  Pos2 mid = (node->low + node->high) / 2;
  REP (i, 4) {
    m_pool.emplace_back();
    AutoKDNode *child = &m_pool.back();
    node->next[i] = child;
    child->depth = node->depth + 1;

    bool highx = (i & 1);
    bool highy = (i & 2);
    child->low =
      Pos2(!highx ? node->low.x : mid.x, !highy ? node->low.y : mid.y);
    child->high =
      Pos2(!highx ? mid.x : node->high.x, !highy ? mid.y : node->high.y);
    child->coord = 2 * node->coord + IPos2(highx, highy);
  }
}

bool is_backface_point(const Pos &cam_pos, const Dir &cam_dir,
                       const Pos &world_pos, const Dir &face_dir) {
  return glm::dot(cam_dir, face_dir) > 0 ||
         glm::dot(world_pos - cam_pos, face_dir) > 0;
}

bool AutoKDTree::should_split(const EnvDesc &env, const AutoKDNode *node,
                              const std::array<Pos, 4> &corners_screen_space,
                              bool &nbackface) const {
  nbackface = true;

  utils::ScopedDebugCtx sctx1(ctx, "should_split");
  std::array<Pos, 4> corners_screen_space_poly;
  REP (i, 4)
    corners_screen_space_poly[i] =
      corners_screen_space[square_coord_to_poly[i]];

  for (auto &tr_id : square_tr_id_from_coord) {
    utils::ScopedDebugCtx sctx2(ctx, OPA_STREAM_STR(tr_id));
    std::array<Pos, 3> tr;
    REP (i, 3)
      tr[i] = corners_screen_space[tr_id[i]];

    PlaneSpec cam_plane = PlaneSpec::FromPoints(tr);
    std::array<Pos2, 3> corners_proj;
    REP (i, 3)
      corners_proj[i] = cam_plane.proj(tr[i]);
    corners_proj = reorient_counterclockise(corners_proj);

    std::array<Pos2, 4> boundaries;
    REP (i, 4)
      boundaries[i] = cam_plane.proj(line_plane_intersection(
        LineSpec{ Pos(screen2d_box.get(i), -1), Pos(screen2d_box.get(i), 1) },
        cam_plane));

    std::array<Pos, 3> local_corners; // should not need to recompute those
    REP (i, 3)
      local_corners[i] = ApplyMat(env.screen_to_obj, tr[i]);
    double orig_area = tr_area(local_corners[1] - local_corners[0],
                               local_corners[2] - local_corners[0]);
    Pos normal =
      -get_plane_normal(local_corners[0], local_corners[1], local_corners[2]);
    OPA_YAML_DUMP(ctx, tr, orig_area, corners_screen_space_poly, boundaries);
    OPA_YAML_DUMP(ctx, local_corners, normal);

    for (auto &tr_boundaries_id : square_tr_id_from_coord) {
      utils::ScopedDebugCtx sctx2(ctx, OPA_STREAM_STR(tr_boundaries_id));
      std::array<Pos2, 3> tr_boundaries;
      REP (i, 3)
        tr_boundaries[i] = boundaries[tr_boundaries_id[i]];
      tr_boundaries = reorient_counterclockise(tr_boundaries);

      Point2Vec tr_inter = tr_tr_intersection(corners_proj, tr_boundaries);
      OPA_YAML_DUMP(ctx, corners_proj, tr_boundaries, tr_inter);
      if (tr_inter.size() < 3) continue;
      // OPA_DISP("intersection >> ", tr_inter, corners_proj, tr_boundaries);

      PointVec lifted_inter = ApplyOp<Pos2, Pos>(
        tr_inter, [&](const Pos2 &x) { return cam_plane.lift(x); });
      OPA_YAML_DUMP(ctx, lifted_inter);

      std::vector<Pos> inter_local_pos;
      for (const Pos &screen_pt : lifted_inter) {
        Pos local_pos = ApplyMat(env.screen_to_obj, screen_pt);
        Pos world_pos = env.obj_pos.get_world_pos(local_pos);
        inter_local_pos.push_back(local_pos);
        const Pos camera_front = env.cam_pos.get_front();
        nbackface &= is_backface_point(env.cam_pos.pos(), camera_front,
                                       world_pos, world_pos);
      }

      double inter_area = polygon_area(inter_local_pos);
      OPA_YAML_DUMP(ctx, inter_local_pos, inter_area, orig_area);

      Point2Vec inter_screen = ApplyOp<Pos, Pos2>(
        lifted_inter, [&](const Pos &x) -> Pos2 { return x.xy; });
      double screen_area = polygon_area(inter_screen);
      double covered_screen =
        screen_area / 4; // Screen space [-1,1]^2 -> screen area = 4

      double tex_visible_fraction = inter_area / orig_area;
      double tex_pixels = tex_visible_fraction / pix_width() / pix_width();
      double screen_tex_pixels =
        covered_screen * env.camera.viewport().x * env.camera.viewport().y;
      // Simple check on the area
      OPA_YAML_DUMP(ctx, inter_screen);
      OPA_YAML_DUMP(ctx, tex_pixels, screen_area, covered_screen,
                    screen_tex_pixels, tex_visible_fraction, tr);

      double target_prec = required_precision(node->actif_id == actif_id);
      double wanted_pix = target_prec * screen_tex_pixels;
      OPA_YAML_DUMP(ctx, target_prec, wanted_pix, tex_pixels < wanted_pix,
                    lifted_inter.size());

      if (tex_pixels < wanted_pix) return true;

      for (const Pos &screen_pt : lifted_inter) {
        Pos tex_pt = ApplyMat(env.screen_to_obj, screen_pt);

        // TODO: change if spheroid
        Pos normal = tex_pt; // sphere centered on 0
        for (auto &delta_vec :
             { env.camera.delta_pix_w(), env.camera.delta_pix_h() }) {
          Pos X = Mat3(env.screen_to_obj) * (screen_pt + delta_vec) +
                  glm::column(env.screen_to_obj, 3).xyz();
          Pos Y = Mat3(env.screen_to_obj) * vec_z;
          Pos K = glm::row(env.screen_to_obj, 3).xyz();
          double b =
            glm::dot(K, screen_pt + delta_vec) + env.screen_to_obj[3][3];
          double a = glm::dot(K, vec_z);
          double kX = glm::dot(X, normal);
          double kY = glm::dot(Y, normal);
          double kk = glm::dot(tex_pt, normal);
          double k = (kX - kk * b) / (kk * a - kY);

          Pos screen_np = screen_pt + delta_vec + k * vec_z;
          Pos np = ApplyMat(env.screen_to_obj, screen_np);
          Pos u = np - tex_pt;
          double norm =
            glm::length(u); // 1 pixel on the screen correspond to a vector
                            // of size norm on the sprite plane

          // if (norm < pix_width() / 2) return true;
        }
      }
    }
  }
  return false;
}

double AutoKDTree::pix_width() const { return 1. / 256; } //

bool AutoKDTree::do_split(const EnvDesc &env,
                          std::vector<const AutoKDNode *> &nodes,
                          AutoKDNode *cur) {
  this->prepare_node(cur);
  if (!cur->is_leaf) {
    for (auto &subnode : cur->next) this->_collect(env, nodes, subnode);
    return true;
  }
  return false;
}

void AutoKDTree::collect(const EnvDesc &env,
                         std::vector<const AutoKDNode *> &nodes_out) {
  actif_id++;
  ss = std::ostringstream();
  ctx.reset();
  return _collect(env, nodes_out, root);
}

void AutoKDTree::prepare_add_node(const EnvDesc &env, AutoKDNode *node) {
  node->actif_id = actif_id;

  struct State {
    Tr2D cur_tr;
    int depth = 0;
  };
  typedef OR::Search_dfs<State> Search_t;
  node->tex_trs.clear();

  utils::ScopedDebugCtx sctx1(ctx, "tr-split");
  Search_t search;
  Search_t::SearchParams params;
  params.func = [&](OR::Search<State> *self, const State &state) {
    utils::MaxFinderPair<double, int> mf;

    OPA_CHECK(state.depth < 15, this->ctx.flush_and_str());
    auto to_screen = [&](const Pos2 &x) -> Pos {
      return env.camera.to_screen_pos(
        ApplyMat(env.to_cam_space, earth::uv_to_pos(node, x)));
    };

    Tr3D screen_space_points;
    Tr2D screen_space_points_2d;

    opa::math::game::BoxAASpec_Gen2D xybox;
    REP (i, 3) {
      OPA_YAML_PUSH(ctx, ApplyMat(env.obj_pos.mat_to_world(),
                                  earth::uv_to_pos(node, state.cur_tr[i])));
      screen_space_points[i] = to_screen(state.cur_tr[i]);
      screen_space_points_2d[i] = screen_space_points[i].xy;
      xybox.update(screen_space_points_2d[i]);
    }
    if (are_aligned(screen_space_points[0], screen_space_points[1],
                    screen_space_points[2]))
      return;
    if (!tr_in_screen(screen_space_points)) return;
    if (!xybox.intersects(screen2d_box)) return;

    utils::ScopedDebugCtx sctx2(ctx);
    OPA_YAML_DUMP(ctx, state.cur_tr, screen_space_points);

    REP (i, 3) {
      int ni = (i + 1) % 3;
      Pos2 np = (state.cur_tr[i] + state.cur_tr[ni]) / 2;
      Pos2 real_cam_np = to_screen(np);
      Pos2 line_cam_np =
        (screen_space_points_2d[i] + screen_space_points_2d[ni]) / 2;
      double l =
        glm::length(screen_space_to_pix_units(env, real_cam_np - line_cam_np));
      OPA_YAML_DUMP(ctx, l, real_cam_np, line_cam_np, np, i);
      mf.update(l, i);
    }

    // TODO: undo
    if (mf.get_cost() > this->params.split_tr_norm) {

      int i = mf.get();
      Pos2 np = (state.cur_tr[i] + state.cur_tr[(i + 1) % 3]) / 2;
      self->add(State{ { state.cur_tr[i], np, state.cur_tr[(i + 2) % 3] },
                       state.depth + 1 });
      self->add(
        State{ { np, state.cur_tr[(i + 1) % 3], state.cur_tr[(i + 2) % 3] },
               state.depth + 1 });
    } else {
      node->tex_trs.push_back(state.cur_tr);
    }
  };

  search.init(params);
  for (auto &face_vids : square_tr_id_from_coord) {
    State s0;
    REP (i, 3)
      s0.cur_tr[i] = Pos2(face_vids[i] & 1, (face_vids[i] >> 1));
    search.add(s0);
  }
  search.start();
}

void AutoKDTree::_collect(const EnvDesc &env,
                          std::vector<const AutoKDNode *> &nodes,
                          AutoKDNode *cur) {
  utils::ScopedDebugCtx sctx1(ctx);
  opa::math::game::BoxAASpec_Gen3D box;
  bool backface = true;
  Pos camera_front = env.cam_pos.get_front();
  env.camera.setup();

  std::array<Pos, 4> corners_screen_space;
  bool force_split = cur->depth < kMinDepth;
  bool want_split = force_split;
  if (!want_split) {
    std::array<Pos, 4> screen_space_corners;
    std::array<Pos, 4> world_corners;

    bool outside = false;
    REP (i, 4) {
      auto corner_pos = earth::get_corner_pos(cur, i);
      auto world_pos = env.obj_pos.get_world_pos(corner_pos);
      world_corners[i] = world_pos;
      OPA_YAML_PUSH(ctx, world_pos);
      OPA_YAML_PUSH(ctx, corner_pos);

      auto corner_cam_pos = ApplyMat(env.to_cam_space, corner_pos);
      corners_screen_space[i] = env.camera.to_screen_pos(corner_cam_pos);
      box.update(corners_screen_space[i]);
      double dot = glm::dot(world_pos, camera_front);
      backface &= is_backface_point(env.cam_pos.pos(), camera_front, world_pos,
                                    world_pos);
    }

    OPA_YAML_DUMP(ctx, camera_front);
    if (backface) return;

    bool in_screen = false;
    for (auto &tr_id : square_tr_id_from_coord) {
      Tr3D tr;
      REP (i, 3)
        tr[i] = corners_screen_space[tr_id[i]];
      in_screen |= tr_in_screen(tr);

      std::vector<HyperPlaneSpec> hps =
        env.camera.hps_cam_space(env.cam_pos); // TODO: compute once
      for (auto &hp : hps) {
        bool all_in = true;
        for (auto &pt : world_corners) {
          if (!hp.is_in(pt)) {
            all_in = false;
            break;
          }
        }
        OPA_YAML_PUSH(ctx, all_in, hp.plane.dir, hp.plane.v);
        if (all_in) return;
      }
    }
    if (!in_screen) return;
  }
  bool contains = (box.high.z >= -1) && (box.low.z <= 1);
  if (!want_split && !contains) return;

  if (!want_split) {
    opa::math::game::BoxAASpec_Gen2D xybox = box.subbox<Pos2, 2>({ { 0, 1 } });

    if (!xybox.intersects(screen2d_box)) return;

    bool nbackface = false;
    want_split = should_split(env, cur, corners_screen_space, nbackface);
    if (nbackface) return;
  }

  bool has_split = want_split && do_split(env, nodes, cur);
  if (!force_split) {
    OPA_YAML_DUMP(ctx, box, corners_screen_space, want_split, has_split,
                  cur->depth);
  }
  if (has_split)
    ;
  else {
    prepare_add_node(env, cur);
    OPA_YAML_DUMP(ctx, cur->tex_trs);
    cur->debug = ctx.cur_str();
    nodes.push_back(cur);
  }
}

void EnvDesc::setup() {
  to_cam_space = cam_pos.mat_to_local() * obj_pos.mat_to_world();
  screen_to_obj =
    obj_pos.mat_to_local() * cam_pos.mat_to_world() * camera.iperspective();
  OPA_DISP("ON KAPPA ", cam_pos.rot(), cam_pos.pos(), screen_to_obj);
}

EnvDesc init_earth_env() {
  EnvDesc desc;
  desc.camera.update_viewport(1280, 1024);
  desc.camera.near_field() = 1e-3;
  desc.camera.far_field() = 1e2;
  desc.camera.setup();
  desc.obj_pos.scale() = kEarthRadiusM;
  desc.obj_pos.pos() = Pos(0, 0, 0);
  desc.obj_pos.rot() = opa::math::game::quat_look_at_safe(vec_x, vec_z);

  desc.cam_pos.pos() =
    ApplyMat(desc.obj_pos.mat_to_world(),
             earth::latlng_to_pos(Pos2(45.7, 4.85) / 360 * 2 * PI)) *
    (kEarthRadiusM + kEarthRadiusM * 0.001) / kEarthRadiusM;

  desc.cam_pos.rot() =
    opa::math::game::quat_look_at_safe(-desc.cam_pos.pos(), vec_z);
  desc.setup();
  return desc;
}

SceneObjPtr create_earth_obj(const TileResourceManager &tm,
                             const std::vector<const AutoKDNode *> &nodes,
                             Game *game, bool obj_per_node) {
  auto obj = std::make_shared<TrianglesSceneObj>();
  obj->init(game);
  game->register_obj(obj, game->scene());

  for (auto &node : nodes) {
    TrianglesSceneObj *cobj = obj.get();
    if (obj_per_node) {
      auto node_obj = std::make_shared<TrianglesSceneObj>();
      node_obj->init(game);
      game->register_obj(node_obj, obj.get());
      cobj = node_obj.get();
      node_obj->set_debug(node->debug);
    }

    for (auto &tex_tr : node->tex_trs) {
      auto tex = tm.get_tex(*node);
      TexturedTriangle &tr = cobj->tr().triangles().get_new2();
      tr.init(tex);
      REP (i, 3) {
        tr.data()[i].pos = earth::uv_to_pos(node, tex_tr[i]);
        tr.data()[i].texpos = tex_tr[i];
      }
      tr.refresh();
    }
  }
  return std::static_pointer_cast<SceneObj, TrianglesSceneObj>(obj);
}

} // namespace engine
} // namespace opa

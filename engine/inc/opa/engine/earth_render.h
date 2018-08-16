#pragma once

#include <opa/engine/Camera.h>
#include <opa/engine/PosObject.h>
#include <opa/engine/common.h>
#include <opa/engine/conf.h>
#include <opa/engine/primitives/Base.h>
#include <opa/math/game/conf.h>
#include <opa/math/game/geo_2d.h>
#include <opa/math/game/intersection.h>
#include <opa/math/game/mesh_util.h>
#include <opa/utils/csv.h>
#include <opa_common.h>
#include <opa/utils/debug.h>

#define _OPA_YAML_DUMP2(dest, a) dest.set(#a, (a))
#define _OPA_YAML_DUMP3(dest, a, b)                                            \
  _OPA_YAML_DUMP2(dest, a), _OPA_YAML_DUMP2(dest, b)
#define _OPA_YAML_DUMP4(dest, a, b, c)                                         \
  _OPA_YAML_DUMP3(dest, a, b), _OPA_YAML_DUMP2(dest, c)
#define _OPA_YAML_DUMP5(dest, a, b, c, d)                                      \
  _OPA_YAML_DUMP4(dest, a, b, c), _OPA_YAML_DUMP2(dest, d)
#define _OPA_YAML_DUMP6(dest, a, b, c, d, e)                                   \
  _OPA_YAML_DUMP5(dest, a, b, c, d), _OPA_YAML_DUMP2(dest, e)
#define _OPA_YAML_DUMP7(dest, a, b, c, d, e, f)                                \
  _OPA_YAML_DUMP6(dest, a, b, c, d, e), _OPA_YAML_DUMP2(dest, f)
#define _OPA_YAML_DUMP8(dest, a, b, c, d, e, f, g)                             \
  _OPA_YAML_DUMP7(dest, a, b, c, d, e, f), _OPA_YAML_DUMP2(dest, g)
#define OPA_YAML_DUMP(...)                                                     \
  OPA_DISPATCH(_OPA_YAML_DUMP, __VA_ARGS__)(__VA_ARGS__)

#define _OPA_YAML_PUSH2(dest, a) dest.push_back(#a, (a))
#define _OPA_YAML_PUSH3(dest, a, b)                                            \
  _OPA_YAML_PUSH2(dest, a), _OPA_YAML_PUSH2(dest, b)
#define _OPA_YAML_PUSH4(dest, a, b, c)                                         \
  _OPA_YAML_PUSH3(dest, a, b), _OPA_YAML_PUSH2(dest, c)
#define _OPA_YAML_PUSH5(dest, a, b, c, d)                                      \
  _OPA_YAML_PUSH4(dest, a, b, c), _OPA_YAML_PUSH2(dest, d)
#define _OPA_YAML_PUSH6(dest, a, b, c, d, e)                                   \
  _OPA_YAML_PUSH5(dest, a, b, c, d), _OPA_YAML_PUSH2(dest, e)
#define _OPA_YAML_PUSH7(dest, a, b, c, d, e, f)                                \
  _OPA_YAML_PUSH6(dest, a, b, c, d, e), _OPA_YAML_PUSH2(dest, f)
#define _OPA_YAML_PUSH8(dest, a, b, c, d, e, f, g)                             \
  _OPA_YAML_PUSH7(dest, a, b, c, d, e, f), _OPA_YAML_PUSH2(dest, g)
#define OPA_YAML_PUSH(...)                                                     \
  OPA_DISPATCH(_OPA_YAML_PUSH, __VA_ARGS__)(__VA_ARGS__)
namespace opa {
namespace engine {
class Game;

struct EnvDesc {
  opa::engine::Camera camera;
  opa::engine::PosObject cam_pos;
  opa::engine::PosObject obj_pos;
  opa::Mat4 to_cam_space;
  opa::Mat4 screen_to_obj;
  void setup();
};

struct AutoKDNode {
  std::array<AutoKDNode *, 4> next{};
  int depth;
  Pos2 low, high;
  IPos2 coord;
  bool is_leaf = false;
  int actif_id = -1;
  std::vector<Tr2D> tex_trs;
  std::string debug;
};

class TileResourceManager {
public:
  TileResourceManager(const std::string &resources_dir, bool mosaic = false)
      : m_dir(resources_dir), m_mosaic(mosaic) {
    m_cp.resize_rand(1000);
  }

  std::string get_tile_path(const AutoKDNode &node) const {
    return glib::strings::Substitute("$0/tile_$1_$2_$3.png", m_dir, node.depth,
                                     node.coord.x, node.coord.y);
  }

  TexturePtr create_tex(const AutoKDNode &node) const;

  TexturePtr get_tex(const AutoKDNode &node) const {
    if (!rmp.count(&node)) rmp[&node] = create_tex(node);
    return rmp[&node];
  }

  mutable std::map<const AutoKDNode *, TexturePtr> rmp;

private:
  bool m_mosaic;
  ColorPalette m_cp;
  std::string m_dir;
};

struct AutoKDParams {
  double precision = 4;
  double split_tr_norm = 1;
};


class AutoKDTree {
public:
  AutoKDTree();

  void prepare_node(AutoKDNode *node);

  bool should_split(const EnvDesc &env, const AutoKDNode *node,
                    const std::array<Pos, 4> &corners, bool &nbackface) const;

  double pix_width() const; //

  bool do_split(const EnvDesc &env, std::vector<const AutoKDNode *> &nodes,
                AutoKDNode *cur);

  void _collect(const EnvDesc &env, std::vector<const AutoKDNode *> &nodes,
                AutoKDNode *cur);
  void collect(const EnvDesc &env, std::vector<const AutoKDNode *> &nodes_out);
  double required_precision(bool is_actif) const;

  void prepare_add_node(const EnvDesc &env, AutoKDNode *node);

  AutoKDNode *root;
  opa::math::game::BoxAASpec_Gen2D screen2d_box;
  AutoKDParams params;

  std::string debug() const { return ss.str(); }

  std::deque<AutoKDNode> m_pool;
  int actif_id = -1;
  mutable utils::DebugCtx ctx;

private:
  mutable std::ostringstream ss;
};

SceneObjPtr create_earth_obj(const TileResourceManager &tm,
                             const std::vector<const AutoKDNode *> &nodes,
                             Game *game, bool obj_per_node);

EnvDesc init_earth_env();
} // namespace engine
} // namespace opa

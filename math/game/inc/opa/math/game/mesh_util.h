#pragma once

#include <opa/math/game/conf.h>
#include <opa/math/game/mesh.h>

OPA_NAMESPACE_DECL3(opa, math, game)

extern const FaceIndexList square_tr_id;
extern const FaceIndexList box_tr_id;
extern const FaceIndexList box_face_id;
extern const FaceIndexList square_face_id;

class KDTree {
public:
  typedef utils::IdType KDId;

  struct KDNode {
    KDId id;
    BoxAASpec box;
    KDId left, right;
    bool is_leaf() const { return left == utils::InvalidId; }

    // only valid for leaves
    int point_id() const { return right; }
  };

  void setup(const PointVec &tb);
  int get_closest(const Pos &pt) const;
  std::vector<int> get_neighborhood(const Pos &pt, double dist) const;

  void reset() {
    m_root = utils::InvalidId;
    m_tb.clear();
    m_kds.clear();
  }

private:
  struct GlobalClosestStatus {
    utils::MinFinder<std::pair<double, int> > best;
  };

  void find_neighborhood(KDId id, const Pos &pt, double dist,
                         std::vector<int> &neighbours) const;

  void find_closest(KDId id, const Pos &pt, GlobalClosestStatus &status) const;

  struct BuildState {
    std::vector<int> lst[3];

    static BuildState FromPoints(const PointVec &points) {
      BuildState res;
      std::vector<std::pair<double, int> > axis_and_id(points.size());
      REP (i, points.size())
        axis_and_id[i].second = i;

      REP (coord, 3) {
        REP (i, points.size())
          axis_and_id[i].first = points[axis_and_id[i].second][coord];
        std::sort(ALL(axis_and_id));
        REP (i, points.size())
          res.lst[coord].push_back(axis_and_id[i].second);
      }
      return res;
    }

    BuildState filter_list(const std::multiset<int> &to_keep) const {
      BuildState new_state;
      REP (coord, 3) {
        for (auto &v : lst[coord])
          if (to_keep.count(v)) new_state.lst[coord].push_back(v);
      }
      return new_state;
    }

    int n() const { return lst[0].size(); }

    void split2(const PointVec &points, int dim, double splitv,
                BuildState &left_state, BuildState &right_state) const {
      std::multiset<int> left;
      std::multiset<int> right;
      split(points, dim, splitv, left, right);
      left_state = this->filter_list(left);
      right_state = this->filter_list(right);
    }

    void split(const PointVec &points, int dim, double splitv,
               std::multiset<int> &left, std::multiset<int> &right) const {
      for (auto &v : lst[dim]) {
        if (points[v][dim] <= splitv)
          left.insert(v);
        else
          right.insert(v);
      }
    }

    Range1D get_range(int dim, const PointVec &points) const {
      OPA_CHECK0(lst[0].size() > 0);
      return Range1D{ points[lst[dim].front()][dim],
                      points[lst[dim].back()][dim] };
    }

    BoxAASpec compute_box(const PointVec &points) const {
      BoxAASpec res;
      REP (i, 3) {
        auto range = get_range(i, points);
        res.low[i] = range.low;
        res.high[i] = range.high;
      }
      return res;
    }
  };

  KDId build(const BuildState &state);
  std::pair<double, int> pick_split_dim(const BuildState &state) const;

  KDNode *get(KDId id) {
    if (id == utils::InvalidId) return nullptr;
    return &m_kds.get(id);
  }

  KDId m_root;
  PointVec m_tb;
  utils::ObjectPool<KDNode> m_kds;
};

class PointMatcher {
public:
  struct PointMatcheParams {
    bool activate = false;
    double eps = 1e-4;
  };

  void reset() {
    m_uj.reset(0);
    m_tree.reset();
  }

  PointMatcher() : m_params(PointMatcheParams()) {}
  PointMatcher(const PointMatcheParams &params) : m_params(params) {}
  bool load_batch(const PointVec &tb);
  int query(const Pos &pos) const;
  OPA_ACCESSOR_R(algo::UnionJoin, m_uj, uj);
  PointVec cur_batch;
  int final_count() const { return uj().count(); }

private:
  PointMatcheParams m_params;
  KDTree m_tree;
  algo::UnionJoin m_uj;
};

static PointMatcher::PointMatcheParams kPointMatcherDefaults =
  PointMatcher::PointMatcheParams();

struct GraphWithPosData {
  algo::FastGraph graph;
  PointVec points;
  std::vector<std::vector<int> > faces;
};

class FaceCollection;

class MeshBuilder {
public:
  MeshBuilder(Mesh *mesh = nullptr,
              PointMatcher::PointMatcheParams &params = kPointMatcherDefaults)
      : m_matcher(params) {
    if (mesh == nullptr) mesh = new Mesh;
    m_mesh = mesh;
  }

  // get or create
  VertexId get_vid(const Pos &pos);

  EdgeId get_edge(VertexId v0, VertexId v1);

  void add_triangle(const Pos &p0, const Pos &p1, const Pos &p2);
  void add_triangle(const Triangle &tr) {
    add_triangle(std::get<0>(tr), std::get<1>(tr), std::get<2>(tr));
  }
  void load_triangles(const std::vector<Triangle> &triangles);
  void add_face(const std::vector<Pos> &pos);

  Mesh *get_mesh();
  void repair_orientation();

  static void PolyhedraFromCloud(Mesh *mesh,
                                 const std::vector<VertexInfoData> &cloud);
  static void PolyhedraFromCloud(Mesh *mesh, const PointVec &cloud) {
    return PolyhedraFromCloud(mesh, vec_point_to_info_data(cloud));
  }
  static void PolygonFromCloud(Mesh *mesh,
                               const std::vector<VertexInfoData> &cloud);
  static void FromTetraedra(Mesh *mesh, const std::vector<VertexInfoData> &tb);

  static bool ConvexPolyhedra_IsInside(const Mesh &mesh, const Pos &pos);
  static void ConvexPolyhedra_Extend(Mesh *mesh,
                                     const std::vector<VertexInfoData> &cloud);
  static void CheckConvexPolyhedra(const Mesh &mesh);
  SPTR(GraphWithPosData) get_graph();
  void prepare();

  void fill(const FaceCollection &fc);

  bool is_manifold_mesh();

private:
  void compute();
  void compute_face(const PointVec &pv);

  std::vector<PointVec> m_faces;
  bool m_computed = false;
  bool m_prepared = false;
  PointMatcher m_matcher;
  Mesh *m_mesh;
  utils::TwoWayRemapper<Pos, VertexId> m_pos_vid;
  EdgeIndex m_edge_index;
  SPTR(GraphWithPosData) m_graph;
};

typedef std::function<bool(const PointVec &face)> mesh_face_filter_t;
class FaceCollection {
public:
  FaceCollection(bool triangulate = false) : m_triangulate(triangulate) {}
  FaceCollection(const std::vector<PointVec> &faces) : m_faces(faces) {}

  FaceCollection &push(const std::vector<PointVec> &faces) {
    for (auto &face : faces) push_face(face);
    return *this;
  }

  FaceCollection &push(const Triangle &tr) {
    return this->push_face(tr_to_vec(tr));
  }
  FaceCollection &push(const PointVec &face) { return this->push_face(face); }

  template <int N> FaceCollection &push(const std::array<Pos, N> &face) {
    PointVec pt(ALL(face));
    return this->push_face(pt);
  }

  FaceCollection &push(const FaceCollection &tc) { return push(tc.faces()); }

  FaceCollection &push_indexed(const PointVec &points,
                               const FaceIndexList &index_list);

  virtual FaceCollection &push_face(const PointVec &face);

  FaceCollection &load_stl(glib::StringPiece filename);
  FaceCollection &load_stl_from_data(glib::StringPiece data);
  FaceCollection &filter(mesh_face_filter_t filter_func);
  FaceCollection &add_plane(const PlaneSpec &plane, const Pos &center,
                            double dist, const Dir &base_dir = vec_0);
  FaceCollection &add_box(const Box2DSpec &box, const PlaneSpec &plane,
                          const Dir &front);
  FaceCollection &add_box(const BoxSpec &box);
  FaceCollection &add_box(const BoxAASpec &box) {
    return add_box(box.to_box());
  }

  template <class T> SPTR(T) convert() const {
    T *res = new T;
    this->to_mesh(res);
    return SPTR(T)(res);
  }
  FaceCollection &clear() {
    m_faces.clear();
    return *this;
  }

  void toggle_orientation();
  SPTR(MeshBuilder) to_mesh_builder(Mesh *mesh = nullptr) const;
  Mesh *to_mesh(Mesh *mesh = nullptr) const;
  OPA_ACCESSOR(std::vector<PointVec>, m_faces, faces);

  void write_stl(utils::BufferWriter &writer) const;
  void write_stl(glib::StringPiece filename) const {
    utils::BufferFileWriter writer{ filename };
    write_stl(writer);
  }
  FaceCollection triangulate() const;

  void to_stl_buf(std::string *res) const {
    utils::BufferWriter tmp;
    write_stl(tmp);
    *res = tmp.get();
  }

  void whiten();

protected:
  std::vector<PointVec> m_faces;
  bool m_triangulate;
};

void compute_new_center(const Mat4 &proj, const PointVec &point_cloud,
                        Pos *new_pos, Rot *new_rot, double ratio = 1);
void compute_new_center(const Mat4 &proj,
                        const std::vector<MeshData> &mesh_data, Pos *new_pos,
                        Rot *new_rot, double ratio = 1);

static inline PointVec to_point_cloud(const std::vector<MeshData> &mesh_data) {
  PointVec point_cloud;
  for (auto &e : mesh_data) {
    e.add_to_point_cloud(&point_cloud);
  }
  return point_cloud;
}

std::vector<SPTR(Mesh)> mesh_cut_plan(const Mesh &mesh,
                                      const HyperPlaneSpec &plane);
std::vector<SPTR(Mesh)> mesh_cut_box(const Mesh &mesh, const BoxSpec &box);
OPA_NAMESPACE_DECL3_END

#pragma once

#include <opa/algo/graph.h>
#include <opa/math/game/base_impl.h>
#include <opa/math/game/conf.h>
#include <opa/math/game/graph.h>
#include <opa/utils/DataStruct.h>

#include <lemon/concepts/graph.h>
#include <lemon/list_graph.h>
#include <lemon/matching.h>
#include <lemon/smart_graph.h>
#include <glib/gtl/optional.h>

OPA_NAMESPACE_DECL3(opa, math, game)
namespace proto {
class MeshFaceGraphList;
}

class MeshExport;
class MeshBuilder;
class FaceCollection;

typedef opa::utils::IdType GraphObjId;
typedef GraphObjId EdgeId;
typedef GraphObjId VertexId;
typedef GraphObjId FaceId;

static const FaceId FACE_NONE = -1;
static const EdgeId EDGE_NONE = -1;
static const VertexId VERTEX_NONE = -1;
static const GraphObjId ATTR_NONE = -1;

struct FaceGraph {
  SPTR(lemon::SmartGraph) graph;
  SPTR(lemon::SmartGraph::ArcMap<double>) edge_map;
};

struct FaceGraph2 {
  SPTR(algo::FastGraph) graph;
  std::vector<double> edge_costs;
};

struct VertexAttribute {
  GraphObjId id;
  Pos pos;
  TexUV uv;
  Dir normal;
};

class Edge : public opa::utils::IdObj {
public:
  VertexId start;
  FaceId face; // left side
  EdgeId opp;
  EdgeId next;
  GraphObjId vertex_attr;
  Edge();
  void reset();
};

class Vertex : public opa::utils::IdObj {
public:
  EdgeId edge; // one edge start at vertex
};

class Face : public opa::utils::IdObj {
public:
  EdgeId edge;
};

class EdgeIndex {
public:
  void add(VertexId a, VertexId b, EdgeId e);
  EdgeId get(VertexId a, VertexId b) const;
  bool has(VertexId a, VertexId b) const {
    return m_index.count(std::make_pair(a, b)) != 0;
  }

private:
  std::map<std::pair<VertexId, VertexId>, EdgeId> m_index;
};

typedef s32 VertexInfoId;
struct VertexInfoData {
  VertexInfoId id;
  Pos pos;
  VertexInfoData() {}
  VertexInfoData(VertexInfoId id, Pos pos) {
    this->id = id;
    this->pos = pos;
  }
};

inline std::vector<VertexInfoData> vec_point_to_info_data(const PointVec &tb) {
  std::vector<VertexInfoData> res(tb.size());
  REP (i, tb.size()) {
    res[i].id = i;
    res[i].pos = tb[i];
  }
  return res;
}

class Mesh : public opa::utils::Initable {

public:
  Mesh() {}
  virtual void init() override;
  virtual ~Mesh();

  double area() const;
  ;
  // double det() const;
  Vertex &get_vertex(VertexId id);
  Edge &get_edge(EdgeId id);
  Face &get_face(FaceId id);
  VertexAttribute &get_attr(GraphObjId id);

  const Vertex &get_vertex(VertexId id) const;
  const Edge &get_edge(EdgeId id) const;
  const Face &get_face(FaceId id) const;
  const VertexAttribute &get_attr(GraphObjId id) const;

  GraphObjId vertex_attr(VertexId vid, FaceId fid) const;
  GraphObjId vertex_attr(EdgeId id) const;
  const VertexAttribute &get_attr_for_edge(EdgeId id) const {
    return get_attr(vertex_attr(id));
  }

  EdgeId face_get_one_edge(FaceId id) const;
  EdgeId vertex_get_one_edge(VertexId id) const;
  FaceId vertex_get_one_face(VertexId id) const;
  GraphObjId vertex_get_one_attr(VertexId) const;

  EdgeId find_vertex_edge(VertexId vid, FaceId fid) const;

  FaceId get_edge_face(EdgeId id) const;
  // std::vector<EdgeId> get_edges(VertexId id) const;
  // std::vector<FaceId> get_faces(VertexId id) const;
  VertexId end(EdgeId id) const;
  VertexId start(EdgeId id) const;
  EdgeId opp(EdgeId id) const;
  EdgeId next(EdgeId id) const;
  FaceId face(EdgeId id) const;
  FaceId opp_face(EdgeId id) const;
  EdgeId rot_edge(EdgeId id) const;
  GraphObjId attr(EdgeId id) const;
  VertexAttribute &get_vertex_attr(VertexId id) {
    return get_attr(vertex_get_one_attr(id));
  }
  const VertexAttribute &get_vertex_attr(VertexId id) const {
    return get_attr(vertex_get_one_attr(id));
  }

  void set_vertex_attr(VertexId vid, GraphObjId aid);

  VertexId get_face_vertex(FaceId id) const;
  void setup_face(FaceId fid, const std::vector<VertexId> &vids,
                  const EdgeIndex &index);
  std::vector<FaceId> list_faces() const;
  std::vector<VertexId> list_vertices() const;
  std::vector<EdgeId> list_edges() const;

  Dir get_normal(FaceId face) const;
  Pos get_pos(VertexId id) const;
  std::vector<Pos> get_face_vertices(FaceId face) const;
  std::vector<VertexId> get_face_vertex_ids(FaceId face) const;
  PlaneSpec get_face_plane(FaceId face) const;

  std::vector<std::pair<FaceId, EdgeId> > get_adj_faces(FaceId face) const;
  Pos get_face_center(FaceId face) const;
  double face_face_dist(const Pos &f1_center, const Pos &f2_center, const EdgeId &edge) const;

  void do_remove_edge(EdgeId id) const;
  void do_remove_vertex(VertexId id) const;
  void do_remove_face(FaceId id) const;

  template <class OutputType> class BaseMeshIterator {
  public:
    typedef BaseMeshIterator<OutputType> SelfType;

    void init(const Mesh *mesh, const EdgeId &edge) {
      m_mesh = mesh;
      m_edge = edge;
      m_start_edge = edge;
      seen.insert(edge);
    }

    OutputType operator*() const { return get(m_edge); }
    const OutputType *operator->() const {
      m_output = get(m_edge);
      return &m_output;
    }

    SelfType &operator++() {
      m_edge = advance(m_edge);
      if (m_edge == m_start_edge) m_edge = EDGE_NONE;
      OPA_CHECK(seen.insert(m_edge).second, m_edge, seen, m_start_edge);
      return *this;
    }
    bool operator!=(const SelfType &it) const { return m_edge != it.m_edge; }

    virtual EdgeId advance(EdgeId cur) = 0;
    virtual OutputType get(EdgeId cur) const = 0;

  protected:
    mutable OutputType m_output;
    const Mesh *m_mesh;

  private:
    EdgeId m_start_edge;
    EdgeId m_edge;
    std::unordered_set<EdgeId> seen;
  };

  class Face_EdgeIterator : public BaseMeshIterator<EdgeId> {
  protected:
    virtual EdgeId advance(EdgeId cur) override { return m_mesh->next(cur); }
    virtual EdgeId get(EdgeId cur) const override { return cur; }
  };

  class Face_VertexIterator
    : public BaseMeshIterator<std::pair<VertexId, GraphObjId> > {
  protected:
    virtual EdgeId advance(EdgeId cur) override { return m_mesh->next(cur); }
    virtual std::pair<VertexId, GraphObjId> get(EdgeId cur) const override {
      return MP(m_mesh->start(cur), m_mesh->attr(cur));
    }
  };

  class Vertex_EdgeIterator : public BaseMeshIterator<EdgeId> {
  protected:
    virtual EdgeId advance(EdgeId cur) override {
      return m_mesh->next(m_mesh->opp(cur));
    }
    virtual VertexId get(EdgeId cur) const override { return cur; }
  };

  class Vertex_FaceIterator : public BaseMeshIterator<FaceId> {
  protected:
    virtual EdgeId advance(EdgeId cur) override {
      return m_mesh->next(m_mesh->opp(cur));
    }
    virtual VertexId get(EdgeId cur) const override {
      return m_mesh->get_edge_face(cur);
    }
  };

  template <class MeshIterator, class InputType> class BaseMeshWalker {
  public:
    BaseMeshWalker(const Mesh &mesh, InputType input) : m_mesh(mesh) {
      m_input = input;
    }
    MeshIterator begin() const {
      auto res = MeshIterator();
      res.init(&m_mesh, get_edge_from_input(m_input));
      return res;
    }
    MeshIterator end() const {
      auto res = MeshIterator();
      res.init(&m_mesh, EDGE_NONE);
      return res;
    }

  protected:
    virtual EdgeId get_edge_from_input(InputType input) const = 0;
    InputType m_input;
    const Mesh &m_mesh;
  };

  template <class MeshIterator>
  class Walker_FromFace : public BaseMeshWalker<MeshIterator, FaceId> {
  public:
    Walker_FromFace(const Mesh &mesh, FaceId input)
        : BaseMeshWalker<MeshIterator, FaceId>(mesh, input) {}
    virtual EdgeId get_edge_from_input(FaceId face) const override {
      return this->m_mesh.get_face(face).edge;
    }
  };

  template <class MeshIterator>
  class Walker_FromVertex : public BaseMeshWalker<MeshIterator, VertexId> {
  public:
    Walker_FromVertex(const Mesh &mesh, VertexId input)
        : BaseMeshWalker<MeshIterator, VertexId>(mesh, input) {}
    virtual EdgeId get_edge_from_input(VertexId vid) const override {
      return this->m_mesh.get_vertex(vid).edge;
    }
  };

  typedef Walker_FromVertex<Vertex_FaceIterator> Vertex_FaceWalker;
  typedef Walker_FromVertex<Vertex_EdgeIterator> Vertex_EdgeWalker;
  typedef Walker_FromFace<Face_VertexIterator> Face_VertexWalker;
  typedef Walker_FromFace<Face_EdgeIterator> Face_EdgeWalker;

  GraphObjId create_attr();
  typedef std::map<VertexId, VertexInfoId> OrigMapType;
  OPA_ACCESSOR(OrigMapType, m_orig_map, orig_map);
  OPA_ACCESSOR_R(std::set<VertexId>, m_vertices.used(), vertices);
  OPA_ACCESSOR_R(std::set<EdgeId>, m_edges.used(), edges);
  OPA_ACCESSOR_R(std::set<FaceId>, m_faces.used(), faces);

  void add_to_point_cloud(PointVec &dest) const {
    for (auto &vid : vertices()) {
      dest.push_back(get_vertex_attr(vid).pos);
    }
  }

  PointVec point_cloud() const {
    PointVec res;
    add_to_point_cloud(res);
    return res;
  }

  void replace_vertices_pos(const PointVec &pv);

  std::string vertex_str(VertexId id) const;
  std::string face_str(EdgeId id) const;
  std::string edge_str(FaceId id) const;
  std::string str() const;
  void whiten();

  // Mesh utils
  std::vector<std::shared_ptr<Mesh> > split_connected() const;
  SPTR(FaceCollection) to_faces() const;
  void to_faces(FaceCollection *fc) const;
  SPTR(Mesh) correct_orientation() const;

  opa::math::game::proto::MeshFaceGraphList &
  get_face_graph(opa::math::game::proto::MeshFaceGraphList *face_graph) const;
  SPTR(Mesh) close_mesh() const;
  SPTR(FaceGraph) compute_face_graph() const;
  SPTR(FaceGraph2) compute_face_graph2() const;
  SPTR(algo::FastGraph) compute_graph() const;
  std::vector<int> get_face_ids(FaceId face) const;
  Pos get_face_normal(FaceId face) const;

  OPA_ACCESSOR_R(utils::ObjectPool<Face>, m_faces, faces_raw);
  OPA_ACCESSOR_R(utils::Remapper<int>, m_rmp, rmp);

  const algo::FastGraph& graph()const;

  void compute_remap() const {
    m_rmp.reset();
    for (auto &v : m_vertices.used()) {
      m_rmp.get(v);
    }
  }
  bool is_correctly_oriented() const;

  bool is_connected() const;

protected:
  void init_edge(EdgeId id, const Pos &pos);
  mutable utils::Remapper<int> m_rmp;
  mutable glib::optional<algo::FastGraph> m_graph;

  // init face utils

  utils::ObjectPool<VertexAttribute> m_attrs;
  utils::ObjectPool<Vertex> m_vertices;
  utils::ObjectPool<Edge> m_edges;
  utils::ObjectPool<Face> m_faces;
  OrigMapType m_orig_map;

  // std::unique_ptr<graph::Graph<graph::VertexBase, graph::EdgeBase>> m_graph;
  friend class MeshExport;
  friend class MeshBuilder;
};
// OPA_DECL_SPTR(Mesh, MeshPtr);

class FaceDisplay {
  FaceDisplay(const Mesh *mesh, FaceId face) {
    this->mesh = mesh;
    this->face = face;
  }
  std::ostream &operator<<(std::ostream &os);

private:
  const Mesh *mesh;
  FaceId face;
};

struct MeshData {
  Mesh *mesh;
  glm::mat4 to_cam;

  void add_to_point_cloud(std::vector<Pos> *point_cloud) const {
    PointVec local_cloud = mesh->point_cloud();
    for (auto &pt : local_cloud) {
      point_cloud->push_back(ApplyMat(to_cam, pt));
    }
  }
};

OPA_NAMESPACE_DECL3_END

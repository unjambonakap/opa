#include <opa/math/common/stats.h>
#include <opa/math/game/geo_2d.h>
#include <opa/math/game/intersection.h>
#include <opa/math/game/mesh.h>
#include <opa/math/game/mesh_export.h>
#include <opa/math/game/mesh_util.h>
#include <opa/math/game/proto/common.pb.h>
#include <opa/utils/string.h>

// TODO: Fix get_one_vertex_attr
const float eps = 1e-6;
using namespace std;
using namespace opa::utils;
using namespace opa::math::game::proto;

OPA_NAMESPACE_DECL3(opa, math, game)
Mesh::~Mesh() { OPA_TRACES(this); }

Edge::Edge() { reset(); }
void Edge::reset() {
  start = VERTEX_NONE;
  face = FACE_NONE;
  opp = EDGE_NONE;
  next = EDGE_NONE;
  vertex_attr = ATTR_NONE;
}

void EdgeIndex::add(VertexId a, VertexId b, EdgeId e) { m_index[MP(a, b)] = e; }

EdgeId EdgeIndex::get(VertexId a, VertexId b) const {
  const auto &it = m_index.find(MP(a, b));
  if (it == m_index.end()) return EDGE_NONE;
  return it->second;
}

GraphObjId Mesh::create_attr() { return m_attrs.get_new(); }

void Mesh::init_edge(EdgeId id, const Pos &pos) {
  GraphObjId aid = get_edge(id).vertex_attr = m_attrs.get_new();
  OPA_CHECK0(aid != ATTR_NONE);
  get_attr(aid).pos = pos;
}

void Mesh::setup_face(FaceId fid, const std::vector<VertexId> &vids,
                      const EdgeIndex &index) {

  int n = vids.size();
  Face &face = get_face(fid);
  OPA_CHECK0(vids.size() > 0);

  REP (i, n) {
    VertexId v0id = vids[i];
    VertexId v1id = vids[(i + 1) % n];
    VertexId v2id = vids[(i + 2) % n];
    Vertex &v0 = get_vertex(v0id);
    Vertex &v1 = get_vertex(v1id);

    EdgeId e0id = index.get(v0id, v1id);
    EdgeId e1id = index.get(v1id, v0id);

    Edge &e0 = get_edge(e0id);
    Edge &e1 = get_edge(e1id);

    e0.start = v0id;
    e0.opp = e1id;
    e0.face = fid;
    e0.next = index.get(v1id, v2id);

    e1.opp = e0id;

    v0.edge = e0id;
    face.edge = e0id;
  }
}

Vertex &Mesh::get_vertex(VertexId id) { return m_vertices.get(id); }
Edge &Mesh::get_edge(EdgeId id) { return m_edges.get(id); }
Face &Mesh::get_face(FaceId id) { return m_faces.get(id); }
VertexAttribute &Mesh::get_attr(GraphObjId id) { return m_attrs.get(id); }

const Vertex &Mesh::get_vertex(VertexId id) const { return m_vertices.get(id); }
const Edge &Mesh::get_edge(EdgeId id) const { return m_edges.get(id); }
const Face &Mesh::get_face(FaceId id) const { return m_faces.get(id); }
const VertexAttribute &Mesh::get_attr(GraphObjId id) const {
  return m_attrs.get(id);
}

FaceId Mesh::face(EdgeId id) const { return get_edge(id).face; }
FaceId Mesh::opp_face(EdgeId id) const { return face(opp(id)); }
VertexId Mesh::end(EdgeId id) const { return start(next(id)); }
EdgeId Mesh::rot_edge(EdgeId id) const { return next(opp(id)); }
FaceId Mesh::get_edge_face(EdgeId id) const { return get_edge(id).face; }
VertexId Mesh::start(EdgeId id) const { return get_edge(id).start; }
EdgeId Mesh::opp(EdgeId id) const { return get_edge(id).opp; }
EdgeId Mesh::next(VertexId id) const { return get_edge(id).next; }
GraphObjId Mesh::attr(EdgeId id) const { return get_edge(id).vertex_attr; }

EdgeId Mesh::face_get_one_edge(FaceId id) const { return get_face(id).edge; }
EdgeId Mesh::vertex_get_one_edge(VertexId id) const {
  return get_vertex(id).edge;
}
FaceId Mesh::vertex_get_one_face(VertexId id) const {
  return get_edge(vertex_get_one_edge(id)).face;
}
GraphObjId Mesh::vertex_get_one_attr(VertexId id) const {
  return vertex_attr(vertex_get_one_edge(id));
}

GraphObjId Mesh::vertex_attr(VertexId vid, FaceId fid) const {
  return vertex_attr(find_vertex_edge(vid, fid));
}
VertexId Mesh::get_face_vertex(FaceId fid) const {
  return start(get_face(fid).edge);
}

GraphObjId Mesh::vertex_attr(EdgeId id) const {
  return get_edge(id).vertex_attr;
}

EdgeId Mesh::find_vertex_edge(VertexId vid, FaceId fid) const {
  for (auto eid : Mesh::Vertex_EdgeWalker(*this, vid))
    if (face(eid) == fid) return eid;
  return EDGE_NONE;
}

Pos Mesh::get_pos(VertexId id) const {
  OPA_DISP("query vert ", id);
  for (auto d : Vertex_EdgeWalker(*this, id)) {
    const Edge &e = get_edge(d);
    if (e.vertex_attr != ATTR_NONE) return get_attr(e.vertex_attr).pos;
  }
  OPA_CHECK0(false);
}

Dir Mesh::get_normal(FaceId face) const {
  auto walker = Mesh::Face_VertexWalker(*this, face);
  auto it = walker.begin();
  auto &a0 = get_attr((*it).ND);
  ++it;

  auto &a1 = get_attr((*it).ND);
  ++it;

  while (it != walker.end()) {
    auto &a2 = get_attr((*it).ND);
    ++it;
    if (!are_aligned(a0.pos, a1.pos, a2.pos))
      return get_plane_normal(a0.pos, a1.pos, a2.pos);
  }
  OPA_CHECK(false, "degenerate face");
}

std::vector<FaceId> Mesh::list_faces() const {
  return std::vector<FaceId>(ALL(m_faces.used()));
}
std::vector<EdgeId> Mesh::list_edges() const {
  return std::vector<EdgeId>(ALL(m_edges.used()));
}
std::vector<VertexId> Mesh::list_vertices() const {
  return std::vector<VertexId>(ALL(m_vertices.used()));
}

void Mesh::set_vertex_attr(VertexId vid, GraphObjId aid) {
  for (auto eid : Mesh::Vertex_EdgeWalker(*this, vid))
    get_edge(eid).vertex_attr = aid;
}

void Mesh::do_remove_edge(EdgeId id) const {}
void Mesh::do_remove_vertex(VertexId id) const {}
void Mesh::do_remove_face(FaceId id) const {}

void Mesh::init() { opa::utils::Initable::init(); }

std::string Mesh::vertex_str(VertexId id) const {
  const VertexAttribute &attr = get_vertex_attr(id);
  std::ostringstream os;
  os << stdsprintf("vid=%d,pos=", id) << attr.pos;
  os << ",el=[";
  for (auto e : Mesh::Vertex_EdgeWalker(*this, id)) os << e << ",";
  os << "],fl=[";
  for (auto f : Mesh::Vertex_FaceWalker(*this, id)) os << f << ",";
  os << "]";
  return os.str();
}

std::string Mesh::edge_str(EdgeId id) const {

  OPA_DISP("On edge ", id, opp(id), start(id), start(opp(id)), next(id));
  std::ostringstream os;
  // OPA_TRACES(get_vertex_attr(start(id)).pos);
  os << stdsprintf("eid=%d,start=%d,next=%d,opp=%d,face=%d,fopp=%d=", id,
                   start(id), next(id), opp(id), face(id), opp_face(id));
  return os.str();
}
std::string Mesh::face_str(FaceId id) const {

  OPA_DISP("on face ", id);
  std::ostringstream os;
  os << stdsprintf("fid=%d", id);
  os << ",el=[";
  for (auto e : Mesh::Face_EdgeWalker(*this, id)) os << e << ",";
  os << "],vl=[";
  for (auto f : Mesh::Face_VertexWalker(*this, id)) os << f.first << ",";
  os << "]";
  return os.str();
}

std::string Mesh::str() const {
  std::ostringstream os;
  os << "<<<< MESH >>>>\n";

  os << "\nEdges\n";
  for (auto e : list_edges()) os << edge_str(e) << endl;

  os << "Faces\n";
  for (auto f : list_faces()) os << face_str(f) << endl;

  os << "\nVertices\n";
  for (auto v : list_vertices()) os << vertex_str(v) << endl;

  os << "\n#### MESH ####\n";
  return os.str();
}

std::vector<SPTR(Mesh)> Mesh::split_connected() const {
  std::unordered_set<FaceId> seen;
  std::vector<SPTR(Mesh)> res;

  for (auto &fid : list_faces()) {
    std::queue<FaceId> q;

    if (seen.count(fid)) continue;
    seen.insert(fid);
    q.push(fid);

    MeshBuilder builder;
    while (!q.empty()) {
      FaceId cur = q.front();
      q.pop();
      std::vector<Pos> pos_list;
      for (auto eid : Mesh::Face_EdgeWalker(*this, cur)) {
        pos_list.push_back(get_attr(eid).pos);
        // OPA_DISP0(eid, pos_list.back());
        FaceId opf = opp_face(eid);
        if (opf != FACE_NONE) {
          if (!seen.count(opf)) {
            seen.insert(opf);
            q.push(opf);
          }
        } else {
        }
      }
      builder.add_face(pos_list);
    }
    res.emplace_back(builder.get_mesh());
  }
  return res;
}

opa::math::game::proto::MeshFaceGraphList &Mesh::get_face_graph(
  opa::math::game::proto::MeshFaceGraphList *face_graph) const {
  for (auto &fid : list_faces()) {
    MeshFaceGraph *face = face_graph->add_face();
    face->set_id(fid);
    for (auto &other : get_adj_faces(fid))
      if (other.first != FACE_NONE) face->add_adj(other.first);
    PosToProto(get_face_center(fid), face->mutable_pos());
  }
  return *face_graph;
}

void Mesh::whiten() {
  math::common::DataWhitener<Pos, 3> whitener;
  whitener.add(this->point_cloud());
  whitener.compute();
  OPA_DISP0(whitener.remap(this->point_cloud()));
  for (auto &attrid : m_attrs.used()) {
    auto &attr = m_attrs.get(attrid);
    attr.pos = whitener.remap(attr.pos);
  }
}

void Mesh::to_faces(FaceCollection *fc) const {
  auto faces = list_faces();

  for (auto &fid : faces) {
    auto tb = this->get_face_vertices(fid);
    fc->push_face(tb);
  }
}

SPTR(FaceCollection) Mesh::to_faces() const {
  SPTR(FaceCollection) res = std::make_unique<FaceCollection>();
  to_faces(res.get());
  return res;
}

std::vector<VertexId> Mesh::get_face_vertex_ids(FaceId face) const {
  std::vector<VertexId> tb;
  for (auto eid : Mesh::Face_EdgeWalker(*this, face)) {
    tb.push_back(start(eid));
  }
  return tb;
}

PointVec Mesh::get_face_vertices(FaceId face) const {
  PointVec tb;
  for (auto eid : Mesh::Face_EdgeWalker(*this, face)) {
    tb.push_back(get_attr(eid).pos);
  }
  return tb;
}

SPTR(Mesh) Mesh::close_mesh() const {
  /*
  auto trs = this->to_triangles();

  auto edges = list_edges();

  std::set<EdgeId> seen_edge;
  for (auto eid : edges) {
    FaceId f = face(eid);
    if (f != FACE_NONE)
      continue;
    if (seen_edge.count(eid)) continue;
    std::vector<EdgeId> lst;

    EdgeId cur = eid;
    while(true){
      cur = get_attr(cur).next;
      lst.push_back(cur);
      if (cur == eid) break;
    }

    OPA_CHECK0(lst.size()>=2);
    lst.push_back(lst[0]);
    lst.push_back(lst[1]);
    REP(i, lst.size()-2){

    }
  }
  return trs->to_mesh();
  */
  return nullptr;
}

std::vector<int> Mesh::get_face_ids(FaceId face) const {
  std::vector<int> res;
  for (auto f : Mesh::Face_VertexWalker(*this, face)) {
    res.push_back(rmp().get_or_die(f.first));
  }
  return res;
}

SPTR(algo::FastGraph) Mesh::compute_graph() const {
  SPTR(algo::FastGraph)
  graph =
    std::make_shared<algo::FastGraph>(m_vertices.size(), algo::Mode::REMAP);

  for (auto face : this->faces()) {
    auto tb = get_face_vertex_ids(face);
    tb.push_back(tb[0]);
    REP (i, tb.size() - 1) { 
      graph->add_bidirectional(tb[i], tb[i + 1]); }
  }

  return graph;
}

Pos Mesh::get_face_normal(FaceId fid) const {
  PointVec tb;
  for (auto &vid : get_face_vertex_ids(fid)) {
    tb.push_back(this->get_vertex_attr(vid).pos);
  }
  return face_normal(tb);
}

PlaneSpec Mesh::get_face_plane(FaceId face) const {
  return PlaneSpec::FromPoints(get_face_vertices(face));
}

bool Mesh::is_correctly_oriented() const {
  utils::MinFinderPair<Pos, VertexId, GlmVecXYCmp<Pos> > minfinder{
    GlmVecXYCmp<Pos>(base_eps)
  };

  for (auto &vid : vertices()) {
    minfinder.update(get_vertex_attr(vid).pos, vid);
  }
  FaceId fid = vertex_get_one_face(minfinder.get());
  Dir normal = get_face_normal(fid);
  OPA_CHECK0(glm::length2(normal.xy()) > base_eps);

  if (std::abs(normal.x) < base_eps) return normal.y < 0;
  return normal.x < 0;
}

SPTR(Mesh) Mesh::correct_orientation() const {
  SPTR(FaceCollection) fc = this->to_faces();
  if (!is_correctly_oriented()) fc->toggle_orientation();
  return SPTR(Mesh)(fc->to_mesh());
}

double Mesh::area() const {
  double res = 0;
  for (auto face : this->faces()) {
    auto tb = get_face_vertices(face);
    auto normal = get_face_normal(face);
    double sgn = OPA_BITSIGN(glm::dot(tb[0], normal) > 0);
    res += tetrahedron_area(tb[0], tb[1], tb[2]);

    OPA_CHECK0(tb.size() == 3);
  }
  return res;
}
std::vector<std::pair<FaceId, EdgeId> > Mesh::get_adj_faces(FaceId face) const {
  std::vector<std::pair<FaceId, EdgeId> > res;
  for (auto eid : Mesh::Face_EdgeWalker(*this, face)) {
    FaceId opf = opp_face(eid);
    res.emplace_back(opf, eid);
  }
  return res;
}

Pos Mesh::get_face_center(FaceId face) const {
  int count = 0;
  Pos center;
  for (auto eid : Mesh::Face_EdgeWalker(*this, face)) {
    center += get_attr(eid).pos;
    ++count;
  }
  center /= count;
  return center;
}

SPTR(FaceGraph2) Mesh::compute_face_graph2() const {
  SPTR(FaceGraph2) res = std::make_shared<FaceGraph2>();
  res->graph = std::make_shared<algo::FastGraph>();

  utils::Remapper<FaceId> face_rmp;
  std::unordered_map<FaceId, Pos> face_to_center;
  for (auto &fid : list_faces()) {
    face_rmp.get(fid);
    face_to_center[fid] = get_face_center(fid);
  }

  res->graph->reset(face_to_center.size(), 0);
  for (const FaceId &fid : list_faces()) {
    for (const auto &other : get_adj_faces(fid)) {
      if (other.first == FACE_NONE) continue;
      if (fid < other.first) {
        int eid = res->graph->add_bidirectional(face_rmp.rget(fid),
                                                face_rmp.rget(other.first));

        double dist = face_face_dist(face_to_center[fid],
                                     face_to_center[other.first], other.second);
        OPA_CHECK(res->edge_costs.size() == eid, res->edge_costs.size(), eid);
        res->edge_costs.push_back(dist);
      }
    }
  }
  return res;
}

SPTR(FaceGraph) Mesh::compute_face_graph() const {
  SPTR(FaceGraph) res = std::make_shared<FaceGraph>();
  res->graph = std::make_shared<lemon::SmartGraph>();
  res->edge_map =
    std::make_shared<lemon::SmartGraph::ArcMap<double> >(*res->graph);

  utils::Remapper<FaceId> face_rmp;
  std::unordered_map<FaceId, Pos> face_to_center;
  for (auto &fid : list_faces()) {
    face_rmp.get(fid);
    face_to_center[fid] = get_face_center(fid);
    res->graph->addNode();
  }

  for (const FaceId &fid : list_faces()) {
    for (const auto &other : get_adj_faces(fid)) {
      if (other.first == FACE_NONE) continue;
      if (fid < other.first) {
        auto edge = res->graph->addEdge(
          res->graph->nodeFromId(face_rmp.rget(fid)),
          res->graph->nodeFromId(face_rmp.rget(other.first)));

        double dist = face_face_dist(face_to_center[fid],
                                     face_to_center[other.first], other.second);

        res->edge_map->set(res->graph->direct(edge, 0), dist);
        res->edge_map->set(res->graph->direct(edge, 1), dist);
      }
    }
  }

  return res;
}

double Mesh::face_face_dist(const Pos &f1_center, const Pos &f2_center,
                            const EdgeId &edge) const {

  Pos edge_p1 = get_attr_for_edge(edge).pos;
  Pos edge_p2 = get_attr_for_edge(opp(edge)).pos;
  Pos norm_f1 = get_plane_normal(f1_center, edge_p1, edge_p2);
  Pos norm_f2 = get_plane_normal(f2_center, edge_p1, edge_p2);
  double angle = -vec_get_angle(norm_f1, norm_f2);
  Pos proj_p2 =
    quat_from_vec_rot_safe(edge_p2 - edge_p1, angle) * (f2_center - edge_p1) +
    edge_p1;
  return glm::distance(proj_p2, f1_center);
}

bool Mesh::is_connected() const { return graph().is_connected(); }

const algo::FastGraph &Mesh::graph() const {
  if (!m_graph.has_value()) {
    m_graph = *compute_graph();
  }
  return m_graph.value();
}

void Mesh::replace_vertices_pos(const PointVec &pv) {
  OPA_CHECK0(pv.size() == vertices().size());
  auto it1 = pv.begin();
  auto it2 = vertices().begin();
  for (; it1 != pv.end(); ++it1, ++it2) {
    get_vertex_attr(*it2).pos = *it1;
  }
}

OPA_NAMESPACE_DECL3_END

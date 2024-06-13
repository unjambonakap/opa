#include <opa/math/game/mesh_util.h>

#include <opa/math/common/stats.h>
#include <opa/math/game/geo_2d.h>
#include <opa/math/game/intersection.h>
#include <opa/math/game/mesh_export.h>
#include <opa/math/game/point_cloud.h>
#include <opa/or/grid_search.h>
#include <opa/utils/buffer_reader.h>

#define SQUARE_IDS(a, b, c, d)                                                 \
  { a, b, c }, { a, c, d }
#define SQUARE_IDS_V(a, b, c, d) SQUARE_IDS(d, c, b, a)

#define BOX_IDS(a, b, c, d, e, f, g, h)                                        \
  SQUARE_IDS_V(a, b, c, d)                                                     \
  , SQUARE_IDS(e, f, g, h), SQUARE_IDS(a, b, f, e), SQUARE_IDS(b, c, g, f),    \
    SQUARE_IDS(c, d, h, g), SQUARE_IDS(d, a, e, h),

#define SQUARE_FACE(a, b, c, d)                                                \
  { a, b, c, d }
#define ISQUARE_FACE(a, b, c, d) SQUARE_FACE(d, c, b, a)
#define BOX_FACES(a, b, c, d, e, f, g, h)                                      \
  ISQUARE_FACE(a, b, c, d)                                                     \
  , SQUARE_FACE(e, f, g, h), SQUARE_FACE(a, b, f, e), SQUARE_FACE(b, c, g, f), \
    SQUARE_FACE(c, d, h, g), SQUARE_FACE(d, a, e, h)

constexpr double eps = 1e-4;
OPA_NAMESPACE_DECL3(opa, math, game)

const FaceIndexList square_tr_id = { SQUARE_IDS(0, 1, 3, 2) };
const FaceIndexList box_tr_id = { BOX_IDS(0, 1, 3, 2, 4, 5, 7, 6) };
const FaceIndexList box_face_id = { BOX_FACES(0, 1, 3, 2, 4, 5, 7, 6) };
const FaceIndexList square_face_id = { SQUARE_FACE(0, 1, 3, 2) };

template <typename T>
std::vector<T> get_single_path(const std::vector<T> &container,
                               const std::function<bool(const T &x)> &pred) {
  std::vector<T> path;
  FEV (it, container) {
    if (!pred(*it)) break;
    path.push_back(*it);
  }
  OPA_CHECK0(path.size() != container.size());
  std::reverse(ALL(path));
  FE (it, container) {
    if (pred(*it))
      path.push_back(*it);
    else if (!path.empty())
      break;
  }
  return path;
}

std::vector<int> KDTree::get_neighborhood(const Pos &pt, double dist) const {
  std::vector<int> res;
  this->find_neighborhood(m_root, pt, dist, res);
  return res;
}

void KDTree::setup(const PointVec &tb) {
  m_tb = tb;
  BuildState state = BuildState::FromPoints(m_tb);
  m_root = this->build(state);
}

void KDTree::find_neighborhood(KDId id, const Pos &pt, double dist,
                               std::vector<int> &neighbours) const {
  const KDNode &cur = m_kds.get(id);
  if (cur.is_leaf()) {
    double cur_dist = glm::distance(m_tb[cur.point_id()], pt);
    if (cur_dist < dist) neighbours.push_back(cur.point_id());
    return;
  } else {
    // OPA_TRACES(pt, dist, cur.box.min_dist2(pt), cur.box.low, cur.box.high);
    if (cur.box.min_dist(pt) > dist) return;
  }

  find_neighborhood(cur.left, pt, dist, neighbours);
  find_neighborhood(cur.right, pt, dist, neighbours);
}
void KDTree::find_closest(KDId id, const Pos &pt,
                          GlobalClosestStatus &status) const {
  const KDNode &cur = m_kds.get(id);
  if (cur.is_leaf()) {
    status.best.update(
      { glm::distance(m_tb[cur.point_id()], pt), cur.point_id() });
    return;
  }

  KDId left = cur.left, right = cur.right;
  if (!status.best.has()) {
    const KDNode &left_node = m_kds.get(cur.left);
    if (!left_node.box.in(pt)) std::swap(left, right);
  } else {
    if (cur.box.min_dist(pt) >= status.best.get().first) return;
  }
  find_closest(left, pt, status);
  find_closest(right, pt, status);
}

int KDTree::get_closest(const Pos &pt) const {
  KDTree::GlobalClosestStatus status;
  this->find_closest(m_root, pt, status);
  return status.best.get().second;
}

KDTree::KDId KDTree::build(const KDTree::BuildState &state) {
  KDId left = utils::InvalidId, right = utils::InvalidId;

  double splitv;
  int splitdim;
  std::tie(splitv, splitdim) = this->pick_split_dim(state);
  // OPA_DISP("Building ", state.lst[0], splitv, splitdim);
  if (splitdim != -1) {
    BuildState left_state, right_state;
    state.split2(m_tb, splitdim, splitv, left_state, right_state);
    // OPA_DISP("Res of split ", left_state.lst[0], right_state.lst[0]);

    left = this->build(left_state);
    right = this->build(right_state);

  } else {
    right = state.lst[0][0];
    for (auto &ix : state.lst[0]) {
      double dist = glm::length(m_tb[right] - m_tb[ix]);
      OPA_CHECK(dist < 2 * eps, dist, m_tb[right], m_tb[ix], eps);
    }
  }

  KDNode &node = m_kds.get_new2();
  node.box = state.compute_box(m_tb);
  /*
  OPA_DISP0(node.box.low, node.box.high);
  REP (j, 3) {
    OPA_DISP("on coord ", j);
    for (auto &e : state.lst[j])
      OPA_DISP0(m_tb[e][j]);
  }
  */
  node.left = left;
  node.right = right;
  return node.id;
}

std::pair<double, int> KDTree::pick_split_dim(const BuildState &state) const {
  if (state.n() == 1) return { 0, -1 };
  utils::MaxFinder<std::pair<double, int> > range_and_dim;
  // for (auto &e : state.lst[0]) {
  //  OPA_DISP0(m_tb[e]);
  //}

  REP (j, 3) {
    // OPA_TRACES(-state.get_range(j, m_tb).get_length(), j);
    range_and_dim.update({ state.get_range(j, m_tb).get_length(), j });
  }
  auto best = range_and_dim.get();
  if (best.first < eps) return { 0, -1 };
  return { state.get_range(best.second, m_tb).get_middle(), best.second };
}

bool PointMatcher::load_batch(const PointVec &tb) {
  cur_batch = tb;
  if (tb.size() == 0) return true;
  // OPA_DISP("Loading batch ", tb.size());
  m_tree = KDTree();
  OPA_TRACEM("start setup");
  m_tree.setup(tb);
  m_uj.reset(tb.size());
  std::multiset<int> done;
  OPA_DISP0(tb.size());

  REP (i, tb.size()) {
    if (done.count(i)) continue;
    std::vector<int> group = m_tree.get_neighborhood(tb[i], 3 * eps);
    OPA_CHECK0(group.size() > 0);
    for (auto &e : group) {
      // if (m_uj.size(e) != 1) return false;
      done.insert(e);
    }

    REP (j, group.size()) { m_uj.merge(i, group[j]); }
  }
  OPA_DISP("after loading ", this->final_count());
  return true;
}

int PointMatcher::query(const Pos &pos) const {
  int pt = m_tree.get_closest(pos);
  return m_uj.root(pt);
}

EdgeId next_boundary_edge(const Mesh &mesh,
                          const std::set<FaceId> &visible_faces, EdgeId cur) {
  cur = mesh.next(cur);

  while (true) {
    if (!visible_faces.count(mesh.opp_face(cur))) return cur;
    cur = mesh.rot_edge(cur);
  }
}

VertexId MeshBuilder::get_vid(const Pos &pos) {
  if (!m_pos_vid.hasu(pos)) {
    VertexId id = m_mesh->m_vertices.get_new();
    m_pos_vid.add(pos, id);
  }
  return m_pos_vid.uv(pos);
}

EdgeId MeshBuilder::get_edge(VertexId v0, VertexId v1) {
  if (!m_edge_index.has(v0, v1)) {
    Edge *e0, *e1;
    m_mesh->m_edges.get_new(2, &e0, &e1);
    m_mesh->init_edge(e0->id, m_pos_vid.vu(v0));
    m_mesh->init_edge(e1->id, m_pos_vid.vu(v1));
    m_edge_index.add(v0, v1, e0->id);
    m_edge_index.add(v1, v0, e1->id);
    // OPA_DISP("Adding edge ", e0->id, v0, m_pos_vid.vu(v0), e0->start);
    // OPA_DISP("Adding edge ", e1->id, v1, m_pos_vid.vu(v1), e1->start);
  }
  return m_edge_index.get(v0, v1);
}

void MeshBuilder::add_triangle(const Pos &p0, const Pos &p1, const Pos &p2) {
  add_face({ p0, p1, p2 });
}

void MeshBuilder::add_face(const std::vector<Pos> &pos) {
  OPA_CHECK0(pos.size() > 0);
  m_faces.push_back(pos);
}

void MeshBuilder::compute_face(const PointVec &pv) {
  std::vector<VertexId> vids;
  for (auto &e : pv) vids.push_back(get_vid(e));
  std::vector<EdgeId> eids;
  for (int i = 0; i < vids.size(); ++i)
    eids.push_back(get_edge(vids[i], vids[(i + 1) % vids.size()]));
  // OPA_TRACES(eids);
  FaceId fid = m_mesh->m_faces.get_new();
  m_mesh->setup_face(fid, vids, m_edge_index);
}

void MeshBuilder::load_triangles(const std::vector<Triangle> &triangles) {
  for (auto &triangle : triangles) {
    add_triangle(triangle);
  }
}

void MeshBuilder::fill(const FaceCollection &fc) {
  for (const auto &face : fc.faces()) {
    add_face(face);
  }
}

void MeshBuilder::compute() {
  prepare();
  for (auto &face : m_faces) {
    PointVec normed;
    for (auto &pt : face)
      normed.push_back(m_matcher.cur_batch[m_matcher.query(pt)]);
    // puts("==============");

    // OPA_DISP0(face);
    // OPA_DISP0(normed);
    // OPA_DISP0(tr_normal(vec_to_tr(normed)));
    // OPA_DISP0(face_normal_convex(normed));
    compute_face(normed);
  }
}

void MeshBuilder::prepare() {
  if (m_prepared) return;
  m_prepared = true;
  std::vector<Pos> pos;
  m_matcher.reset();
  for (auto &face : m_faces) {
    pos.insert(pos.end(), ALL(face));
  }
  OPA_CHECK0(m_matcher.load_batch(pos));
}

void MeshBuilder::repair_orientation() {
  // prepare();
  // std::vector<PointVec> remapped_faces;
  // for (auto &face : m_faces) {
  //  PointVec normed;
  //  for (auto &pt : face)
  //    normed.push_back(m_matcher.cur_batch[m_matcher.query(pt)]);
  //  remapped_faces.push_back(normed);
  //}
  // utils::MinFinderPair<Pos, int, GlmVecXYCmp<Pos> > minfinder{
  //  GlmVecXYCmp<Pos>(base_eps)
  //};

  // REP(i, remapped_faces.size()){
  //  for (auto &pt : remapped_faces[i]){
  //  minfinder.update(pt, i);
  //}
  // Dir normal = face_normal_convex(remapped_faces[minfinder.get()]);
  // OPA_CHECK0(glm::length2(normal.xy()) > base_eps);

  // if (std::abs(normal.x) < base_eps) return normal.y < 0;
  // return normal.x < 0;
}

bool MeshBuilder::is_manifold_mesh() {
  auto graph = get_graph();
  for (const auto &edge : graph->graph.list_edges()) {
    if (graph->graph.get_edge_count(edge.from, edge.to) != 1) return false;
  }
  return true;
}

SPTR(GraphWithPosData) MeshBuilder::get_graph() {
  if (m_graph) return m_graph;
  prepare();
  SPTR(GraphWithPosData) res = std::make_shared<GraphWithPosData>();
  utils::Remapper<int> reprs = m_matcher.uj().get_repr();

  for (auto &repr : reprs.mp()) {
    res->points.push_back(m_matcher.cur_batch[repr]);
  }
  res->graph.reset(reprs.size(), algo::Mode::DIGRAPH | algo::Mode::MULTIGRAPH);

  OPA_DISP0(m_faces.size());
  for (auto &face : m_faces) {
    std::vector<int> pos;
    for (auto &pt : face) pos.push_back(reprs.get(m_matcher.query(pt)));
    res->faces.push_back(pos);

    REP (i, pos.size()) {
      int a = pos[i];
      int b = pos[(i + 1) % pos.size()];
      // OPA_DISP("adding edge ", a, b, face[i], face[(i + 1) % pos.size()]);
      res->graph.adde(a, b);
    }
  }
  m_graph = res;
  return res;
}

Mesh *MeshBuilder::get_mesh() {
  if (!m_computed) compute();
  m_computed = true;
  // OPA_DISP0(m_mesh->str());
  return m_mesh;
}

void FaceCollection::whiten() {
  math::common::DataWhitener<Pos, 3> whitener;
  for (auto &face : m_faces) whitener.add(face);
  whitener.compute();
  for (auto &face : m_faces) face = whitener.remap(face);
}

void FaceCollection::toggle_orientation() {
  for (auto &face : m_faces) std::reverse(ALL(face));
}

FaceCollection &FaceCollection::add_plane(const PlaneSpec &plane,
                                          const Pos &center, double dist,
                                          const Dir &base_dir) {
  Dir u = base_dir;
  if (base_dir == vec_0) u = plane.get_rand_plane_vec();
  Dir v = glm::cross(u, plane.dir);
  u *= dist;
  v *= dist;
  Pos center2 = center - u * 0.5 - v * 0.5;
  PointVec face;

  REP (i, 2) {
    REP (j, 2) { face.push_back(center2 + u * i + v * j); }
  }
  return push_indexed(face, square_face_id);
}

FaceCollection &FaceCollection::push_indexed(const PointVec &points,
                                             const FaceIndexList &index_list) {
  for (auto &ids : index_list) {
    PointVec face;
    for (auto &id : ids) face.push_back(points[id]);
    this->push_face(face);
  }
  return *this;
}

FaceCollection &FaceCollection::add_box(const Box2DSpec &box,
                                        const PlaneSpec &plane,
                                        const Dir &front) {
  PointVec vertices;
  REP (i, 4)
    vertices.push_back(plane.lift(box.get(i), front));
  push_indexed(vertices, square_face_id);
  return *this;
}

FaceCollection &FaceCollection::add_box(const BoxSpec &box) {
  PointVec vertices;
  REP (i, 8)
    vertices.push_back(box.get(i));
  push_indexed(vertices, box_face_id);
  return *this;
}

FaceCollection &FaceCollection::load_stl_from_data(std::string_view data) {
  utils::BufferReader reader(data);
  reader.get(80);
  u32 ntr = reader.read_u32();
  OPA_CHECK(!data.starts_with("solid"), "Only binary stl are supported");

  MeshBuilder res;
  REP (i, ntr) {
    reader.readn<f32>(3);
    PointVec face;
    REP (j, 3) {
      std::vector<float> tmp = reader.readn<f32>(3);
      face.emplace_back(tmp[0], tmp[1], tmp[2]);
    }
    this->push_face(face);

    reader.read<u16>();
  }
  return *this;
}

FaceCollection &FaceCollection::load_stl(std::string_view filename) {
  std::string content = utils::read_file(filename);
  return load_stl_from_data(content);
}

FaceCollection &FaceCollection::filter(mesh_face_filter_t filter_func) {
  REP (i, m_faces.size()) {
    if (!filter_func(m_faces[i])) {
      std::swap(m_faces[i], m_faces.back());
      m_faces.pop_back();
    }
  }
  return *this;
}

SPTR(MeshBuilder) FaceCollection::to_mesh_builder(Mesh *mesh) const {
  SPTR(MeshBuilder) builder(new MeshBuilder(mesh));
  builder->fill(*this);
  return builder;
}

Mesh *FaceCollection::to_mesh(Mesh *mesh) const {
  return to_mesh_builder(mesh)->get_mesh();
}

FaceCollection &FaceCollection::push_face(const PointVec &face) {
  OPA_CHECK0(face.size() >= 3);
  if (!m_triangulate || face.size() == 3) {
    m_faces.push_back(face);
    return *this;
  } else {
    return FaceCollection::push(triangulate_polygon(face));
  }
}

FaceCollection FaceCollection::triangulate() const {
  FaceCollection res(true);
  for (auto &face : m_faces) {
    res.push_face(face);
  }
  return res;
}

void FaceCollection::write_stl(utils::BufferWriter &writer) const {
  writer.write_repeated<u8>(0, 80);
  writer.write_u32(m_faces.size());

  for (const auto &tr : m_faces) {
    OPA_CHECK0(tr.size() == 3);
    writer.write_repeated<f32>(0, 3);
    int count = 0;
    writer << tr[0] << tr[1] << tr[2];
    writer.write<u16>(0);
  }
}

void MeshBuilder::FromTetraedra(Mesh *mesh,
                                const std::vector<VertexInfoData> &tb) {
  OPA_CHECK0(tb.size() == 4);

  std::vector<VertexId> vertices = mesh->m_vertices.get_new_vec(4);
  EdgeIndex index;

  REP (i, 4)
    mesh->orig_map()[vertices[i]] = tb[i].id;

  REP (i, 4)
    REP (j, i) {
      Edge *e0, *e1;
      mesh->m_edges.get_new(2, &e0, &e1);
      mesh->init_edge(e0->id, tb[i].pos);
      mesh->init_edge(e1->id, tb[j].pos);

      index.add(vertices[i], vertices[j], e0->id);
      index.add(vertices[j], vertices[i], e1->id);
    }

  REP (i, 4) {
    int i0 = (i + 0) % 4;
    int i1 = (i + 1) % 4;
    int i2 = (i + 2) % 4;
    int i3 = (i + 3) % 4;

    Dir normal = get_plane_normal(tb[i0].pos, tb[i1].pos, tb[i2].pos);

    FaceId fid = mesh->m_faces.get_new();

    std::vector<VertexId> vl = { vertices[i0], vertices[i1], vertices[i2] };

    if (point_signed_dist_to_plane(tb[i3].pos - tb[i0].pos, normal) > 0)
      std::reverse(ALL(vl));
    mesh->setup_face(fid, vl, index);
  }
}

void MeshBuilder::PolyhedraFromCloud(Mesh *mesh,
                                     const std::vector<VertexInfoData> &cloud) {
  int max_z = 0, min_z = 0;
  int n = cloud.size();

  std::vector<int> used(n, 0);
  REP (i, n) {
    OPA_DISP0(cloud[i].pos);
    if (cloud[i].pos.z < cloud[min_z].pos.z) min_z = i;
    if (cloud[i].pos.z > cloud[max_z].pos.z) max_z = i;
  }
  used[min_z] = 1;
  used[max_z] = 1;
  OPA_CHECK0(min_z != max_z);

  int best1 = -1;
  REP (i, n)
    if (!used[i]) {
      if (best1 == -1 ||
          dist_line(cloud[i].pos, cloud[min_z].pos, cloud[max_z].pos) >
            dist_line(cloud[best1].pos, cloud[min_z].pos, cloud[max_z].pos))
        best1 = i;
    }
  OPA_CHECK0(best1 != -1);

  used[best1] = 1;
  const Pos &p0 = cloud[min_z].pos;
  const Pos &p1 = cloud[max_z].pos;
  const Pos &p2 = cloud[best1].pos;
  const Dir plane_normal = get_plane_normal(p0, p1, p2);

  int best2 = 0;
  REP (i, n)
    if (!used[i]) {
      if (point_dist_to_plane(cloud[i].pos - p0, plane_normal) >
          point_dist_to_plane(cloud[best2].pos - p0, plane_normal))
        best2 = i;
    }
  OPA_CHECK0(!used[best2]);
  used[best2] = 1;

  std::vector<VertexId> tmp1, tmp2;
  std::vector<int> remap1 = { min_z, max_z, best1, best2 };

  std::vector<VertexInfoData> cloud1;
  for (auto x : remap1) cloud1.pb(cloud[x]);
  FromTetraedra(mesh, cloud1);

  std::vector<VertexInfoData> cloud2;
  REP (i, n)
    if (!used[i]) cloud2.pb(cloud[i]);

  ConvexPolyhedra_Extend(mesh, cloud2);
}

bool MeshBuilder::ConvexPolyhedra_IsInside(const Mesh &mesh, const Pos &pos) {
  for (const auto &x : mesh.faces()) {
    Dir normal = mesh.get_normal(x);
    auto &vpos = mesh.get_vertex_attr(mesh.get_face_vertex(x)).pos;
    if (OPA_FLOAT_GT(point_signed_dist_to_plane(pos - vpos, normal), 0, eps))
      return false;
  }
  return true;
}

void MeshBuilder::PolygonFromCloud(Mesh *mesh,
                                   const std::vector<VertexInfoData> &cloud) {
  int n = cloud.size();
  OPA_ASSERT0(n >= 3);
  int p0, p1, p2;
  p0 = 0;
  p1 = 1;
  p2 = -1;

  FOR (i, 2, n)
    if (!are_aligned(cloud[p1].pos, cloud[p2].pos, cloud[i].pos)) {
      p2 = i;
      break;
    }
  OPA_ASSERT0(p2 != -1);

  std::vector<Pos2> proj(n);
  Dir front = glm::normalize(cloud[p1].pos - cloud[p0].pos);
  Dir normal = get_plane_normal(cloud[p0].pos, cloud[p1].pos, cloud[p2].pos);

  FOR (i, 3, cloud.size()) {
    OPA_ASSERT0(OPA_FLOAT_EQ(
      point_dist_to_plane(cloud[i].pos - cloud[p0].pos, normal), 0, eps));
  }

  REP (i, n)
    proj[i] = get_plane_coord(cloud[i].pos - cloud[p0].pos, normal, front);

  auto hull = compute_convex_hull(proj);

  int nx = hull.size();
  hull.pb(hull[0]);

  std::vector<VertexId> vertices = mesh->m_vertices.get_new_vec(nx);
  std::vector<EdgeId> edges = mesh->m_edges.get_new_vec(nx * 2);
  vertices.pb(vertices[0]);

  FaceId nf = mesh->m_faces.get_new();

  EdgeIndex idx;

  REP (i, nx) {
    mesh->orig_map()[vertices[i]] = hull[i];
    mesh->init_edge(edges[2 * i], cloud[hull[i]].pos);
    mesh->init_edge(edges[2 * i + 1], cloud[hull[i + 1]].pos);

    idx.add(vertices[i], vertices[i + 1], edges[2 * i]);
    idx.add(vertices[i + 1], vertices[i], edges[2 * i + 1]);
  }
  vertices.pop_back();
  mesh->setup_face(nf, vertices, idx);
}

void MeshBuilder::ConvexPolyhedra_Extend(
  Mesh *mesh, const std::vector<VertexInfoData> &cloud) {
  OPA_DISP("Status >>> ", mesh->str());

  struct PointInfo {
    PointInfo() {}
    std::set<FaceId> visible_faces;
    std::set<FaceId> aligned_faces;
  };
  struct FaceInfo {
    std::set<int> bad_points;
    FaceId id;
    bool operator<(const FaceInfo &a) const {
      if (bad_points.size() != a.bad_points.size())
        return bad_points.size() > a.bad_points.size();
      return id < a.id;
    }
  };
  OPA_DECL_SPTR(FaceInfo, FaceInfoPtr);

  int n = cloud.size();

  std::vector<PointInfo> info(n);
  std::map<FaceId, FaceInfoPtr> cur_faces;
  std::map<VertexId, int> m_point_map;
  std::set<int> rem_points;
  std::set<FaceInfoPtr, utils::PtrComparator<FaceInfoPtr> > faces_rank;
  REP (i, n)
    rem_points.insert(i);

  std::map<VertexId, Pos> mp;

  for (auto face_id : mesh->list_faces()) {
    OPA_DISP("QUERY FACE ", face_id);
    auto pt = mesh->get_pos(mesh->get_face_vertex(face_id));
    auto normal = mesh->get_normal(face_id);
    FaceInfoPtr ptr(new FaceInfo);
    ptr->id = face_id;

    for (auto x : rem_points) {
      if (OPA_FLOAT_GT(point_signed_dist_to_plane(cloud[x].pos - pt, normal), 0,
                       eps)) {
        info[x].visible_faces.insert(face_id);
        ptr->bad_points.insert(x);
      } else if (OPA_FLOAT_EQ(
                   point_signed_dist_to_plane(cloud[x].pos - pt, normal), 0,
                   eps)) {
        info[x].aligned_faces.insert(face_id);
      }
    }
    cur_faces[face_id] = ptr;
    faces_rank.insert(ptr);
  }

  REP (i, n)
    if (info[i].visible_faces.size() == 0) rem_points.erase(i);

  while (rem_points.size()) {
    OPA_CHECK0(faces_rank.size() > 0);
    FaceInfoPtr face_pick = *faces_rank.begin();
    printf("picking ID=%d\n", face_pick->id);
    auto pt = mesh->get_pos(mesh->get_face_vertex(face_pick->id));
    auto normal = mesh->get_normal(face_pick->id);
    const auto &cnd = face_pick->bad_points;
    int best =
      *std::max_element(ALL(cnd), [&cloud, &pt, &normal](int a, int b) {
        return point_dist_to_plane(cloud[a].pos - pt, normal) <
               point_dist_to_plane(cloud[b].pos - pt, normal);
      });

    EdgeId boundary_edge = EDGE_NONE;
    auto &visible_faces = info[best].visible_faces;

    for (auto see : visible_faces) {
      for (auto eid : Mesh::Face_EdgeWalker(*mesh, see)) {
        FaceId other = mesh->opp_face(eid);

        // not boundary edge -> go on
        if (!visible_faces.count(other)) {
          boundary_edge = eid;
          goto out;
        }
      }
    }
  out:
    OPA_CHECK0(boundary_edge != EDGE_NONE);
    std::vector<EdgeId> boundary_edges;
    boundary_edges.pb(boundary_edge);

    puts("\n");
    puts("\n");
    puts("========================");
    puts("========================");
    puts("========================");
    puts("========================");
    OPA_DISP("picking >> ", cloud[best].pos);
    for (auto tmp : visible_faces) OPA_DISP("can see >> ", tmp);
    OPA_DISP0(visible_faces);
    OPA_DISP0(info[best].aligned_faces);

    while (true) {
      OPA_DISP("update face ", mesh->face(boundary_edge), boundary_edge);
      mesh->get_face(mesh->face(boundary_edge)).edge = boundary_edge;
      mesh->get_vertex(mesh->start(boundary_edge)).edge = boundary_edge;

      boundary_edge = next_boundary_edge(*mesh, visible_faces, boundary_edge);

      if (boundary_edge == boundary_edges[0]) break;
      boundary_edges.pb(boundary_edge);
    }
    for (auto x : boundary_edges)
      printf("got boudnary edge >> %d %d\n", mesh->start(x), mesh->end(x));

    std::map<VertexId, EdgeId> vertex_to_edge;
    for (auto &x : boundary_edges) vertex_to_edge[mesh->start(x)] = x;

    std::vector<EdgeId> aligned_path =
      get_single_path<EdgeId>(boundary_edges, [&](const EdgeId &eid) -> bool {
        return info[best].aligned_faces.count(
          mesh->get_edge(mesh->opp(eid)).face);
      });

    std::set<VertexId> removed_boundary_vertices;
    for (int i = 1; i < aligned_path.size(); ++i)
      removed_boundary_vertices.emplace(mesh->start(aligned_path[i]));

    std::vector<VertexId> boundary_vertices;

    for (const auto &x : boundary_edges) {
      VertexId v = mesh->start(x);
      if (removed_boundary_vertices.count(v)) continue;
      boundary_vertices.push_back(v);
    }

    std::set<EdgeId> boundary_edges_set(ALL(boundary_edges));
    std::set<VertexId> boundary_vertices_set(ALL(boundary_vertices));

    std::set<FaceId> aligned_faces;
    for (auto &aligned_edge : aligned_path) {
      FaceId fid = mesh->face(mesh->opp(aligned_edge));
      aligned_faces.insert(fid);
    }

    OPA_CHECK0(aligned_faces.size() <= 2);
    std::set<std::pair<VertexId, VertexId> > blacklist_boundary_edge;
    if (aligned_faces.size() >= 1) {
      boundary_vertices.erase(
        std::remove_if(ALL(boundary_vertices),
                       [&](const VertexId &vid) {
                         return removed_boundary_vertices.count(vid);
                       }),
        boundary_vertices.end());
      int vstart = mesh->start(aligned_path[0]);
      int vend = mesh->end(aligned_path.back());
      blacklist_boundary_edge.emplace(vstart, vend);
      blacklist_boundary_edge.emplace(vend, vstart);
    }

    int nx = boundary_vertices.size();
    VertexId new_vertex = mesh->m_vertices.get_new();
    mesh->orig_map()[new_vertex] = cloud[best].id;

    EdgeIndex idx;
    REP (i, nx) {
      EdgeId ea = mesh->m_edges.get_new();
      EdgeId eb = mesh->m_edges.get_new();
      OPA_DISP0("Adding edge ", ea, eb, boundary_vertices[i], new_vertex);
      mesh->init_edge(ea, cloud[best].pos);
      mesh->init_edge(eb, mesh->get_pos(boundary_vertices[i]));
      idx.add(new_vertex, boundary_vertices[i], ea);
      idx.add(boundary_vertices[i], new_vertex, eb);
    }

    for (auto &x : boundary_edges) {
      VertexId va = mesh->start(x);
      VertexId vb = mesh->end(x);
      if (removed_boundary_vertices.count(va)) continue;
      if (removed_boundary_vertices.count(vb)) continue;
      if (blacklist_boundary_edge.count(MP(va, vb))) continue;
      idx.add(va, vb, x);
      idx.add(vb, va, mesh->opp(x));
    }

    std::set<EdgeId> drop_edge;
    std::set<VertexId> drop_vertex;
    for (auto face : visible_faces) {
      for (auto eid : Mesh::Face_EdgeWalker(*mesh, face)) {
        if (boundary_edges_set.count(eid)) continue;
        VertexId vid = mesh->start(eid);
        drop_edge.insert(eid);
        if (!boundary_vertices_set.count(vid))
          drop_vertex.insert(mesh->start(eid));
      }
    }
    for (auto &x : aligned_path){
      drop_edge.insert(x);
      drop_edge.insert(mesh->opp(x));
    }
    drop_vertex.insert(ALL(removed_boundary_vertices));

    if (aligned_faces.size() == 1) {
      FaceId curid = *aligned_faces.begin();
      std::vector<VertexId> nface_vertices;
      EdgeId cur = mesh->next(mesh->opp(aligned_path[0]));
      EdgeId last = mesh->opp(aligned_path.back());

      nface_vertices.push_back(new_vertex);
      while (cur != last) {
        int v1 = mesh->start(cur);
        int v2 = mesh->end(cur);
        int opp = mesh->opp(cur);
        idx.add(v1, v2, cur);
        idx.add(v2, v1, opp);
        nface_vertices.push_back(mesh->start(cur));
        cur = mesh->next(cur);
      }

      nface_vertices.push_back(mesh->start(last));
      //for (auto &xid : { nface_vertices[1], nface_vertices.back() }) {
      //  EdgeId ea = mesh->m_edges.get_new();
      //  EdgeId eb = mesh->m_edges.get_new();
      //  mesh->init_edge(ea, cloud[best].pos);
      //  mesh->init_edge(eb, mesh->get_pos(xid));
      //  idx.add(new_vertex, xid, ea);
      //  idx.add(xid, new_vertex, eb);
      //}

      OPA_DISP("Adding aigned ", nface_vertices, aligned_path,
               mesh->get_face_vertex_ids(curid));

      mesh->setup_face(curid, nface_vertices, idx);

    } else if (aligned_faces.size() == 2) {
      OPA_CHECK0(false);
    }

    struct NewFaceInfo {
      FaceId adj_face;
      FaceId new_id;
    };

    std::vector<NewFaceInfo> new_face_infos;
    std::vector<FaceId> new_faces;
    auto circular_boundary_vertices =
      utils::CreateCircularView(boundary_vertices);
    REP (i, nx) {
      VertexId va = circular_boundary_vertices[i];
      VertexId vb = circular_boundary_vertices[i + 1];
      if (blacklist_boundary_edge.count(MP(va, vb))) continue;

      FaceId new_face = mesh->m_faces.get_new();
      std::vector<VertexId> tb = { va, vb, new_vertex };

      mesh->setup_face(new_face, tb, idx);
      auto tmp = new FaceInfo;
      tmp->id = new_face;
      cur_faces[new_face].reset(tmp);
      OPA_DISP("new2 ", tb);

      new_face_infos.push_back(
        NewFaceInfo{ mesh->face(idx.get(vb, va)), new_face });
    }

    for (auto x : drop_edge) {
      OPA_DISP("drop edge ", x);
      mesh->do_remove_edge(x);
      mesh->m_attrs.remove(mesh->get_edge(x).vertex_attr);
      mesh->m_edges.remove(x);
    }

    for (auto x : drop_vertex) {
      OPA_DISP("drop vert >> ", x);
      mesh->orig_map().erase(x);
      mesh->m_vertices.remove(x);
    }

    for (auto face : visible_faces) {
      OPA_DISP("REMOVE FACE ", face);
      mesh->m_faces.remove(face);
    }

    // Preparing for next round
    for (auto &new_face_info : new_face_infos) {
      auto &f1 = cur_faces[new_face_info.adj_face];
      auto &f2 = cur_faces[face_pick->id];
      FaceId fid = new_face_info.new_id;
      auto face_info = cur_faces[fid];

      auto face_normal = mesh->get_normal(fid);
      auto pt = mesh->get_pos(mesh->get_face_vertex(fid));

      for (auto p : utils::MergeWalker<std::set<int> >(
             { &f1->bad_points, &f2->bad_points })) {
        if (OPA_FLOAT_GT(
              point_signed_dist_to_plane(cloud[p].pos - pt, face_normal), 0,
              eps)) {
          info[p].visible_faces.insert(fid);
          face_info->bad_points.insert(p);
        }
      }
      faces_rank.insert(face_info);
    }

    std::vector<int> to_remove;
    for (auto x : visible_faces) {
      for (auto p : cur_faces[x]->bad_points) {
        auto &pinfo = info[p].visible_faces;
        OPA_CHECK0(pinfo.count(x));
        pinfo.erase(x);
        if (pinfo.size() == 0) to_remove.pb(p);
      }
      OPA_CHECK0(faces_rank.count(cur_faces[x]));
      faces_rank.erase(cur_faces[x]);
      OPA_CHECK0(faces_rank.count(cur_faces[x]) == 0);
      OPA_CHECK0(cur_faces.count(x));
      cur_faces.erase(x);
    }

    for (auto x : to_remove) {
      rem_points.erase(x);
      info[x].visible_faces.clear();
    }
    CheckConvexPolyhedra(*mesh);
  }
}

void MeshBuilder::CheckConvexPolyhedra(const Mesh &mesh) {
  OPA_DISP("CUR STATE >> ", mesh.str());
  for (auto f : mesh.list_faces()) {
    auto normal = mesh.get_normal(f);
    auto pt = mesh.get_pos(mesh.get_face_vertex(f));
    for (auto x : mesh.list_vertices()) {
      auto pos = mesh.get_vertex_attr(x).pos;
      float v = point_signed_dist_to_plane(pos - pt, normal);
      OPA_DISP("at >> ", f, x, pos, pt, normal, v, eps);
      OPA_CHECK0(OPA_FLOAT_LE(v, 0, eps));
    }
  }
}

class CenteringHelper {
public:
  // no roll allowed
  // x is front, z up
  CenteringHelper(const Mat4 &proj, const PointVec &point_cloud, double ratio)
      : m_proj(proj), m_ratio(ratio) {
    m_iproj = glm::inverse(m_proj);
    m_fovy = compute_perspective_fovy(m_iproj) * m_ratio;
    m_fovx = compute_perspective_fovx(m_iproj) * m_ratio;
    Pos center = opa::math::game::gravity_center(point_cloud);
    m_center_rot = quat_look_at_safe(center, vec_z);
    Rot irot = glm::inverse(m_center_rot);

    for (auto &pt : point_cloud) {
      m_point_cloud.push_back(irot * pt);
    }
  }

  bool checkit(const Pos &pos, Rot *new_rot) {
    // only need to check for pitch
    double mx_high = -inf;
    double mx_low = inf;
    double my_high = -inf;
    double my_low = inf;
    for (auto &pt : m_point_cloud) {
      Pos v = pt - pos;
      double ang_x = atan2(v.y, v.x);
      mx_high = std::max(mx_high, ang_x);
      mx_low = std::min(mx_low, ang_x);
      double ang_y = atan2(v.z, v.x);
      my_high = std::max(my_high, ang_y);
      my_low = std::min(my_low, ang_y);
    }
    if (mx_high - mx_low > m_fovx) return false;
    if (my_high - my_low > m_fovy) return false;

    int nstep = 20;
    OR::GridSearchFunc tmp;
    tmp.func = { [&](const std::vector<double> &state) {
      Mat4 rot = glm::inverse(
        glm::mat4_cast(quat_from_euler(Pos(state[0], state[1], 0))));
      for (auto &pt : m_point_cloud) {
        if (!InFOV(m_proj, rot * Pos4(pt - pos, 1.), 0.1)) return 1;
      }
      return 0;
    } };

    auto result = OR::DoGridSearch(
      OR::GridSearchBounds::Build(
        { std::tuple<double, double, int>{ mx_low, mx_high, nstep },
          std::tuple<double, double, int>{ my_low, my_high, nstep } }),

      tmp);
    if (result.score == 1) return false;

    if (new_rot != nullptr) {
      OPA_DISP("got score ", result.score);
      *new_rot = quat_from_euler(Pos(result.state[0], result.state[1], 0));
    }
    return true;
  }

  void compute(Pos *new_pos, Rot *new_rot) {
    Pos cur_pos{ 0, 0, 0 };
    double min_pos = inf;
    for (auto &pt : m_point_cloud)
      min_pos = std::min<double>(min_pos, glm::dot(base_look_dir, pt));

    double step = 1.;
    cur_pos += base_look_dir * min_pos;
    while (true) {
      cur_pos -= base_look_dir * step;
      OPA_DISP0(step);
      if (checkit(cur_pos, nullptr)) break;
      step *= 2;
      OPA_CHECK0(step < inf);
    }
    cur_pos += base_look_dir * step;
    double L = 0.;
    double H = step;
    REP (ntry, 30) {
      double M = (L + H) / 2;
      if (checkit(cur_pos - base_look_dir * M, nullptr))
        H = M;
      else
        L = M;
    }
    *new_pos = cur_pos - base_look_dir * H;
    checkit(*new_pos, new_rot);
    *new_pos = m_center_rot * *new_pos;
    *new_rot = m_center_rot * *new_rot;
  }

  const Pos base_look_dir = { 1, 0, 0 };

  PointVec m_point_cloud;
  Rot m_center_rot;
  const Mat4 &m_proj;
  Mat4 m_iproj;
  double m_fovy;
  double m_fovx;
  double m_ratio;
};

void compute_new_center(const Mat4 &proj, const PointVec &point_cloud,
                        Pos *new_pos, Rot *new_rot, double ratio) {
  CenteringHelper helper(proj, point_cloud, ratio);
  helper.compute(new_pos, new_rot);
  if (0) {
    for (auto &pt : point_cloud) {
      Pos cur = pt - *new_pos;
      Pos tsf =
        ApplyMat(proj, ApplyMat(glm::inverse(glm::mat4_cast(*new_rot)), cur));
      OPA_DISP0(tsf);
    }
  }
}

void compute_new_center(const Mat4 &proj,
                        const std::vector<MeshData> &mesh_data, Pos *new_pos,
                        Rot *new_rot, double ratio) {
  compute_new_center(proj, to_point_cloud(mesh_data), new_pos, new_rot, ratio);
}

std::vector<SPTR(Mesh)> mesh_cut_plan(const Mesh &mesh,
                                      const HyperPlaneSpec &plane) {
  // must be triangle faces, otherwise use
  // SPTR(Mesh) proc_mesh = mesh.to_triangles()->convert<Mesh>();
  FaceCollection tc{ true };
  std::vector<std::pair<Pos, FaceId> > all_new_points;
  mesh.compute_remap();
  OPA_CHECK0(mesh.compute_graph()->is_connected());

  for (const auto &face : mesh.faces()) {
    std::vector<Pos> face_pos = mesh.get_face_vertices(face);
    Triangle tr = vec_to_tr(face_pos);
    std::vector<Pos> in_id, out_id, on_id;
    std::for_each(ALL(face_pos), [&](const Pos &pos) {
      if (plane.is_on(pos))
        on_id.push_back(pos);
      else if (plane.is_in(pos))
        in_id.push_back(pos);
      else
        out_id.push_back(pos);
    });
    OPA_DISP0(face_pos, on_id, in_id, out_id);

    for (auto &pt : on_id) all_new_points.emplace_back(pt, face);

    if (out_id.empty()) {
      tc.push(vec_to_tr(face_pos));
      OPA_TRACES(tc.faces().back());
    } else if (in_id.empty()) {
      continue;
    } else {
      std::vector<Pos> new_points;
      REP (i, in_id.size()) {
        REP (j, out_id.size()) {
          new_points.push_back(
            line_plane_intersection({ in_id[i], out_id[j] }, plane.plane));
        }
      }
      if (in_id.size() == 2) {
        tc.push(realign_tr(Triangle{ in_id[0], in_id[1], new_points[0] }, tr));
        OPA_TRACES(tc.faces().back());
        tc.push(
          realign_tr(Triangle{ in_id[1], new_points[0], new_points[1] }, tr));
        OPA_TRACES(tc.faces().back());
      } else {
        std::vector<Pos> new_tr = on_id;
        new_tr.insert(new_tr.end(), ALL(in_id));
        new_tr.insert(new_tr.end(), ALL(new_points));
        auto new_tr2 = vec_to_tr(new_tr);
        tc.push(realign_tr(new_tr2, tr));
        OPA_TRACES(tc.faces().back(), tr_normal(tr));
      }
      for (auto &new_point : new_points)
        all_new_points.emplace_back(new_point, face);
    }
  }

  if (all_new_points.size() > 0) {
    // All new points stores all the points on the new face.
    // We have to order them to triangularize.
    // First we do some points matching.
    // Then we know on which face each point was created.
    // A graph is created with points as vertices and face as edges.
    // If all goes well, should have a biconnected graph.
    // Each cycle corresponds to a face.
    std::vector<Pos> to_match;
    for (auto &pt : all_new_points) to_match.push_back(pt.first);
    PointMatcher matcher;
    OPA_CHECK0(matcher.load_batch(to_match));
    std::vector<std::pair<int, FaceId> > remapped_points;
    for (auto &pt : all_new_points) {
      remapped_points.emplace_back(matcher.query(pt.first), pt.second);
    }

    std::unordered_map<int, std::vector<FaceId> > vertex_to_face;
    std::unordered_map<FaceId, std::vector<int> > face_to_vertex;
    for (auto &pt : remapped_points) {
      vertex_to_face[pt.first].push_back(pt.second);
      face_to_vertex[pt.second].push_back(pt.first);
    }

    algo::FastGraph graph(vertex_to_face.size(), algo::Mode::REMAP);
    for (const auto &face : face_to_vertex) {
      OPA_CHECK(face.second.size() <= 2, face);
      if (face.second.size() == 2) {
        graph.add_bidirectional(face.second[0], face.second[1]);
      }
    }

    auto cc = graph.list_connected_components();
    OPA_TRACES(remapped_points, cc);
    OPA_TRACES(to_match);
    OPA_TRACES(cc.size());
    for (const auto &c : cc) {
      OPA_TRACES(c);
      if (c.size() == 1) continue;
      OPA_CHECK0(c.size() >= 3);
      for (auto &v : c) {
        for (auto &x : graph.get_edges(v)) OPA_DISP0(x.from, x.to);
        puts("");
        OPA_CHECK(graph.deg(v) == 2, graph.deg(v), to_match[v]);
      }
      PointVec points;
      for (int e : c) points.push_back(to_match[e]);
      OPA_TRACES(points);
      PointVec realigned = realign_face(points, -plane.plane.dir);
      OPA_CHECK(glm::dot(face_normal(realigned), -plane.plane.dir) > 0,
                plane.plane.dir, face_normal(realigned));
      tc.push(realigned);
    }
    for (auto &tr : tc.faces()) {
      OPA_TRACES(tr, face_normal(tr));
    }
    OPA_DISP("PPLANE DIR ", plane.plane.dir);
  }

  {
    proto::MeshList meshes_proto;
    Exporter(ExportType::STL)
      .to_mesh(tc, meshes_proto.add_mesh())
      ->set_name("DEBUG");
    Exporter(ExportType::STL)
      .to_mesh(mesh, meshes_proto.add_mesh())
      ->set_name("ORIG");

    std::ofstream tmpfile("/tmp/debug.out", std::ios::binary);
    meshes_proto.SerializeToOstream(&tmpfile);
  }

  return tc.convert<Mesh>()->split_connected();
}

std::vector<SPTR(Mesh)> mesh_cut_box(const Mesh &mesh, const BoxSpec &box) {

  std::vector<SPTR(const Mesh)> cur = { make_dummy_sptr(&mesh) };
  std::vector<SPTR(Mesh)> res;
  std::vector<HyperPlaneSpec> hps = get_box_hyperplans(box);

  REP (i, hps.size()) {
    if (i != 0) {
      for (auto &x : res) cur.push_back(std::const_pointer_cast<Mesh>(x));
      res.clear();
    }
    for (auto &x : cur) {
      std::vector<SPTR(Mesh)> tb = mesh_cut_plan(*x, hps[i]);
      res.insert(res.end(), ALL(tb));
    }
    cur.clear();
  }
  return res;
}

OPA_NAMESPACE_DECL3_END

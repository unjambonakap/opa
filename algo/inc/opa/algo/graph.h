#pragma once

#include <opa/utils/map_util.h>
#include <opa/algo/base.h>
#include <opa_common.h>
#include <type_traits>

OPA_NAMESPACE_DECL2(opa, algo)

#define OPA_FOREACH_EDGE_OBJ(e, v, obj)                                                            \
  for (int e = (obj).last[(v)]; e != -1; e = (obj).edges[e].prev)

struct GraphDistData {
  std::vector<std::vector<int> > dists;
  std::vector<std::vector<int> > next;
  std::vector<int> compute_path(int a, int b) const;
};

struct EdgeData {
  EdgeData() {}
  EdgeData(int prev, int next, int from, int to, int id, bool forward)
      : prev(prev), next(next), from(from), to(to), id(id), forward(forward) {}
  int prev, next;
  int from;
  int to;
  int id;
  int data;
  bool forward = true;
  bool removed = false;
};

struct DijkstraParams {
  std::unordered_map<int, int> edge_cost;
  int s;
  int t;
  int max_dist = -1;
};

struct DijkstraResult {
  std::vector<int> path;
  int cost;
};

class Graph_EdgeIterator;
// Keep in sync with mode union.struct bitset
enum Mode {
  NONE = 0,
  REMAP = 1,
  MULTIGRAPH = 2,
  DIGRAPH = 4,
  DYNAMIC_SIZE = 8,
};
const int INVALID = -1;

enum AttrType {
  ATTR_INT,
  ATTR_STR,
};

struct AttributeContainer {
  std::unordered_map<int, int> int_map;
  std::unordered_map<int, std::string> str_map;
  AttrType attr;

  std::string str() const {

    std::stringstream ss;
    ss << (attr == AttrType::ATTR_INT ? "int" : "str") << ","
       << std::max(int_map.size(), str_map.size()) << std::endl;
    for (auto &kv : int_map) ss << kv << std::endl;
    for (auto &kv : str_map) ss << kv << std::endl;

    return ss.str();
  }
  template <typename T>
  void extract_map(const std::unordered_map<int, T> &src, std::unordered_map<int, T> &dest,
                   const std::vector<int> &vals) const {
    for (auto &v : vals) {
      const T *x = glib::gtl::FindOrNull(src, v);
      if (x != nullptr) dest[v] = *x;
    }
  }

  SPTR(AttributeContainer) extract(const std::vector<int> &vals) const {
    SPTR(AttributeContainer) res = std::make_shared<AttributeContainer>();

    res->attr = attr;
    if (attr == AttrType::ATTR_INT)
      extract_map(int_map, res->int_map, vals);
    else if (attr == AttrType::ATTR_STR)
      extract_map(str_map, res->str_map, vals);
    return res;
  }

  template <typename T> T get(int node);
  template <typename T> void set(int node, const T &val);
  template <typename T> static AttributeContainer *New();
};

template <> inline AttributeContainer *AttributeContainer::New<int>() {
  auto res = new AttributeContainer();
  res->attr = AttrType::ATTR_INT;
  return res;
}

template <> inline AttributeContainer *AttributeContainer::New<std::string>() {
  auto res = new AttributeContainer();
  res->attr = AttrType::ATTR_STR;
  return res;
}

template <> inline std::string AttributeContainer::get(int node) {
  OPA_CHECK0(attr == AttrType::ATTR_STR);
  return glib::gtl::FindWithDefault(str_map, node, "");
}

template <> inline int AttributeContainer::get(int node) {
  OPA_CHECK0(attr == AttrType::ATTR_INT);
  return glib::gtl::FindWithDefault(int_map, node, 0);
}

template <> inline void AttributeContainer::set(int node, const std::string &val) {
  OPA_CHECK0(attr == AttrType::ATTR_STR);
  str_map[node] = val;
}

template <> inline void AttributeContainer::set(int node, const int &val) {
  OPA_CHECK0(attr == AttrType::ATTR_INT);
  int_map[node] = val;
}

struct AttributeMap {

  template <typename T> void set_attr(const std::string &attr, int node, const T &val) {
    attribute_map.insert(attr, make_sptr(AttributeContainer::New<T>())).first->set(node, val);
  }
  template <typename T> T get_attr(const std::string &attr, int node) {
    return attribute_map[attr]->get<T>(node);
  }

  const AttributeContainer &get_attribute_container(const std::string &name) const {
    auto v = glib::gtl::FindOrNull(attribute_map, name);
    OPA_CHECK0(v != nullptr);
    return **v;
  }

  void extract(const std::vector<int> &nodes, AttributeMap *dest) const {
    for (auto &kv : attribute_map) dest->attribute_map[kv.first] = kv.second->extract(nodes);
  }

  std::string str() const {
    std::stringstream ss;
    for (auto &kv : attribute_map) {
      std::string attr = kv.first;
      ss << RAW_OPA_DISP_VARS(attr);
      ss << kv.second->str();
    }
    return ss.str();
  }

  std::unordered_map<std::string, SPTR(AttributeContainer)> attribute_map;
};

class FastGraph {

public:
  union {
    struct {
      u8 remap : 1;
      u8 multigraph : 1;
      u8 digraph : 1;
      u8 dynamic_size : 1;
    };
    u32 raw;
  } mode;

  struct VertexData {
    int deg = 0;
    int deg_in = 0;
    int deg_out = 0;
  };

  struct DijkstraCtx {
    std::vector<double> edge_costs;
    std::vector<int> prev;
    std::vector<double> dists;
    std::vector<int> sources;
  };

  struct BfsCtx {
    std::unordered_set<int> vis;
    std::vector<int> sources;
    std::vector<int> waiting;
  };

  int get_edge_id(int eid) const { return eid; }

  FastGraph(int n = 0, int mode = Mode::MULTIGRAPH) { reset(n, mode); }

  void reset(int n, int mode = Mode::MULTIGRAPH) {
    this->mode.raw = mode;
    this->n = n;
    this->edge_count = 0;
    last = std::vector<int>(n, -1);
    rlast = std::vector<int>(n, -1);
    vertices = std::vector<VertexData>(n);
    m_seen_edges.clear();
    if (this->mode.dynamic_size)
      m_rmp.reset(-1);
    else
      m_rmp.reset(n);
  }

  bool hase(int a, int b, bool norm = true) const;

  int add_bidirectional(int a, int b, bool norm = true) {
    if (a < b) std::swap(a, b);
    if (!mode.multigraph) {
      if (!m_seen_edges.insert(std::make_pair(a, b)).second) return INVALID;
    }
    a = normv(a, norm);
    b = normv(b, norm);

    int id = adde_internal(a, b);
    adde_internal(b, a);
    OPA_CHECK0(edges.size() % 2 == 0);
    return get_edge_id(id);
  }

  int adde(int a, int b, bool norm = true) {
    a = normv(a, norm);
    b = normv(b, norm);
    int id = adde_internal(a, b, true);
    adde_internal(b, a, false);
    OPA_CHECK0(edges.size() % 2 == 0);
    return get_edge_id(id);
  }

  int normv(int a, bool norm = true) const {
    if (!norm || !mode.remap) return a;
    return m_rmp.get_or_die(a);
  }

  int normv(int a, bool norm = true) {
    int res;
    if (!norm || !mode.remap)
      res = a;
    else
      res = m_rmp.get(a);
    if (mode.dynamic_size) {
      while (n <= res) {
        vertices.emplace_back();
        last.push_back(-1);
        rlast.push_back(-1);
        ++n;
      }
    }
    return res;
  }

  int inormv(int a, bool norm = true) const {
    if (!norm || !mode.remap) return a;
    return m_rmp.rget(a);
  }

  int deg(int a, bool norm = true) const {
    auto &d = vertices[normv(a, norm)];
    return mode.digraph ? d.deg_in + d.deg_out : d.deg_in;
  }

  int deg_in(int a, bool norm = true) const {
    auto &d = vertices[normv(a, norm)];
    return d.deg_in;
  }

  int deg_out(int a, bool norm = true) const {
    auto &d = vertices[normv(a, norm)];
    return d.deg_out;
  }

  int get_edge_count(int a, int b, bool norm = true) const;

  void dfs_cc(int a, std::unordered_set<int> &seen, std::vector<int> &lst) const {
    if (seen.count(a)) return;
    seen.insert(a);
    lst.push_back(inormv(a));

    for (int e = last[a]; e != -1; e = edges[e].prev) {
      dfs_cc(edges[e].to, seen, lst);
    }
  }

  Graph_EdgeIterator edge_iter(int vertex, bool norm = true) const;
  void dijkstra(DijkstraCtx &ctx);
  DijkstraResult dijkstra_adv(DijkstraParams &params);
  void bfs_one(BfsCtx &ctx);

  bool is_connected() const;
  std::string str() const;

  GraphDistData &all_shortest_paths(GraphDistData *data) const;

  std::vector<std::vector<int> > list_connected_components(bool want_single_vertices = true) const;
  std::vector<std::vector<int> > get_path_cover();
  std::vector<std::vector<int> > get_path_cover_dumb();

  std::vector<SPTR(FastGraph)> split_cc() const;
  OPA_ACCESSOR_R(opa::utils::Remapper<int>, m_rmp, rmp);

  // Must be connected, chinese postman
  std::vector<int> get_cover_walk();
  std::vector<int> get_cover_walk_dumb(bool need_dup = true);
  std::vector<int> get_eulerian_cycle(bool norm = false);
  std::vector<int> get_eulerian_path(int start, bool norm = true);
  int n_edges() const { return mode.digraph ? edge_count : edge_count / 2; }

  int n;
  int edge_count;
  int node_count() const { return mode.remap ? m_rmp.size() : n; }
  void add_node(int v) { normv(v); }
  AttributeMap attr_map;

  std::vector<VertexData> vertices;
  std::vector<int> last, rlast;
  std::vector<EdgeData> get_edges(int a, bool norm = true) const;
  std::vector<EdgeData> get_redges(int a, bool norm = true) const;
  std::vector<EdgeData> get_all_edges(int a, bool norm = true) const;
  std::vector<EdgeData> get_edges2(int a, int b, bool norm = true) const;
  std::vector<std::pair<pii, pii> > find_set_neighbours(const std::vector<std::vector<int> > &lst,
                                                        bool norm = true) const;
  std::vector<int> get_successors(int a, bool norm = true) const;
  std::vector<int> get_neighbours(int a, bool norm = true) const;
  std::vector<int> get_predecessors(int a, bool norm = true) const;
  std::unordered_map<pii, int> edge_map;
  std::unordered_map<pii, int> redge_map;
  std::vector<EdgeData> edges;

  SPTR(FastGraph) make_bidirectional() const;
  SPTR(FastGraph) clone() const;
  SPTR(FastGraph)
  subgraph(const std::vector<int> &nodes, bool norm = true, bool neighborhood = false) const;

  std::vector<EdgeData> list_edges() const {
    std::vector<EdgeData> res;
    REP (i, node_count()) {
      for (auto ed : get_edges(i, false)) {
        ed.from = inormv(ed.from);
        ed.to = inormv(ed.to);
        res.pb(ed);
      }
    }
    return res;
  }

  EdgeData &get_edge_data(int eid) { return edges[eid]; }
  const EdgeData &get_edge_data(int eid) const { return edges[eid]; }

  void remove_bidirectional_edge(int edge_id) { remove_edge(edge_id); }

  void remove_edge(int a, int b, bool norm = true) {
    a = normv(a, norm);
    b = normv(b, norm);
    int ix = edge_map[MP(a, b)];
    OPA_CHECK0(edges.size() % 2 == 0);
    remove_edge(ix);
  }
  void remove_edge(int edge_id) {
    remove_edge_internal(edge_id);
    remove_edge_internal(edge_id ^ 1);
    OPA_CHECK0(edges.size() % 2 == 0);
  }
  void remove_edge_internal(int edge_id);
  void disconnect_node(int node, bool norm = true);

  std::vector<int> list_nodes(bool norm = true) const {
    std::vector<int> res;
    int nc = this->node_count();
    REP (i, nc)
      res.pb(inormv(i, norm));
    return res;
  }

private:
  int adde_internal(int a, int b, bool forward = true) {
    OPA_CHECK(a < n && b < n, a, b, n);
    int cur_id = edges.size();

    int real_id;
    if (!mode.multigraph) {
      if (forward)
        real_id = glib::gtl::LookupOrInsert(&edge_map, MP(a, b), cur_id);
      else
        real_id = glib::gtl::LookupOrInsert(&redge_map, MP(a, b), cur_id);

      if (real_id != cur_id) return real_id;
    }

    if (forward) {
      ++edge_count;
      ++vertices[a].deg_out;
      ++vertices[b].deg_in;
    }

    auto &lastp = forward ? last : rlast;

    edges.push_back(EdgeData(lastp[a], -1, a, b, cur_id, forward));

    if (lastp[a] != -1) edges[lastp[a]].next = cur_id;
    lastp[a] = cur_id;
    return cur_id;
  }

  struct EulerianRecState {
    std::vector<int> order;
  };

  void eulerian_rec(int pos, EulerianRecState &state);
  struct EulerianRecStateRobust {
    std::vector<int> order;
    std::vector<std::vector<int> > paths;
    std::vector<pii> stack_and_order;
  };

  void eulerian_rec_robust(int pos, EulerianRecStateRobust &state);
  void dfs_edge(int pos, std::vector<int> &vis);

  opa::utils::Remapper<int> m_rmp;
  std::set<std::pair<int, int> > m_seen_edges;
};

class Graph_EdgeIterator {
public:
  const FastGraph *m_graph;
  int m_eid;
  Graph_EdgeIterator() {}

  Graph_EdgeIterator(const FastGraph *graph, int vertex) {
    m_graph = graph;
    m_eid = m_graph->last[vertex];
  }
  Graph_EdgeIterator begin() { return *this; }
  Graph_EdgeIterator end() {
    Graph_EdgeIterator res;
    res.m_eid = -1;
    res.m_graph = m_graph;
    return res;
  }

  const EdgeData &operator*() const {
    OPA_CHECK0(m_eid != -1);
    return m_graph->get_edge_data(m_eid);
  }
  const EdgeData *operator->() const { return &**this; }

  Graph_EdgeIterator &operator++() {
    m_eid = (*this)->prev;
    return *this;
  }
  bool operator!=(const Graph_EdgeIterator &it) const { return m_eid != it.m_eid; }
};

inline Graph_EdgeIterator FastGraph::edge_iter(int vertex, bool norm) const {
  vertex = normv(vertex, norm);
  return Graph_EdgeIterator(this, vertex);
}

OPA_NAMESPACE_DECL2_END

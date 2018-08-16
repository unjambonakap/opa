#pragma once

#include <opa/algo/base.h>
#include <opa_common.h>
#include <type_traits>

OPA_NAMESPACE_DECL2(opa, algo)

#define OPA_FOREACH_EDGE_OBJ(e, v, obj)                                        \
  for (int e = (obj).last[(v)]; e != -1; e = (obj).edges[e].prev)

struct GraphDistData {
  std::vector<std::vector<int> > dists;
  std::vector<std::vector<int> > next;
  std::vector<int> compute_path(int a, int b) const;
};

struct EdgeData {
  int prev, next;
  int from;
  int to;
  int id;
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

  int get_edge_id(int eid) const { return eid >> (mode.digraph == 0); }

  FastGraph(int n = 0, int mode = Mode::MULTIGRAPH) { reset(n, mode); }

  void reset(int n, int mode = Mode::MULTIGRAPH) {
    this->mode.raw = mode;
    this->n = n;
    this->edge_count = 0;
    last = std::vector<int>(n, -1);
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
    int id = adde(a, b, norm);
    adde(b, a, norm);
    return get_edge_id(id);
  }

  int adde(int a, int b, bool norm = true) {
    a = normv(a, norm);
    b = normv(b, norm);
    OPA_CHECK(a < n && b < n, a, b, n);
    ++edge_count;
    ++vertices[a].deg;
    ++vertices[b].deg;
    int cur_id = edges.size();
    edges.push_back(EdgeData{ last[a], -1, a, b, cur_id });
    if (last[a] != -1) edges[last[a]].next = cur_id;
    last[a] = cur_id;
    return cur_id;
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
        ++n;
      }
    }
    return res;
  }

  int inormv(int a) const {
    if (!mode.remap) return a;
    return m_rmp.rget(a);
  }

  int deg(int a, bool norm = true) const {
    int x = vertices[normv(a, norm)].deg;
    return mode.digraph ? x : x / 2;
  }

  int get_edge_count(int a, int b, bool norm = true) const;

  void dfs_cc(int a, std::unordered_set<int> &seen,
              std::vector<int> &lst) const {
    if (seen.count(a)) return;
    seen.insert(a);
    lst.push_back(inormv(a));

    for (int e = last[a]; e != -1; e = edges[e].prev) {
      dfs_cc(edges[e].to, seen, lst);
    }
  }

  Graph_EdgeIterator edge_iter(int vertex, bool norm = true) const;
  void dijkstra(DijkstraCtx &ctx);
  void bfs_one(BfsCtx &ctx);

  bool is_connected() const;
  std::string str() const;

  GraphDistData &all_shortest_paths(GraphDistData *data) const;

  std::vector<std::vector<int> >
  list_connected_components(bool want_single_vertices = true) const;
  std::vector<std::vector<int> > get_path_cover();
  std::vector<std::vector<int> > get_path_cover_dumb();

  std::vector<SPTR(FastGraph)> split_cc() const;
  OPA_ACCESSOR_R(opa::utils::Remapper<int>, m_rmp, rmp);

  // Must be connected, chinese postman
  std::vector<int> get_cover_walk();
  std::vector<int> get_cover_walk_dumb();
  std::vector<int> get_eulerian_cycle(bool norm=false);
  std::vector<int> get_eulerian_path(int start, bool norm = true);
  int n_edges() const { return mode.digraph? edge_count : edge_count / 2; }

  int n;
  std::vector<VertexData> vertices;
  std::vector<EdgeData> edges;
  std::vector<int> last;
  std::vector<EdgeData> get_edges(int a, bool norm = true) const;
  std::vector<EdgeData> get_edges2(int a, int b, bool norm = true) const;
  std::vector<std::pair<pii, pii> >
  find_set_neighbours(const std::vector<std::vector<int> > &lst,
                      bool norm = true) const;
  int edge_count;

  SPTR(FastGraph) make_bidirectional() const;
  SPTR(FastGraph) clone() const;

private:
  void remove_bidirectional_edge(int edge_id) {
    remove_edge(edge_id);
    remove_edge(edge_id ^ 1);
  }

  void remove_edge(int edge_id) {
    auto &edge_data = edges[edge_id];
    --vertices[edge_data.from].deg;
    --vertices[edge_data.to].deg;
    --edge_count;

    if (edge_data.next == -1) {
      last[edge_data.from] = edge_data.prev;
    } else {
      edges[edge_data.next].prev = edge_data.prev;
    }
    if (edge_data.prev != -1) edges[edge_data.prev].next = edge_data.next;
  }

  struct EulerianRecState {
    std::vector<int> order;
  };

  void eulerian_rec(int pos, EulerianRecState &state);
  struct EulerianRecStateRobust {
    std::vector<int> order;
    std::vector<std::vector<int>> paths;
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
    return m_graph->edges[m_eid];
  }
  const EdgeData *operator->() const { return &**this; }

  Graph_EdgeIterator &operator++() {
    m_eid = (*this)->prev;
    return *this;
  }
  bool operator!=(const Graph_EdgeIterator &it) const {
    return m_eid != it.m_eid;
  }
};

inline Graph_EdgeIterator FastGraph::edge_iter(int vertex, bool norm) const {
  vertex = normv(vertex, norm);
  return Graph_EdgeIterator(this, vertex);
}

OPA_NAMESPACE_DECL2_END

#include <opa/algo/graph.h>

#include <glib/strings/substitute.h>
#include <lemon/concepts/graph.h>
#include <lemon/list_graph.h>
#include <lemon/matching.h>
#include <lemon/smart_graph.h>
#include <opa/utils/DataStruct.h>
#include <opa/utils/string.h>

OPA_NAMESPACE_DECL2(opa, algo)
#define OPA_FOREACH_EDGE(e, v)                                                 \
  for (int e = last[(v)]; e != -1; e = edges[e].prev)

std::vector<int> GraphDistData::compute_path(int a, int b) const {
  std::vector<int> res;
  while (a != b) {
    res.push_back(a);
    a = this->next[a][b];
  }
  res.push_back(b);
  return res;
}

std::vector<std::vector<int> >
FastGraph::list_connected_components(bool want_single_vertices) const {
  std::vector<std::vector<int> > res;
  std::unordered_set<int> seen;

  REP (i, n) {
    if (seen.count(i)) continue;
    res.emplace_back();
    dfs_cc(i, seen, res.back());
    if (!want_single_vertices && res.back().size() == 1) res.pop_back();
  }
  return res;
}

SPTR(FastGraph) FastGraph::make_bidirectional() const {
  SPTR(FastGraph) res(new FastGraph(n, Mode::REMAP));
  REP (i, n) {
    OPA_FOREACH_EDGE(e, i) {
      res->add_bidirectional(inormv(i), inormv(edges[e].to));
    }
  }
  return res;
}
SPTR(FastGraph) FastGraph::clone() const {
  return std::make_shared<FastGraph>(*this);
}

std::vector<SPTR(FastGraph)> FastGraph::split_cc() const {
  auto cc_list = list_connected_components();
  std::vector<SPTR(FastGraph)> res(cc_list.size());
  REP (i, cc_list.size()) {
    res[i] = std::make_shared<FastGraph>(cc_list[i].size(), true);
    for (auto &v : cc_list[i]) {
      int v2 = normv(v);
      OPA_FOREACH_EDGE(e, v2) {
        if (e & 1) continue;
        res[i]->add_bidirectional(v, inormv(edges[e].to));
      }
    }
  }
  return res;
}
std::vector<int> FastGraph::get_eulerian_cycle(bool norm) {
  REP (i, n)
    if (deg(i, false) > 0) return get_eulerian_path(norm ? inormv(i) : i, norm);
  OPA_CHECK0(false);
}

void FastGraph::eulerian_rec(int pos, EulerianRecState &state) {
  std::vector<int> stack;
  stack.push_back(pos);
  while (!stack.empty()) {
    pos = stack.back();

    int eid = last[pos];
    if (eid == -1) {
      state.order.push_back(pos);
      stack.pop_back();
    } else {
      auto &e = edges[eid];
      remove_bidirectional_edge(eid);
      stack.push_back(e.to);
    }
  }
}

void FastGraph::dfs_edge(int pos, std::vector<int> &vis) {
  vis.push_back(pos);
  while (true) {
    pos = vis.back();

    int eid = last[pos];
    if (eid == -1) {
      return;
    } else {
      auto &e = edges[eid];
      remove_bidirectional_edge(eid);
      vis.push_back(e.to);
    }
  }
}

void FastGraph::eulerian_rec_robust(int pos, EulerianRecStateRobust &state) {
  std::vector<int> stack;
  stack.push_back(pos);

  while (!stack.empty()) {
    pos = stack.back();
    while (state.stack_and_order.size() >= 2) {
      if (state.stack_and_order[state.stack_and_order.size() - 2].ST !=
          state.stack_and_order.back().ST) {
        break;
      }

      auto a = state.stack_and_order.back();
      state.stack_and_order.pop_back();
      auto b = state.stack_and_order.back();
      state.stack_and_order.pop_back();

      std::vector<int> npath;
      FOR (i, a.ND, state.order.size())
        npath.push_back(state.order[i]);
      npath.push_back(inormv(pos));
      OPA_CHECK(state.order.size() > a.ND, a, state.order);
      state.order.resize(a.ND);
      while (state.order.size() > b.ND)
        npath.push_back(state.order.back()), state.order.pop_back();

      state.paths.emplace_back(std::move(npath));
    }

    int eid = last[pos];
    if (eid == -1) {

      stack.pop_back();
      if (state.stack_and_order.size() == 0 ||
          state.stack_and_order.back().first <= stack.size()) {
        state.stack_and_order.emplace_back(stack.size(), state.order.size());
      } else {
        state.stack_and_order.back().ST = stack.size();
      }
      state.order.push_back(inormv(pos));

    } else {
      auto &e = edges[eid];
      remove_bidirectional_edge(eid);
      stack.push_back(e.to);
    }
  }
  if (state.order.size() > 1) {
    state.paths.emplace_back(std::move(state.order));
  }
}

std::vector<int> FastGraph::get_eulerian_path(int start, bool norm) {
  start = normv(start, norm);
  FastGraph::EulerianRecState state;
  eulerian_rec(start, state);
  std::vector<int> res;
  for (auto &v : state.order) {
    res.push_back(norm ? inormv(v) : v);
  }
  std::reverse(ALL(res));
  return res;
}

std::vector<std::vector<int> > FastGraph::get_path_cover_dumb() {
  std::vector<std::vector<int> > paths;
  REP (i, n)
    if (deg(i, false) & 1) {
      paths.emplace_back();
      dfs_edge(i, paths.back());
      if (paths.back().size() <= 1) paths.pop_back();
    }

  std::vector<std::vector<int> > cycles;
  for (auto &path : paths) {
    std::vector<int> npath;
    REPV (i, path.size()) {
      npath.push_back(i);
      int v = normv(path[i]);
      if (!deg(v, false)) continue;
      std::vector<int> cycle = get_eulerian_path(v, false);
      if (cycle.size() > 1) npath.insert(npath.end(), ALL(cycle));
    }
  }

  REP (i, n)
    if (deg(i, false)) {
      paths.emplace_back(std::move(get_eulerian_path(i, false)));
    }
  return paths;
}

std::vector<std::vector<int> > FastGraph::get_path_cover() {
  std::vector<std::vector<int> > res;
  REP (i, n)
    if (deg(i, false) & 1) {
      EulerianRecStateRobust tmp;
      eulerian_rec_robust(i, tmp);
      for (auto &p : tmp.paths) res.push_back(p);
    }

  if (n_edges() > 0) {
    REP (i, n)
      if (deg(i, false) > 0) {
        EulerianRecStateRobust tmp;
        eulerian_rec_robust(i, tmp);
        for (auto &p : tmp.paths) res.push_back(p);
      }
  }
  OPA_CHECK0(n_edges() == 0);
  return res;
}

std::vector<int> FastGraph::get_cover_walk_dumb() {
  int ne = edges.size();
  // duplicate all edges, all vertices are now odd.
  REP (i, ne)
    adde(edges[i].from, edges[i].to, false);
  return get_eulerian_cycle();
}

std::vector<int> FastGraph::get_cover_walk() {
  utils::Remapper<int> odd_remap;
  REP (i, n)
    if (deg(i, false) & 1) odd_remap.get(i);

  if (odd_remap.size() == 0) {
    return get_eulerian_cycle();
  }
  // OPA_TRACES(odd_remap.mp());

  GraphDistData dist_data;
  all_shortest_paths(&dist_data);

  lemon::SmartGraph graph;
  REP (i, odd_remap.size())
    graph.addNode();

  lemon::SmartGraph::EdgeMap<int> edge_map(graph);
  REP (i, odd_remap.size())
    REP (j, i) {
      int oi = odd_remap.rget(i);
      int oj = odd_remap.rget(j);
      int c = dist_data.dists[oi][oj];
      OPA_CHECK(c < INT_MAX / 2, oi, oj, c);
      edge_map.set(graph.addEdge(graph.nodeFromId(i), graph.nodeFromId(j)), c);
    }

  lemon::MaxWeightedPerfectMatching<lemon::SmartGraph> max_weight_matching(
    graph, edge_map);
  max_weight_matching.run();
  std::vector<std::tuple<int, int, int> > data; // cost, i, j
  REP (i, odd_remap.size()) {
    auto node = graph.nodeFromId(i);
    auto edge = max_weight_matching.matching(node);
    int target = graph.id(graph.target(edge));
    if (target < i) continue;
    int ri = odd_remap.rget(i);
    int rj = odd_remap.rget(target);
    data.emplace_back(dist_data.dists[ri][rj], ri, rj);
  }

  int maxe = std::max_element(ALL(data)) - data.begin();
  FastGraph work_graph = *this;
  work_graph.mode.multigraph = true;
  // Adding all except most costly edge
  REP (i, data.size()) {
    if (i == maxe) continue;
    int a = std::get<1>(data[i]);
    int b = std::get<2>(data[i]);
    std::vector<int> res = dist_data.compute_path(a, b);
    // OPA_TRACES(a, b, res);
    REP (j, res.size() - 1) {
      OPA_CHECK(this->hase(res[j], res[j + 1], false), res[j], res[j + 1],
                dist_data.dists[res[j]][res[j + 1]]);
      work_graph.add_bidirectional(res[j], res[j + 1], false);
    }
  }

  REP (i, work_graph.n)
    OPA_DISP0(work_graph.deg(i, false) & 1);

  std::vector<int> res =
    work_graph.get_eulerian_path(std::get<1>(data[maxe]), false);
  OPA_DISP0(res.size(), this->n_edges());
  OPA_DISP0(res);
  std::set<std::pair<int, int> > lst;
  REP (i, res.size() - 1) {
    lst.insert(MP(normv(res[i]), normv(res[i + 1])));
    OPA_CHECK(this->hase(res[i], res[i + 1]), res[i], res[i + 1]);
  }

  for (auto &e : edges) {
    OPA_CHECK0(lst.count(MP(e.to, e.from)) || lst.count(MP(e.from, e.to)));
  }
  return res;
}
bool FastGraph::hase(int a, int b, bool norm) const {
  a = normv(a, norm);
  b = normv(b, norm);
  OPA_FOREACH_EDGE(e, a) {
    if (edges[e].to == b) return true;
  }
  return false;
}

bool FastGraph::is_connected() const {
  algo::UnionJoin uj(n);
  for (auto &edge : edges) uj.merge(edge.from, edge.to);
  return uj.count() == 1;
}

GraphDistData &FastGraph::all_shortest_paths(GraphDistData *data) const {
  data->dists.resize(n);
  data->next.resize(n);
  int inf = INT_MAX / 2;
  OPA_DISP0(n);
  REP (i, n) {
    data->dists[i] = std::vector<int>(n, inf);
    data->next[i] = std::vector<int>(n, -1);
    data->dists[i][i] = 0;
  }
  for (auto &e : edges) {
    data->dists[e.from][e.to] = 1;
    data->next[e.from][e.to] = e.to;
  }

  REP (k, n) {
    REP (i, n)
      REP (j, n) {
        int nc = data->dists[i][k] + data->dists[k][j];
        if (nc < data->dists[i][j]) {
          data->dists[i][j] = nc;
          data->next[i][j] = data->next[i][k];
        }
      }
  }
  return *data;
}

std::vector<EdgeData> FastGraph::get_edges(int a, bool norm) const {
  a = normv(a, norm);
  std::vector<EdgeData> res;
  OPA_FOREACH_EDGE(e, a) { res.push_back(edges[e]); }
  return res;
}

std::vector<EdgeData> FastGraph::get_edges2(int a, int b, bool norm) const {
  a = normv(a, norm);
  b = normv(b, norm);
  std::vector<EdgeData> res;
  OPA_FOREACH_EDGE(e, a) if (edges[e].to == b) { res.push_back(edges[e]); }
  return res;
}

void FastGraph::dijkstra(DijkstraCtx &ctx) {
  if (ctx.prev.empty()) ctx.prev.resize(n, -1);
  if (ctx.dists.empty()) ctx.dists.resize(n, 1e100);

  std::vector<int> vis(n, 0);
  SmallPQ<std::pair<double, int> > q;
  for (auto &s : ctx.sources) {
    ctx.dists[s] = 0;
    ctx.prev[s] = -1;
    q.push(MP(0, s));
  }

  while (q.size() > 0) {
    double cost = q.top().first;
    int cur = q.top().second;
    q.pop();
    if (vis[cur]) continue;
    vis[cur] = 1;

    OPA_FOREACH_EDGE(e, cur) {
      int nxt = edges[e].to;
      OPA_CHECK0(get_edge_id(e) <= ctx.edge_costs.size());
      double edge_cost = ctx.edge_costs[get_edge_id(e)];
      double nc = edge_cost + cost;
      if (ctx.dists[nxt] > nc) {
        ctx.dists[nxt] = nc;
        ctx.prev[nxt] = cur;
        q.push(MP(nc, nxt));
      }
    }
  }
}

void FastGraph::bfs_one(BfsCtx &ctx) {
  for (auto &x : ctx.sources) {
    OPA_FOREACH_EDGE(e, x) {
      int next = edges[e].to;
      if (ctx.vis.insert(next).second) {
        ctx.waiting.push_back(next);
      }
    }
  }
}

int FastGraph::get_edge_count(int a, int b, bool norm) const {
  a = normv(a, norm);
  b = normv(b, norm);
  int cnt = 0;
  OPA_FOREACH_EDGE(e, a) {
    if (edges[e].to == b) ++cnt;
  }
  return cnt;
}

std::string FastGraph::str() const {
  std::stringstream ss;
  ss << glib::strings::Substitute("n_edges=$0, n_vertices=$1\n", n_edges(),
                                  vertices.size());
  REP (i, vertices.size()) {
    ss << glib::strings::Substitute("Vertex=$0, real_id=$1, deg=$2\n", i,
                                    inormv(i), vertices[i].deg);
    for (auto &edge : get_edges(i, false)) {
      ss << RAW_OPA_DISP_VARS(inormv(edge.to), edge.prev, edge.next)
         << std::endl;
    }
  }

  return ss.str();
}

std::vector<std::pair<pii, pii> >
FastGraph::find_set_neighbours(const std::vector<std::vector<int> > &lst,
                               bool norm) const {
  std::unordered_map<int, int> node_to_set;
  REP (i, lst.size()) {
    for (auto &v : lst[i]) {
      node_to_set[normv(v, norm)] = i;
    }
  }

  std::vector<std::pair<pii, pii> > res;
  for (auto &k : node_to_set) {
    int self_set = k.second;
    for (auto &e : edge_iter(k.first, false)) {
      int other = FindWithDefault(node_to_set, e.to, -1);
      if (other == -1 || other == self_set) continue;
      if (other <= self_set)
        continue; // superset of all previous conditions, clarity
      res.emplace_back(MP(inormv(k.first), self_set), MP(inormv(e.to), other));
    }
  }
  return res;
}

OPA_NAMESPACE_DECL2_END

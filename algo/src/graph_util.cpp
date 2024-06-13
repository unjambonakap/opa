#include <opa/algo/graph_util.h>
#include <opa/utils/DataStruct.h>

OPA_NAMESPACE_DECL2(opa, algo)

char deser_buf[2048];
class DigraphConnectedComponentsHelper {
public:
  DigraphConnectedComponentsHelper(const FastGraph &graph) : graph(graph) {}
  void compute() {
    n = graph.n;
    depth.resize(n, -1);
    cc.resize(n, -1);

    REP (i, n) {
      dfs(i, n - i);
      OPA_CHECK0(opened_vertices.size() == 0);
      OPA_CHECK0(opened_comps.size() == 0);
    }
  }

  int dfs(int a, int d) {
    if (depth[a] != -1) return cc[a] == -1 ? depth[a] : d;
    depth[a] = d;
    int pos_v = opened_vertices.size();
    int pos_c = opened_comps.size();
    int md = d;

    for (auto &nv : graph.get_successors(a, false)) {
      int nd = dfs(nv, d + 1);
      if (nv > d) {
        opened_comps.push_back(cc[nv]);
      }
      md = std::min(md, nd);
    }

    opened_vertices.push_back(a);
    if (md == d) {
      CComponent cur;
      int ccid = res.size();
      for (; opened_vertices.size() > pos_v; opened_vertices.pop_back()) {
        int v = opened_vertices.back();
        cc[v] = ccid;
        cur.nodes.push_back(graph.inormv(v));
      }
      for (; opened_comps.size() > pos_c; opened_comps.pop_back()) {
        cur.next.push_back(opened_comps.back());
      }
      res.push_back(cur);
    }
    return md;
  }

  std::vector<int> opened_vertices;
  std::vector<int> opened_comps;
  std::vector<int> depth;
  std::vector<int> cc;
  int n;
  std::vector<CComponent> res;
  const FastGraph &graph;
};

std::vector<CComponent> compute_digraph_connected_components(const FastGraph &graph) {

  DigraphConnectedComponentsHelper helper(graph);
  helper.compute();
  return helper.res;
}

struct TreeCtx {
  const FastGraph &tree;
  std::vector<int> par;
  std::vector<int> depth;
  std::vector<int> rdepth;
  std::vector<std::vector<int> > children;
  int n;
  int root;
  std::vector<int> dfs_order;
  std::vector<int> tree_id;
  std::vector<std::vector<int> > automorphisms;

  utils::MaxFinderPair<int, std::pair<int, int> > longuest_path;

  TreeCtx(const FastGraph &tree, int root = 0) : tree(tree) {
    this->root = root;
    n = tree.n;
    children.resize(n);
    par.resize(n);
    depth.resize(n);
    rdepth.resize(n);
    dfs_init(root, -1);
  }

  int dfs_init(int a, int p) {
    par[a] = p;
    depth[a] = (p == -1 ? 0 : depth[p] + 1);
    dfs_order.push_back(a);
    int mx = 0;
    for (auto &b : tree.get_successors(a)) {
      if (b == p) continue;
      children[a].push_back(b);
      mx = std::max(mx, dfs_init(b, a) + 1);
    }
    rdepth[a] = mx;
    return mx;
  }

  std::pair<int, int> get_longuest_path() {
    build_longuest_path(root);
    return longuest_path.get();
  }

  std::pair<int, int> build_longuest_path(int a) {
    std::pair<int, int> md(a, 0);
    utils::KMaxFinderPair<2, int, int> kmf;

    for (auto &nx : tree.get_successors(a)) {
      if (nx == par[a]) continue;
      auto cur = build_longuest_path(nx);
      kmf.update(cur.second, cur.first);
    }

    if (kmf.has(0)) {
      md = MP(kmf.get(0), kmf.get_cost(0) + 1);
      longuest_path.update(md.second + 1, MP(a, md.first));
    }
    if (kmf.has(1)) {
      longuest_path.update(kmf.get_cost(0) + kmf.get_cost(1), MP(kmf.get(0), kmf.get(1)));
    }
    return md;
  }

  void compute_automorphism(const std::vector<int> &nodes) {
    std::unordered_map<int, std::vector<int> > mp;
    for (auto &node : nodes) mp[tree_id[node]].push_back(node);

    for (auto &s : mp) {
      automorphisms.push_back(s.second);
      std::vector<int> next;
      for (auto &x : s.second) next.insert(next.end(), ALL(children[x]));
      compute_automorphism(next);
    }
  }
};

std::pair<int, int> get_tree_center(const FastGraph &tree) {
  std::pair<int, int> res;

  TreeCtx ctx(tree);
  auto lp = ctx.get_longuest_path();
  int a = lp.first;
  int b = lp.second;

  while (ctx.depth[a] > ctx.depth[b]) a = ctx.par[a];
  while (ctx.depth[b] > ctx.depth[a]) b = ctx.par[b];
  while (a != b) a = ctx.par[a], b = ctx.par[b];
  int mid = a;
  a = lp.first;
  b = lp.second;

  int da = ctx.depth[a] - ctx.depth[mid];
  int db = ctx.depth[b] - ctx.depth[mid];
  int len = da + db;

  int midp = len / 2;
  int x = (da > midp ? a : b);
  REP (i, midp) x = ctx.par[x];
  if (len % 2 == 1) return { x, ctx.par[x] };
  return { x, -1 };
}

std::vector<std::vector<int> > get_tree_automorphism_partitions(const FastGraph &tree) {
  auto center = get_tree_center(tree);
  OPA_CHECK0(center.second == -1); // atm only support vertex centered trees

  TreeCtx ctx(tree, center.first);
  std::vector<std::vector<int> > rdepth_to_nodes(ctx.n);
  REP (i, ctx.n) rdepth_to_nodes[ctx.rdepth[i]].push_back(i);

  int nid = 0;
  std::vector<int> tree_id(ctx.n);
  for (auto &lst : rdepth_to_nodes) {
    std::vector<std::pair<std::vector<int>, int> > node_to_children;
    for (auto &x : lst) {
      std::vector<int> sig;
      for (auto &c : ctx.children[x]) sig.push_back(tree_id[c]);
      std::sort(ALL(sig));
      node_to_children.emplace_back(sig, x);
    }
    std::sort(ALL(node_to_children));
    REP (i, node_to_children.size()) {
      int cid =
        i == 0 || (node_to_children[i].first != node_to_children[i - 1].first) ? nid++ : nid - 1;
      tree_id[node_to_children[i].second] = cid;
    }
  }

  ctx.tree_id = tree_id;
  ctx.compute_automorphism({ ctx.root });
  return ctx.automorphisms;
}

std::vector<int> topological_ordering(const FastGraph &graph) {
  std::vector<int> res;
  std::unordered_map<int, int> degree_rem;
  std::queue<int> q;
  REP (i, graph.node_count()) {
    degree_rem[i] = graph.vertices[i].deg_in;
    if (degree_rem[i] == 0) q.push(i);
  }

  while (!q.empty()) {
    int v = q.front();
    q.pop();
    res.pb(graph.inormv(v));
    for (auto &edge : graph.get_edges(v, false)) {
      if (--degree_rem[edge.to] == 0) q.push(edge.to);
    }
  }
  if (res.size() != graph.node_count()) return {};
  OPA_CHECK0(res.size() == graph.node_count());

  return res;
}

UPTR(FastGraph)
compress_digraph(const FastGraph &graph, const std::unordered_set<int> &keep_nodes) {

  FastGraph graph_nokeep(graph.n, Mode::DIGRAPH);
  for (auto &edge : graph.list_edges()) {
    if (keep_nodes.count(edge.from) || keep_nodes.count(edge.to)) continue;
    graph_nokeep.adde(edge.from, edge.to);
  }

  auto ordering = topological_ordering(graph_nokeep);
  std::reverse(ALL(ordering));

  std::unordered_map<int, std::unordered_map<int, utils::MinFinder<int> > > node_to_next;
  for (auto &v : ordering) {
    std::unordered_map<int, utils::MinFinder<int> > entry;
    for (auto &e : graph.get_successors(v)) {
      if (keep_nodes.count(e))
        entry[e].update(1);
      else {
        for (auto &u : node_to_next[e]) entry[u.first].update(u.second.get() + 1);
        entry.insert(ALL(node_to_next[e]));
      }
    }
    node_to_next[v] = entry;
  }

  auto res = std::make_unique<FastGraph>(0, Mode::DIGRAPH | Mode::REMAP | Mode::DYNAMIC_SIZE);
  REP (v, graph.n) {
    if (!keep_nodes.count(v)) continue;
    std::unordered_map<int, utils::MinFinder<int> > nexts;
    for (auto &e : graph.get_successors(v)) {
      nexts.insert(ALL(node_to_next[e]));
    }
    res->normv(v);
    for (auto &to : nexts) {
      int eid = res->adde(v, to.first);
      res->get_edge_data(eid).data = to.second.get();
    }
  }
  graph.attr_map.extract(utils::to_vector(keep_nodes), &res->attr_map);
  OPA_CHECK0(res->node_count() == keep_nodes.size());

  return res;
}

void contract_edge(FastGraph &graph, int a, int b) {
  auto b_pred = graph.get_predecessors(b);
  auto b_suc = graph.get_successors(b);
  for (auto &nx : b_pred) {
    graph.remove_edge(nx, b);
    if (nx == a) continue;
    OPA_CHECK0(nx != b);
    graph.adde(nx, a);
  }

  for (auto &nx : b_suc) {
    graph.remove_edge(b, nx);
    if (nx == a) continue;
    OPA_CHECK0(nx != b);
    graph.adde(a, nx);
  }
}

SPTR(algo::UnionJoin)
compress_paths(FastGraph &graph, std::unordered_set<int> blacklist) {
  std::unordered_set<pii> targets;
  auto uj = std::make_shared<algo::UnionJoin>(graph.node_count());

  REP (i, graph.node_count()) {
    if (graph.vertices[i].deg_in == 1) targets.insert(MP(graph.get_predecessors(i, false)[0], i));
    if (graph.vertices[i].deg_out == 1) targets.insert(MP(i, graph.get_successors(i, false)[0]));
  }

  for (auto &target : targets) {
    if (blacklist.count(graph.inormv(target.first)) || blacklist.count(graph.inormv(target.second)))
      continue;

    int a = uj->root(target.first);
    int b = uj->root(target.second);
    if (a == b) continue;
    if (graph.deg_out(a, false) != 1 && graph.deg_in(b, false)) continue;
    if (graph.deg(a, false) < graph.deg(b, false)) std::swap(a, b);

    contract_edge(graph, graph.inormv(a), graph.inormv(b));
    OPA_CHECK0(graph.deg(b, false) == 0);
    OPA_CHECK0(graph.deg(b, false) ==
               graph.get_successors(b, false).size() + graph.get_predecessors(b, false).size());
    OPA_CHECK0(graph.deg(a, false) ==
               graph.get_successors(a, false).size() + graph.get_predecessors(a, false).size());
    uj->set_repr(a, b);
  }
  return uj;
}

void serialize_attr_t(FILE *f, const std::unordered_map<int, int> &mp) {
  fprintf(f, "%lu\n", mp.size());
  for (auto &kv : mp) {
    fprintf(f, "%d,%d\n", kv.first, kv.second);
  }
}

void serialize_attr_t(FILE *f, const std::unordered_map<int, std::string> &mp) {
  fprintf(f, "%lu\n", mp.size());
  for (auto &kv : mp) {
    fprintf(f, "%d,%s\n", kv.first, kv.second.c_str());
  }
}

void serialize_attr(FILE *f, const AttributeContainer &container) {
  fprintf(f, "%d\n", container.attr);
  if (container.attr == AttrType::ATTR_INT)
    serialize_attr_t(f, container.int_map);
  else if (container.attr == AttrType::ATTR_STR)
    serialize_attr_t(f, container.str_map);
}

void serialize_attrmap(FILE *f, const AttributeMap &attr_map) {
  fprintf(f, "%lu\n", attr_map.attribute_map.size());
  for (const auto &kv : attr_map.attribute_map) {
    fprintf(f, "%s\n", kv.first.c_str());
    serialize_attr(f, *kv.second);
  }
}

void serialize_graph(FILE *f, const FastGraph &graph) {
  fprintf(f, "%d,%d\n", graph.node_count(), graph.n_edges());
  REP (i, graph.node_count()) {
    fprintf(f, "%d\n", graph.inormv(i));
  }
  for (auto &edge : graph.list_edges()) {
    if (!edge.forward) continue;
    fprintf(f, "%d,%d,%d\n", edge.from, edge.to, edge.data);
  }
  serialize_attrmap(f, graph.attr_map);
}

void deserialize_attr_t(FILE *f, std::unordered_map<int, int> &mp) {
  int count;
  fscanf(f, "%d", &count);
  REP (i, count) {
    int a, b;
    fscanf(f, "%d,%d", &a, &b);
    mp[a] = b;
  }
}

void deserialize_attr_t(FILE *f, std::unordered_map<int, std::string> &mp) {
  int count;
  fscanf(f, "%d", &count);
  REP (i, count) {
    int node;
    fscanf(f, "%d,", &node);
    fgets(deser_buf, sizeof(deser_buf), f);
    deser_buf[strlen(deser_buf) - 1] = 0;

    mp[node] = deser_buf;
  }
}

void deserialize_attr(FILE *f, AttributeContainer &container) {
  fscanf(f, "%d", &container.attr);
  if (container.attr == AttrType::ATTR_INT)
    deserialize_attr_t(f, container.int_map);
  else if (container.attr == AttrType::ATTR_STR)
    deserialize_attr_t(f, container.str_map);
}

void deserialize_attrmap(FILE *f, AttributeMap *attr_map) {
  int attrcount;
  fscanf(f, "%d", &attrcount);
  REP (i, attrcount) {
    fscanf(f, "%s\n", deser_buf);
    auto attrc = std::make_shared<AttributeContainer>();
    attr_map->attribute_map[deser_buf] = attrc;
    deserialize_attr(f, *attrc);
  }
}

FastGraph deserialize_graph(FILE *f, int mode) {
  int n, m;
  fscanf(f, "%d,%d", &n, &m);
  FastGraph res(n, mode);
  REP (i, n) {
    int a;
    fscanf(f, "%d,", &a);
    res.add_node(a);
  }

  REP (i, m) {
    int a, b, c;
    OPA_DISP0(a, b, c);
    fscanf(f, "%d,%d,%d", &a, &b, &c);
    res.adde(a, b);
  }
  deserialize_attrmap(f, &res.attr_map);

  return res;
}

OPA_NAMESPACE_DECL2_END

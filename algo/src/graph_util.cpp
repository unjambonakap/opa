#include <opa/algo/graph_util.h>
#include <opa/utils/DataStruct.h>

OPA_NAMESPACE_DECL2(opa, algo)

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
    if (depth[a] != -1) return depth[a];
    depth[a] = d;
    int pos_v = opened_vertices.size();
    int pos_c = opened_comps.size();
    int md = d;

    OPA_FOREACH_EDGE_OBJ(e, a, graph) {
      int nv = graph.edges[e].to;
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
        cur.nodes.push_back(v);
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

std::vector<CComponent>
compute_digraph_connected_components(const FastGraph &graph) {

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
  std::vector<std::vector<int>> automorphisms;

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
    OPA_FOREACH_EDGE_OBJ(e, a, tree) {
      int b = tree.edges[e].to;
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

    OPA_FOREACH_EDGE_OBJ(e, a, tree) {
      int nx = tree.edges[e].to;
      if (nx == par[a]) continue;
      auto cur = build_longuest_path(nx);
      kmf.update(cur.second, cur.first);
    }

    if (kmf.has(0)) {
      md = MP(kmf.get(0), kmf.get_cost(0) + 1);
      longuest_path.update(md.second + 1, MP(a, md.first));
    }
    if (kmf.has(1)) {
      longuest_path.update(kmf.get_cost(0) + kmf.get_cost(1),
                           MP(kmf.get(0), kmf.get(1)));
    }
    return md;
  }

  void compute_automorphism(const std::vector<int> &nodes){
    std::unordered_map<int, std::vector<int>> mp;
    for (auto & node : nodes)
      mp[tree_id[node]].push_back(node);

    for (auto &s : mp){
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
  REP (i, midp)
    x = ctx.par[x];
  if (len % 2 == 1) return { x, ctx.par[x] };
  return { x, -1 };
}

std::vector<std::vector<int> >
get_tree_automorphism_partitions(const FastGraph &tree) {
  auto center = get_tree_center(tree);
  OPA_CHECK0(center.second == -1); // atm only support vertex centered trees

  TreeCtx ctx(tree, center.first);
  std::vector<std::vector<int> > rdepth_to_nodes(ctx.n);
  REP (i, ctx.n)
    rdepth_to_nodes[ctx.rdepth[i]].push_back(i);

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
        i == 0 || (node_to_children[i].first != node_to_children[i - 1].first)
          ? nid++
          : nid-1;
      tree_id[node_to_children[i].second] = cid;
    }
  }

  ctx.tree_id = tree_id;
  ctx.compute_automorphism({ctx.root});
  return ctx.automorphisms;
}

OPA_NAMESPACE_DECL2_END

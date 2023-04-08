#pragma once

#include <opa/algo/graph.h>
#include <opa/utils/hash.h>

OPA_NAMESPACE_DECL2(opa, algo)

template <class T> class DigraphComputation {
public:
  typedef std::function<T(int a, const FastGraph &graph,
                          const std::function<T(int)>)>
    FuncType;
  DigraphComputation(const FastGraph &graph, const FuncType &func)
      : m_graph(graph), m_func(func) {
    comp_mem.init(this);
  }

  T comp_(int node) const {
    return m_func(node, m_graph, std::bind(&DigraphComputation<T>::get, this,
                                           std::placeholders::_1));
  }

  T get(int node) const { return CALL1(DigraphComputation<T>, comp, node); }

  utils::MemoizeHelperClass<T, DigraphComputation<T>, int> comp_mem;
  FuncType m_func;
  FastGraph m_graph;
};

static int dist_comp(int x, const FastGraph &graph,
              const std::function<int(int)> &getter) {
  utils::MinFinder<int> mf;

  for (auto &e : graph.edge_iter(x)) {
    mf.update(getter(e.to) + 1);
  }
  return mf.get_or(0);
}

static int grundy_comp(int x, const FastGraph &graph,
                const std::function<int(int)> &getter) {
  std::set<int> mex;

  for (auto &e : graph.edge_iter(x)) {
    mex.insert(getter(e.to));
  }
  for (int i = 0;; ++i) {
    if (!mex.count(i)) return i;
  }
}


struct CComponent {
  std::vector<int> nodes;
  std::vector<int> next;
};


std::pair<int,int> get_tree_center(const FastGraph &tree);
std::vector<std::vector<int>> get_tree_automorphism_partitions(const FastGraph &tree);
void contract_edge(FastGraph &graph, int a, int b);


std::vector<CComponent>
compute_digraph_connected_components(const FastGraph &graph);

std::vector<int> topological_ordering(const FastGraph &graph);
UPTR(FastGraph) compress_digraph(const FastGraph &graph,
    const std::unordered_set<int> &keep_nodes);
SPTR(algo::UnionJoin) compress_paths(FastGraph &graph, std::unordered_set<int> blacklist) ;

void serialize_graph(FILE *f, const FastGraph &graph) ;
FastGraph deserialize_graph(FILE *f, int mode) ;



OPA_NAMESPACE_DECL2_END

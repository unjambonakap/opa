#include <gtest/gtest.h>

#include <opa/algo/base.h>
#include <opa/algo/graph.h>
#include <opa/algo/graph_util.h>
#include <opa/algo/sat2.h>
#include <opa/utils/range.h>
#include <unordered_set>

using namespace opa::utils;
using namespace opa::algo;
using namespace std;

void init() {}

TEST(Graph, Test1) {
  FastGraph graph(10, Mode::REMAP);

  graph.add_bidirectional(1, 20);
  graph.add_bidirectional(20, 5);
  graph.add_bidirectional(1, 5);
  graph.add_bidirectional(2, 5);
  graph.add_bidirectional(2, 3);
  graph.add_bidirectional(3, 4);
  graph.add_bidirectional(3, 6);
  REP (i, graph.n) OPA_DISP0(graph.inormv(i), graph.deg(i, false));

  auto res = graph.get_cover_walk();
  OPA_TRACES(res);
}

TEST(CC, Test1) {
  FastGraph graph(7, Mode::DIGRAPH);

  graph.adde(0, 1);
  graph.adde(1, 2);
  graph.adde(2, 3);
  graph.adde(3, 4);
  graph.adde(4, 5);
  graph.adde(5, 3);
  graph.adde(2, 1);
  graph.adde(4, 6);
  auto res = compute_digraph_connected_components(graph);
  for (auto &cc : res) {
    OPA_DISP0(cc.nodes, cc.next);
  }
}

TEST(SAT2, Test1) {
  Sat2<int> sat;
  sat += sat.term(0) == false;
  sat += (sat.term(0) || 1) == true;
  sat += (sat.term(1) + 2) == true;
  sat += (sat.term(1) ^ 3) == true;
  sat += (sat.term(0) + 4) == true;
  sat += (sat.term(0) ^ 5) == true;
  OPA_CHECK0(sat.compute());
  REP (i, 6) OPA_DISP0(i, sat.get(i));
}

TEST(SAT2, Test2) {
  Sat2<pii> sat;
  sat += sat.term(0, 1) == false;
  sat += (sat.term(0, 3) || pii(1, 3)) == true;
  sat += sat.term(0, 1) == sat.term(2, 3);
  OPA_CHECK0(sat.compute());
  OPA_DISP0(sat.get2(0, 1));
  OPA_DISP0(sat.get2(0, 3));
  OPA_DISP0(sat.get2(1, 3));
  OPA_DISP0(sat.get2(2, 3));
}

TEST(SAT2, Test3) {
  Sat2<pii> sat;
  sat += sat.term(0, 1) * (s8)-1 + sat.term(0, 0) * (s8)-1 == (s8)-2;
  OPA_CHECK0(sat.compute());
  OPA_DISP0(sat.get2(0, 0));
  OPA_DISP0(sat.get2(0, 1));
}

TEST(SAT2, TestFail) {
  Sat2<int> sat;
  sat += sat.term(0) == false;
  sat += (sat.term(0)) == true;
  OPA_CHECK0(!sat.compute());
}

TEST(PathCover, Test1) {
  FastGraph graph(0, Mode::DYNAMIC_SIZE);
  graph.add_bidirectional(0, 1);
  graph.add_bidirectional(1, 2);
  graph.add_bidirectional(2, 3);
  graph.add_bidirectional(3, 4);
  graph.add_bidirectional(5, 1);
  graph.add_bidirectional(6, 3);
  auto cover = graph.get_path_cover();
  OPA_DISP0(cover);
}

TEST(PathCover, Test2) {
  FastGraph graph(0, Mode::DYNAMIC_SIZE);
  graph.add_bidirectional(0, 1);
  graph.add_bidirectional(1, 2);
  graph.add_bidirectional(2, 3);
  graph.add_bidirectional(3, 4);
  graph.add_bidirectional(5, 1);
  graph.add_bidirectional(6, 3);
  graph.add_bidirectional(2, 7);
  graph.add_bidirectional(7, 8);
  graph.add_bidirectional(8, 3);
  auto graph2 = graph;
  auto cover = graph.get_path_cover();
  OPA_DISP0(cover);
  auto cover2 = graph2.get_path_cover_dumb();
  OPA_DISP0(cover2);
}

TEST(TreeAutomorphism, Test1) {
  FastGraph tree(4, Mode::NONE);
  tree.add_bidirectional(0, 1);
  tree.add_bidirectional(1, 2);
  tree.add_bidirectional(1, 3);
  OPA_DISP0(get_tree_center(tree));

  FastGraph tree2(6, Mode::NONE);
  tree2.add_bidirectional(0, 1);
  tree2.add_bidirectional(0, 2);
  tree2.add_bidirectional(2, 3);
  tree2.add_bidirectional(3, 4);
  tree2.add_bidirectional(3, 5);
  OPA_DISP0(get_tree_center(tree2));
}

TEST(TreeAutomorphism, Test2) {
  FastGraph tree(4, Mode::NONE);
  tree.add_bidirectional(0, 1);
  tree.add_bidirectional(1, 2);
  tree.add_bidirectional(1, 3);
  OPA_DISP0(get_tree_center(tree));
  OPA_DISP0(get_tree_automorphism_partitions(tree));
}

TEST(DigraphCompress, Test1) {
  FastGraph graph(10, Mode::DIGRAPH);
  graph.adde(0, 1);
  graph.adde(1, 2);
  graph.adde(2, 3);
  graph.adde(3, 4);
  graph.adde(4, 5);
  graph.adde(4, 6);
  graph.adde(1, 7);
  std::unordered_set<int> keep = { 0, 5, 6, 7, 3 };
  auto ng = compress_digraph(graph, keep);
  OPA_DISP0(ng->str());
}

TEST(TopologicalOrdering, Test1) {
  FastGraph graph(10, Mode::DIGRAPH);
  graph.adde(0, 1);
  graph.adde(1, 2);
  graph.adde(2, 3);
  graph.adde(3, 4);
  graph.adde(4, 5);
  graph.adde(4, 6);
  graph.adde(1, 7);
  OPA_DISP0(graph.get_successors(4));
  OPA_DISP0(graph.get_predecessors(4));
  auto tp = topological_ordering(graph);
  OPA_DISP0(tp);
}

TEST(ContractEdge, Test1) {
  FastGraph graph(10, Mode::DIGRAPH);
  graph.adde(0, 1);
  graph.adde(1, 2);
  graph.adde(2, 3);
  graph.adde(3, 4);
  graph.adde(4, 5);
  graph.adde(4, 6);
  graph.adde(1, 7);

  contract_edge(graph, 0, 1);

  contract_edge(graph, 3, 4);
  OPA_DISP0(graph.str());
}

TEST(CompressPaths, Test1) {
  FastGraph graph(10, Mode::DIGRAPH);
  graph.adde(0, 1);
  graph.adde(1, 2);
  graph.adde(2, 3);
  graph.adde(3, 4);
  graph.adde(4, 5);
  graph.adde(4, 6);
  graph.adde(1, 7);
  graph.adde(5, 6);

  OPA_DISP0(graph.str());
  compress_paths(graph, { 1 });
  OPA_DISP0(graph.str());
}

TEST(CompressPaths, Test2) {
  FastGraph graph(4, Mode::DIGRAPH);
  graph.adde(0, 1);
  graph.adde(0, 2);
  graph.adde(1, 3);
  graph.adde(2, 3);

  OPA_DISP0(graph.str());
  compress_paths(graph, { 1 });
  OPA_DISP0(graph.str());
}

TEST(Subgraph, Test1) {
  FastGraph graph(10, Mode::DIGRAPH);
  graph.adde(0, 1);
  graph.adde(1, 2);
  graph.adde(2, 3);
  graph.adde(3, 4);
  graph.adde(4, 5);
  graph.adde(4, 6);
  graph.adde(1, 7);
  graph.adde(5, 6);

  OPA_DISP0(graph.str());
  auto sg = graph.subgraph({ 1, 2, 3, 6 }, true, true);
  OPA_DISP0(sg->str());
  auto sg2 = sg->subgraph({ 1, 3 }, true, true);
  OPA_DISP0(sg2->str());
}

struct BrushSolver {
  VV_t(bool) src, target;
  int n() const { return src.size(); }
  int m() const { return src[0].size(); }

  int key(int x, bool is_row) const { return 1 + x + (is_row ? 0 : n()); }

  std::vector<std::pair<int, bool> > solve() const {
    FastGraph graph(1 + n() + m(), Mode::DIGRAPH);
    constexpr int root = 0;
    constexpr bool row_action_col = false;
    REP (i, n())
      REP (j, m()) {
        int row_action = key(i, true);
        int col_action = key(j, false);
        bool tc = target[i][j];
        bool sc = src[i][j];
        if (tc == row_action_col)
          graph.adde(col_action, row_action);
        else
          graph.adde(row_action, col_action);

        if (sc != tc) {
          graph.adde(root, tc == row_action_col ? row_action : col_action);
        }
      }

    std::unordered_set<int> seen;
    V_t(int) lst;
    graph.dfs_cc(root, seen, lst);
    auto sg = graph.subgraph(lst);
    auto ord = topological_ordering(*sg);
    if (ord.empty()) return { { -1, false } };
    ord.erase(ord.begin());
    return ord |
           STD_TSFX(x < key(0, false) ? MP(x - key(0, true), true) : MP(x - key(0, false), false)) |
           STD_VEC;
  }
};

TEST(Brush, Test1) {
  BrushSolver solver{ .src = { { 1, 1, 0 }, { 1, 0, 0 }, { 0, 1, 1 } },
                      .target = { { 1, 0, 1 }, { 1, 0, 0 }, { 1, 1, 1 } } };

  auto res = solver.solve();
  OPA_DISP0(res);
}

GTEST_API_ int main(int argc, char **argv) {
  printf("Running main() from gtest_main.cc\n");
  testing::InitGoogleTest(&argc, argv);
  puts("1");
  opa::init::opa_init(argc, argv);
  puts("2");
  return RUN_ALL_TESTS();
}

#pragma once

#include <opa/algo/graph.h>
#include <opa/math/game/conf.h>
#include <opa/utils/DataStruct.h>

OPA_NAMESPACE_DECL3(opa, math, game)

struct EdgeBase : public opa::utils::IdObj {
public:
  GraphObjId attr_id;
  VertexId start;
  VertexId end;
};

class VertexBase : public opa::utils::IdObj {
public:
  GraphObjId attr_id;
};

template <class Vertex, class Edge> class Graph : public opa::utils::Initable {
public:
  virtual void init() override { opa::utils::Initable::init(); }

  virtual Vertex &create_vertex() {
    Vertex &v = m_vertices.get_new2();
    return v;
  }

  virtual Edge &add_edge(VertexId a, VertexId b) {
    Edge &edge = m_edges.get_new2();
    edge.start = a;
    edge.end = b;
    return edge;
  }

  Vertex &vertex(VertexId id) { return m_vertices.get(id); }

  Edge &edge(EdgeId id) { return m_edges.get(id); }
  /*
  SPTR(algo::FastGraph) compute_graph() const {
    SPTR(algo::FastGraph)
    graph = std::make_shared<algo::FastGraph>(m_vertices.size(), true);

    for (auto edge : m_edges.size()) {
      graph->add_bidirectional(edge.start, edge.end);
    }

    return graph;
  }
  */

private:
  utils::ObjectPool<Edge> m_edges;
  utils::ObjectPool<Vertex> m_vertices;
};

class VertexWithPos : public VertexBase {
public:
  Pos pos;
};

using GraphWithPos = Graph<VertexWithPos, EdgeBase>;

OPA_NAMESPACE_DECL3_END

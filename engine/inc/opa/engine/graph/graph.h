#pragma once

#include <opa/engine/conf.h>
#include <opa/engine/Drawer.h>
#include <opa/engine/Program.h>
#include <opa/engine/Scene.h>
#include <opa/math/game/graph.h>
#include <opa/utils/DataStruct.h>
#include <type_traits>

OPA_NAMESPACE_DECL3(opa, engine, graph)

extern Program point_prog;
extern Program line_prog;

extern Attribute point_unif_cam_left;
extern Attribute point_unif_cam_up;
extern Attribute point_unif_view;
extern Attribute point_attr_sq_pos;
extern Attribute point_attr_pos;
extern Attribute point_attr_color;

Program *get_program_point();

extern Attribute line_attr_sq_pos;
extern Attribute line_attr_pos1;
extern Attribute line_attr_pos2;
extern Attribute line_attr_color;
extern Attribute line_unif_cam_left;
extern Attribute line_unif_cam_up;
extern Attribute line_unif_view;

Program *get_program_line();

struct GraphEdgeAttribute {
  OPA_IDOBJ_FIELD();
  Pos4 pos1;
  Pos4 pos2;
  Col col;
  GraphEdgeAttribute() {}
};

struct GraphVertexAttribute {
  OPA_IDOBJ_FIELD();
  Pos4 pos;
  Col col;
  GraphVertexAttribute() {}
};

struct RawVertex2d {
  Pos2 pos;
};

template <class Vertex, class Edge>
class DisplayGraph : public opa::math::game::Graph<Vertex, Edge> {
public:
  enum Visibility : int {
    Hidden = 0,
    Visible = 1,
  };

  virtual void init() override;

  void draw(const DrawData &data, SceneObj *cur);
  virtual Vertex &create_vertex() override;
  virtual Edge &add_edge(VertexId a, VertexId b) override;

  GraphEdgeAttribute &edge_attr(EdgeId id) {
    return m_attr_edges.get(this->edge(id).attr_id);
  }
  GraphVertexAttribute &vertex_attr(VertexId id){
    return m_attr_vertices.get(this->vertex(id).attr_id);
  }

  void set_vertex_visibility(VertexId id, Visibility vis);
  void set_edge_visibility(EdgeId id, Visibility vis);

  void update_edge(EdgeId id);

  SceneObj *create_scene_obj() {
    SceneObj *res = new SceneObj;
    auto drawer = new FuncDrawer;
    drawer->init([this, res](const DrawData &data) { this->draw(data, res); });
    res->drawer().register_drawer(DrawerPtr(drawer));
    return res;
  }

private:
  utils::BucketContBufferContainer<GraphEdgeAttribute> m_attr_edges;
  utils::BucketContBufferContainer<GraphVertexAttribute> m_attr_vertices;
  utils::ObjectPool<Edge> m_edges;
  utils::ObjectPool<Vertex> m_vertices;

  Buffer m_vertex_buffer;
  Buffer m_pos_buffer;
  Buffer m_index_buffer;
  std::vector<u32> m_index;
  Program *m_point_prog;
  Program *m_line_prog;

  std::vector<RawVertex2d> m_vertex_poly;
};

template <class Vertex, class Edge> void DisplayGraph<Vertex, Edge>::init() {
  opa::math::game::Graph<Vertex, Edge>::init();
  m_point_prog = get_program_point();
  m_line_prog = get_program_line();

  m_attr_edges.create_bucket(Visibility::Visible);
  m_attr_edges.create_bucket(Visibility::Hidden);

  m_attr_vertices.create_bucket(Visibility::Visible);
  m_attr_vertices.create_bucket(Visibility::Hidden);

  m_index_buffer.init(GL_ELEMENT_ARRAY_BUFFER);
  m_vertex_buffer.init(GL_ARRAY_BUFFER);
  m_pos_buffer.init(GL_ARRAY_BUFFER);

  m_vertex_poly.resize(4);
  REP (i, 4)
    m_vertex_poly[i].pos = Pos2(i / 2, i % 2);

  m_pos_buffer.bind();
  ScopedBuffer scope1(m_pos_buffer);
  glBufferData(GL_ARRAY_BUFFER, sizeof(m_vertex_poly[0]) * m_vertex_poly.size(),
               m_vertex_poly.data(), GL_STATIC_DRAW);
}

template <class Vertex, class Edge>
void DisplayGraph<Vertex, Edge>::draw(const DrawData &data, SceneObj *cur) {
  if (m_index.size() == 0)
    return;
  {
    m_point_prog->bind();

    point_unif_view.bind();
    point_unif_cam_up.bind();
    point_unif_cam_left.bind();
    point_attr_pos.bind();
    point_attr_color.bind();
    point_attr_sq_pos.bind();
    GL_CHECK(glUniformMatrix4fv(point_unif_view.id(), 1, GL_FALSE,
                                glm::value_ptr(data.proj())));
    Dir cam_up = data.base()->pos()->get_up(cur);
    Dir cam_left = data.base()->pos()->get_left(cur);

    GL_CHECK(glUniform3fv(point_unif_cam_up.id(), 1, glm::value_ptr(cam_up)));
    GL_CHECK(
      glUniform3fv(point_unif_cam_left.id(), 1, glm::value_ptr(cam_left)));

    m_vertex_buffer.bind();
    glBufferData(
      GL_ARRAY_BUFFER, m_attr_vertices.get_bucket_size(Visibility::Visible) *
                         sizeof(GraphVertexAttribute),
      m_attr_vertices.get_start(Visibility::Visible), GL_DYNAMIC_DRAW);
    m_vertex_buffer.unbind();

    int id_attr = 0;
    m_vertex_buffer.vertex_bind(id_attr, 0, sizeof(GraphVertexAttribute));
    int id_pos = 1;
    m_pos_buffer.vertex_bind(id_pos, 0, sizeof(m_vertex_poly[0]));

    glVertexAttribFormat(point_attr_sq_pos.id(), 2, GL_FLOAT, GL_FALSE,
                         offsetof(RawVertex2d, pos));
    glVertexAttribBinding(point_attr_sq_pos.id(), id_pos);

    glVertexAttribFormat(point_attr_pos.id(), 4, GL_FLOAT, GL_FALSE,
                         offsetof(GraphVertexAttribute, pos));
    GL_CHECK(glVertexAttribBinding(point_attr_pos.id(), id_attr));

    glVertexAttribFormat(point_attr_color.id(), 3, GL_FLOAT, GL_FALSE,
                         offsetof(GraphVertexAttribute, col));
    glVertexAttribBinding(point_attr_color.id(), id_attr);

    glVertexBindingDivisor(id_pos, 0);
    glVertexBindingDivisor(id_attr, 1);

    glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, m_index.size());
    glVertexBindingDivisor(id_attr, 0);

    point_attr_pos.unbind();
    point_attr_color.unbind();
    point_attr_sq_pos.unbind();

    m_pos_buffer.vertex_unbind(id_attr);
    m_vertex_buffer.vertex_unbind(id_pos);

    point_unif_view.unbind();
    point_unif_cam_up.unbind();
    point_unif_cam_left.unbind();
  }

  {
    m_line_prog->bind();

    line_unif_view.bind();
    line_unif_cam_up.bind();
    line_unif_cam_left.bind();
    line_attr_pos1.bind();
    line_attr_pos2.bind();
    line_attr_color.bind();
    line_attr_sq_pos.bind();
    GL_CHECK(glUniformMatrix4fv(line_unif_view.id(), 1, GL_FALSE,
                                glm::value_ptr(data.proj())));
    Dir cam_up = data.base()->pos()->get_up(cur);
    Dir cam_left = data.base()->pos()->get_left(cur);

    GL_CHECK(glUniform3fv(line_unif_cam_up.id(), 1, glm::value_ptr(cam_up)));
    GL_CHECK(
      glUniform3fv(line_unif_cam_left.id(), 1, glm::value_ptr(cam_left)));

    m_vertex_buffer.bind();
    glBufferData(GL_ARRAY_BUFFER,
                 m_attr_edges.get_bucket_size(Visibility::Visible) *
                   sizeof(GraphEdgeAttribute),
                 m_attr_edges.get_start(Visibility::Visible), GL_DYNAMIC_DRAW);
    m_vertex_buffer.unbind();

    int id_attr = 0;
    m_vertex_buffer.vertex_bind(id_attr, 0, sizeof(GraphEdgeAttribute));
    int id_pos = 1;
    m_pos_buffer.vertex_bind(id_pos, 0, sizeof(m_vertex_poly[0]));

    glVertexAttribFormat(line_attr_sq_pos.id(), 2, GL_FLOAT, GL_FALSE,
                         offsetof(RawVertex2d, pos));
    glVertexAttribBinding(line_attr_sq_pos.id(), id_pos);

    glVertexAttribFormat(line_attr_pos1.id(), 4, GL_FLOAT, GL_FALSE,
                         offsetof(GraphEdgeAttribute, pos1));
    GL_CHECK(glVertexAttribBinding(line_attr_pos1.id(), id_attr));

    glVertexAttribFormat(line_attr_pos2.id(), 4, GL_FLOAT, GL_FALSE,
                         offsetof(GraphEdgeAttribute, pos2));
    GL_CHECK(glVertexAttribBinding(line_attr_pos2.id(), id_attr));

    glVertexAttribFormat(line_attr_color.id(), 3, GL_FLOAT, GL_FALSE,
                         offsetof(GraphEdgeAttribute, col));
    glVertexAttribBinding(line_attr_color.id(), id_attr);

    glVertexBindingDivisor(id_pos, 0);
    glVertexBindingDivisor(id_attr, 1);

    glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, m_index.size());
    glVertexBindingDivisor(id_attr, 0);

    line_attr_pos1.unbind();
    line_attr_pos2.unbind();
    line_attr_color.unbind();
    line_attr_sq_pos.unbind();

    m_pos_buffer.vertex_unbind(id_attr);
    m_vertex_buffer.vertex_unbind(id_pos);

    line_unif_view.unbind();
    line_unif_cam_up.unbind();
    line_unif_cam_left.unbind();
  }
}

template <class Vertex, class Edge>
void DisplayGraph<Vertex, Edge>::update_edge(EdgeId id) {
  Edge &e = this->edge(id);
  GraphEdgeAttribute &attr_e = edge_attr(e.attr_id);
  attr_e.pos1 = vertex_attr(this->vertex(e.start).attr_id).pos;
  attr_e.pos2 = vertex_attr(this->vertex(e.end).attr_id).pos;
}

template <class Vertex, class Edge>
void DisplayGraph<Vertex, Edge>::set_vertex_visibility(VertexId id,
                                                       Visibility visibility) {
  m_attr_vertices.set_bucket(this->vertex(id).attr_id, visibility);
}

template <class Vertex, class Edge>
Vertex &DisplayGraph<Vertex, Edge>::create_vertex() {
  Vertex &v = opa::math::game::Graph<Vertex, Edge>::create_vertex();
  v.attr_id = m_attr_vertices.get_new(Visibility::Visible);
  return v;
}

template <class Vertex, class Edge>
Edge &DisplayGraph<Vertex, Edge>::add_edge(VertexId a, VertexId b) {
  Edge &edge = opa::math::game::Graph<Vertex, Edge>::add_edge(a, b);
  edge.attr_id = m_attr_edges.get_new(Visibility::Visible);
  return edge;
}

OPA_NAMESPACE_DECL3_END

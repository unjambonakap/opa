#include "Points.h"
#include "Drawer.h"
#include "Program.h"
#include "resources.h"

OPA_NAMESPACE_DECL2(opa, engine)

Program point_prog;
Program line_prog;

// Attribute point_unif_cam_left(Attribute::UNIFORM, "cam_left");
// Attribute point_unif_cam_up(Attribute::UNIFORM, "cam_up");
Attribute point_unif_view(Attribute::UNIFORM, "view");
Attribute point_unif_r_fact(Attribute::UNIFORM, "r_fact");
// Attribute point_attr_sq_pos(Attribute::ATTRIB, "sq_pos");
Attribute point_attr_pos(Attribute::ATTRIB, "pos");
Attribute point_attr_color(Attribute::ATTRIB, "color");

Program *get_program_point() {
  if (!point_prog.is_init()) {
    point_prog.init();
    point_prog.add_shader(GL_FRAGMENT_SHADER, point_fs);
    point_prog.add_shader(GL_GEOMETRY_SHADER, point_gs);
    point_prog.add_shader(GL_VERTEX_SHADER, point_vs);

    point_prog.register_attribute(&point_attr_pos);
    point_prog.register_attribute(&point_attr_color);
    // point_prog.register_attribute(&point_attr_sq_pos);
    // point_prog.register_attribute(&point_unif_cam_left);
    // point_prog.register_attribute(&point_unif_cam_up);
    point_prog.register_attribute(&point_unif_view);
    point_prog.register_attribute(&point_unif_r_fact);
    point_prog.setup();
  }
  return &point_prog;
}

PointCloud::PointCloud() {
  m_program = 0;
  m_n = 0;
}

void PointCloud::init(int n) {
  Buffer::init(GL_ARRAY_BUFFER);
  m_data.resize(n);
  m_n = n;

  glGenVertexArrays(1, &vao);
  m_program = get_program_point();
}

void PointCloud::fini() {
  glDeleteVertexArrays(1, &vao);
  Buffer::fini();
}

void PointCloud::refresh() {
  bind();
  glBufferData(GL_ARRAY_BUFFER, sizeof(PointData) * m_n, m_data.data(),
               GL_STATIC_DRAW);
  unbind();
}

void PointCloud::draw(const DrawData &draw_data) {
  check_init();
  glBindVertexArray(vao);

  m_program->bind();
  ScopedAttribute a1(point_attr_pos);
  ScopedAttribute a2(point_attr_color);
  // refresh();

  bind();

  GL_CHECK(glVertexAttribPointer(point_attr_pos.id(), 3, GL_FLOAT, GL_FALSE,
                                 sizeof(PointData),
                                 (void *)offsetof(PointData, pos)));

  GL_CHECK(glVertexAttribPointer(point_attr_color.id(), 2, GL_FLOAT, GL_FALSE,
                                 sizeof(PointData),
                                 (void *)offsetof(PointData, col)));
  unbind();
  GL_CHECK(;);

  GL_CHECK(glUniformMatrix4fv(point_unif_view.id(), 1, GL_FALSE,
                              glm::value_ptr(draw_data.proj())));
  GL_CHECK(glUniform1f(point_unif_r_fact.id(), m_r_fact));
  GL_CHECK(glDrawArrays(GL_POINTS, 0, m_n));
  glBindVertexArray(0);
}

OPA_NAMESPACE_DECL2_END

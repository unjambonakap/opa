#include "Strip.h"
#include "Drawer.h"
#include "Program.h"
#include "resources.h"

using namespace std;

OPA_NAMESPACE_DECL2(opa, engine)

Attribute s1_attrib_pos(Attribute::ATTRIB, "pos");
Attribute s1_attrib_texpos(Attribute::ATTRIB, "texpos");
Attribute s1_attrib_view(Attribute::UNIFORM, "view");
Attribute s1_attrib_texture(Attribute::UNIFORM, "texture1");

Program s1_prog;

Program *get_program_s1() {
  if (!s1_prog.is_init()) {
    s1_prog.init();
    s1_prog.add_shader(GL_FRAGMENT_SHADER, s1_fs);
    s1_prog.add_shader(GL_VERTEX_SHADER, s1_vs);

    s1_prog.register_attribute(&s1_attrib_pos);
    s1_prog.register_attribute(&s1_attrib_texpos);
    s1_prog.register_attribute(&s1_attrib_view);
    s1_prog.register_attribute(&s1_attrib_texture);
    s1_prog.setup();
  }
  return &s1_prog;
}

Attribute ground_attrib_pos(Attribute::ATTRIB, "pos");
Attribute ground_attrib_x_color(Attribute::ATTRIB, "x_color");
Attribute ground_attrib_view(Attribute::UNIFORM, "view");

Program ground_prog;

Program *get_program_ground() {
  if (!ground_prog.is_init()) {
    ground_prog.init();
    ground_prog.add_shader(GL_FRAGMENT_SHADER, ground_fs);
    ground_prog.add_shader(GL_VERTEX_SHADER, ground_vs);

    ground_prog.register_attribute(&ground_attrib_pos);
    ground_prog.register_attribute(&ground_attrib_x_color);
    ground_prog.register_attribute(&ground_attrib_view);
    ground_prog.setup();
  }
  return &ground_prog;
}

TexturedStrip::TexturedStrip() {
  m_texture = 0;
  m_program = 0;
  m_data = 0;
  m_n = 0;
}

void TriangleContainer::draw(const DrawData &draw_data) {
  for (auto &x : triangles().used()) triangles().get(x).draw(draw_data);
  for (auto &x : strips().used()) strips().get(x).draw(draw_data);
}
void TexturedStrip::fill_tr(int id, Pos const *&p1, Pos const *&p2,
                            Pos const *&p3) const {
  p1 = &m_data[id].pos;
  p2 = &m_data[id + 1].pos;
  p3 = &m_data[id + 2].pos;
}

void TexturedStrip::init(u32 n, TexturePtr texture) {
  Buffer::init(GL_ARRAY_BUFFER);
  m_data = new TextureData[n + 2];
  m_texture = texture;
  m_n = n;

  glGenVertexArrays(1, &vao);
  m_program = get_program_s1();
}

void TexturedStrip::fini() {
  delete m_data;
  m_data = 0;
  glDeleteVertexArrays(1, &vao);
  Buffer::fini();
}

void TexturedStrip::refresh() {
  bind();
  glBufferData(GL_ARRAY_BUFFER, sizeof(TextureData) * (m_n + 2), m_data,
               GL_STATIC_DRAW);
  unbind();
}

void TexturedStrip::draw(const DrawData &draw_data) {
  check_init();
  glBindVertexArray(vao);

  m_program->bind();
  ScopedAttribute a1(s1_attrib_pos);
  ScopedAttribute a2(s1_attrib_texpos);

  refresh();

  vertex_bind(s1_attrib_pos.id(), 0, sizeof(TextureData));
  vertex_bind(s1_attrib_texpos.id(), 0, sizeof(TextureData));

  glVertexAttribFormat(s1_attrib_pos.id(), 3, GL_FLOAT, GL_FALSE,
                       offsetof(TextureData, pos));

  glVertexAttribFormat(s1_attrib_texpos.id(), 2, GL_FLOAT, GL_FALSE,
                       offsetof(TextureData, texpos));

  glVertexAttribBinding(s1_attrib_pos.id(), s1_attrib_pos.id());
  glVertexAttribBinding(s1_attrib_texpos.id(), s1_attrib_texpos.id());

  glActiveTexture(GL_TEXTURE0);
  m_texture->bind();
  GL_CHECK(glUniform1i(s1_attrib_texture.id(), 0));

  GL_CHECK(glUniformMatrix4fv(s1_attrib_view.id(), 1, GL_FALSE,
                              glm::value_ptr(draw_data.proj())));

  GL_CHECK(glDrawArrays(GL_TRIANGLE_STRIP, 0, m_n + 2));
  glBindVertexArray(0);
}

// TexturedTriangle
void TexturedTriangle::fill_tr(Pos const *&p1, Pos const *&p2,
                               Pos const *&p3) const {
  p1 = &m_data[0].pos;
  p2 = &m_data[1].pos;
  p3 = &m_data[2].pos;
}

TexturedTriangle::TexturedTriangle() {
  m_texture = 0;
  m_program = 0;
}

void TexturedTriangle::init(TexturePtr texture) {
  Buffer::init(GL_ARRAY_BUFFER);
  m_program = get_program_s1();
  m_texture = texture;
  glGenVertexArrays(1, &vao);
}

void TexturedTriangle::refresh() {
  bind();
  glBufferData(GL_ARRAY_BUFFER, sizeof(TextureData) * 3, m_data,
               GL_STATIC_DRAW);
  unbind();
}

void TexturedTriangle::fini() {
  Buffer::fini();
  glDeleteVertexArrays(1, &vao);
}

void TexturedTriangle::draw(const DrawData &draw_data) {
  check_init();

  glBindVertexArray(vao);
  m_program->bind();
  s1_attrib_pos.bind();
  s1_attrib_texpos.bind();

  //REP(i,3) OPA_DISP0(ApplyMat(draw_data.proj(), m_data[i].pos));

  bind();
  GL_CHECK(glUniformMatrix4fv(s1_attrib_view.id(), 1, GL_FALSE,
                              glm::value_ptr(draw_data.proj())));

  GL_CHECK(glVertexAttribPointer(s1_attrib_pos.id(), 3, GL_FLOAT, GL_FALSE,
                                 sizeof(TextureData),
                                 (void *)offsetof(TextureData, pos)));

  GL_CHECK(glVertexAttribPointer(s1_attrib_texpos.id(), 2, GL_FLOAT, GL_FALSE,
                                 sizeof(TextureData),
                                 (void *)offsetof(TextureData, texpos)));
  unbind();
  GL_CHECK(;);

  glActiveTexture(GL_TEXTURE0);
  GL_CHECK(;);
  m_texture->bind();
  GL_CHECK(;);
  GL_CHECK(glUniform1i(s1_attrib_texture.id(), 0));
  GL_CHECK(glUniformMatrix4fv(s1_attrib_view.id(), 1, GL_FALSE,
                              glm::value_ptr(draw_data.proj())));
  GL_CHECK(glDrawArrays(GL_TRIANGLES, 0, 3));
  s1_attrib_pos.unbind();
  s1_attrib_texpos.unbind();
  glBindVertexArray(0);
}

// XColorTriangle

XColorTriangle::XColorTriangle() { m_program = 0; }

void XColorTriangle::init() {
  Buffer::init(GL_ARRAY_BUFFER);
  t_b2.init(GL_ARRAY_BUFFER);
  m_program = get_program_ground();
}
void XColorTriangle::fini() { Buffer::fini(); }

void XColorTriangle::refresh() {
  if (0) {
    bind();
    glBufferData(GL_ARRAY_BUFFER, sizeof(m_data), m_data, GL_STATIC_DRAW);

    m_program->bind();
    ground_attrib_pos.bind();
    ground_attrib_x_color.bind();

    glVertexAttribPointer(ground_attrib_pos.id(), 3, GL_FLOAT, GL_FALSE,
                          sizeof(XColorData),
                          (void *)offsetof(XColorData, pos));
    glVertexAttribPointer(ground_attrib_x_color.id(), 3, GL_FLOAT, GL_FALSE,
                          sizeof(XColorData),
                          (void *)offsetof(XColorData, color));

    ground_attrib_pos.unbind();
    ground_attrib_x_color.unbind();
  } else {
    REP (i, 3)
      t_pos[i] = m_data[i].pos;
    REP (i, 3)
      t_col[i] = m_data[i].color;
    m_program->bind();
    ground_attrib_pos.bind();
    ground_attrib_x_color.bind();

    bind();
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 9, t_pos, GL_STATIC_DRAW);
    glVertexAttribPointer(ground_attrib_pos.id(), 3, GL_FLOAT, GL_FALSE,
                          3 * sizeof(float), 0);
    unbind();
    t_b2.bind();
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 9, t_col, GL_STATIC_DRAW);
    glVertexAttribPointer(ground_attrib_x_color.id(), 3, GL_FLOAT, GL_FALSE,
                          3 * sizeof(float), 0);
    t_b2.unbind();
    ground_attrib_pos.unbind();
    ground_attrib_x_color.unbind();
  }
}

void XColorTriangle::draw(const DrawData &draw_data) {

  m_program->bind();
  ground_attrib_pos.bind();
  ground_attrib_x_color.bind();
  auto proj = draw_data.proj();

  GL_CHECK(glUniformMatrix4fv(ground_attrib_view.id(), 1, GL_FALSE,
                              glm::value_ptr(proj)));

  bind();
  glVertexAttribPointer(ground_attrib_pos.id(), 3, GL_FLOAT, GL_FALSE,
                        3 * sizeof(float), 0);
  unbind();
  t_b2.bind();
  glVertexAttribPointer(ground_attrib_x_color.id(), 3, GL_FLOAT, GL_FALSE,
                        3 * sizeof(float), 0);
  t_b2.unbind();
  GL_CHECK(glDrawArrays(GL_TRIANGLES, 0, 3));
  ground_attrib_pos.unbind();
  ground_attrib_x_color.unbind();
}

OPA_NAMESPACE_DECL2_END

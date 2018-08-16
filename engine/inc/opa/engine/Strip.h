#pragma once
#include <opa/engine/conf.h>

OPA_NAMESPACE_DECL2(opa, engine)
class DrawData;
class Program;

Program *get_program_ground();
Program *get_program_s1();

struct TextureData {
  Pos pos;
  glm::vec2 texpos;
};

class TexturedStrip : public Buffer, public opa::utils::IdObj {
public:
  TexturedStrip();

  virtual void init(u32 n, TexturePtr texture);
  virtual void fini();
  void refresh();

  void draw(const DrawData &draw_data);

  OPA_ACCESSOR_PTR(TextureData, m_data, data)
  OPA_ACCESSOR_R(u32, m_n, n)
  void fill_tr(int id, Pos const *&p1, Pos const *&p2, Pos const *&p3) const;

private:
  virtual void init() {}

  u32 m_n;
  TextureData *m_data; // owned
  TexturePtr m_texture;
  Program *m_program;
  GLuint vao;
};
OPA_DECL_SPTR(TexturedStrip, TexturedStripPtr)

class TexturedTriangle : public Buffer, public opa::utils::IdObj {
public:
  TexturedTriangle();
  virtual void init(TexturePtr texture);
  virtual void fini();
  void refresh();

  void draw(const DrawData &draw_data);

  OPA_ACCESSOR_PTR(TextureData, m_data, data)
  void fill_tr(Pos const *&p1, Pos const *&p2, Pos const *&p3) const;

private:
  virtual void init() {}

  TextureData m_data[3];
  TexturePtr m_texture;
  Program *m_program;
  GLuint vao;
};
OPA_DECL_SPTR(TexturedTriangle, TexturedTrianglePtr)

struct XColorData {
  glm::vec3 pos;
  glm::vec3 color;
};

class XColorTriangle : public Buffer {
public:
  XColorTriangle();
  virtual void init();
  virtual void fini();
  void refresh();

  void draw(const DrawData &draw_data);

  OPA_ACCESSOR_PTR(XColorData, m_data, data)
private:
  XColorData m_data[3];
  Program *m_program;
  Buffer t_b2;
  glm::vec3 t_pos[3];
  glm::vec3 t_col[3];
};

class TriangleContainer {
public:
  void draw(const DrawData &draw_data);
  OPA_ACCESSOR(opa::utils::ObjectPool<TexturedTriangle>, m_triangles,
               triangles);
  OPA_ACCESSOR(opa::utils::ObjectPool<TexturedStrip>, m_strips, strips);

  opa::math::game::FaceCollection to_tr_collection() const {
    opa::math::game::FaceCollection res;
    for (auto &x : triangles().used()) {
      auto &tr = triangles().get(x);
      res.push(opa::PointVec { tr.data()[0].pos, tr.data()[1].pos, tr.data()[2].pos });
    }
    for (auto &x : strips().used()) {
      const auto &strip = strips().get(x);
      REP (i, strip.n()) {
        res.push(opa::PointVec{ strip.data()[i].pos, strip.data()[i + 1].pos,
                   strip.data()[i + 2].pos });
      }
    }
    return res;
  }

private:
  opa::utils::ObjectPool<TexturedTriangle> m_triangles;
  opa::utils::ObjectPool<TexturedStrip> m_strips;
};
OPA_DECL_SPTR(TriangleContainer, TriangleContainerPtr);

OPA_NAMESPACE_DECL2_END

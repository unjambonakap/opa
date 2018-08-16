#pragma once
#include <opa/engine/conf.h>

OPA_NAMESPACE_DECL2(opa, engine)

class Program;

struct Color {
  int col[3];

  OPA_ACCESSOR(int, col[0], r);
  OPA_ACCESSOR(int, col[1], g);
  OPA_ACCESSOR(int, col[2], b);
  Color() { r() = b() = g() = 0; }
  Color(int r, int g, int b) {
    this->r() = r;
    this->g() = g;
    this->b() = b;
  }
  Col to_col() const { return Col(r(), g(), b()) / 255; }

  static Color from_col(const Col &c) {
    Col tmp = c * 255.;
    Color res(tmp.x, tmp.y, tmp.z);
    return res;
  }

  bool operator<(const Color &x) const {
    REP (i, 3)
      if (col[i] != x.col[i]) return col[i] < x.col[i];
    return 0;
  }
};
static std::ostream &operator<<(std::ostream &os, const Color &col) {
  return os << col.to_col();
}

static void CheckOpenGLError(const char *stmt, const char *fname, int line) {
  GLenum err = glGetError();
  if (err != GL_NO_ERROR) {
    printf("OpenGL error %08x, at %s:%i - for %s\n", err, fname, line, stmt);
    OPA_PRINT_STACKTRACE();
    abort();
  }
}

struct Ray {
  Pos pos;
  Dir dir;
};

#if OPA_DEBUG
#define GL_CHECK(stmt)                                                         \
  do {                                                                         \
    stmt;                                                                      \
    CheckOpenGLError(#stmt, __FILE__, __LINE__);                               \
  } while (0)
#else
#define GL_CHECK(stmt) stmt
#endif

class Time {
public:
  Time() { m_us = 0; }
  Time(u64 v) : m_us(v) {}
  OPA_ACCESSOR(u64, m_us, us)
private:
  u64 m_us;
};

class Buffer : public opa::utils::Initable {
public:
  virtual void init(GLenum mode) {
    m_mode = mode;
    opa::utils::Initable::init();
    glGenBuffers(1, &m_bid);
  }

  virtual void fini() override {
    glDeleteBuffers(1, &m_bid);
    opa::utils::Initable::fini();
  }

  void bind() {
    check_init();
    glBindBuffer(m_mode, m_bid);
  }

  void vertex_bind(GLuint binding_id, u32 offset, u32 stride) {
    glBindVertexBuffer(binding_id, m_bid, offset, stride);
  }

  void vertex_unbind(GLuint id) { glBindVertexBuffer(id, 0, 0, 0); }

  void unbind() {
    check_init();
    glBindBuffer(m_mode, 0);
  }
  OPA_ACCESSOR_R(GLuint, m_bid, bid)

private:
  virtual void init() override {}
  GLenum m_mode;
  GLuint m_bid;
};

class ScopedBuffer {
public:
  ScopedBuffer(Buffer &buf) : m_buf(buf) { m_buf.bind(); }
  ~ScopedBuffer() { m_buf.unbind(); }

private:
  Buffer &m_buf;
};

class Attribute : public opa::utils::Initable {
public:
  enum AttributeType {
    ATTRIB,
    UNIFORM,
  };

  Attribute(AttributeType type, const std::string &name) {
    m_type = type;
    m_name = name;
  }

  virtual void init(Program *program);

  void bind();
  void unbind();

  OPA_ACCESSOR_R(GLint, m_id, id)

private:
  virtual void init() override {}
  AttributeType m_type;
  std::string m_name;
  GLint m_id;
};

class ScopedAttribute {
public:
  ScopedAttribute(Attribute &attr) : m_attr(attr) { m_attr.bind(); }
  ~ScopedAttribute() { m_attr.unbind(); }

private:
  Attribute &m_attr;
};

class Image;
OPA_DECL_SPTR(Image, ImagePtr);
// todo: template
class Image : public opa::utils::Initable {
public:
  enum ImageFormat : int {
    RGB,
    Grey,
  };
  virtual void init(int w, int h, ImageFormat format);
  virtual void fini() override;

  Color get_col(int x, int y) const;
  TexUV get_uv(int x, int y) const;
  void set_col(int x, int y, Color col);

  OPA_ACCESSOR_PTR(u8, m_buf.get(), buf);
  OPA_ACCESSOR_R(int, m_w, w);
  OPA_ACCESSOR_R(int, m_h, h);
  OPA_ACCESSOR_R(ImageFormat, m_format, format);

  ImagePtr extract(int xi, int yi, int w, int h) const;
  void dump(std::ostream &os) const;
  void load(std::istream &is);
  bool load_png(const std::string &fname);
  GLenum get_gl_enum() const;
  int bytesize() const { return m_stride * m_h; }

private:
  virtual void init() override {}

  u8 *get_ptr(int x, int y);
  const u8 *get_ptr(int x, int y) const;
  ImageFormat m_format;
  int m_pix_size;
  int m_stride;
  std::shared_ptr<u8> m_buf;
  int m_w;
  int m_h;
};
OPA_DECL_SPTR(Image, ImagePtr);

std::ostream &operator<<(std::ostream &os, const Image &img);

class Texture : public opa::utils::Initable {
public:
  virtual void init(const Image &image);
  virtual void fini() override;
  OPA_ACCESSOR_R(GLuint, m_id, id);

  void bind();
  void unbind();

protected:
  GLuint m_id;
};
OPA_DECL_SPTR(Texture, TexturePtr)

class ImageResampler {

public:
  struct ImageEntry {
    Image *img;
    int x, y;
    ImageEntry(Image *img, int x, int y) {
      this->img = img;
      this->x = x;
      this->y = y;
    }
  };
  struct Params {
    std::vector<ImageEntry> entries;
    Params() {}
    Params(const std::vector<ImageEntry> &tb) { entries = tb; }
  };
  ImagePtr resample(const Params &params, int dw, int dh);

private:
  Col get(int x, int y);
  int mx, my;
  const Params *m_params;
};

typedef u32 TextureId;
class TexturePool {
public:
  TextureId add(TexturePtr texture);
  void remove(TextureId id);
  TexturePtr get(TextureId id);

private:
  TextureId next_id();
  std::map<TextureId, TexturePtr> m_pool;
  TextureId m_cur;
};

class ColorPalette {
public:
  TexturePtr get_rand() const;
  void resize_rand(int num);
  TexturePtr get_color(const Color &color);

private:
  TexturePtr gen_color_texture(const Color &color) const;
  std::vector<TexturePtr> m_rand_pool;
  std::map<Color, TexturePtr> m_color_cache;
};

class MultiColorTexture : public Texture {
public:
  virtual void init(int ncolors);
  TexUV get(int col);
  virtual void fini() override;
  int size() const { return m_ncolors; }

private:
  virtual void init(const Image &image) override {}
  int m_ncolors;
  Image m_img;
};

enum VBOType : int { VBO_TriangleStrip, VBO_TriangleFan, VBO_Triangle };

class VBO : public Buffer {
public:
  virtual void init(VBOType type);
  void draw();
  std::vector<int> &get();
  void refresh();

private:
  virtual void init(GLenum mode) override {}
  std::vector<int> m_index;
  VBOType m_type;
  GLenum m_native_type;
};

class TrVBO : public VBO {
public:
  virtual void init() override;
  void add_tr(int a, int b, int c);

private:
  virtual void init(VBOType type) override{};
};

class FanVBO : public VBO {
public:
  virtual void init() override;
  void set_center(int a);
  void add_tr(int a, int b);

private:
  virtual void init(VBOType type) override{};
};

OPA_NAMESPACE_DECL2_END

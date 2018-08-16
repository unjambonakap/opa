#include "common.h"

#include "Program.h"
#include "procgen/Texture.h"
#include <lodepng.h>

using namespace opa::math::game;
using namespace std;

OPA_NAMESPACE_DECL2(opa, engine)
void Attribute::init(Program *program) {
  opa::utils::Initable::init();

  switch (m_type) {
  case ATTRIB:
    m_id = glGetAttribLocation(program->id(), m_name.c_str());
    break;

  case UNIFORM:
    m_id = glGetUniformLocation(program->id(), m_name.c_str());
    break;

  default:
    OPA_CHECK(false, "bad attrib");
    break;
  }
  OPA_CHECK(m_id != -1, "bad shader");
}

void Attribute::bind() { GL_CHECK(glEnableVertexAttribArray(m_id)); }
void Attribute::unbind() { glDisableVertexAttribArray(m_id); }

void Texture::init(const Image &image) {
  opa::utils::Initable::init();

  glGenTextures(1, &m_id);
  bind();
  // glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  GLenum format = image.get_gl_enum();

  GL_CHECK(glTexImage2D(GL_TEXTURE_2D, 0, format, image.w(), image.h(), 0,
                        format, GL_UNSIGNED_BYTE, image.buf()));
  unbind();
}
void Texture::bind() {
  check_init();
  glBindTexture(GL_TEXTURE_2D, m_id);
}

void Texture::unbind() {
  glBindTexture(GL_TEXTURE_2D, 0);
}

void Texture::fini() {
  if (is_init()) {
    glDeleteTextures(1, &m_id);
  }
  opa::utils::Initable::fini();
}
GLenum Image::get_gl_enum() const {
  switch (m_format) {
    OPA_CASE(ImageFormat::RGB, return GL_RGB;)
    OPA_CASE(ImageFormat::Grey, return GL_LUMINANCE;)
  default:
    OPA_ABORT0(true);
    break;
  }
}

void Image::init(int w, int h, ImageFormat format) {
  opa::utils::Initable::init();
  m_w = w;
  m_h = h;
  m_format = format;

  switch (format) {
    OPA_CASE(ImageFormat::RGB, m_pix_size = 3;)
    OPA_CASE(ImageFormat::Grey, m_pix_size = 1;)
  default:
    OPA_ABORT0(true);
    break;
  };
  m_stride = m_pix_size * m_w;

  m_buf.reset(new u8[bytesize()]);
}

ImagePtr Image::extract(int xi, int yi, int w, int h) const {
  ImagePtr res(new Image);
  res->init(w, h, m_format);

  REP (i, h) {
    u8 *dest = res->get_ptr(0, i);
    const u8 *src = get_ptr(xi, yi + i);
    memcpy(dest, src, w * m_pix_size);
  }
  return res;
}

void Image::fini() { opa::utils::Initable::fini(); }

bool Image::load_png(const std::string &fname) {

  std::vector<unsigned char> image; // the raw pixels
  unsigned width, height;
  unsigned error =
    lodepng::decode(image, width, height, fname, LodePNGColorType::LCT_RGB);

  if (error) {
    OPA_DISP("Got err >> ", lodepng_error_text(error));
    return false;
  }

  OPA_CHECK_EQ0(width, m_w);
  OPA_CHECK_EQ0(height, m_h);

  memcpy(m_buf.get(), image.data(), bytesize());
  return true;
}

const u8 *Image::get_ptr(int x, int y) const {
  return m_buf.get() + y * m_stride + x * m_pix_size;
}

u8 *Image::get_ptr(int x, int y) {
  return m_buf.get() + y * m_stride + x * m_pix_size;
}

void Image::dump(std::ostream &os) const {
  check_init();
  os.write((const char *)get_ptr(0, 0), m_stride * m_h);
}
void Image::load(std::istream &is) {
  check_init();
  is.read((char *)get_ptr(0, 0), m_stride * m_h);
}

Color Image::get_col(int x, int y) const {
  const u8 *pos = get_ptr(x, y);
  switch (m_format) {
    OPA_CASE(ImageFormat::RGB, return Color(pos[0], pos[1], pos[2]);)
    OPA_CASE(ImageFormat::Grey, return Color(pos[0], 0, 0);)
  }
  OPA_ABORT0(true);
}

TexUV Image::get_uv(int x, int y) const {
  return TexUV((x + 0.5) / m_w, (m_h - 0.5 - y) / m_h);
}
void Image::set_col(int x, int y, Color col) {
  u8 *pos = get_ptr(x, y);
  switch (m_format) {
    OPA_CASE(ImageFormat::RGB, pos[0] = col.r(); pos[1] = col.g();
             pos[2] = col.b();)
    OPA_CASE(ImageFormat::Grey, pos[0] = col.r();)
  }
}

std::ostream &operator<<(std::ostream &os, const Image &img) {
  REP (j, img.h()) {
    REP (i, img.w()) {
      if (img.format() == Image::ImageFormat::Grey)
        os << img.get_col(i, j).r() << " ";
      else
        os << img.get_col(i, j) << " ";
    }
    os << "\n";
  }
  return os;
}

ImagePtr ImageResampler::resample(const Params &params, int dw, int dh) {
  m_params = &params;
  Image *res = new Image();
  res->init(dw, dh, params.entries[0].img->format());

  int ow = 0;
  int oh = 0;
  for (auto &e : m_params->entries) {
    ow = max(ow, e.img->w() + e.x);
    oh = max(oh, e.img->h() + e.y);
  }
  mx = ow;
  my = oh;

  float cw = 1. * ow / dw;
  float ch = 1. * oh / dh;
  Pos2 scale = Pos2(cw, ch);

  REP (i, dw)
    REP (j, dh) {
      // todo: halfpoints
      float sum = 0;
      Col col;
      REP (k, 4) {
        Pos2 p = Pos2(i + k % 2, j + k / 2) * scale;
        float coeff = max(0.f, 1 - glm::length2(Pos2(int(p.x), int(p.y)) - p));
        // OPA_DISP0(k, p, coeff, get(p.x, p.y));
        sum += coeff;
        col += coeff * get(p.x, p.y);
      }
      col /= sum;
      res->set_col(i, j, Color::from_col(col));
    }

  return ImagePtr(res);
}

Col ImageResampler::get(int x, int y) {
  x = min(x, mx - 1);
  y = min(y, my - 1);
  for (auto &e : m_params->entries) {
    if (x < e.x || e.x + e.img->w() <= x) continue;
    if (y < e.y || e.y + e.img->h() <= y) continue;
    return e.img->get_col(x - e.x, y - e.y).to_col();
  }
  // OPA_CHECK0(false);
  return Col();
}

void MultiColorTexture::init(int ncolors) {
  m_ncolors = ncolors;
  m_img.init(ncolors, 1, Image::ImageFormat::RGB);
  REP (i, ncolors)
    m_img.set_col(i, 0, Color::from_col(vec_rand_uni()));

  Texture::init(m_img);
}

TexUV MultiColorTexture::get(int col) { return m_img.get_uv(col, 0); }

void MultiColorTexture::fini() {
  Texture::fini();
  m_img.fini();
}

TextureId TexturePool::next_id() { return m_cur++; }
TextureId TexturePool::add(TexturePtr texture) {
  TextureId id = next_id();
  m_pool[id] = texture;
  return id;
}
void TexturePool::remove(TextureId id) { m_pool.erase(id); }

TexturePtr TexturePool::get(TextureId id) {
  OPA_ASSERT0(m_pool.count(id));
  return m_pool[id];
}

TexturePtr ColorPalette::get_rand() const {
  OPA_CHECK0(m_rand_pool.size() > 1);
  int id = opa::math::common::rng() % m_rand_pool.size();
  return m_rand_pool[id];
}

void ColorPalette::resize_rand(int num) {
  if (num < m_rand_pool.size()) m_rand_pool.resize(num);

  while (m_rand_pool.size() < num) {
    Color col = Color::from_col(vec_rand_uni());
    m_rand_pool.push_back(gen_color_texture(col));
  }
}
TexturePtr ColorPalette::get_color(const Color &color) {
  if (!m_color_cache.count(color))
    m_color_cache[color] = gen_color_texture(color);
  return m_color_cache[color];
}

TexturePtr ColorPalette::gen_color_texture(const Color &color) const {
  auto texture = new UniformTexture();
  texture->init(UniformTexture::UniformTextureParameters(color));
  return TexturePtr(texture);
}

void VBO::init(VBOType type) {
  Buffer::init(GL_ELEMENT_ARRAY_BUFFER);
  m_type = type;
  switch (m_type) {
    OPA_CASE(VBOType::VBO_Triangle, m_native_type = GL_TRIANGLES;);
    OPA_CASE(VBOType::VBO_TriangleStrip, m_native_type = GL_TRIANGLE_STRIP;);
    OPA_CASE(VBOType::VBO_TriangleFan, m_native_type = GL_TRIANGLE_FAN;);
  default:
    OPA_CHECK0(false);
  }
}

void VBO::draw() {
  ScopedBuffer aa(*this);
  glDrawElements(m_native_type, m_index.size(), GL_UNSIGNED_INT, 0);
}

void VBO::refresh() {
  ScopedBuffer aa(*this);
  GL_CHECK(glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                        sizeof(m_index[0]) * m_index.size(), m_index.data(),
                        GL_STATIC_DRAW));
}

std::vector<int> &VBO::get() { return m_index; }

void TrVBO::init() { VBO::init(VBOType::VBO_Triangle); }

void TrVBO::add_tr(int a, int b, int c) {
  get().pb(a);
  get().pb(b);
  get().pb(c);
}

void FanVBO::init() {
  VBO::init(VBOType::VBO_TriangleFan);
  get().resize(1);
}
void FanVBO::set_center(int a) { get()[0] = a; }
void FanVBO::add_tr(int a, int b) {
  get().pb(a);
  get().pb(b);
}

OPA_NAMESPACE_DECL2_END

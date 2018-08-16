#include "RenderStreamer.h"
#include "Renderer.h"
#include <EGL/egl.h>
#include <GL/glext.h>

DEFINE_bool(stream_rendering, false, "");
DEFINE_string(stream_fname, "", "");
DEFINE_int32(stream_nqueue, 3, "");

using namespace opa::math::game;
OPA_NAMESPACE_DECL2(opa, engine)

void RenderStreamer::init(const Renderer *renderer) {
  opa::utils::Initable::init();

  this->renderer = renderer;
  OPA_CHECK0(renderer->is_init());
  if (!FLAGS_stream_rendering) return;

  if (!FLAGS_stream_fname.empty()) {
    ofs = std::ofstream(FLAGS_stream_fname, std::ofstream::binary);
  }
  img.init(renderer->w(), renderer->h(), Image::ImageFormat::RGB);

  const int ndelay = FLAGS_stream_nqueue;
  this->m_container.init(ndelay);

  glGenBuffers(m_container.n, m_container.objs.data());
  REP (i, m_container.n) {
    glBindBuffer(GL_PIXEL_PACK_BUFFER, m_container.objs[i]);
    glBufferData(GL_PIXEL_PACK_BUFFER, img.w() * img.h() * 4, 0,
                 GL_STATIC_READ);
  }
  glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
  GL_CHECK(;);
}

void RenderStreamer::fini() {
  opa::utils::Initable::fini();
  glDeleteBuffers(m_container.n, m_container.objs.data());
}

void RenderStreamer::update() {
  this->check_init();
  if (!FLAGS_stream_rendering) return;
  GL_CHECK(;);

  if (0) {
    // use get fbo
    GL_CHECK(;);
    glReadPixels(0, 0, img.w(), img.h(), img.get_gl_enum(), GL_UNSIGNED_BYTE,
                 img.buf());
  } else {
    glBindBuffer(GL_PIXEL_PACK_BUFFER, m_container.next_write());

    // glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glReadPixels(0, 0, img.w(), img.h(), GL_RGBA, GL_UNSIGNED_BYTE, 0);
    GL_CHECK(;);

    if (m_container.can_read()) {
      auto tmp = m_container.next_read();
      glBindBuffer(GL_PIXEL_PACK_BUFFER, tmp);
      GL_CHECK(;);
      GLubyte *src = (GLubyte *)glMapBuffer(GL_PIXEL_PACK_BUFFER, GL_READ_ONLY);
      if (src != nullptr) {
        OPA_CHECK0(src);
        glUnmapBuffer(GL_PIXEL_PACK_BUFFER); // release pointer to the mapped
                                             // buffer
      } else {
        OPA_TRACEM("Fuu src null");
      }
      if (!FLAGS_stream_fname.empty()) {
        ofs.write((const char *)src, img.w() * img.h() * 4);
      }
    }
    glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
    GL_CHECK(;);
  }
}

OPA_NAMESPACE_DECL2_END

#include "procgen/Texture.h"

OPA_NAMESPACE_DECL2(opa, engine)

UniformTexture::UniformTexture() { m_buf = 0; }

void UniformTexture::init(const UniformTextureParameters &params) {
  Image img;
  img.init(1, 1, Image::ImageFormat::RGB);
  img.set_col(0, 0, params.color);

  Texture::init(img);
}

void UniformTexture::fini() {
  delete m_buf;
  m_buf = 0;
  Texture::fini();
}

OPA_NAMESPACE_DECL2_END

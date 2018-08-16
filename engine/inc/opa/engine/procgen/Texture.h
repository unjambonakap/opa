#pragma once

#include <opa/engine/conf.h>

OPA_NAMESPACE_DECL2(opa, engine)

class UniformTexture : public Texture {
public:
  struct UniformTextureParameters {
    UniformTextureParameters(const Color &color) { this->color = color; }
    UniformTextureParameters() {}

    Color color;
  };

  UniformTexture();
  virtual void init(const UniformTextureParameters &params);
  virtual void fini() override;

protected:
  u8 *m_buf;
};
OPA_DECL_SPTR(UniformTexture, UniformTexturePtr)

// class GaussianTexture : public Texture {
//  public:
//    struct GaussianTextureParameters {
//        int w, h;
//        struct Entry {
//            Pos2 center;
//            glm::mat2 cvar;
//            Color A;
//            Color B;
//        };
//
//        std::vector<Entry> entries;
//    };
//
//
//
//
//    virtual void init(const GaussianTextureParameters &params);
//
//  private:
//    virtual void init(const UniformTextureParameters &params);
//};

OPA_NAMESPACE_DECL2_END

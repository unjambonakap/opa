#pragma once

#include <opa/engine/Intersection.h>
#include <opa/engine/Scene.h>
#include <opa/engine/Strip.h>
#include <opa/engine/conf.h>
#include <opa/engine/primitives/Base.h>

OPA_NAMESPACE_DECL2(opa, engine)

class Sphere : public PrimitiveSceneObj {
public:
  Sphere();
  virtual void init(Game *game, float radius);

  void build(u32 npoly);

  void build_strip(TexturedStrip *strip, float zl, float zh, int ntr);

  virtual opa::math::game::Mesh *
  get_mesh(double precision = 0.) const override {
    return m_mesh_cache.get_or_insert(
      precision, [&]() { return this->build_mesh(precision); });
  }

private:
  opa::math::game::Mesh *build_mesh(double precision = 0.) const;

  virtual void init(Game *game) override {}
  float m_r;
  TexturePtr m_texture; // in the meantime
  u8 *t_buf;

  utils::LRUCache<double, opa::math::game::Mesh> m_mesh_cache;
};
OPA_DECL_SPTR(Sphere, SpherePtr)

OPA_NAMESPACE_DECL2_END

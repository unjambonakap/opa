#pragma once

#include <opa/engine/Intersection.h>
#include <opa/engine/Scene.h>
#include <opa/engine/Strip.h>
#include <opa/engine/conf.h>

OPA_NAMESPACE_DECL2(opa, engine)

class PrimitiveSceneObj : public SceneObj {
public:
  virtual opa::math::game::Mesh *
  get_mesh(double precision = 0.) const override {
    if (!m_mesh_cache) {
      m_mesh_cache.reset(tr().to_tr_collection().to_mesh());
    }
    return m_mesh_cache.get();
  }

  OPA_ACCESSOR(TriangleContainer, m_tr, tr);

protected:
  mutable UPTR(opa::math::game::Mesh) m_mesh_cache;
  TriangleContainer m_tr;
};

class Square : public PrimitiveSceneObj {
public:
  Square();
  virtual void init(Game *game, float side, TexturePtr texture,
                    const std::vector<TexUV> &texpos);
};

class Cube : public PrimitiveSceneObj {
public:
  virtual void init(Game *game, float side, TexturePtr *textures);
};

class Icosahedron : public PrimitiveSceneObj {
public:
  virtual void init(Game *game, float side, TexturePtr *textures);
};

class TrianglesSceneObj : public PrimitiveSceneObj {
public:
  virtual void init(Game *game);
};


OPA_NAMESPACE_DECL2_END

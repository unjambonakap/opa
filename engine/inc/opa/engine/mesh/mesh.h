#pragma once

#include <opa/engine/common.h>
#include <opa/engine/conf.h>
#include <opa/engine/graph/graph.h>
#include <opa/math/game/mesh.h>
#include <opa/utils/DataStruct.h>

OPA_NAMESPACE_DECL2(opa, engine)

class DisplayMesh : public opa::math::game::Mesh {

public:
  virtual void init() override;
  void prepare_draw();

  OPA_ACCESSOR(TrVBO, m_vbo, vbo);
  OPA_ACCESSOR(Buffer, m_buffer, buffer);
  OPA_ACCESSOR_PTR(Texture, m_texture.get(), texture);
  void set_texture(TexturePtr texture) { m_texture = texture; }

  OPA_ACCESSOR(OrigMapType, m_orig_map, orig_map);

protected:
  // init face utils

  TexturePtr m_texture;
  TrVBO m_vbo;
  Buffer m_buffer;
};
OPA_DECL_SPTR(DisplayMesh, DisplayMeshPtr);

class MeshObj : public SceneObj {
public:
  virtual void init(Game *game, DisplayMeshPtr dmesh);

  virtual opa::math::game::Mesh *
  get_mesh(double precision = 0.) const override {
    return m_dmesh.get();
  };
  DisplayMeshPtr m_dmesh;
};

OPA_NAMESPACE_DECL2_END

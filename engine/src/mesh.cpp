#include "mesh/mesh.h"
#include "mesh/MeshDrawer.h"
#include "mesh/MeshIsec.h"
#include <opa/math/game/intersection.h>
#include <opa/utils/string.h>

// TODO: Fix get_one_vertex_attr
const float eps = 1e-6;
using namespace std;
using namespace opa::math::game;
using namespace opa::utils;

OPA_NAMESPACE_DECL2(opa, engine)

void MeshObj::init(Game *game, DisplayMeshPtr dmesh) {
  SceneObj::init(game);
  m_dmesh = dmesh;
  auto drawer = new mesh::MeshDrawer;
  drawer->init(mesh::MeshDrawer::Params(dmesh));
  auto isec = new mesh::MeshIsec;
  isec->init(dmesh);

  this->drawer().register_drawer(DrawerPtr(drawer));
  this->isec().register_isec(IsecPtr(isec));
}

void DisplayMesh::prepare_draw() {
  m_vbo.get().clear();
  GL_CHECK(;);
  for (auto &x : m_faces.used()) {
    vector<VertexId> tb;
    for (auto y : Mesh::Face_VertexWalker(*this, x))
      tb.pb(y.ND);
    FOR (i, 1, tb.size() - 1) { m_vbo.add_tr(tb[0], tb[i], tb[i + 1]); }
  }
  m_vbo.refresh();

  ScopedBuffer aa(m_buffer);
  GL_CHECK(glBufferData(GL_ARRAY_BUFFER,
                        sizeof(VertexAttribute) * m_attrs.tot_size(),
                        m_attrs.get_start(), GL_DYNAMIC_DRAW));
}

void DisplayMesh::init() {
  opa::math::game::Mesh::init();
  m_buffer.init(GL_ARRAY_BUFFER);
  m_vbo.init();
}

OPA_NAMESPACE_DECL2_END

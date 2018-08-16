#pragma ocne

#include <opa/engine/Drawer.h>
#include <opa/engine/mesh/mesh.h>

OPA_NAMESPACE_DECL3(opa, engine, mesh)

class MeshDrawer : public Drawer {
public:
  struct Params {
    Params(DisplayMeshPtr mesh) { this->mesh = mesh; }
    Params() {}
    DisplayMeshPtr mesh;
  };

  virtual void init(const Params &params);
  virtual void draw(const DrawData &data);

private:
  Params m_params;
  Program *m_program;
};
OPA_NAMESPACE_DECL3_END

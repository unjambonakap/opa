#pragma once

#include <opa/engine/conf.h>
#include <opa/engine/common.h>
#include <opa/engine/Intersection.h>
#include <opa/engine/mesh/mesh.h>

OPA_NAMESPACE_DECL3(opa, engine, mesh)

class MeshIsec : public Isec {

public:
  struct Params {
    DisplayMeshPtr mesh;
    Params() {}
    Params(DisplayMeshPtr mesh) { this->mesh = mesh; }
  };
  virtual void init(Params params);
  virtual bool find_intersection(const Ray &ray,
                                 IntersectionResult &res) override;

private:
  Params m_params;
};

OPA_NAMESPACE_DECL3_END

#include "mesh/MeshIsec.h"
#include <opa/math/game/intersection.h>

// TODO: Fix get_one_vertex_attr
const float eps = 1e-7;
using namespace std;
using namespace opa::math::game;

OPA_NAMESPACE_DECL3(opa, engine, mesh)

void MeshIsec::init(Params params) {
  Isec::init();
  m_params = params;
}

bool MeshIsec::find_intersection(const Ray &ray, IntersectionResult &res) {
  auto faces = m_params.mesh->list_faces();
  OPA_DISP_VARS(faces.size());

  for (auto x : faces) {
    std::vector<Pos> tb;

    for (auto d : Mesh::Face_VertexWalker(*m_params.mesh, x)) {
      const VertexAttribute &attr = m_params.mesh->get_attr(d.ND);
      tb.pb(attr.pos);
    }

    Pos center;
    for (auto &p : tb) center += p;
    center /= tb.size();

    std::vector<Pos2> proj_tb;
    Pos front = glm::normalize(tb[0] - center);
    Dir normal = get_plane_normal(tb[0] - center, tb[1] - center);

    for (auto &p : tb) proj_tb.pb(get_plane_coord(p - center, normal, front));
    Pos proj_point;

    if (!proj_point_plane(ray.pos - center, ray.dir, normal, proj_point))
      continue;

    Pos2 proj_p2 = get_plane_coord(proj_point, normal, front);
    if (inside_polygon(proj_tb, proj_p2)) {
      OPA_DISPL("FOUND AT >>> ", (proj_point + center), obj());
      out(tb);
      res.update_success(proj_point + center, ray, obj());
    }
  }

  return res.found;
}

OPA_NAMESPACE_DECL3_END

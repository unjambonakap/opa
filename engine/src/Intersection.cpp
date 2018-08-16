#include "Intersection.h"
#include "SceneTreeHelper.h"
#include <opa/engine/Scene.h>
#include <opa/engine/Strip.h>
#include <opa/math/game/intersection.h>

using namespace opa::math::game;
using namespace std;

OPA_NAMESPACE_DECL2(opa, engine)

Pos IntersectionResult::get_world_pos() { return obj->pos()->get_pos_to(pos); }

void Isec::update_intersection(const Ray &ray, IntersectionResult &best) {
  IntersectionResult cur;
  find_intersection(ray, cur);
  best.update(cur);
}

void IntersectionHandler::update_intersection(const Ray &ray,
                                              IntersectionResult &best) {
  IntersectionResult cur;
  find_intersection(ray, cur);
  best.update(cur);
}

void IntersectionHandler::register_isec(IsecPtr isec) {
  m_isec.insert(isec);
  isec->set_obj(m_obj);
}

bool IntersectionHandler::find_intersection(const Ray &ray,
                                            IntersectionResult &res,
                                            bool from_cam) {
  check_init();
  OPA_ASSERT0(res.found == false);
  res.reset();

  Ray local_ray;
  if (from_cam) {
    local_ray = ray;

  } else {
    local_ray.pos = obj()->pos_obj().get_local_pos(ray.pos);
    local_ray.dir = obj()->pos_obj().get_local_dir(ray.dir);
  }

  for (auto &isec : m_isec) isec->update_intersection(local_ray, res);

  for (const auto &elem : obj()->tree_struct()->children) {
    IntersectionResult res2;
    elem->obj->isec().update_intersection(local_ray, res2);
    res.update(res2);
  }

  return res.found;
}

static bool find_tr_intersection(const Ray &ray, const glm::vec3 &t1,
                                 const glm::vec3 &t2, const glm::vec3 &t3,
                                 glm::vec3 &res) {
  if (!line_parallelogram_intersection(ray.pos, ray.dir, t1, t2, t3, res, true))
    return false;
  return true;
}

bool SphereIsec::find_intersection(const Ray &ray, IntersectionResult &res) {
  check_init();
  res.found = false;

  if (!line_sphere_intersection(ray.pos, glm::normalize(ray.dir), m_params.r,
                                res.pos))
    return false;
  res.update_success(res.pos, ray, obj());
  return true;
}

bool TriangleContainerIsec::find_intersection(const Ray &ray,
                                              IntersectionResult &res) {
  res.reset();
  const Pos *p1, *p2, *p3;

  IntersectionResult cur;
  for (auto xid : m_params.tc->strips().used()) {
    const TexturedStrip &x = m_params.tc->strips().get(xid);
    REP (i, x.n()) {
      x.fill_tr(i, p1, p2, p3);
      if (find_tr_intersection(ray, *p1, *p2, *p3, cur.pos)) {
        cur.update_success(cur.pos, ray, obj());
        if (res.update(cur)) {
          res.debug = OPA_DISP_STR(*p1, *p2, *p3, ray.pos, ray.dir);
        }
      }
    }
  }

  for (auto xid : m_params.tc->triangles().used()) {
    const TexturedTriangle &x = m_params.tc->triangles().get(xid);
    x.fill_tr(p1, p2, p3);
    if (find_tr_intersection(ray, *p1, *p2, *p3, cur.pos)) {
      cur.update_success(cur.pos, ray, obj());
      if (res.update(cur)) {
        res.debug = OPA_DISP_STR(*p1, *p2, *p3, ray.pos, ray.dir);
      }
    }
  }
  return res.found;
}
OPA_NAMESPACE_DECL2_END

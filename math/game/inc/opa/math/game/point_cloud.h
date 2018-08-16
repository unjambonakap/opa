#pragma once

#include <opa/math/game/base_impl.h>
#include <opa/math/game/conf.h>
#include <opa/math/game/quat.h>
#include <opa/or/grid_search.h>

OPA_NAMESPACE_DECL3(opa, math, game)

double get_sphere_radius(const PointVec &points, const Pos &sphere_pos);

SphereSpec get_min_enclosing_sphere(const PointVec &points);
PointVec get_border_points(const PointVec &points, const SphereSpec &sphere);

PointVec box_cloud(const Pos &p0, const Pos &p1, const Pos &p2, const Pos &q0);

BoxSpec compute_best_box_vnaze(const PointVec &points);
BoxSpec compute_best_box(const PointVec &points);
BoxSpec compute_best_box_dumb(const PointVec &points);
BoxAASpec compute_aabb_box(const PointVec &points, const Rot &rot);
inline BoxAASpec compute_aabb_box_norot(const PointVec &points) {
  return compute_aabb_box(points, kIdentRot);
}

Box2DAASpec compute_aabb_box2d(const Point2Vec &points,
                               const Rot2 &rot = kIdentRot2);
Box2DSpec compute_best_box2d(const Point2Vec &points);
Box2DSpec compute_best_box2d_dumb(const Point2Vec &points);

Box2DAASpec compute_plan_aabb(const PointVec &points, const PlaneSpec &plane,
                              const Dir &front);

PointVec whiten(const PointVec &pv);
OPA_NAMESPACE_DECL3_END

#pragma once

#include <opa/math/game/base.h>
#include <opa/math/game/conf.h>
#include <opa/math/game/quat.h>

OPA_NAMESPACE_DECL3(opa, math, game)

// the two must intersect
Pos line_plane_intersection(const LineSpec &line, const PlaneSpec &plane);

bool proj_point_plane(const glm::vec3 &p, const glm::vec3 &d,
                      const glm::vec3 &n, glm::vec3 &res);
bool line_sphere_intersection(const glm::vec3 &p, const glm::vec3 &d, float r,
                              glm::vec3 &res);

bool line_parallelogram_intersection(const glm::vec3 &p, const glm::vec3 &d,
                                     const glm::vec3 &t1, const glm::vec3 &t2,
                                     const glm::vec3 &t3, glm::vec3 &res,
                                     bool is_tr);
bool line_parallelogram_intersection2(const glm::vec3 &p, const glm::vec3 &d,
                                      const glm::vec3 &t2, const glm::vec3 &t3,
                                      glm::vec3 &res, bool is_tr);

bool point_in_parallelogram(const glm::vec3 &p, const glm::vec3 &t2,
                            const glm::vec3 &t3, bool is_tr);

bool are_boxes_intersecting(const BoxSpec &a, const BoxSpec &b);
double get_boxes_intersection_area(const BoxSpec &a, const BoxSpec &b);
double get_boxes_intersection_area_grid(const BoxSpec &a, const BoxSpec &b, int grid_size = 20);

OPA_NAMESPACE_DECL3_END

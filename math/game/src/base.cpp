#include "opa/math/game/base.h"

#include <opa/math/common/Utils.h>
#include <opa/math/game/quat.h>
#include <opa/utils/string.h>

using namespace std;
static const float eps = 1e-4;

OPA_NAMESPACE_DECL3(opa, math, game)

std::uniform_real_distribution<float> unif_distrib(0., 1.);

glm::vec3 vec_rand_uni() {
  glm::vec3 res;
  REP (i, 3)
    res[i] = unif_distrib(opa::math::common::rng);
  return glm::normalize(res);
}

Pos2 vec2_rand_uni() {
  Pos2 res;
  // who said inefficient
  REP (i, 2)
    res[i] = unif_distrib(opa::math::common::rng);
  return glm::normalize(res);
}

float vec2_cross(const glm::vec2 &a, const glm::vec2 &b) {
  return a.x * b.y - a.y * b.x;
}

double positive_angle(double ang) {
  if (ang < 0) ang += 2 * OPA_PI;
  return ang;
}
float vec2_ang(const glm::vec2 &a) {
  float ang = atan2(a.y, a.x);
  return positive_angle(ang);
}

float tr_area(const glm::vec3 &a, const glm::vec3 &b) {
  return glm::length(glm::cross(a, b)) / 2;
}

bool are_aligned(const glm::vec3 &a, const glm::vec3 &b, const glm::vec3 &c) {
  return glm::isNull(glm::cross(glm::normalize(b - a), glm::normalize(c - a)),
                     eps);
}

bool inside_polygon_unsafe(const std::vector<glm::vec2> &tb,
                           const glm::vec2 &x) {
  Point2Vec centered;
  Pos2 g = gravity_center(tb);
  REP (i, tb.size()) {
    Pos2 pt = tb[i] - g;
    centered.push_back(pt);
  }
  Pos2 x2 = x - g;
  return inside_polygon(centered, x2, -1);
}

bool inside_polygon(const std::vector<glm::vec2> &tb, const glm::vec2 &x,
                    int front) {
  int n = tb.size();

  float want = vec2_ang(x);
  float cur_ang = vec2_ang(tb[0]);
  OPA_DISP0(tb, x);
  REP (i, n) {
    const Pos2 &p1 = tb[(i) % n];
    const Pos2 &p2 = tb[(i + 1) % n];
    float next_ang = vec2_ang(p2);

    if (next_ang < cur_ang) next_ang += 2 * OPA_PI;
    if (want < cur_ang) want += 2 * OPA_PI;
    OPA_DISP0(cur_ang, next_ang, want);

    if (want <= next_ang) {
      float tmp = vec2_cross(p2 - p1, x - p1);
      // OPA_DISP0(i, front, p2, p1, x, tmp, tb);
      return OPA_FLOAT_GE(tmp, 0, eps);
    }
    cur_ang = next_ang;
  }
  OPA_CHECK0(false);
  return false;
}

glm::vec3 get_plane_normal(const glm::vec3 &a, const glm::vec3 &b) {
  return glm::normalize(glm::cross(a, b));
}
glm::vec3 get_plane_normal(const glm::vec3 &a, const glm::vec3 &b,
                           const glm::vec3 &c) {
  OPA_CHECK(!are_aligned(a, b, c), a, b, c);
  return get_plane_normal(b - a, c - a);
}

float point_signed_dist_to_plane(const glm::vec3 &x, const glm::vec3 &n) {
  OPA_ASSERT0(glm::isNormalized(n, eps));
  return glm::dot(x, n);
}

float point_dist_to_plane(const glm::vec3 &x, const glm::vec3 &n) {
  return abs(point_signed_dist_to_plane(x, n));
}

glm::vec3 get_plane_orth_proj(const glm::vec3 &x, const glm::vec3 &n) {
  glm::vec3 proj = x - point_signed_dist_to_plane(x, n) * n;
  return proj;
}

glm::vec2 get_plane_coord(const glm::vec3 &x, const glm::vec3 &n,
                          const glm::vec3 &front, const Pos &left) {
  glm::vec3 proj = get_plane_orth_proj(x, n);
  return glm::vec2(glm::dot(proj, front), glm::dot(proj, left));
}

Pos2 solve_sys(const Mat2 &a, const Pos2 &b) {
  Pos2 res;
  double d = a[0][0] * a[1][1] - a[0][1] * a[1][0];
  return Pos2(b[0] * a[1][1] - b[1] * a[0][1], a[0][0] * b[1] - a[1][0] * b[0]) / d;
}

glm::vec2 get_plane_coord_safe(const glm::vec3 &x, const glm::vec3 &front,
                               const Pos &left) {
  Mat2 m;
  m[0][0] = glm::dot(front, front);
  m[0][1] = glm::dot(left, front);
  m[1][1] = glm::dot(left, left);
  m[1][0] = m[0][1];
  return solve_sys(m, Pos2(glm::dot(front, x), glm::dot(left, x)));
}

glm::vec2 get_plane_coord(const glm::vec3 &x, const glm::vec3 &n,
                          const glm::vec3 &front) {
  glm::vec3 left = glm::cross(n, front);
  return get_plane_coord(x, n, front, left);
}

float dist_line(const glm::vec3 &x, const glm::vec3 &a, const glm::vec3 &b) {
  return glm::dot(x - a, glm::normalize(b - a));
}

bool proj_point_plane(const glm::vec3 &p, const glm::vec3 &d,
                      const glm::vec3 &n, glm::vec3 &res) {
  float dot_d = glm::dot(n, d);
  float dot_p = glm::dot(n, p);
  res = p;

  if (!OPA_FLOAT_EQ(dot_d, 0, eps))
    res += -dot_p / dot_d * d;
  else if (!OPA_FLOAT_EQ(dot_p, 0, eps))
    return false;
  return true;
}

std::vector<HyperPlaneSpec> get_box_hyperplans(const BoxSpec &box) {
  std::vector<HyperPlaneSpec> res;
  REP (dir, 3) {
    Pos pos = box.get(0);
    pos += box.vec((dir + 1) % 3) / 2;
    pos += box.vec((dir + 2) % 3) / 2;
    Dir plane_dir = box.vec(dir);

    res.push_back(HyperPlaneSpec{ PlaneSpec::FromPointAndVec(pos, plane_dir) });
    res.push_back(HyperPlaneSpec{
      PlaneSpec::FromPointAndVec(pos + plane_dir, -plane_dir) });
  }
  return res;
}

BoxSpec get_in_box_space(const BoxSpec &src, const BoxSpec &target) {
  BoxSpec res;
  Rot r = rot_to_box_space(src);
  res.corner = r * (target.corner - src.corner);
  REP (i, 3)
    res.v[i] = r * target.v[i];
  return res;
}

OPA_NAMESPACE_DECL3_END

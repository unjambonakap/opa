#include "opa/math/game/intersection.h"

#include "opa/math/game/base_impl.h"

using namespace std;
static const float eps = 1e-6;

OPA_NAMESPACE_DECL3(opa, math, game)

bool line_sphere_intersection(const glm::vec3 &p, const glm::vec3 &d, float r,
                              glm::vec3 &res) {
  glm::vec3 proj;

  OPA_CHECK0(proj_point_plane(p, d, d, proj));

  if (glm::length(proj) > r) return false;
  float dk = sqrt(r * r - glm::length2(proj));

  res = proj - dk * d;
  return true;
}

bool line_parallelogram_intersection(const glm::vec3 &p, const glm::vec3 &d,
                                     const glm::vec3 &t1, const glm::vec3 &t2,
                                     const glm::vec3 &t3, glm::vec3 &res,
                                     bool is_tr) {
  bool has = line_parallelogram_intersection2(p - t1, glm::normalize(d),
                                              t2 - t1, t3 - t1, res, is_tr);
  res += t1;
  return has;
}

bool line_parallelogram_intersection2(const glm::vec3 &p, const glm::vec3 &d,
                                      const glm::vec3 &t2, const glm::vec3 &t3,
                                      glm::vec3 &res, bool is_tr) {
  glm::vec3 n = get_plane_normal(t2, t3);
  if (!proj_point_plane(p, d, n, res)) return false;
  Pos2 coord = get_plane_coord_safe(res, t2, t3);
  if (coord.x < 0 || coord.y < 0 || coord.x+coord.y > 1) return false;
  res = coord.x * t2 + coord.y * t3;
  //return point_in_parallelogram(res, t2, t3, is_tr);
  return true;
}

bool point_in_parallelogram(const glm::vec3 &p, const glm::vec3 &t2,
                            const glm::vec3 &t3, bool is_tr) {

  float l2 = glm::length(t2);
  float l3 = glm::length(t3);
  glm::vec3 u = t2 / l2;
  glm::vec3 v = t3 / l3;
  glm::vec3 v2 = v - u * glm::dot(u, v);

  float b = glm::dot(p, v2);
  if (!OPA_FLOAT_BETWEEN(b, 0, l3, eps)) return false;
  float a = glm::dot(p - v * b, u);
  if (!OPA_FLOAT_BETWEEN(a, 0, l2, eps)) return false;

  if (is_tr && !OPA_FLOAT_LE(a / l2 + b / l3, 1, eps)) return false;
  return true;
}

Pos line_plane_intersection(const LineSpec &line, const PlaneSpec &plane) {
  double da = plane.get_signed_dist(line.a);
  double db = plane.get_signed_dist(line.b);

  double ca = db / (db - da);
  return line.a * ca + line.b * (1 - ca);
}

bool are_boxes_intersecting(const BoxSpec &a, const BoxSpec &b) {
  if (a.area() < b.area()) return are_boxes_intersecting(b, a);
  BoxSpec b2 = get_in_box_space(a, b);
  BoxAASpec a2 = BoxAASpec::FromBoxSpace(a);

  auto edges = b.edges();
  PointVec points = b2.corners();
  ;
  for (const auto &pt : points)
    if (a2.in(pt)) return true;

  for (auto &edge : edges) {
    Pos pa, pb;
    pa = points[edge.first];
    pb = points[edge.second];
    REP (j, 3) {
      Range1D r1 = Range1D::MakeRange(pa[j], pb[j]);
      int oc1 = (j + 1) % 3;
      int oc2 = (j + 2) % 3;

      if (r1.is_almost_empty()) continue;
      for (double v : { (double)0., (double)a2.high[j] }) {
        if (r1.is_in(v)) {
          double k = (v - pa[j]) / (pb[j] - pa[j]);
          Pos proj_p = pa * (1 - k) + k * pb;
          if (proj_p[oc1] < -base_eps || proj_p[oc1] > a2.high[oc1] + base_eps)
            continue;
          if (proj_p[oc2] < -base_eps || proj_p[oc2] > a2.high[oc1] + base_eps)
            continue;
          return true;
        }
      }
    }
  }

  return false;
}

double get_boxes_intersection_area(const BoxSpec &a, const BoxSpec &b) {
  if (a.area() < b.area()) return get_boxes_intersection_area(b, a);
  BoxSpec b2 = get_in_box_space(a, b);
  BoxAASpec a2 = BoxAASpec::FromBoxSpace(a);

  std::array<std::vector<double>, 3> axis_points;
  REP (j, 3) {
    axis_points[j].push_back(0);
    axis_points[j].push_back(a2.high[j]);
  }

  for (auto &pt : b2.corners()) {
    REP (j, 3) { axis_points[j].push_back(pt[j]); }
  }

  bool has_out = false;
  REP (j, 3) {
    auto &cur = axis_points[j];
    std::sort(ALL(cur));
    if (cur[0] < -base_eps || cur.back() > a2.high[j] + base_eps) {
      has_out = true;
      break;
    }
  }
  if (!has_out) return b.area();
  OPA_CHECK(false, "not implemented");
}

double get_boxes_intersection_area_grid(const BoxSpec &a, const BoxSpec &b,
                                        int grid_size) {
  if (a.area() < b.area())
    return get_boxes_intersection_area_grid(b, a, grid_size);
  BoxSpec b2 = get_in_box_space(a, b);
  BoxAASpec a2 = BoxAASpec::FromBoxSpace(a);
  Rot r = rot_to_box_space(a);
  REP (i, 8) {
    OPA_DISP0(b2.get(i), r * (b.get(i) - a.corner), r * (a.get(i) - a.corner),
              a2.get(i));
  }

  std::vector<std::vector<double> > grid_def;
  REP (i, 3)
    grid_def.push_back(utils::Range<double>::Build_range(0, 1, grid_size).tb());

  u32 sum = 0;
  OR::CrossProdGen<double> cp(grid_def, [&](const std::vector<double> &tb) {
    Pos p = b2.get(tb);
    if (a2.in(p)) ++sum;
  });

  double step = std::pow((1. / grid_size), 3);
  double ratio = sum * step;
  return ratio * b2.area();
}

OPA_NAMESPACE_DECL3_END

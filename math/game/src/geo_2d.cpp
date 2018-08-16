#include <opa/math/game/geo_2d.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/partition_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Partition_traits_2<K> Traits;
typedef Traits::Point_2 Point_2;
typedef Traits::Polygon_2 Polygon_2;
typedef Polygon_2::Vertex_iterator Vertex_iterator;
typedef std::list<Polygon_2> Polygon_list;

OPA_NM_MATH_GAME

std::vector<int> compute_convex_hull(const Point2Vec &cloud) {

  int n = cloud.size();

  std::vector<int> hull;
  std::vector<std::pair<Pos2, int> > points;
  REP (i, n)
    points.pb(MP(cloud[i], i));

  std::sort(ALL(points), vec_cmp(base_eps));
  utils::make_unique(points, vec_cmp_eq(base_eps, 1));
  n = points.size();

  if (n == 1) return { 0 };

  REP (hull_dir, 2) {

    std::vector<int> stack;
    REP (i, n) {
      if (points[i].second == -1) continue;

      while (stack.size() >= 2) {
        int last = stack.size() - 1;
        int b = stack[last];
        int a = stack[last - 1];
        int c = i;
        if (OPA_FLOAT_GE(vec2_cross(points[b].ST - points[a].ST,
                                    points[c].ST - points[a].ST),
                         0, base_eps))
          stack.pop_back();
        else
          break;
      }
      stack.push_back(i);
    }

    if (hull_dir) {
      std::reverse(ALL(stack));
    }
    REP (i, stack.size() - 1)
      hull.pb(points[stack[i]].second);

    if (hull_dir == 0) {
      FOR (i, 1, stack.size() - 1)
        points[stack[i]].second = -1;
      REP (i, n)
        points[i].ST.y *= -1;
    }
  }
  // trigo walk
  std::reverse(ALL(hull));
  return hull;
}

FaceIndexList triangulate_polygon_indices(const Point2Vec &polygon) {
  FaceIndexList res;
  Polygon_2 pg2;
  Polygon_list partition_polys;
  utils::Remapper<Pos2> rmp;
  Point2Vec tb = reorient_counterclockise(polygon);

  for (auto &pt : tb) {
    Point_2 pt2(pt.x, pt.y);
    pg2.push_back(pt2);
    rmp.get(pt);
  }
  CGAL::approx_convex_partition_2(pg2.vertices_begin(), pg2.vertices_end(),
                                  std::back_inserter(partition_polys));

  for (const Polygon_2 &conv_poly : partition_polys) {
    std::vector<int> ids;
    for (auto &pt : conv_poly.container()) {
      Pos2 pt2 = Pos2(pt.x(), pt.y());
      ids.push_back(rmp.get_or_die(pt2));
    }

    auto lst = triangulate_convex_polygon_indices(ids);
    res.insert(res.end(), ALL(lst));
  }

  return res;
}

FaceIndexList triangulate_polygon_indices(const PointVec &polygon) {
  Point2Vec proj_points;
  PlaneSpec plane = PlaneSpec::FromPoints(polygon);
  auto front = plane.get_rand_plane_vec();
  for (auto &pt : polygon) proj_points.push_back(plane.proj(pt, front));
  return triangulate_polygon_indices(proj_points);
}

int line_line_intersection_2d(const LineSpec2D &l1, const LineSpec2D &l2,
                              Pos2 &p1, Pos2 &p2) {
  Pos2 v = l1.b - l1.a;
  Pos2 w = l2.b - l2.a;
  double lv2 = glm::length2(v);
  double pa = glm::dot(l2.a - l1.a, v) / lv2;
  double pb = glm::dot(l2.b - l1.a, v) / lv2;
  if (std::min(pa, pb) > 1) return 0;
  if (std::max(pa, pb) < 0) return 0;

  Pos2 orth_norm = glm::normalize(vec2_orth(v));
  double dist = glm::dot(l1.a - l2.a, orth_norm);
  double k = glm::dot(w, orth_norm);

  if (std::abs(dist) < base_eps) {
    if (std::abs(k) > base_eps) return 0;
    double x1 = std::max<double>(0, std::min(pa, pb));
    double x2 = std::min<double>(1, std::max(pa, pb));

    p1 = from_line_coord(l1, x1);
    if (x1 + base_eps < x2) {
      p2 = from_line_coord(l1, x2);
      return 2;
    }
    return 1;
  }


  dist /= k;
  if (dist < 0 || dist > 1) return 0;
  double d1 = pa + (pb - pa) * dist;
  if (d1 < 0 || d1 > 1) return 0;
  p1 = from_line_coord(l2, dist);

  OPA_FLOAT_EQ(glm::dot(p1 - l1.a, v) / lv2, d1, base_eps);
  return 1;
}

int line_line_intersection_2d_v2(const LineSpec2D &l1, const LineSpec2D &l2,
                                 Pos2 &p1, Pos2 &p2) {
  Pos2 v1 = l1.b - l1.a;
  Pos2 v2 = l2.b - l2.a;
  double det = vec2_cross(v1, -v2);
  Pos2 b = l2.a - l1.a;

  if (std::abs(det) < base_eps) {

    if (std::abs(vec2_cross(b, v1)) > base_eps) return 0;
    Pos2 v1norm = glm::normalize(v1);
    double k1 = glm::dot(v1norm, b);
    double k2 = glm::dot(v1norm, l2.b - l1.a);

    if (k1 > k2) std::swap(k1, k2);
    k1 = std::max<double>(k1, 0);
    k2 = std::min<double>(k2, 1);
    if (k1 > k2 + base_eps) return 0;

    p1 = from_line_coord(l1, k1);
    if (k1 + base_eps > k2) return 1;
    p2 = from_line_coord(l1, k2);
    return 2;
  }

  double k1 = vec2_cross(b, -v2) / det;
  double k2 = vec2_cross(v1, b) / det;

  if (k1 < 0 || k1 > 1 || k2 < 0 || k2 > 1) return 0;
  p1 = from_line_coord(l1, k1);
  return 1;
}

Point2Vec tr_tr_intersection(const Tr2D &tr1, const Tr2D &tr2) {
  Point2Vec res;

  REP (i, 3) {
    if (inside_polygon_v2(tr2, tr1[i])) res.push_back(tr1[i]);
    if (inside_polygon_v2(tr1, tr2[i])) res.push_back(tr2[i]);
  }
  REP (i, 3) {
    REP (j, 3) {
      Pos2 p1, p2;
      int cnt = math::game::line_line_intersection_2d(
        LineSpec2D{ tr1[i], tr1[(i + 1) % 3] },
        LineSpec2D{ tr2[j], tr2[(j + 1) % 3] }, p1, p2);
      if (cnt > 0) res.push_back(p1);
      if (cnt > 1) res.push_back(p2);
    }
  }
  auto ids = math::game::compute_convex_hull(res);
  return utils::select_ids(res, ids);
}
OPA_NM_MATH_GAME_END

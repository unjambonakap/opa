#pragma once

#include <opa/math/game/base.h>
#include <opa/math/game/conf.h>

OPA_NM_MATH_GAME

inline Pos2 from_line_coord(const LineSpec2D &l, double x) {
  return l.a + (l.b - l.a) * x;
}

std::vector<int> compute_convex_hull(const Point2Vec &cloud);

inline FaceIndexList
triangulate_convex_polygon_indices(const std::vector<int> &ids = {}) {
  int n = ids.size();
  FaceIndexList res;
  REP (i, n - 2)
    res.push_back({ ids[0], ids[i + 1], ids[i + 2] });
  return res;
}

inline FaceIndexList triangulate_convex_polygon_indices(int n) {
  FaceIndexList res;
  REP (i, n - 2)
    res.push_back({ 0, i + 1, i + 2 });
  return res;
}

FaceIndexList triangulate_polygon_indices(const Point2Vec &polygon);
FaceIndexList triangulate_polygon_indices(const PointVec &polygon);

template <class T>
std::vector<std::vector<T> >
triangulate_polygon(const std::vector<T> &polygon) {
  auto normal = face_normal(polygon);
  auto tr_list = triangulate_polygon_indices(polygon);
  std::vector<std::vector<T> > res;
  for (auto &tr : tr_list) {
    std::vector<T> cur;
    for (auto &e : tr) cur.push_back(polygon[e]);
    res.push_back(realign_face(cur, normal));
  }
  return res;
}

template <class ContainerT>
bool inside_polygon_v2(const ContainerT &tb, const Pos2 &x) {
  REP (i, tb.size()) {
    if (vec2_cross(tb[(i + 1) % tb.size()] - tb[i], x - tb[i]) < -base_eps)
      return false;
  }
  return true;
}

inline float tr_area(const Pos2 &a, const Pos2 &b) {
  return std::abs(a[0] * b[1] - a[1] * b[0]) / 2;
}

inline double polygon_area(const Point2Vec &polygon) {
  double area = 0;
  for (auto &face : triangulate_convex_polygon_indices(polygon.size())) {
    area += tr_area(polygon[face[1]] - polygon[face[0]],
                    polygon[face[2]] - polygon[face[0]]);
  }
  return area;
}


inline double polygon_area(const PointVec &polygon) {
  double area = 0;
  for (auto &face : triangulate_convex_polygon_indices(polygon.size())) {
    area += tr_area(polygon[face[1]] - polygon[face[0]],
                    polygon[face[2]] - polygon[face[0]]);
  }
  return area;
}

Point2Vec tr_tr_intersection(const Tr2D &tr1, const Tr2D &tr2);
int line_line_intersection_2d_v2(const LineSpec2D &l1, const LineSpec2D &l2,
                                 Pos2 &p1, Pos2 &p2);

int line_line_intersection_2d(const LineSpec2D &l1, const LineSpec2D &l2,
                              Pos2 &p1, Pos2 &p2);
OPA_NM_MATH_GAME_END

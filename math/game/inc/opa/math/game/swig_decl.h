#pragma once

#include <opa_common.h>
#include <opa/math/game/base.h>
#include <opa/math/game/geo_2d.h>

// trnasition to swig::from?
OPA_NM_MATH_GAME
inline bool inside_polygon_v2_tr2d(const Tr2D &tb, const Pos2 &x) {
  return inside_polygon_v2(tb, x);
}

inline bool inside_polygon_v2_point2vec(const Point2Vec &tb, const Pos2 &x) {
  return inside_polygon_v2(tb, x);
}

inline bool is_counterclockwise_tr2d(const Tr2D &face) {
  return is_counterclockwise(face);
}
inline bool is_counterclockwise_point2vec(const Point2Vec &face) {
  return is_counterclockwise(face);
}

inline Tr2D reorient_counterclockise_tr2d(const Tr2D &face) {
  return reorient_counterclockise(face);
}
inline Point2Vec reorient_counterclockise_point2vec(const Point2Vec &face) {
  return reorient_counterclockise(face);
}

OPA_NM_MATH_GAME_END

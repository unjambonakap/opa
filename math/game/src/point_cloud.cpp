#include <opa/math/game/geo_2d.h>
#include <opa/math/game/point_cloud.h>
#include <opa/math/game/quat.h>
#include <opa/or/adaptative_search.h>
#include <opencv2/opencv.hpp>
#include <opa/math/common/stats.h>

template <class T> using PointAndScore = opa::OR::PointAndScore<T>;
using opa::OR::C1Manifold;
using namespace cv;

OPA_NAMESPACE_DECL3(opa, math, game)
const double eps = 1e-6;

double get_sphere_radius(const PointVec &points, const Pos &sphere_pos) {
  double best = 0;
  for (auto &p : points) {
    best = std::max<double>(best, glm::length(sphere_pos - p));
  }
  return best;
}

PointVec get_border_points(const PointVec &points, const SphereSpec &sphere) {
  std::vector<Pos> res;
  for (auto &p : points) {
    if (std::abs(glm::length(sphere.center - p) - sphere.radius) < eps)
      res.push_back(p);
  }
  return res;
}

SphereSpec get_min_enclosing_sphere(const PointVec &points) {
  SphereSpec res;
  /*
     1. go to center of gravity + find appropriate radius
     2. while all border points directions are strictly contained in a half
     plane
        3. ternary search in improving direction

     */

  res.center = gravity_center(points);
  res.radius = get_sphere_radius(points, res.center);

  REP (nstep, 100) {
    auto border_points = get_border_points(points, res);
    Pos improve_dir =
      glm::normalize(gravity_center(border_points) - res.center);

    double L = 0, H = 1e100;
    for (auto &p : border_points) {
      H = std::min<double>(H, glm::dot(p - res.center, improve_dir));
    }
    if (H < eps) break;

    double score_l = res.radius;
    double score_h = get_sphere_radius(points, res.center + improve_dir * H);
    REP (step2, 20) {
      double T1 = (2 * L + H) / 3;
      double T2 = (L + 2 * H) / 3;
      double score_t1 =
        get_sphere_radius(points, res.center + improve_dir * T1);
      double score_t2 =
        get_sphere_radius(points, res.center + improve_dir * T2);
      if (score_t1 < score_t2)
        H = T2, score_h = score_t2;
      else
        L = T1, score_l = score_t1;
    }
    res.center = res.center + improve_dir * L;
    res.radius = score_l;
  }

  return res;
}

PointVec box_cloud(const Pos &p0, const Pos &p1, const Pos &p2, const Pos &q0) {
  PointVec res = { p0, p1, p2, p1 + p2 - p0 };
  REP (i, 4)
    res.push_back(res[i] - (q0 - p0));
  return res;
}

BoxAASpec compute_aabb_box(const PointVec &points, const Rot &rot) {
  BoxAASpec res;
  Rot irot = glm::inverse(rot);
  for (auto &pt : points) {
    res.update(irot * pt);
  }
  return res;
}

Box2DAASpec compute_aabb_box2d(const Point2Vec &points, const Rot2 &rot) {
  Box2DAASpec res;
  Rot2 irot = glm::inverse(rot);
  for (auto &pt : points) {
    res.update(irot * pt);
  }
  return res;
}

Rot find_best_box_init_rot(const PointVec &points) {

  int n = points.size();
  Mat data_pts = Mat(n, 3, CV_64FC1);
  REP (i, n) {
    data_pts.at<double>(i, 0) = points[i].x;
    data_pts.at<double>(i, 1) = points[i].y;
    data_pts.at<double>(i, 2) = points[i].z;
  }

  // Perform PCA analysis
  PCA pca_analysis(data_pts, Mat(), cv::PCA::DATA_AS_ROW);
  std::vector<Pos> eigen_vecs(3);
  REP (i, 3) {
    eigen_vecs[i] = Pos(pca_analysis.eigenvectors.at<double>(i, 0),
                        pca_analysis.eigenvectors.at<double>(i, 1),
                        pca_analysis.eigenvectors.at<double>(i, 2));
  }
  if (glm::length(eigen_vecs[0]) < base_eps) {
    eigen_vecs[0] = vec_x;
    eigen_vecs[1] = vec_y;
  } else if (glm::length(eigen_vecs[1]) < base_eps) {
    eigen_vecs[1] = vec_ortho_rand(eigen_vecs[0]);
  }
  //OPA_TRACES(eigen_vecs, eigen_vecs[0], eigen_vecs[1], quat_look_at_safe(eigen_vecs[0], eigen_vecs[1]));
  return quat_look_at_safe(eigen_vecs[0], eigen_vecs[1]);
}

BoxSpec compute_best_box_dumb(const PointVec &points) {
  int n = points.size();
  Rot rot = find_best_box_init_rot(points);
  OR::AggregateStopCriteria criteria(
    { OR::RoundStoppingCriteria(40), OR::MinimaStoppingCriteria(5) });

  RotManifold manifold;
  typedef OR::ConstNDimManagedDesc<Pos, RotManifold> CurManagedDesc;
  CurManagedDesc managed_desc(manifold, points.size(), 50,
                              RotManifold::from_rot(rot));

  typedef OR::ManagedSimpleAdaptativeStrategy<Pos, RotManifold, CurManagedDesc>
    Strategy;

  auto strategy = Strategy(managed_desc, criteria, manifold);
  auto search = OR::AdaptativeSearcher<Pos, Strategy>(
    [&](const Pos &a) {
      auto res = compute_aabb_box(points, RotManifold::to_rot(a));
      return res.area();
    },
    strategy);
  auto result = search.find_minimum();
  rot = RotManifold::to_rot(result.point);

  BoxSpec res = compute_aabb_box(points, rot).to_box();
  res.corner = rot * res.corner;
  REP (i, 3)
    res.v[i] = rot * res.v[i];
  return res;
}

BoxSpec compute_best_box_vnaze(const PointVec &points) {
  int n = points.size();
  Rot rot = find_best_box_init_rot(points);
  OR::AggregateStopCriteria criteria(
    { OR::RoundStoppingCriteria(40), OR::MinimaStoppingCriteria(5) });

  RotManifold manifold;
  typedef OR::IterativeSearcher<Rot> Searcher;
  auto search = Searcher(
    [&](const Rot &a, const Searcher::IterationParams &params) {
      Pos cnd = RotManifold::from_rot(a);
      typedef OR::ConstNDimManagedDesc<Pos, RotManifold> CurManagedDesc;
      CurManagedDesc managed_desc(manifold, points.size(), 50, cnd,
                                  { rand() % 3 });

      typedef OR::ManagedSimpleAdaptativeStrategy<Pos, RotManifold,
                                                  CurManagedDesc>
        Strategy;
      auto strategy = Strategy(managed_desc, criteria, manifold);
      auto search = OR::AdaptativeSearcher<Pos, Strategy>(
        [&](const Pos &a) {
          auto res = compute_aabb_box(points, RotManifold::to_rot(a));
          return res.area();
        },
        strategy);
      double score;
      auto result = search.find_minimum();
      return PointAndScore<Rot>{ RotManifold::to_rot(cnd), result.score };
    },
    criteria, quat_dist);
  auto result = search.find_minimum(rot);
  rot = result.point;

  BoxSpec res = compute_aabb_box(points, rot).to_box();
  res.corner = rot * res.corner;
  REP (i, 3)
    res.v[i] = rot * res.v[i];
  return res;
}

BoxSpec compute_best_box_planar(const PointVec &points) {
  PlaneSpec plane = PlaneSpec::FromPoints(points);
  Point2Vec proj_points;
  auto front = plane.get_rand_plane_vec();
  for (auto &pt : points) proj_points.push_back(plane.proj(pt, front));

  Box2DSpec box2d = compute_best_box2d(proj_points);

  BoxSpec box;
  Dir dz = plane.dir * base_eps / 2;
  box.corner = plane.lift(box2d.corner, front) - dz;
  Pos p1 = plane.lift(box2d.get(1), front) - dz - box.corner;
  Pos p2 = plane.lift(box2d.get(2), front) - dz - box.corner;
  box.v[0] = p1;
  box.v[1] = p2;
  box.v[2] = dz * 2;
  return box;
}

BoxSpec compute_best_box(const PointVec &points) {
  int n = points.size();
  if (find_non_planar_point(points).empty()) {
    auto non_aligned = find_non_aligned_point(points);
    if (non_aligned.empty()) {
      int other_id = find_other_point(points);
      if (other_id == -1) return BoxSpec({ points[0], vec_mat * base_eps });
      Dir vec_dir = glm::normalize(points[other_id] - points[0]);

      utils::MinFinderPair<double, Pos> minf, maxf;
      for (auto &pt : points) {
        minf.update(glm::dot(pt - points[0], vec_dir), pt);
        maxf.update(-glm::dot(pt - points[0], vec_dir), pt);
      }
      Dir veca = maxf.get() - minf.get();
      Dir vecb = vec_ortho_rand(veca);
      Dir vecc = get_plane_normal(veca, vecb);
      return BoxSpec{ minf.get(), { veca, vecb * base_eps, vecc * base_eps } };
    }
    return compute_best_box_planar(points);
  }

  Rot rot = find_best_box_init_rot(points);
  OR::AggregateStopCriteria criteria(
    { OR::RoundStoppingCriteria(40), OR::MinimaStoppingCriteria(5) });

  RotManifold manifold;
  typedef OR::IterativeSearcher<Rot> Searcher;
  auto search = Searcher(
    [&](const Rot &a, const Searcher::IterationParams &params) {
      int dir = rand() % 3;
      Pos up = vec_tb[(dir + 2) % 3];
      Pos z_plane = rot * up;
      Pos front = rot * vec_tb[dir];
      Point2Vec proj_points;
      for (auto &e : points) {
        proj_points.push_back(get_plane_coord(e, z_plane, front));
      }
      Box2DSpec res = compute_best_box2d(proj_points);
      double rot_angle = d2_rot_to_angle(res.v[0]);
      rot = rot * quat_from_vec_rot_safe(up, rot_angle);
      BoxAASpec res2 = compute_aabb_box(points, rot);
      return PointAndScore<Rot>{ rot, res2.area() };
    },
    criteria, quat_dist);
  rot = search.find_minimum(rot).point;

  BoxSpec res = compute_aabb_box(points, rot).to_box();
  res.corner = rot * res.corner;
  REP (i, 3)
    res.v[i] = rot * res.v[i];

#if OPA_DEBUG
  BoxSpec tmp = res;
  tmp.expand(1 + 1e-3);
  for (auto &pt: points){
    OPA_CHECK(tmp.in(pt), tmp.str(), pt, tmp.dist2(pt), points);
  }
#endif
  return res;
}

Box2DSpec compute_best_box2d(const Point2Vec &points) {
  Box2DSpec res;
  std::vector<int> hull_ids = compute_convex_hull(points);
  int n = hull_ids.size();

  if (hull_ids.size() == 1) {
    res.corner = points[hull_ids[0]];
    res.v[0] = vec2_x * 1e-3;
    res.v[1] = vec2_y * 1e-3;
    return res;
  }
  REP (i, n)
    hull_ids.push_back(hull_ids[i]);
  REP (i, 2)
    hull_ids.push_back(hull_ids[i % n]);

  std::vector<double> angles;
  REP (i, 2 * n + 1)
    angles.push_back(vec2_ang(points[hull_ids[i + 1]] - points[hull_ids[i]]));

  int pright = 0, pup = 0, pleft = 0;
  utils::MinFinderPair<double, Box2DSpec> best;

  REP (pdown, n) {
    for (; positive_angle(angles[pright] - angles[pdown]) < PI / 2; ++pright)
      ;
    pup = std::max(pup, pright);
    for (; pup != pdown + n && positive_angle(angles[pup] - angles[pdown]) < PI;
         ++pup)
      ;
    pleft = std::max(pleft, pup);
    for (; pleft != pdown + n &&
           positive_angle(angles[pleft] - angles[pdown]) < 3 * PI / 2;
         ++pleft)
      ;
    Box2DSpec cnd;
    Pos2 bottom = points[hull_ids[pdown]];
    Pos2 right = points[hull_ids[pright]];
    Pos2 up = points[hull_ids[pup]];
    Pos2 left = points[hull_ids[pleft]];
    Pos2 l1 =
      glm::normalize(points[hull_ids[pdown + 1]] - points[hull_ids[pdown]]);
    Pos2 l2 = vec2_orth(l1);
    cnd.corner = bottom + l1 * glm::dot(left - bottom, l1);
    cnd.v[0] = l1 * glm::dot(right - left, l1);
    cnd.v[1] = l2 * glm::dot(up - bottom, l2);
    best.update(cnd.area(), cnd);
  }

  return best.get();
}

Box2DSpec compute_best_box2d_dumb(const Point2Vec &points) {
  typedef OR::ConstRManagedDesc<double, C1Manifold> CurManagedDesc;
  CurManagedDesc managed_desc(CircleManifold, points.size() * 2, 50);
  OR::RoundStoppingCriteria stop_criteria(40);

  typedef OR::ManagedSimpleAdaptativeStrategy<double, C1Manifold,
                                              CurManagedDesc>
    Strategy;
  auto strategy = Strategy(managed_desc, stop_criteria, CircleManifold);
  auto search = OR::AdaptativeSearcher<double, Strategy>(
    [&](const double &a) {
      auto res = compute_aabb_box2d(points, d2_rot(a * 2 * PI));
      return res.area();
    },
    strategy);
  auto result = search.find_minimum();
  Rot2 rot = d2_rot(result.point * 2 * PI);
  Box2DSpec res = compute_aabb_box2d(points, rot).to_box();
  res.corner = rot * res.corner;
  res.v = rot * res.v;
  return res;
}

Box2DAASpec compute_plan_aabb(const PointVec &points, const PlaneSpec &plane,
                              const Dir &front) {
  Point2Vec tb;
  for (auto &pt : points) {
    tb.push_back(plane.proj(pt, front));
  }
  return compute_aabb_box2d(tb);
}

PointVec whiten(const PointVec &pv) {
  math::common::DataWhitener<Pos, 3> whitener;
  whitener.add(pv);
  whitener.compute();
  OPA_DISP0(pv);
  OPA_DISP0(whitener.remap(pv));
  return whitener.remap(pv);
}

OPA_NAMESPACE_DECL3_END

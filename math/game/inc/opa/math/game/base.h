#pragma once

#include <opa/math/common/rng.h>
#include <opa/math/game/conf.h>
#include <opa/math/game/consts.h>
#include <opa/utils/base.h>

namespace opa {
static inline Pos Pos4ToPos(const Pos4 &pos) { return pos.xyz() / pos.w; }
static inline Pos ApplyMat(const Mat4 &mat, const Pos &pos) {
  return Pos4ToPos(mat * Pos4(pos, 1.));
}
} // namespace opa

OPA_NM_MATH_GAME

using opa::utils::Range1D;

extern std::uniform_real_distribution<float> unif_distrib;

float vec2_cross(const glm::vec2 &a, const glm::vec2 &b);
float vec2_ang(const glm::vec2 &a);
static inline Pos2 vec2_orth(const Pos2 &a) { return Pos2(-a.y, a.x); }
static inline Mat2 d2_rot(double ang) {
  double c = std::cos(ang);
  double s = std::sin(ang);
  // x1, y1, x2, y2, column major
  return Mat2(c, -s, s, c);
}
static inline double d2_rot_to_angle(const Pos2 &rot) {
  return atan2(rot.y, rot.x);
}

float tr_area(const glm::vec3 &a, const glm::vec3 &b);

inline double tetrahedron_area_signed(const Pos &a, const Pos &b,
                                      const Pos &c) {
  return glm::dot(a, glm::cross(b, c)) / 6;
}
inline double tetrahedron_area(const Pos &a, const Pos &b, const Pos &c) {
  return std::abs(tetrahedron_area_signed(a, b, c));
}

inline double tetrahedron_area(const Pos &a, const Pos &b, const Pos &c,
                               const Pos &d) {
  return tetrahedron_area(b - a, c - a, d - a);
}

bool are_aligned(const glm::vec3 &a, const glm::vec3 &b, const glm::vec3 &c);
inline bool are_aligned2(const Pos2 &a, const Pos2 &b, const Pos2 &c) {
  return std::abs(glm::dot(b - a, c - a)) < base_eps;
}
inline bool are_aligned(const Pos2 &a, const Pos2 &b, const Pos2 &c) {
  return are_aligned2(a, b, c);
}

double positive_angle(double ang);

glm::vec3 vec_rand_uni();
Pos2 vec2_rand_uni();

template <typename T, int N> T vec_rand_uni_t() {
  T res;
  REP (i, N)
    res[i] = unif_distrib(opa::math::common::rng);
  return glm::normalize(res);
}

inline double acos_safe(double v) {
  if (v < -1.) v = -1.;
  if (v > 1) v = 1;
  return std::acos(v);
}

static inline float vec_get_angle(const glm::vec3 &a, const glm::vec3 &b) {
  return acos_safe(glm::dot(a, b)); // suppose unit vec
}
static inline float vec_get_angle_safe(const glm::vec3 &a, const glm::vec3 &b) {
  return acos_safe(glm::dot(glm::normalize(a), glm::normalize(b))); // suppose
                                                                    // unit vec
}
static inline float vec_get_angle_safe(const glm::vec2 &a, const glm::vec2 &b) {
  return glm::angle(glm::normalize(a), glm::normalize(b));
}

template <typename VecType>
VecType vec_make_ortho(const VecType &a, const VecType &b) {
  return glm::normalize(a - glm::dot(a, b) * b);
}
template <typename VecType>
VecType vec_make_ortho_safe(const VecType &a, const VecType &b) {
  return vec_make_ortho(glm::normalize(a), glm::normalize(b));
}

inline glm::vec3 vec_ortho_rand(const glm::vec3 &a) {
  return vec_make_ortho(vec_rand_uni(), a);
}

static inline bool InFOV(const Mat4 &proj_mat, const Pos4 &pos,
                         double pad = 0.) {
  Pos tmp = Pos4ToPos(proj_mat * pos);
  return std::abs(tmp.x) < 1 - pad && std::abs(tmp.y) < 1 - pad;
}

static inline double compute_perspective_fovx(const Mat4 &imat) {
  Pos p1 = ApplyMat(imat, Pos(1, 0, 1));
  Pos p2 = ApplyMat(imat, Pos(-1, 0, 1));
  return vec_get_angle_safe(p1.xy(), p2.xy());
}

static inline double compute_perspective_fovy(const Mat4 &imat) {
  Pos p1 = ApplyMat(imat, Pos(0, 1, 1));
  Pos p2 = ApplyMat(imat, Pos(0, -1, 1));
  return vec_get_angle_safe(p1.xz(), p2.xz());
}

template <class T> T gravity_center(const std::vector<T> &points) {
  T res;
  for (auto &e : points) res += e;
  res /= points.size();
  return res;
}

float dist_line(const glm::vec3 &x, const glm::vec3 &a, const glm::vec3 &b);

float point_signed_dist_to_plane(const glm::vec3 &x, const glm::vec3 &n);
float point_dist_to_plane(const glm::vec3 &x, const glm::vec3 &n);
glm::vec3 get_plane_orth_proj(const glm::vec3 &x, const glm::vec3 &n);
// project on plane with normal n, through O. Front vector gives x vec
glm::vec2 get_plane_coord(const glm::vec3 &x, const glm::vec3 &n,
                          const glm::vec3 &front);
glm::vec2 get_plane_coord(const glm::vec3 &x, const glm::vec3 &n,
                          const glm::vec3 &front, const glm::vec3 &left);

glm::vec2 get_plane_coord_safe(const glm::vec3 &x, const glm::vec3 &front,
                               const glm::vec3 &left);
glm::vec3 get_plane_normal(const glm::vec3 &a, const glm::vec3 &b,
                           const glm::vec3 &c);
glm::vec3 get_plane_normal(const glm::vec3 &a, const glm::vec3 &b);

inline double get_plane_dir(const Pos2 &a, const Pos2 &b) {
  return a.x * b.y - a.y * b.x;
}
inline double get_plane_dir(const Pos2 &a, const Pos2 &b, const Pos2 &c) {
  return get_plane_dir(b - a, c - a);
}
inline double get_plane_normal(const Pos2 &a, const Pos2 &b) {
  return get_plane_dir(a, b);
}
inline double get_plane_normal(const Pos2 &a, const Pos2 &b, const Pos2 &c) {
  return get_plane_dir(b - a, c - a);
}

bool inside_polygon_unsafe(const std::vector<glm::vec2> &tb,
                           const glm::vec2 &x);
// convex, centered around 0, front is the first vertex above y=0 line, with x>0
bool inside_polygon(const std::vector<glm::vec2> &tb, const glm::vec2 &x,
                    int front = 0);

inline int find_other_point(const PointVec &tb, const double eps = base_eps) {
  FOR (i, 1, tb.size())
    if (glm::length(tb[i] - tb[0]) > eps) return i;
  return -1;
}
inline std::vector<int> find_non_aligned_point(const PointVec &tb,
                                               const double eps = base_eps) {
  int other_pt = find_other_point(tb, eps);
  if (other_pt == -1) return {};
  FOR (j, other_pt + 1, tb.size()) {
    if (!are_aligned(tb[0], tb[other_pt], tb[j])) return { other_pt, j };
  }
  return {};
}

inline std::vector<int> find_non_planar_point(const PointVec &tb) {
  std::vector<int> naligned = find_non_aligned_point(tb);
  if (naligned.empty()) return {};
  Dir plane_dir = get_plane_normal(tb[0], tb[naligned[0]], tb[naligned[1]]);
  FOR (j, naligned.back() + 1, tb.size()) {
    if (std::abs(glm::dot(tb[j] - tb[0], plane_dir)) > base_eps)
      return { naligned[0], naligned[1], j };
  }
  return {};
}

template <class ContainerT>
inline double face_normal_convex2(const ContainerT &face) {
  FOR (i, 1, face.size()) {
    int c = (i + 1) % face.size();
    if (!are_aligned2(face[i - 1], face[i], face[c])) {
      return get_plane_dir(face[c] - face[i], face[i - 1] - face[i]);
    }
  }
  OPA_CHECK(false, face);
  return 0;
}

template <class ContainerT>
inline Pos face_normal_convex(const ContainerT &face) {
  FOR (i, 1, face.size()) {
    int c = (i + 1) % face.size();
    if (!are_aligned(face[i - 1], face[i], face[c])) {
      Pos ndir = get_plane_normal(face[c] - face[i], face[i - 1] - face[i]);
      return ndir;
    }
  }
  OPA_CHECK(false, face);
  return Pos();
}

template <class ContainerT> inline double face_normal2(const ContainerT &face) {

  int mine = std::min_element(ALL(face), GlmVecXCmp<Pos2>()) - face.begin();
  int next = (mine + 1) % face.size();
  int prev = (mine + face.size() - 1) % face.size();
  return (d2_rot_to_angle(face[next] - face[mine]) <
          d2_rot_to_angle(face[prev] - face[mine]))
           ? 1
           : -1;
}

template <class ContainerT> inline Dir face_normal(const ContainerT &face) {
  Dir normal = face_normal_convex(face);
  Point2Vec tb;
  Pos front = face[1] - face[0];
  for (auto &pt : face) tb.push_back(get_plane_coord(pt, normal, front));
  return normal * face_normal2(tb);
}

inline Pos tr_normal(const Triangle &tr) {
  Pos t1, t2, t3;
  std::tie(t1, t2, t3) = tr;
  t2 -= t1;
  t3 -= t1;
  return glm::cross(t2, t3);
}

inline std::vector<Pos> tr_to_vec(const Triangle &tr) {
  return { std::get<0>(tr), std::get<1>(tr), std::get<2>(tr) };
}
inline Triangle vec_to_tr(const std::vector<Pos> &tr) {
  return Triangle{ tr[0], tr[1], tr[2] };
}

inline Triangle realign_tr(const Triangle &dest_tr, const Dir &dir) {
  Dir n1 = tr_normal(dest_tr);
  if (glm::dot(n1, dir) > 0) return dest_tr;
  return Triangle{ std::get<0>(dest_tr), std::get<2>(dest_tr),
                   std::get<1>(dest_tr) };
}

inline Triangle realign_tr(const Triangle &dest_tr, const Triangle &model_tr) {
  return realign_tr(dest_tr, tr_normal(model_tr));
}

inline Point2Vec realign_face(const Point2Vec &face, double dir) {
  double n1 = face_normal2(face);
  if (std::signbit(n1) == std::signbit(dir)) return face;
  Point2Vec res = face;
  std::reverse(ALL(res));
  return res;
}

template <typename ContainerT>
inline bool is_counterclockwise(const ContainerT &face) {
  return face_normal2(face) > 0;
}

template <class ContainerT>
inline ContainerT reorient_counterclockise(const ContainerT &face) {
  if (is_counterclockwise(face)) return face;
  ContainerT res = face;
  std::reverse(ALL(res));
  return res;
}

inline PointVec realign_face(const PointVec &face, const Dir &dir) {
  if (glm::dot(face_normal(face), dir) > 0) return face;
  PointVec res = face;
  std::reverse(ALL(res));
  return res;
}

struct PlaneSpec {
  Pos dir;
  double v;
  Pos default_front;

  PlaneSpec() {}

  PlaneSpec(const Pos &dir, double v) : dir(dir), v(v) {
    default_front = vec_ortho_rand(dir);
  }

  template <class ContainerT>
  static PlaneSpec FromPoints(const ContainerT &points) {
    return PlaneSpec::FromPointAndVec(points[0], face_normal(points));
  }
  static PlaneSpec FromPointAndVec(const Pos &pt, const Pos &dir) {
    Pos ndir = glm::normalize(dir);
    return PlaneSpec(ndir, glm::dot(ndir, pt));
  }
  double get_dist(const Pos &pt) const { return std::abs(get_signed_dist(pt)); }
  double get_signed_dist(const Pos &pt) const { return glm::dot(dir, pt) - v; }
  Dir get_rand_plane_vec() const { return vec_ortho_rand(dir); }

  Pos2 proj(const Pos &pos, const Dir &front) const {
    return get_plane_coord(pos, dir, front);
  }

  Pos lift(const Pos2 &pos, const Dir &front) const {
    Pos y = glm::cross(dir, front);
    return dir * v + pos[0] * front + pos[1] * y;
  }

  Pos2 proj(const Pos &pos) const { return proj(pos, default_front); }

  Pos lift(const Pos2 &pos) const { return lift(pos, default_front); }
  bool coplanar(const PlaneSpec &p, double eps = base_eps) const {
    return coplanar(p.dir, eps) &&
           OPA_FLOAT_EQ(p.v, v, eps);
  }

  bool coplanar(const Dir &p, double eps = base_eps) const {
    return OPA_FLOAT_EQ(glm::dot(p, dir), 0, eps);
  }
};

struct LineSpec {
  Pos a, b;
};
struct LineSpec2D {
  Pos2 a, b;
};

struct HyperPlaneSpec {
  PlaneSpec plane;

  bool is_in_strict(const Pos &pt, const double eps = base_eps) const {
    return plane.get_signed_dist(pt) >= eps;
  }
  bool is_in(const Pos &pt, const double eps = base_eps) const {
    return plane.get_signed_dist(pt) >= -eps;
  }
  bool is_on(const Pos &pt, const double eps = base_eps) const {
    return plane.get_dist(pt) < eps;
  }
};

template <typename VecType> std::vector<VecType> generate_basis();

template <> inline std::vector<Pos2> generate_basis<Pos2>() {
  std::vector<Pos2> res(2);
  res[0] = vec_rand_uni_t<Pos2, 2>();
  res[1] = vec_make_ortho(vec_rand_uni_t<Pos2, 2>(), res[0]);
  return res;
}

template <> inline std::vector<Pos> generate_basis<Pos>() {
  std::vector<Pos> res(3);
  res[0] = vec_rand_uni_t<Pos, 3>();
  res[1] = vec_make_ortho(vec_rand_uni_t<Pos, 3>(), res[0]);
  res[2] = glm::cross(res[0], res[1]);
  return res;
}

template <int N, typename VecType, typename MatType> struct BoxSpec_Gen {

  VecType corner;
  // column = direction
  MatType v;
  typedef BoxSpec_Gen<N, VecType, MatType> SelfType;
  bool in(const VecType &pos) const {
    VecType tmp = (pos - corner) * v;

    REP (i, N) {
      double val = tmp[i] / glm::length2(v[i]);
      if (val < -base_eps || val > 1 + base_eps) return false;
    }
    return true;
  }

  bool in_strict(const VecType &pos) const {
    VecType tmp = (pos - corner) * v;

    REP (i, N) {
      double val = tmp[i] / glm::length2(v[i]);
      if (val < 0 || val > 1) return false;
    }
    return true;
  }

  bool in(const std::tuple<VecType, VecType, VecType> &tup) const {
    return in(std::get<0>(tup)) && in(std::get<1>(tup)) && in(std::get<2>(tup));
  }
  bool in(const std::vector<VecType> &tb) const {
    for (auto &e : tb)
      if (!in(e)) return false;
    return true;
  }
  bool in_strict(const std::vector<VecType> &tb) const {
    for (auto &e : tb)
      if (!in_strict(e)) return false;
    return true;
  }

  bool in_tb(const std::vector<VecType> &tb) const { return in(tb); }

  void set(int dir, const VecType &pos) { v[dir] = pos; }

  double dist2(const VecType &pos) const {
    double d = 0;
    VecType p2 = pos - corner;
    REP (i, N) {
      double l = glm::length(v[i]);
      double x = glm::dot(p2, v[i]) / l;
      if (x < 0)
        ;
      else if (x < l)
        x = 0;
      else
        x -= l;
      d += x * x;
    }

    return d;
  }

  std::vector<VecType> corners() const {
    std::vector<VecType> res;
    REP (i, 1 << N)
      res.push_back(get(i));
    return res;
  }

  std::vector<std::pair<int, int> > edges() const {
    std::vector<std::pair<int, int> > res;
    REP (a, 1 << N)
      REP (toggle, N) {
        int b = a ^ MASK(toggle);
        if (a < b) res.emplace_back(a, b);
      }
    return res;
  }

  VecType center() const { return corner + v * VecType(1. / 2); }

  double dist_to_center(const VecType &pos) const {
    return glm::length2(pos - center());
  }

  double dist_to_face(const VecType &pos) const {
    utils::MinFinder<double> mf;
    VecType p2 = pos - corner;
    REP (i, N) {
      double l = glm::length(v[i]);
      double x = glm::dot(pos, v[i]) / l;
      mf.update(std::abs(x * x));
      mf.update(std::abs((l - x) * (l - x)));
    }
    return mf.get();
  }

  double area() const {
    double res = 1.;
    REP (i, N)
      res *= glm::length(v[i]);
    return res;
  }

  bool almost_empty() const {
    REP (i, 3)
      if (glm::length(v[i]) < base_eps * 10) return true;
    return false;
  }

  SelfType &make_safe(double len) {
    REP (i, N) {
      double length = glm::length(v[i]);
      if (length < len) v[i] *= len / length;
    }
    return *this;
  }

  SelfType &expand(double scale) {
    corner = corner - v * VecType((scale - 1) / 2);
    v = v * scale;
    return *this;
  }

  VecType vec(int dir) const { return v[dir]; }
  double side_len(int dir) const { return glm::length(vec(dir)); }

  VecType get(int coord) const {
    VecType res = corner;
    REP (i, N)
      res += (coord >> i & 1 ? v[i] : VecType(0));
    return res;
  }

  VecType get(const std::vector<double> coords) const {
    VecType res = corner;
    REP (i, N)
      res += coords[i] * v[i];
    return res;
  }

  static SelfType Rand() {
    SelfType res;
    res.corner = vec_rand_uni_t<VecType, N>();
    auto basis = generate_basis<VecType>();
    REP (i, N)
      res.v[i] = basis[i] * unif_distrib(math::common::rng);
    return res;
  }

  double main_axis_dist_normed(const VecType &p, int axis) const {
    VecType p2 = p - this->center();
    return glm::dot(vec(axis), p2) / glm::length2(vec(axis)) * 2;
  }

  std::string str() const {
    if (N == 2) {
      return absl::Substitute(
        "Box2D=(corner=$0, vx=$1, vy=$2, area=$3)", toStr(corner), toStr(v[0]),
        toStr(v[1]), area());
    } else {
      return absl::Substitute(
        "Box=(corner=$0, vx=$1, vy=$2, vz=$3, volume=$4)", toStr(corner),
        toStr(v[0]), toStr(v[1]), toStr(v[2]), area());
    }
  }
  OPA_DECL_COUT_OPERATOR(SelfType);
};

using BoxSpec = BoxSpec_Gen<3, Pos, Mat3>;
using Box2DSpec = BoxSpec_Gen<2, Pos2, Mat2>;

template <int N, typename VecType, int NINF = 100> struct BoxAASpec_Gen {
  typedef BoxAASpec_Gen<N, VecType, NINF> SelfType;
  typedef typename VecType::value_type Type;
  typedef BoxAASpec_Gen<1, glm::tvec1<Type, glm::highp>, NINF> IntervalType;

  VecType low = VecType(std::pow(10, NINF));
  VecType high = VecType(-std::pow(10, NINF));

  Range1D get_range(int dim) const { return Range1D{ low[dim], high[dim] }; }
  IntervalType get_range2(int dim) const {
    IntervalType res;
    res.low.x = low[dim];
    res.high.x = high[dim];
    return res;
  }

  template <typename NVecType, int M>
  BoxAASpec_Gen<M, NVecType, NINF> subbox(const std::array<int, M> &dims) {
    BoxAASpec_Gen<M, NVecType, NINF> res;
    REP (i, M) {
      res.low[i] = low[dims[i]];
      res.high[i] = high[dims[i]];
    }
    return res;
  }

  double min_dist(const VecType &pt) const {
    VecType v;
    REP (j, N)
      v[j] = get_range(j).dist(pt[j]);
    return glm::length(v);
  }

  double min_dist2(const VecType &pt) const {
    VecType v;
    REP (j, N)
      v[j] = get_range(j).dist(pt[j]);
    return glm::length2(v);
  }
  bool empty() const { return low[0] > high[0]; }

  double area() const {
    if (this->degenerate()) return 0.;
    double res = 1.;
    REP (i, N)
      res *= (high[i] - low[i]);
    return res;
  }

  bool degenerate() const { return glm::any(glm::lessThanEqual(high, low)); }

  bool in_fast(const VecType &pos) const {
    return glm::all(glm::greaterThanEqual(pos, low)) &&
           glm::all(glm::lessThanEqual(pos, high));
  }

  bool in(const VecType &pos, double eps = base_eps) const {
    return !glm::any(glm::greaterThanEqual(low - pos, VecType(eps))) &&
           !glm::any(glm::greaterThanEqual(pos - high, VecType(eps)));
  }

  bool almost_empty() const {
    REP (i, N)
      if (glm::length(high[i] - low[i]) <= base_eps * 10) return true;
    return false;
  }

  bool contains(const SelfType &peer) const {
    return glm::all(lessThanEqual(low, peer.low)) &&
           glm::all(greaterThanEqual(high, peer.high));
  }

  bool intersects(const SelfType &peer) const {
    SelfType res;
    VecType nmin = glm::max(low, peer.low);
    VecType nmax = glm::min(high, peer.high);
    return glm::all(glm::lessThanEqual(nmin, nmax));
  }

  SelfType intersection(const SelfType &peer) const {
    SelfType res;
    VecType nmin = glm::max(low, peer.low);
    VecType nmax = glm::min(high, peer.high);
    nmax = glm::max(nmax, nmin);
    return SelfType().update(nmin).update(nmax);
  }

  SelfType &update(const VecType &pos) {
    low = glm::min(low, pos);
    high = glm::max(high, pos);
    return *this;
  }

  SelfType &update(const std::vector<VecType> &v) {
    for (auto &pos : v) {
      this->update(pos);
    }
    return *this;
  }

  SelfType get_union(const SelfType &peer) const {
    SelfType res;
    res.low = glm::min(low, peer.low);
    res.high = glm::max(high, peer.high);
    return res;
  }

  SelfType cut_toward_small(int dim, const Type &cutv) const {
    SelfType res;
    res.low = low;
    res.high = high;
    res.high[dim] = cutv;
    return res;
  }

  SelfType cut_toward_large(int dim, const Type &cutv) const {
    SelfType res;
    res.low = low;
    res.high = high;
    res.low[dim] = cutv;
    return res;
  }

  VecType get(int coord) const {
    VecType res;
    REP (i, N)
      res[i] = (coord >> i & 1 ? high[i] : low[i]);
    return res;
  }

  VecType center() const { return (low + high) / 2; }
  VecType dim() const { return (high - low); }

  BoxAASpec_Gen &from_center_and_dim(const VecType &c, const VecType &d) {

    low = c - d / 2;
    high = low + d;
    return *this;
  }

  BoxAASpec_Gen<N, VecType> &expand(double ratio) {
    from_center_and_dim(center(), dim() * ratio);
    return *this;
  }

  std::string str() const {
    return OPA_STREAM_STR("BoxAA(" << RAW_OPA_DISP_VARS(low, high) << ")");
  }

  OPA_DECL_EQ_OPERATOR(SelfType, low, high);
  OPA_DECL_COUT_OPERATOR(SelfType);
};

using BoxAASpec_Gen3D = BoxAASpec_Gen<3, Pos>;
using BoxAASpec_Gen2D = BoxAASpec_Gen<2, Pos2>;

std::vector<HyperPlaneSpec> get_box_hyperplans(const BoxSpec &box);
BoxSpec get_in_box_space(const BoxSpec &src, const BoxSpec &target);

OPA_NM_MATH_GAME_END

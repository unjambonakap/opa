#pragma once

#include <opa_common.h>
#define GLM_FORCE_SWIZZLE
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/ext.hpp>
#include <glm/glm.hpp>
#include <glm/geometric.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/norm.hpp>

#include <opa/utils/hash.h>
#include <opa/algo/base.h>
#include <opa/algo/graph.h>
#include <opa/utils/DataStruct.h>
#include <opa/utils/buffer_reader.h>
#include <opa/utils/buffer_writer.h>
#include <opa/utils/misc.h>

#define OPA_NM_MATH_GAME                                                       \
  namespace opa {                                                              \
  namespace math {                                                             \
  namespace game {
#define OPA_NM_MATH_GAME_END                                                   \
  }                                                                            \
  }                                                                            \
  }

namespace glm {
std::ostream &operator<<(std::ostream &os, const glm::vec2 &vec);
std::ostream &operator<<(std::ostream &os, const glm::ivec2 &vec);
std::ostream &operator<<(std::ostream &os, const glm::ivec3 &vec);
std::ostream &operator<<(std::ostream &os, const glm::vec4 &vec);
std::ostream &operator<<(std::ostream &os, const glm::vec3 &vec);
std::ostream &operator<<(std::ostream &os, const glm::quat &quat);
std::ostream &operator<<(std::ostream &os, const glm::mat2 &mat);
std::ostream &operator<<(std::ostream &os, const glm::mat3 &mat);
std::ostream &operator<<(std::ostream &os, const glm::mat4 &mat);
} // namespace glm
#define OPA_FLOAT_EQ(a, b, eps) (fabs((a) - (b)) <= eps)
#define OPA_FLOAT_LE(a, b, eps) ((a)-eps <= (b))
#define OPA_FLOAT_GE(a, b, eps) ((b)-eps <= (a))
#define OPA_FLOAT_LT(a, b, eps) ((a) + eps < (b))
#define OPA_FLOAT_GT(a, b, eps) ((b) + eps < (a))

#define OPA_FLOAT_BETWEEN(a, l, h, eps)                                        \
  (OPA_FLOAT_LE(l, a, eps) && OPA_FLOAT_LE(a, h, eps))
#define OPA_VEC_EQ(a, b, eps)                                                  \
  OPA_FLOAT_EQ(glm::distance((a), (b)), 0, (a).length() * eps)

#define OPA_DRANGE_REMAP(L, H, v) ((atan(v) / OPA_PI + 1. / 2) * ((H) - (L)) + (L))
#define OPA_DRANGE_REMAP_CENTERED(D, v) OPA_DRANGE_REMAP(-D, D, v)

template <class T> T deg_to_rad(T a) { return a * 2 * OPA_PI / 360; }
template <class T> T rad_to_deg(T a) { return a * 360 / 2 / OPA_PI; }
static const std::vector<glm::vec2> UnitSquare = {
  glm::vec2(0, 0), glm::vec2(1, 0), glm::vec2(0, 1), glm::vec2(1, 1)
};

template <class T> struct GlmVecXYCmp {
  GlmVecXYCmp(float eps = 0) { this->eps = eps; }
  float eps;
  bool operator()(const T &a, const T &b) const {
    if (OPA_FLOAT_EQ(a.x - b.x, 0, eps)) return OPA_FLOAT_LT(a.y, b.y, eps);
    return OPA_FLOAT_LT(a.x, b.x, eps);
  }
};

template <class T> struct GlmVecXCmp {
  GlmVecXCmp(float eps = 0) { this->eps = eps; }
  float eps;
  bool operator()(const T &a, const T &b) const {
    return OPA_FLOAT_LT(a.x, b.x, eps);
  }
};

template <class T> struct GlmVecYCmp {
  GlmVecYCmp(float eps = 0) { this->eps = eps; }
  float eps;
  bool operator()(const T &a, const T &b) const {
    return OPA_FLOAT_LT(a.y, b.y, eps);
  }
};

template <class T> struct GlmVecZCmp {
  GlmVecZCmp(float eps = 0) { this->eps = eps; }
  float eps;
  bool operator()(const T &a, const T &b) const {
    return OPA_FLOAT_LT(a.z, b.z, eps);
  }
};

struct vec_cmp {

  vec_cmp(float eps, int pair_cmp_type = 0) {
    m_eps = eps;
    m_pair_cmp_type = pair_cmp_type;
  }

  template <class A> bool operator()(const A &a, const A &b) const {
    return a < b;
  }
  template <class A, class B>
  bool operator()(const std::pair<A, B> &a, const std::pair<A, B> &b) const {
    if (m_pair_cmp_type == 0) {
      bool tmp1 = this->operator()(a.ST, b.ST);
      if (tmp1) return true;
      tmp1 = this->operator()(b.ST, a.ST);
      if (tmp1) return false;
      return this->operator()(a.ND, b.ND);
    } else if (m_pair_cmp_type == 1) {
      return this->operator()(a.ST, b.ST);
    } else {

      return this->operator()(a.ND, b.ND);
    }
  }

private:
  float m_eps;
  int m_pair_cmp_type;
};

struct vec_cmp_eq {

  vec_cmp_eq(float eps, int pair_cmp_type = 0) {
    m_eps = eps;
    m_pair_cmp_type = pair_cmp_type;
  }

  template <class A> bool operator()(const A &a, const A &b) const {
    return a == b;
  }
  template <class A, class B>
  bool operator()(const std::pair<A, B> &a, const std::pair<A, B> &b) const {
    if (m_pair_cmp_type == 0) {
      return this->operator()(a.ST, b.ST) && this->operator()(a.ND, b.ND);
    } else if (m_pair_cmp_type == 1) {
      return this->operator()(a.ST, b.ST);
    } else {

      return this->operator()(a.ND, b.ND);
    }
  }

private:
  float m_eps;
  int m_pair_cmp_type;
};

#if !defined(SWIG)

template <>
inline bool vec_cmp::operator()<float>(const float &a, const float &b) const {
  return OPA_FLOAT_LT(a, b, m_eps);
}

#define DO_ONE_CMP(coord)                                                      \
  if (!OPA_FLOAT_EQ(a.coord, b.coord, m_eps)) return a.coord < b.coord;
template <>
inline bool vec_cmp::operator()<glm::vec2>(const glm::vec2 &a,
                                           const glm::vec2 &b) const {
  DO_ONE_CMP(x);
  DO_ONE_CMP(y);
  return 0;
}

template <>
inline bool vec_cmp::operator()<glm::vec3>(const glm::vec3 &a,
                                           const glm::vec3 &b) const {
  DO_ONE_CMP(x);
  DO_ONE_CMP(y);
  DO_ONE_CMP(z);
  return 0;
}

template <>
inline bool vec_cmp_eq::operator()<float>(const float &a,
                                          const float &b) const {
  return OPA_FLOAT_EQ(a, b, m_eps);
}

template <>
inline bool vec_cmp_eq::operator()<glm::vec2>(const glm::vec2 &a,
                                              const glm::vec2 &b) const {
  return OPA_FLOAT_EQ(a.x, b.x, m_eps) && OPA_FLOAT_EQ(a.y, b.y, m_eps);
  ;
}

template <>
inline bool vec_cmp_eq::operator()<glm::vec3>(const glm::vec3 &a,
                                              const glm::vec3 &b) const {
  return OPA_FLOAT_EQ(a.x, b.x, m_eps) && this->operator()(a.yz(), b.yz());
}

#endif

#undef DO_ONE_CMP

OPA_NAMESPACE_DECL1(opa)
const double inf = 1e100;

typedef glm::ivec2 IPos2;
typedef glm::ivec3 IPos3;

typedef glm::vec3 Pos;
typedef glm::vec3 Dir;
typedef glm::vec4 Pos4;
typedef glm::quat Rot;
typedef glm::vec2 Pos2;
typedef glm::vec2 TexUV;
typedef glm::mat4 Mat4;
typedef glm::mat3 Mat3;
typedef glm::mat2 Mat2;
typedef glm::mat2 Rot2;
typedef std::vector<Pos> PointVec;
typedef std::vector<Pos2> Point2Vec;

typedef opa::utils::IdType GraphObjId;
typedef GraphObjId EdgeId;
typedef GraphObjId VertexId;
typedef GraphObjId FaceId;
typedef std::tuple<Pos, Pos, Pos> Triangle;
typedef std::array<Pos2, 3> Tr2D;
typedef std::array<Pos, 3> Tr3D;
typedef std::vector<int> FaceIndex;
typedef std::vector<FaceIndex> FaceIndexList;

const Rot kIdentRot = Rot(1., 0, 0, 0);
const Rot2 kIdentRot2 = Rot2(1.);

struct SphereSpec {
  Pos center;
  double radius;
};

inline double clamp(double x, double a, double b) {
  if (x < a) return a;
  if (x > b) return b;
  return x;
}

template <typename T, typename U = T>
std::vector<U> ApplyOp(const std::vector<T> &tb,
                       const std::function<U(const T &t)> func) {
  std::vector<U> res;
  for (const auto &e : tb) res.push_back(func(e));
  return res;
}

template <class T, class U>
std::vector<U> ApplyOp(const std::vector<T> &tb, U (*func)(const T &t)) {
  std::vector<U> res;
  for (const auto &e : tb) res.push_back(func(e));
  return res;
}

OPA_NAMESPACE_DECL1_END

namespace std {

#if !defined(SWIG)
inline opa::utils::BufferWriter &operator<<(opa::utils::BufferWriter &writer,
                                            const opa::Pos &pos) {
  writer.write<float>(pos.x);
  writer.write<float>(pos.y);
  writer.write<float>(pos.x);
  return writer;
}

template <> struct hash<opa::Pos4> {
  typedef opa::Pos4 argument_type;
  typedef std::size_t result_type;

  result_type operator()(argument_type const &s) const {
    return opa::utils::do_hash(std::string_view((const char *)&s[0], 4 * sizeof(argument_type::value_type)));
  }
};

template <> struct hash<opa::Pos2> {
  typedef opa::Pos2 argument_type;
  typedef std::size_t result_type;

  result_type operator()(argument_type const &s) const {
    return opa::utils::do_hash(std::string_view((const char *)&s[0], 2 * sizeof(argument_type::value_type)));
  }
};

template <> struct hash<opa::Pos> {
  typedef opa::Pos argument_type;
  typedef std::size_t result_type;

  result_type operator()(argument_type const &s) const {
    return opa::utils::do_hash(std::string_view((const char *)&s[0], 3 * sizeof(argument_type::value_type)));
  }
};
#endif
} // namespace std

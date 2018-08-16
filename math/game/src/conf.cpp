#include "opa/math/game/conf.h"

#include <opa/math/common/Utils.h>
#include <opa/utils/string.h>

using namespace std;
static const float eps = 1e-6;

namespace glm {
std::ostream &operator<<(std::ostream &os, const glm::ivec2 &vec) {
  os << opa::utils::stdsprintf("(%d,%d)", vec.x, vec.y);
  return os;
}

std::ostream &operator<<(std::ostream &os, const glm::ivec3 &vec) {
  os << opa::utils::stdsprintf("(%d,%d,%d)", vec.x, vec.y, vec.z);
  return os;
}

std::ostream &operator<<(std::ostream &os, const glm::vec2 &vec) {
  os << opa::utils::stdsprintf("(%f,%f)", vec.x, vec.y);
  return os;
}

std::ostream &operator<<(std::ostream &os, const glm::vec3 &vec) {
  os << opa::utils::stdsprintf("(%f,%f,%f)", vec.x, vec.y, vec.z);
  return os;
}

std::ostream &operator<<(std::ostream &os, const glm::vec4 &vec) {
  os << opa::utils::stdsprintf("(%f,%f,%f,%f)", vec.x, vec.y, vec.z, vec.w);
  return os;
}

std::ostream &operator<<(std::ostream &os, const glm::quat &quat) {
  os << opa::utils::stdsprintf("(%f,%f,%f,%f)", quat.w, quat.x, quat.y, quat.z);
  return os;
}

std::ostream &operator<<(std::ostream &os, const glm::mat2 &mat) {
  std::string res;
  REP (i, mat.length()) {
    REP (j, mat[0].length())
      res += toStr(mat[j][i]) + "\t";
    res += "\n";
  }
  os << res;
  return os;
}

std::ostream &operator<<(std::ostream &os, const glm::mat3 &mat) {
  std::string res;
  REP (i, mat.length()) {
    REP (j, mat[0].length())
      res += toStr(mat[j][i]) + "\t";
    res += "\n";
  }
  os << res;
  return os;
}

std::ostream &operator<<(std::ostream &os, const glm::mat4 &mat) {
  std::string res;
  REP (i, mat.length()) {
    REP (j, mat[0].length())
      res += toStr(mat[j][i]) + "\t";
    res += "\n";
  }
  os << res;
  return os;
}
}

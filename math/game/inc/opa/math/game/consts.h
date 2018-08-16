#pragma once

#include <opa/math/game/conf.h>
#include <opa/or/manifold.h>

namespace opa {
constexpr double base_eps = 1e-6;
}

OPA_NM_MATH_GAME

const glm::vec2 vec2_x(1, 0);
const glm::vec2 vec2_y(0, 1);
const Mat2 vec2_mat = { vec2_x, vec2_y };
const glm::vec2 vec2_1(1, 1);

const Pos vec_eps_base(base_eps, base_eps, base_eps);
const glm::vec3 vec_x(1, 0, 0);
const glm::vec3 vec_y(0, 1, 0);
const glm::vec3 vec_z(0, 0, 1);
const glm::vec3 vec_0(0, 0, 0);
const glm::vec3 vec_1(1, 1, 1);
const std::array<glm::vec3, 3> vec_tb = { vec_x, vec_y, vec_z };
const Mat3 vec_mat = { vec_x, vec_y, vec_z };


OPA_NM_MATH_GAME_END

#pragma once

#include <opa/math/game/base.h>
#include <opa/math/game/quat.h>

OPA_NM_MATH_GAME

class BoxAASpec : public BoxAASpec_Gen3D {
public:
  BoxSpec to_box() const {
    if (empty()) return BoxSpec();
    Pos diff = high - low;
    return BoxSpec{ low, glm::diagonal3x3(diff) };
  }

  static BoxAASpec FromBoxSpace(const BoxSpec &box) {
    BoxAASpec res;
    Rot rot = rot_to_box_space(box);
    res.low = Pos(0);
    res.high = rot * (box.get(7) - box.corner);
    return res;
  }
};

class Box2DAASpec : public BoxAASpec_Gen2D {
public:
  Box2DSpec to_box() const {
    if (empty()) return Box2DSpec();
    Pos2 diff = high - low;
    return Box2DSpec{ low, glm::diagonal2x2(diff) };
  }
};
OPA_NM_MATH_GAME_END

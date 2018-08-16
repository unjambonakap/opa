#pragma once

namespace opa {
namespace math {
namespace common {
class Float;
}
}
}

namespace std {
opa::math::common::Float ceil(const opa::math::common::Float &a);
opa::math::common::Float floor(const opa::math::common::Float &a);
opa::math::common::Float round(const opa::math::common::Float &a);

opa::math::common::Float sqrt(const opa::math::common::Float &a);
opa::math::common::Float cos(const opa::math::common::Float &a);
opa::math::common::Float sin(const opa::math::common::Float &a);
opa::math::common::Float tan(const opa::math::common::Float &a);

opa::math::common::Float atan2(const opa::math::common::Float &a,
                               const opa::math::common::Float &b);
opa::math::common::Float acos(const opa::math::common::Float &a);
opa::math::common::Float asin(const opa::math::common::Float &a);
opa::math::common::Float atan(const opa::math::common::Float &a);
opa::math::common::Float abs(const opa::math::common::Float &a);
opa::math::common::Float pow(const opa::math::common::Float &a,
                             const opa::math::common::Float &b);
}

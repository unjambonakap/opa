#pragma once

#include <opa_common.h>
#include <opa/math/common/bignum.h>
#include <opa/math/common/Utils.h>

OPA_NAMESPACE_DECL2(opa, crypto)

class Dlp {
  public:
    void setup(const opa::math::common::bignum &xn,
               const opa::math::common::bignum &xorder,
               const opa::math::common::bignum &xg,
               const opa::math::common::bignum &xy);

    void setup_main(const opa::math::common::bignum &xn,
                    const opa::math::common::bignum &xg,
                    const opa::math::common::bignum &xy);

    opa::math::common::bignum do_one(const opa::math::common::bignum &suborder,
                                     u32 suborder_pw);

    void compute_result(const std::vector<opa::math::common::bignum> &subres);
    bool check(const OPA_BG &res) const;

    opa::math::common::bignum solve();

    opa::math::common::bignum n;
    opa::math::common::bignum order;
    opa::math::common::bignum g;
    opa::math::common::bignum ig;
    opa::math::common::bignum y;
    opa::math::common::BGFactors factors;

    opa::math::common::bignum ans;
};

OPA_NAMESPACE_DECL2_END

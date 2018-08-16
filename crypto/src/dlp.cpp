#include "dlp.h"

#include <opa/threading/runner.h>
#include <opa/math/common/GF_pBN.h>

using namespace opa::math::common;
using namespace std;

OPA_NAMESPACE_DECL2(opa, crypto)
void Dlp::setup_main(const opa::math::common::bignum &xn,
                     const opa::math::common::bignum &xg,
                     const opa::math::common::bignum &xy) {

    n = xn;
    g = xg;
    ig = g.inv(n);
    y = xy;
    BGFactors orig_factors = factor_large(xn - 1);

    GF_pBN f(xn);
    order = f.compute_order(g, orig_factors);
    puts("GOT ORDER >> ");
    order.disp();
    bignum tmp = order;

    factors.clear();
    for (const auto &x : orig_factors) {
        u32 now = x.ND;
        REP(i, x.ND) {
            if (tmp % x.ST != 0) {
                now = i;
                break;
            }
            tmp /= x.ST;
        }
        if (now > 0)
            factors.pb(MP(x.ST, now));
    }
    puts("FACTORS> >>");
    out(factors);

    cout << "mod=" << n << endl;
    cout << "g=" << g << endl;
    cout << "y=" << y << endl;
    cout << "order=" << order << endl;
}

void Dlp::setup(const opa::math::common::bignum &xn,
                const opa::math::common::bignum &xorder,
                const opa::math::common::bignum &xg,
                const opa::math::common::bignum &xy) {
    n = xn;
    g = xg;
    ig = g.inv(n);
    y = xy;
    order = xorder;

    cout << "mod=" << n << endl;
    cout << "g=" << g << endl;
    cout << "y=" << y << endl;
    cout << "order=" << order << endl;

}

opa::math::common::bignum Dlp::do_one(const opa::math::common::bignum &suborder,
                                      u32 suborder_pw) {
    usleep(1e5);
    u32 ub = suborder.sqrt().getu32(); // otherwise, you aint there boy
    u32 rem = ((suborder + ub - 1) / ub).getu32();
    bignum res = 0;
    bignum cur = y;
    bignum sorder_q = 1;
    bignum igpw = ig;
    bignum mul = g.faste(order / suborder, n);
    bignum imul = ig.faste(order / suborder, n);

    unordered_map<bignum, u32> h;
    bignum tmp = 1;
    bignum step = mul.faste(rem, n);
    REP(j, ub) {
        h[tmp] = j;
        tmp = tmp * step % n;
    }

    printf("DO ONE ub=%d, rem=%d\n", ub, rem);
    g.disp();
    suborder.disp();
    REP(i, suborder_pw) {
        bignum next_sq = sorder_q * suborder;

        puts("");
        printf("ON ROUND %d >> ", i);
        cur.disp();
        bignum target = cur.faste(order / next_sq, n);
        puts("TARGET >");
        target.disp();
        sorder_q.disp();

        bool fd = false;
        REP(j, rem) {
            if (h.count(target)) {
                fd = true;
                bignum add = rem;
                bignum pw = h[target];
                printf("FIND AT %d >> pw=%d\n", j, pw.getu32());
                res.disp();
                add *= pw * sorder_q;
                cur = cur * igpw.faste(pw * rem, n) % n;
                res += add;
                add.disp();
                res.disp();
                break;
            }
            cur = cur * igpw % n;
            res += sorder_q;
            target = target * imul % n;
        }
        OPA_ASSERT(fd, "bad dlp");
        sorder_q = next_sq;
        igpw = igpw.faste(suborder, n);
    }
    return res;
}

void Dlp::compute_result(const std::vector<opa::math::common::bignum> &subres) {
    vector<pair<bignum, bignum> > lst;

    REP(i, subres.size()) {
        bignum tmp = factors[i].ST.pow(factors[i].ND);
        lst.pb(MP(subres[i], tmp));
    }
    ans = crt_coprime(lst);
    ans.disp();
}

opa::math::common::bignum Dlp::solve() {
    vector<bignum> tb;
    REP(i, factors.size())
    tb.pb(do_one(factors[i].ST, factors[i].ND));
    compute_result(tb);
    return ans;
}

bool Dlp::check(const OPA_BG &res) const { return t_faste(g, res, n) == y; }

OPA_NAMESPACE_DECL2_END

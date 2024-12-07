#include <range/v3/all.hpp>


#define STD_FUNC1(a) [&](auto x) { return a(x); }
#define STD_FUNC2(a) [&](auto x0, auto x1) { return a(x0, x1); }
#define STD_TSFX(content) (LQ::transform([&](const auto &x) { return (content); }))
#define STD_TSFY(content) (LQ::transform([&](const auto &y) { return (content); }))
#define STD_TSFXNORET(content) (LQ::transform([&](const auto &x) { (content); }))
#define STD_FOREACHX(entries, content) (QQ::for_each(entries, [&](const auto &x) { (content); }))
#define STD_FOREACHTUPLE2(entry,content)                                                                 \
  (QQ::for_each(entry, [&](const auto &tx) {                                                                     \
    auto &p0 = std::get<0>(tx);                                                                    \
    auto &p1 = std::get<1>(tx);                                                                    \
    content;                                                                                       \
    ;                                                                                              \
  }))
#define STD_FUNCXY(content) ([&](const auto &x, const auto &y) { return (content); })
#define STD_FUNCAB(content) ([&](const auto &a, const auto &b) { return (content); })
#define STD_FUNCX(content) ([&](const auto &x) { return (content); })
#define STD_FUNCY(content) ([&](const auto &y) { return (content); })

#define STD_TSFTUPLE(content)                                                                      \
  (LQ::transform([&](auto tx) {                                                                    \
    auto &p0 = std::get<0>(tx);                                                                    \
    auto &p1 = std::get<1>(tx);                                                                    \
    return (content);                                                                              \
  }))
#define STD_TSFTUPLE3(content)                                                                     \
  (LQ::transform([&](auto tx) {                                                                    \
    auto &p0 = std::get<0>(tx);                                                                    \
    auto &p1 = std::get<1>(tx);                                                                    \
    auto &p2 = std::get<2>(tx);                                                                    \
    return (content);                                                                              \
  }))

#define V_t(t) std::vector<t>
#define VV_t(t) V_t(std::vector<t>)
#define VVV_t(t) VV_t(std::vector<t>)
#define STD_VEC QQ::to<std::vector>
#define STD_SET QQ::to<std::set>
#define STD_SETT(t) QQ::to<std::set<t>>
#define STD_VECT(t) QQ::to<std::vector<t> >
#define STD_REPEAT(v, cnt) (LQ::repeat(v) | LQ::take(cnt))
#define STD_SINGLEV(v) STD_REPEAT(v, 1)
#define STD_SINGLEVT(v,t) STD_REPEAT(v, 1) | STD_VECT(t)
#define STD_RANGE(a, b) LQ::iota(s32(a), s32(b))
namespace LQ = ranges::views;
namespace QQ = ranges;

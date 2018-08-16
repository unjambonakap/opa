#pragma once
#include <opa/math/common/GF_p.h>
#include <opa/math/common/Matrix.h>
#include <opa/math/common/Types.h>
#include <opa/utils/misc.h>
#include <opa_common.h>

OPA_NAMESPACE(opa, math, common)

class RangeCoverage : public opa::utils::Initable,
                      public opa::utils::ProtobufParams {
public:
  RangeCoverage(int high) { init(0, high); }
  RangeCoverage(int low, int high) { init(low, high); }
  RangeCoverage() { init(); }
  RangeCoverage &operator=(const RangeCoverage &x) {
    Initable::operator=(x);
    OPA_COPY_OP(x, m_low, m_high, m_cover, m_unique, m_ignore_unique,
                m_dynamic);
    return *this;
  }
  template <class T> RangeCoverage(std::initializer_list<T> lst) {
    init();
    add_range(lst);
  }

  static RangeCoverage From_binary_vec(const std::vector<u32> &tb) {
    RangeCoverage res;
    REP (i, tb.size())
      if (tb[i])
        res.add(i);
    return res;
  }

  static RangeCoverage From(u64 val, int sz) {
    RangeCoverage res(0, sz);
    REP (i, sz)
      if (val >> i & 1)
        res.add(i);
    return res;
  }

  void from(u64 val) {
    int pos = 0;
    for (int pos = 0; val; val >>= 1, ++pos)
      if (val & 1)
        add(pos);
  }

  template <class T> T to() const {
    T res = 0;
    for (auto v : all())
      res |= 1ull << v;
    return res;
  }

  virtual void init(int low, int high) {
    Initable::init();
    m_low = low;
    m_high = high;
    m_cover.clear();
    m_unique = true;
    m_dynamic = false;
  }

  int dot(const RangeCoverage &peer) const { return (*this & peer).size() & 1; }

  virtual void init() override {
    Initable::init();
    m_low = -1;
    m_high = -1;
    m_unique = false;
    m_cover.clear();
    m_dynamic = true;
  }
  RangeCoverage &clear() {
    m_cover.clear();
    m_unique = true;
    return *this;
  }

  RangeCoverage &filter_less(int v) {
    while (m_cover.size() && *m_cover.begin() < v)
      m_cover.erase(m_cover.begin());
    return *this;
  }

  RangeCoverage &filter_greater(int v) {
    while (m_cover.size() && *m_cover.rbegin() > v)
      m_cover.erase(--m_cover.end());

    return *this;
  }

  RangeCoverage &add(int v) {
    check_init();
    OPA_CHECK0(m_dynamic || (v >= m_low && v < m_high));

    if (m_cover.count(v))
      m_unique = false;
    else
      m_cover.insert(v);
    return *this;
  }

  RangeCoverage &add(const RangeCoverage &r) {
    for (auto x : r.m_cover)
      add(x);
    return *this;
  }

  template <class T> RangeCoverage &add_lst(const T &lst) {
    for (auto x : lst)
      add(x);
    return *this;
  }

  template <class T> RangeCoverage &add_range(std::initializer_list<T> lst) {
    for (auto i : lst)
      add(i);
    return *this;
  }

  RangeCoverage sub(int v) {
    RangeCoverage res = *this;
    return res.ssub(v);
  }

  RangeCoverage sub(const RangeCoverage &r) {
    RangeCoverage res = *this;
    return res.ssub(r);
  }

  RangeCoverage &ssub(int v) {
    m_cover.erase(v);
    return *this;
  }

  RangeCoverage &ssub(const RangeCoverage &r) {
    for (auto x : r.m_cover)
      ssub(x);
    return *this;
  }

  RangeCoverage &do_xor(const RangeCoverage &r) {
    std::set<int> tmp;
    auto walker =
      opa::utils::MergeWalker<std::set<int> >({ &m_cover, &r.m_cover });
    FE (it, walker) {
      if (it.is_active(0) ^ it.is_active(1))
        tmp.insert(*it);
    }
    m_cover = tmp;
    return *this;
  }

  RangeCoverage &do_and(const RangeCoverage &r) {
    std::set<int> tmp;
    auto walker =
      opa::utils::MergeWalker<std::set<int> >({ &m_cover, &r.m_cover });
    FE (it, walker) {
      if (it.is_active(0) && it.is_active(1))
        tmp.insert(*it);
    }
    m_cover.swap(tmp);
    return *this;
  }

  RangeCoverage operator&(const RangeCoverage &peer) const {
    RangeCoverage res = *this;
    res.do_and(peer);
    return res;
  }

  RangeCoverage operator^(const RangeCoverage &peer) const {
    RangeCoverage res = *this;
    res.do_xor(peer);
    return res;
  }

  RangeCoverage operator+(const RangeCoverage &peer) const {
    RangeCoverage res = *this;
    res.add(peer);
    return res;
  }

  RangeCoverage operator-(const RangeCoverage &peer) const {
    RangeCoverage res = *this;
    res.sub(peer);
    return res;
  }

  bool in_range(int v) const { return m_dynamic || (v >= m_low && v < m_high); }
  bool in(int v) const { return m_cover.count(v); }
  bool contains(const RangeCoverage &peer) const {
    for (auto x : peer.all())
      if (!in(x))
        return false;
    return true;
  }

  bool has_one(const RangeCoverage &peer) const {
    for (auto x : peer.all())
      if (in(x))
        return true;
    return false;
  }

  RangeCoverage &sshift(int v) {
    std::set<int> tmp;
    for (auto x : m_cover)
      if (in_range(x + v))
        tmp.insert(x + v);
    m_cover.swap(tmp);
    return *this;
  }

  RangeCoverage shift(int v) const {
    RangeCoverage res = *this;
    res.m_dynamic = true;
    res.sshift(v);
    return res;
  }

  RangeCoverage &add_interval(int a, int n) {
    REP (i, n)
      add(a + i);
    return *this;
  }

  bool covered() const {
    OPA_CHECK0(!m_dynamic);
    return m_cover.size() == m_high - m_low;
  }
  OPA_ACCESSOR_R(bool, m_unique, unique);
  OPA_ACCESSOR_R2(u32, m_cover.size(), size);

  friend std::ostream &operator<<(std::ostream &os, const RangeCoverage &r) {
    os << RAW_OPA_DISP_VARS(r.m_cover);
    return os;
  }

  std::vector<int> to_vec() const { return std::vector<int>(ALL(m_cover)); }
  std::vector<u32> to_binary_vec(int sz) const {
    std::vector<u32> res(sz);
    REP (i, sz)
      res[i] = in(i);
    return res;
  }

  Matrix<u32> to_mat(int sz) const {
    Matrix<u32> res;
    res.initialize(&GF2, 1, sz);
    return res.fromvec(to_binary_vec(sz), opa::math::common::VecType::Row);
  }

  bool operator==(const RangeCoverage &r) const {
    return m_low == r.m_low && m_high == r.m_high && m_unique == r.m_unique &&
           m_cover == r.m_cover;
  }
  bool operator!=(const RangeCoverage &r) const { return !(*this == r); }

  bool operator<(const RangeCoverage &r) const {
    if (m_low != r.m_low)
      return m_low < r.m_low;
    if (m_high != r.m_high)
      return m_high < r.m_high;
    if (!m_ignore_unique && m_unique != r.m_unique)
      return m_unique < r.m_unique;
    if (m_cover.size() != r.m_cover.size())
      return m_cover.size() < r.m_cover.size();
    return m_cover < r.m_cover;
  }

  OPA_ACCESSOR_R(std::set<int>, m_cover, cover);
  OPA_ACCESSOR_R(std::set<int>, m_cover, all);

  OPA_TGEN_IMPL_INIT2(m_low, m_high, m_cover, m_unique, m_ignore_unique,
                      m_dynamic);

private:
  int m_low, m_high;
  std::set<int> m_cover;
  bool m_unique;
  bool m_ignore_unique = false;
  bool m_dynamic;
};

OPA_NAMESPACE_END(opa, math, common)

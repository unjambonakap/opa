#pragma once
#include <opa_common.h>
#include <opa/utils/misc.h>
#include <opa/utils/string.h>
#include <opa/math/common/Utils.h>

// TODO: replace this with MulMatrixF2

OPA_NAMESPACE(opa, math, common)
class RangeCoverage;

class Value : public opa::utils::Initable {
public:
  static Value From(u64 v, s32 sz) {
    Value res(sz);
    res.from(v);
    return res;
  }
  static Value FromBytes(const RangeCoverage &range, opa::stolen::StringRef s);

  template <class T> static Value From2(const T *ptr, int n) {
    int sz = sizeof(T) * 8;
    Value res(n * sz);
    REP (i, n)
      res.from(ptr[i], sz * i, sz);
    return res;
  }

  Value() {}
  Value(s32 sz) { init(sz); }

  Value &operator=(const Value &v) {
    Initable::operator=(v);
    m_bits = v.m_bits;
    m_size = v.m_size;
    return *this;
  }
  bool operator==(const Value &v) const { OPA_EQ_OP(v, m_size, m_bits); }

  static Value Rand(s32 sz) {
    Value res(sz);
    while (sz) {
      int take = std::min(sz, 32);
      u32 cur = rng();
      REP (j, take)
        res.setb(--sz, cur >> j & 1);
    }
    return res;
  }

  int dot(const std::vector<u32> &x) const {
    int res = 0;
    for (auto i : x)
      res ^= getb(i);
    return res;
  }

  int dot(const Value &x) const {
    int res = 0;
    REP (i, m_size)
      res ^= getb(i) & x.getb(i);
    return res;
  }

  void sxor(const RangeCoverage &x);
  void sxor(const Value &x) {
    REP (i, m_size)
      if (x.getb(i))
        toggleb(i);
  }

  bool isz() const {
    REP (i, m_size)
      if (getb(i))
        return false;
    return true;
  }

  virtual void init(int size) {
    Initable::init();
    m_size = size;
    m_bits.resize(size);
  }

  void toggleb(int pos) { setb(pos, getb(pos) ^ 1); }

  void setb(int pos, int v) {
    check_range(pos);
    OPA_CHECK0(v == 0 || v == 1);
    getb(pos) = v;
  }

  template <class T> void from(T a, u32 off = 0, s32 sz = -1) {
    norm_size(off, sz);
    OPA_CHECK0(off + sz <= m_size);
    REP (i, sz) {
      setb(off + i, a & 1);
      a >>= 1;
    }
  }

  template <class T> T to(u32 off = 0, s32 sz = -1) const {
    T res = 0;
    norm_size(off, sz);
    REP (i, sz)
      res = res << 1 | getb(off + sz - 1 - i);
    return res;
  }

  Value get(u32 off, s32 sz = -1) const {
    norm_size(off, sz);
    Value res;
    res.init(sz);
    REP (i, sz)
      res.setb(i, getb(off + i));
    return res;
  }

  int dot(const RangeCoverage &x) const;
  Value and_(const RangeCoverage &x) const;

  void set(const Value &src, u32 off = 0, s32 sz = -1) {
    norm_size(off, sz);
    sz = std::min(sz, src.size());
    REP (i, sz)
      setb(off + i, src.getb(i));
  }

  void norm_size(u32 off, s32 &sz) const {
    if (sz == -1)
      sz = size() - off;
  }

  void check_range(int pos) const {
    OPA_CHECK0(pos >= 0 && pos < m_bits.size());
  }
  u8 getb(int pos) const {
    check_range(pos);
    return m_bits[pos];
  }
  u8 &getb(int pos) {
    check_range(pos);
    return m_bits[pos];
  }

  s32 size() const { return m_size; }

  std::string str() const { return utils::SPrintf("%Lx", to<u64>()); }

  friend std::ostream &operator<<(std::ostream &os, const Value &v) {
    os << v.str();
    return os;
  }

  bool operator<(const Value &v) const {
    if (m_size != v.m_size)
      return m_size < v.m_size;
    // sort by lsb first!
    return m_bits < v.m_bits;
  }

  std::string to_bytes() const {
    std::string res((m_size + 7) >> 3, 0);
    REP (i, m_size)
      res[i >> 3] |= m_bits[i] << (i & 7);
    return res;
  }

  s32 m_size;
  std::vector<u8> m_bits;
};
OPA_NAMESPACE_END(opa, math, common)

#pragma once

#include <opa/math/common/GF_p.h>
#include <opa/math/common/GF_q.h>
#include <opa/math/common/Matrix.h>
#include <opa/math/common/Poly.h>
#include <opa/math/common/PolyRing.h>
#include <opa/math/common/Types.h>
#include <opa/math/common/Utils.h>
#include <opa/math/common/base.h>
#include <opa/math/common/base/range_coverage.h>

OPA_NM_MATH_COMMON

typedef std::array<u64, 2> Gf128_t2;
typedef opa::math::common::Poly<u32> Gf128_t;
extern opa::math::common::GF_q<u32> gf128;
extern Gf128_t gcm_poly;

template <bool Left> class MulMatrixF2 {
public:
  MulMatrixF2() {}
  void init(int N, int M) {
    this->N = N;
    this->M = M;
    N2 = (N + 63) / 64;
    M2 = (M + 63) / 64;
    if (!Left)
      tb.resize(N2 * M);
    else
      tb.resize(M2 * N);
    setz();
  }

  void init_from_gf128_t2(const Gf128_t2 &v) {
    init(128, 1);
    REP (i, N)
      set(i, 0, v[i / 64] >> (i & 0x3f) & 1);
  }

  inline int get_id(int a, int b) const {
    // OPA_CHECK(a >= 0 && a < N && b >= 0 && b < M, a, b, N, M);
    if (Left)
      return a * M2 + (b >> 6);
    else
      return b * N2 + (a >> 6);
  }
  inline int get_bitpos(int a, int b) const {
    if (Left)
      return (b & 0x3f);
    else
      return (a & 0x3f);
  }

  bool is_zero() const {
    for (const auto &x : tb)
      if (x != 0) return false;
    return true;
  }

  void setz() { memset(&tb[0], 0, sizeof(tb[0]) * tb.size()); }

  int count_set() const {
    int res = 0;
    for (const auto &e : tb) {
      res += count_bit(e);
    }
    return res;
  }

  inline u64 get_mask(int a, int b) const { return 1ull << get_bitpos(a, b); }

  inline int get(int a, int b) const {
    return tb[get_id(a, b)] >> get_bitpos(a, b) & 1;
  }

  inline void set(int a, int b, int v) {
    int e = get_id(a, b);
    u64 r = tb[e];
    r = (r & ~get_mask(a, b)) | ((u64)v << get_bitpos(a, b));
    tb[e] = r;
  }

  inline void toggle(int a, int b) {
    int e = get_id(a, b);
    u64 r = tb[e];
    r ^= (u64)1 << get_bitpos(a, b);
    tb[e] = r;
  }

  template <class X> void set_row(int i, const std::vector<X> &data) {
    REP (j, data.size())
      set(i, j, data[j]);
  }

  template <class X> void set_col(int i, const std::vector<X> &data) {
    REP (j, data.size())
      set(j, i, data[j]);
  }

  template <bool V> MulMatrixF2<V> copy() const {
    MulMatrixF2<V> res;
    res.init(N, M);
    REP (i, N)
      REP (j, M)
        res.set(i, j, get(i, j));
    return res;
  }

  int reduce(int lastc) {
    int pos = 0;
    if (lastc == -1) lastc = M;

    for (int i = 0; i < N;) {
      if (pos == lastc) return i;

      int pivot = -1;
      FOR (j, i, N)
        if (get(j, pos) == 1) {
          pivot = j;
          break;
        }

      if (pivot != -1) {
        int idc = pos >> 6;
        FOR (c, idc, M2) { std::swap(tb[pivot * M2 + c], tb[i * M2 + c]); }
        REP (j, N)
          if (j != i && get(j, pos)) {
            FOR (c, idc, M2)
              tb[j * M2 + c] ^= tb[i * M2 + c];
          }
        ++i;
      }
      ++pos;
    }
    return -1;
  }

  template <class X>
  MulMatrixF2<Left> extract_rows(const std::vector<X> &rows) const {
    MulMatrixF2<Left> res;
    res.init(rows.size(), M);
    REP (i, rows.size())
      REP (j, M2)
        res.tb[res.get_id(i, j * 64)] = tb[get_id(rows[i], j * 64)];
    return res;
  }

  // eval only on 128x128 mat
  inline void eval(const Gf128_t2 &in, Gf128_t2 &o) const {
    static_assert(Left == false, "Need left false");
    // OPA_CHECK0(in.size() == 2);

    REP (k, 2) {
      u64 cur = 0;
      REP (i, M) {
        // if (GETB_BLOCK(in, i)) cur ^= tb[get_id(k * 64, i)];
        cur ^= tb[get_id(k * 64, i)] & -GETB_BLOCK(in, i);
      }
      o[k] = cur;
    }
  }

  void set_num_rows(int nrows) {
    static_assert(Left == true, "need true");
    OPA_CHECK0(nrows <= N);
    N = nrows;
    tb.resize(N * M2);
  }

  template <bool V> MulMatrixF2<V> eval_all(const MulMatrixF2<V> &v) const {
    OPA_CHECK0(v.M == 1 && v.N == M);
    MulMatrixF2<V> res;
    res.init(N, 1);
    REP (i, N) {
      int a = 0;
      REP (j, M)
        a ^= v.get(j, 0) & get(i, j);
      res.set(i, 0, a);
    }
    return res;
  }

  void disp() const {
    REP (i, N) {
      REP (j, M)
        printf("%d", get(i, j));
      puts("");
    }
  }

  std::string str(bool header = true) const {
    std::ostringstream oss;
    if (header) {
      oss << "MulMatrix N=" << N << " M=" << M << "\n";
    }
    REP (i, N) {
      REP (j, M)
        oss << get(i, j);
      if (i != N - 1) oss << "\n";
    }
    return oss.str();
  }

  void sadd(const MulMatrixF2<Left> &a) {
    OPA_CHECK0(a.N == N && a.M == M);
    REP (i, tb.size())
      tb[i] ^= a.tb[i];
  }

  void sneg() {
    REP (i, tb.size())
      tb[i] = ~tb[i];
  }
  void sand(const MulMatrixF2<Left> &a) {
    OPA_CHECK0(a.N == N && a.M == M);
    REP (i, tb.size())
      tb[i] &= a.tb[i];
  }
  void sor(const MulMatrixF2<Left> &a) {
    OPA_CHECK0(a.N == N && a.M == M);
    REP (i, tb.size())
      tb[i] |= a.tb[i];
  }

  MulMatrixF2<Left> transpose() const {
    MulMatrixF2<Left> res;
    res.init(M, N);
    REP (i, N)
      REP (j, M)
        res.set(j, i, get(i, j));
    return res;
  }

  MulMatrixF2<false> get_kernel_basis() {
    int lastr = reduce(-1);
    std::vector<std::vector<int> > res;

    if (lastr == -1) lastr = N;

    int j = 0;
    std::vector<int> fixed_vars;
    REP (i, lastr + 1) {
      while (j < M && (i == lastr || !get(i, j))) {
        std::vector<int> cur;
        REP (k, i)
          if (get(k, j)) cur.pb(fixed_vars[k]);
        cur.pb(j);
        res.pb(cur);
        ++j;
      }
      fixed_vars.pb(j);
      ++j;
    }

    MulMatrixF2<false> m;
    m.init(M, res.size());

    REP (i, res.size()) {
      for (auto x : res[i]) m.set(x, i, 1);
    }
    return m;
  }

  std::vector<int> get_col(int a) const {
    std::vector<int> res;
    REP (i, N)
      res.pb(get(i, a));
    return res;
  }

  bool operator==(const MulMatrixF2<Left> &other) const {
    return N == other.N && M == other.M && tb == other.tb;
  }

  friend std::ostream &operator<<(std::ostream &os,
                                  const MulMatrixF2<Left> &v) {
    os << v.str();
    return os;
  }

  Matrix<u32> to_mat() const {
    Matrix<u32> res;
    res.initialize(&GF2, N, M);
    REP (i, N)
      REP (j, M)
        res(i, j) = get(i, j);
    return res;
  }

  std::vector<u64> tb;
  int N2=0, M2=0;
  int N=0, M=0;
};

template <typename TArray>
std::vector<int> bits_to_vec(const TArray &tb, int n) {
  std::vector<int> res;
  REP (i, n) { res.pb(GETB_BLOCK(tb, i)); }
  return res;
}

inline void add_gf128_t2(Gf128_t2 &i1, const Gf128_t2 &i2) {
  i1[0] ^= i2[0];
  i1[1] ^= i2[1];
}

inline std::string to_bytes_gf128_t2(const Gf128_t2 &a) {
  return std::string((const char *)&a[0], 16);
}

inline Gf128_t to_gf128_t(const Gf128_t2 &a) {
  return gf128.import_base(bignum::fromrbytes(to_bytes_gf128_t2(a)));
}

inline Gf128_t2 to_gf128_t2(StringRef s) {
  OPA_CHECK0(s.size() == 16);
  Gf128_t2 res;
  res[0] = *(u64 *)(s.data());
  res[1] = *(u64 *)(s.data() + 8);
  return res;
}

inline Gf128_t to_gf128_t(StringRef s) { return to_gf128_t(to_gf128_t2(s)); }

inline Gf128_t2 to_gf128_t2(const Gf128_t &a) {
  return to_gf128_t2(gf128.export_base(a).getrbytes(16));
}

const int maxb = 1 << 16;
extern u8 cnt_bit[maxb];

inline int get_cnt_bit64(u64 a) {
  return cnt_bit[a & 0xffff] ^ cnt_bit[a >> 16 & 0xffff] ^
         cnt_bit[a >> 32 & 0xffff] ^ cnt_bit[a >> 48 & 0xffff];
}

inline MulMatrixF2<true> mul(const MulMatrixF2<true> &a,
                             const MulMatrixF2<false> &b) {
  MulMatrixF2<true> res;
  OPA_CHECK0(a.M == b.N);

  res.init(a.N, b.M);
  const int N2 = (a.M + 63) / 64;
  REP (i, a.N)
    REP (j, b.M) {
      int x = 0;
      REP (k, N2) {
        x ^=
          get_cnt_bit64(a.tb[a.get_id(i, k * 64)] & b.tb[b.get_id(k * 64, j)]);
      }
      res.set(i, j, x);
    }
  return res;
}

inline MulMatrixF2<false> mul2(const MulMatrixF2<true> &a,
                               const MulMatrixF2<false> &b) {
  MulMatrixF2<false> res;
  OPA_CHECK0(a.M == b.N);

  res.init(a.N, b.M);
  const int N2 = (a.M + 63) / 64;
  REP (i, a.N)
    REP (j, b.M) {
      int x = 0;
      REP (k, N2) {
        x ^=
          get_cnt_bit64(a.tb[a.get_id(i, k * 64)] & b.tb[b.get_id(k * 64, j)]);
      }
      res.set(i, j, x);
    }
  return res;
}

template <bool V1, bool V2>
inline MulMatrixF2<true> mul_slow(const MulMatrixF2<V1> &a,
                                  const MulMatrixF2<V2> &b) {
  MulMatrixF2<true> res;
  OPA_CHECK0(a.M == b.N);

  res.init(a.N, b.M);
  REP (i, a.N)
    REP (j, b.M) {
      int x = 0;
      REP (k, a.M)
        x ^= a.get(i, k) & b.get(k, j);
      res.set(i, j, x);
    }
  return res;
}

template <bool L1, bool L2>
inline MulMatrixF2<true> add(const MulMatrixF2<L1> &a,
                             const MulMatrixF2<L2> &b) {
  MulMatrixF2<true> res;
  res = a.template copy<true>();
  res.sadd(b);
  return res;
}

inline void and_gf128_t2(Gf128_t2 &a, const Gf128_t2 &b) {
  a[0] &= b[0];
  a[1] &= b[1];
}

inline void shiftr_gf128_t2(Gf128_t2 &a) {
  a[0] >>= 1;
  a[0] |= a[1] << 63;
  a[1] >>= 1;
}

inline void shiftl_gf128_t2(Gf128_t2 &a) {
  a[1] <<= 1;
  a[1] |= a[0] >> 63;
  a[0] <<= 1;
}

class Gf128_Fast : public opa::math::common::Field<Gf128_t2> {
public:
  std::vector<Gf128_t2> mul_cache;
  bignum m_order;
  bignum m_inv_pw;

  Gf128_Fast() {
    Field<Gf128_t2>::init(bignum(2).pow(128), 2);
    mul_cache.pb({ 0, 0 });
    REP (i, 2 * 128) { mul_cache.pb(to_gf128_t2(gf128.x().powm(i))); }
    m_order = bignum(2).pow(128) - 1;
    m_inv_pw = m_order - 1;
  }

  virtual Gf128_t2 inv(const Gf128_t2 &a) const {
    auto tmp = this->faste(a, m_inv_pw);
    OPA_CHECK0(this->mul(tmp, a) == getE());
    return tmp;
  }

  virtual Gf128_t2 mul(const Gf128_t2 &a, const Gf128_t2 &b) const {
    // return to_gf128_t2(to_gf128_t(a)*to_gf128_t(b));
    Gf128_t2 res = { 0, 0 };
    u64 f1 = a[0], f2 = a[1];
    int hb = 0;
    REP (i, 128)
      if (GETB_BLOCK(b, i)) hb = i;

    int ha = 0;
    REP (i, 128)
      if (GETB_BLOCK(a, i)) ha = i;

    Gf128_t2 rb = { 0, 0 };
    REP (i, hb + 1)
      if (GETB_BLOCK(b, i)) BIT_BLOCK(rb, hb - i) |= BIT_B(rb, hb - i);

    {
      Gf128_t2 tmp = rb;
      REP (i, hb + 1) {
        Gf128_t2 tmp2 = tmp;
        and_gf128_t2(tmp2, a);
        int x = get_cnt_bit64(tmp2[0]) ^ get_cnt_bit64(tmp2[1]);
        add_gf128_t2(res, mul_cache[x * (1 + hb - i)]);

        shiftr_gf128_t2(tmp);
      }
    }

    {
      Gf128_t2 tmp = rb;
      shiftl_gf128_t2(tmp);

      REP (i, ha) {
        Gf128_t2 tmp2 = tmp;
        and_gf128_t2(tmp2, a);
        int x = get_cnt_bit64(tmp2[0]) ^ get_cnt_bit64(tmp2[1]);
        add_gf128_t2(res, mul_cache[x * (2 + hb + i)]);

        shiftl_gf128_t2(tmp);
      }
    }

    return res;
  }
  virtual Gf128_t2 add(const Gf128_t2 &a, const Gf128_t2 &b) const {
    Gf128_t2 res = a;
    add_gf128_t2(res, b);
    return res;
  }

  virtual Gf128_t2 neg(const Gf128_t2 &a) const { return a; }

  virtual bool isZ(const Gf128_t2 &a) const { return a == this->getZ(); }
  virtual bool isE(const Gf128_t2 &a) const { return a == this->getE(); }
  virtual Gf128_t2 import(const Gf128_t2 &a) const { return a; }

  virtual Gf128_t2 getZ() const { return Gf128_t2({ 0, 0 }); }
  virtual Gf128_t2 getE() const { return Gf128_t2({ 1, 0 }); }
  virtual Gf128_t2 getRandRaw() const { return { rng64(), rng64() }; }
};

class BitVecRepr;
class BitVec : public opa::utils::Initable {
public:
  BitVec() {}
  BitVec(int n) { init(n); }
  virtual void init(int n) {
    opa::utils::Initable::init();
    this->n = n;
    m_mat.init(1, n);
  }
  void setz() { m_mat.setz(); }

  BitVec concat(const BitVec &peer) const {
    BitVec res(n + peer.n);
    res.copy(*this);
    res.copy(peer, n);
    return res;
  }

  void copy(const BitVec &src, u32 off = 0, s32 sz = -1) {
    norm_size(off, sz);
    sz = std::min(sz, src.size());
    REP (i, sz)
      set(off + i, src.get(i));
  }

  template <class T> void from(T a, u32 off = 0, s32 sz = -1) {
    norm_size(off, sz);
    OPA_CHECK0(off + sz <= size());
    REP (i, sz) {
      set(off + i, a & 1);
      a >>= 1;
    }
  }

  int dot(const RangeCoverage &x) const {
    int res = 0;
    for (auto i : x.all()) res ^= get(i);
    return res;
  }

  int dot(const std::vector<u32> &x) const {
    int res = 0;
    for (auto i : x) res ^= get(i);
    return res;
  }

  int dot(const BitVec &x) const { return this->andz(x).count_set() & 1; }

  bool is_zero() const { return m_mat.is_zero(); }
  int count_set() const { return m_mat.count_set(); }

  static BitVec From(u64 v, s32 sz) {
    BitVec res(sz);
    res.from(v);
    return res;
  }
  static BitVec FromBytes(const RangeCoverage &range, opa::stolen::StringRef s);
  template <class T> static BitVec FromVec(const std::vector<T> &vec, int n) {
    BitVec res(n);
    for (auto &x : vec) res.toggle(x);
    return res;
  }

  static BitVec FromRange(const RangeCoverage &range, int n) {
    return BitVec::FromVec(range.to_vec(), n);
  }

  template <class T> static BitVec FromBinaryVec(const std::vector<T> &vec) {
    BitVec res(vec.size());
    REP (i, vec.size())
      res.set(i, vec[i]);
    return res;
  }

  static BitVec FromStr(const StringRef &s) {
    BitVec res(s.size());
    REP (i, s.size())
      res.set(i, s[i] == '1');
    return res;
  }

  template <class T> static BitVec From2(const T *ptr, int n) {
    int sz = sizeof(T) * 8;
    BitVec res(n * sz);
    REP (i, n)
      res.from(ptr[i], sz * i, sz);
    return res;
  }

  BitVec extract(u32 off, s32 sz = -1) const {
    norm_size(off, sz);
    BitVec res;
    res.init(sz);
    REP (i, sz)
      res.set(i, get(off + i));
    return res;
  }

  template <class T> T to(u32 off = 0, s32 sz = -1) const {
    T res = 0;
    norm_size(off, sz);
    REP (i, sz)
      res = res << 1 | get(off + sz - 1 - i);
    return res;
  }

  void norm_size(u32 off, s32 &sz) const {
    if (sz == -1) sz = size() - off;
  }

  int size() const { return this->n; }

  u32 get(int a) const { return m_mat.get(0, a); }
  void set(int a, int v) { m_mat.set(0, a, v); }
  void toggle(int a) { m_mat.toggle(0, a); }
  void xorv(int a, int v) {
    if (v) toggle(a);
  }
  void minus1() {
    for (auto &v : m_mat.tb) {
      bool stop = v != 0;
      v -= 1;
      if (stop) break;
    }
  }

  // set i to perm[i]
  template <class T> BitVec shuffle(const std::vector<T> &perm) const {
    OPA_CHECK0(perm.size() == m_mat.M);
    BitVec res = *this;
    REP (i, perm.size())
      res.set(i, get(perm[i]));
    return res;
  }

  bool isz() const {
    REP (k, m_mat.tb.size())
      if (m_mat.tb[k]) return false;
    return true;
  }

  BitVec xorz(const BitVec &peer) const {
    BitVec res = *this;
    return res.sxorz(peer);
  }
  void xorz(u64 *target) const {
    REP (i, m_mat.tb.size())
      target[i] ^= this->m_mat.tb[i];
  }

  BitVec xorz(const RangeCoverage &x) const { return BitVec(*this).sxorz(x); }

  BitVec &sxorz(const BitVec &peer) {
    this->m_mat.sadd(peer.m_mat);
    return *this;
  }
  BitVec &sxorz(const RangeCoverage &x) {
    for (auto &v : x.all()) {
      toggle(v);
    }
    return *this;
  }

  BitVec orz(const BitVec &peer) const {
    BitVec res = *this;
    res.m_mat.sor(peer.m_mat);
    return res;
  }

  BitVec andz(const BitVec &peer) const {
    BitVec res = *this;
    res.m_mat.sand(peer.m_mat);
    return res;
  }

  BitVec swap(int mid) const {
    OPA_CHECK0(mid * 2 == n);
    BitVec res(n);
    res.copy(extract(0, mid), mid);
    res.copy(extract(mid), 0);
    return res;
  }

  bool operator==(const BitVec &a) const { return m_mat == a.m_mat; }
  bool operator!=(const BitVec &a) const { return !(m_mat == a.m_mat); }

  bool operator<(const BitVec &x) const { return m_mat.tb < x.m_mat.tb; }
  void disp() const { m_mat.disp(); }
  std::string str() const { return m_mat.str(); }

  static BitVec rand(int n);

  friend std::ostream &operator<<(std::ostream &os, const BitVec &v) {
    os << "BitVec: " << v.m_mat.str(false);
    return os;
  }
  Matrix<u32> to_mat() const { return m_mat.to_mat(); }

  MulMatrixF2<false> m_mat;
  int n = 0;

  BitVecRepr to_repr() const;

private:
};

struct BitVecRepr : opa::utils::ProtobufParams {
  std::vector<u64> v;
  int n;

  BitVec to_bitvec() const {

    BitVec res;
    res.init(n);
    res.m_mat.tb = v;
    return res;
  }

  OPA_TGEN_IMPL(v, n);
};
inline BitVecRepr BitVec::to_repr() const {
  BitVecRepr res;
  res.v = m_mat.tb;
  res.n = n;
  return res;
}

OPA_NM_MATH_COMMON_END

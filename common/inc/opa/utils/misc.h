#pragma once

#include <opa/stolen/StringRef.h>
#include <opa/utils/serialize.h>
#include <opa/utils/string.h>

#define OPA_NAMED_PAIR(clname, a, b, t1, t2)                                   \
  struct clname : public std::pair<t1, t2> {                                   \
    clname(const t1 &x, const t2 &y) : std::pair<t1, t2>(x, y) {}              \
    clname(const std::pair<t1, t2> &x) : std::pair<t1, t2>(x) {}               \
    clname() {}                                                                \
    inline t1 &a() { return this->ST; }                                        \
    inline t2 &b() { return this->ND; }                                        \
    inline const t1 &a() const { return this->ST; }                            \
    inline const t2 &b() const { return this->ND; }                            \
  };

namespace opa {
namespace utils {

template <class T> class ScopedHeapAlloc {
public:
  ScopedHeapAlloc(int n) {
    m_buf = new T[n];
    m_size = n;
  }
  ~ScopedHeapAlloc() { delete[] m_buf; }
  T &operator[](int i) { return m_buf[i]; }
  const T &operator[](int i) const { return m_buf[i]; }

  OPA_ACCESSOR_PTR(T, m_buf, buf);
  OPA_ACCESSOR_R(int, m_size, size);
  OPA_ACCESSOR_R2(int, m_size * sizeof(T), bytesize);

private:
  T *m_buf;
  int m_size;
};

template <class T> struct PtrComparator {
  bool operator()(const T &a, const T &b) const { return *a < *b; }
};

class Base {
public:
  virtual ~Base() {}
};
OPA_DECL_SPTR(Base, BasePtr);

class Initable : public Base {
public:
  Initable() { m_init = false; }
  virtual void init();

  virtual void fini();
  bool is_init() const { return m_init; }

  void check_init() const { OPA_CHECK0(m_init == true); }
  void check_not_init() const { OPA_CHECK0(m_init == false); }

  virtual ~Initable();

protected:
  Initable &operator=(const Initable &a) {
    m_init = a.m_init;
    return *this;
  }

private:
  bool m_init;
};

std::string get_hostname();
std::string get_process_fingerprint();

template <class T> std::vector<T> rotatel(const std::vector<T> &a, int cnt) {
  int n = a.size();
  cnt %= n;
  if (cnt < 0) cnt += n;
  std::vector<T> res(n);
  REP (i, n)
    res[i] = a[(cnt + i) % n];
  return res;
}

template <class T> void make_unique(std::vector<T> &x) {
  x.resize(std::unique(ALL(x)) - x.begin());
}

template <class T, class Y> void make_unique(std::vector<T> &x, const Y &cmp) {
  x.resize(std::unique(ALL(x), cmp) - x.begin());
}

template <class T, class U> void container_extend(T &t, const U &u) {
  t.insert(t.end(), ALL(u));
}

template <class Container> class MergeWalker {
public:
  typedef typename Container::const_iterator container_const_iterator;
  typedef typename Container::value_type value_type;

  class ConstMergeIterator {
  public:
    ConstMergeIterator(const std::vector<container_const_iterator> &its,
                       const std::vector<container_const_iterator> &ends) {
      m_its = its;
      m_ends = ends;
      m_ids.resize(its.size());
      choose_cur();
    }
    void choose_cur() {
      m_cur = -1;
      m_selected_id = 0;
      REP (i, m_its.size()) {
        m_ids[i] = -1;
        if (good(i) && (m_cur == -1 || *m_its[i] < *m_its[m_cur])) {
          m_ids[i] = ++m_selected_id;
          m_cur = i;
        } else if (good(i) && *m_its[i] == *m_its[m_cur]) {
          m_ids[i] = m_selected_id;
        } else
          m_ids[i] = -1;
      }
    }

    const value_type &operator*() const {
      OPA_CHECK0(m_cur != -1);
      return *m_its[m_cur];
    }
    ConstMergeIterator &operator++() {
      REPV (i, m_its.size()) {
        if (good(i) && *m_its[i] == *m_its[m_cur]) ++m_its[i];
        if (i == m_cur) break;
      }
      choose_cur();
      return *this;
    }

    bool operator==(const ConstMergeIterator &it) const {
      REP (i, m_its.size()) {
        if (m_its[i] != it.m_its[i]) return false;
      }
      return true;
    }

    bool operator!=(const ConstMergeIterator &it) const {
      return !(*this == it);
    }
    bool is_active(int i) const { return m_ids[i] == m_selected_id; }

  private:
    bool good(int i) const { return m_its[i] != m_ends[i]; }
    int m_cur;
    int m_selected_id;
    std::vector<int> m_ids;
    std::vector<container_const_iterator> m_its;
    std::vector<container_const_iterator> m_ends;
  };

  ConstMergeIterator begin() const {
    std::vector<container_const_iterator> a;
    for (const auto &x : m_lst) a.pb(x->begin());
    return ConstMergeIterator(a, get_end());
  }

  ConstMergeIterator end() const {
    return ConstMergeIterator(get_end(), get_end());
  }

  MergeWalker(const std::vector<const Container *> &lst) { m_lst = lst; }

private:
  std::vector<container_const_iterator> get_end() const {
    std::vector<container_const_iterator> a;
    for (const auto &x : m_lst) a.pb(x->end());
    return a;
  }

  std::vector<const Container *> m_lst;
};

template <class Key, class Item, class KeyComp = std::less<Key> >
class KBestContainer : public Initable {
public:
  virtual void init(int lim) {
    Initable::init();
    m_lim = lim;
  }

  bool should_add(const Key &k) const {
    if (m_data.size() < m_lim) return true;
    // we are better than the worst in the current set.
    if (comparator(m_data.begin()->ST, k)) return true;
    return false;
  }

  void maybe_add(const Key &k, const Item &item) {
    if (!should_add(k)) return;
    add(k, item);
  }

  void add(const Key &k, const Item &item) {
    if (m_seen.count(item)) return;
    m_seen.insert(item);

    m_data.insert(MP(k, item));
    if (m_data.size() > m_lim) {
      m_seen.erase(m_data.begin()->ND);
      m_data.erase(m_data.begin());
    }
  }

  void get_items(std::vector<Item> &res) const {
    res.reserve(m_data.size());
    for (auto &x : m_data) {
      res.pb(x.ND);
    }
  }

  int size() const { return m_data.size(); }

  void clear() {
    m_data.clear();
    m_seen.clear();
  }

private:
  int m_lim;
  std::multimap<Key, Item, KeyComp> m_data;
  std::set<Item> m_seen;
  KeyComp comparator;
};

template <int N> inline u64 rotl(u64 a, int b) {
  if (b < 0) b += N;
  return ((a << b) | (a >> (N - b))) & ((1ull << N) - 1);
}

template <int N> inline u64 rotr(u64 a, int b) {
  if (b < 0) b += N;
  return ((a >> b) | (a << (N - b))) & ((1 << N) - 1);
}

template <class T> inline u32 dot(T a, T b) {
  u32 r = 0;
  while (a && b) {
    r ^= a & b & 1;
    a >>= 1;
    b >>= 1;
  }
  return r;
}

template <class T, class U> struct RemapCompare {
  RemapCompare(const std::map<T, int> &rmp) : rmp(rmp) {}
  bool operator()(const T &a, const T &b) const {
    return rmp.find(a)->ND < rmp.find(b)->ND;
  }

  const std::map<T, int> &rmp;
};

template <class T> struct PairCompareFirst {
  bool operator()(const T &a, const T &b) { return a.ST < b.ST; }
};

template <class To, class From>
std::vector<To> convert(const std::vector<From> &a) {
  std::vector<To> res;
  res.reserve(a.size());
  REP (i, a.size())
    res.emplace_back(a[i]);
  return res;
}

template <class T>
std::vector<std::pair<T, int> > list_to_count(const std::vector<T> &lst) {
  std::map<T, int> rmp;
  for (auto &x : lst) {
    ++rmp[x];
  }
  std::vector<std::pair<T, int> > res(ALL(rmp));
  return res;
}

template <class ContainerT, class ContainerIds>
ContainerT select_ids(const ContainerT &ts, const ContainerIds &ids) {
  ContainerT res;
  for (auto &id : ids) res.push_back(ts[id]);
  return res;
}

template <class ContainerT> class CircularView {

public:
  typedef typename ContainerT::value_type T;

  CircularView(ContainerT &container) : m_container(container) {
    n = container.size();
  }

  T &operator[](int i) { return m_container[i % n]; }
  const T &operator[](int i) const { return m_container[i % n]; }


  int n;
  ContainerT &m_container;
};
template<class ContainerT>
static CircularView<ContainerT>  CreateCircularView(ContainerT &container){
  return CircularView<ContainerT>(container);
}

template <typename T> using BinaryOp = std::function<T(const T &, const T &)>;

} // namespace utils
} // namespace opa

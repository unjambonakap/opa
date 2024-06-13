#pragma once

#include <opa/utils/hash.h>
#include <opa/utils/map_util.h>
#include <opa_common.h>
#include <type_traits>

OPA_NAMESPACE_DECL2(opa, utils)

typedef s32 IdType;
static const IdType InvalidId = -1;

#define OPA_IDOBJ_FIELD() ::opa::utils::IdType id;

struct IdObj {
  IdType id;
};

template<class T> 
std::vector<T> to_vector(const std::unordered_set<T> &set){ return std::vector<T>(ALL(set)); }
template <class T> class ObjContainer {
public:
  template <class U> U *add(U *a) {
    objs.emplace_back((T *)a);
    return a;
  }
  std::deque<UPTR(T)> objs;
};

class IdGenerator {
public:
  IdGenerator() { m_tot = 0; }

  IdType get_new() {
    while (m_free.empty()) {
      int inc = std::max(16, (int)m_tot);
      REP (i, inc) {
        int cur = m_tot + i;
        if (m_blocked.count(cur)) continue;
        m_free.push(cur);
      }
      m_tot += inc;
    }

    IdType res = m_free.front();
    m_free.pop();
    return res;
  }

  void block(IdType id) { m_blocked.insert(id); }

  void free(IdType id) { m_free.push(id); }

  std::set<int> m_blocked;
  std::queue<int> m_free;
  u32 m_tot;
};

template <class Obj> class ObjectPool {
public:
  class ObjectPoolLoader {
  public:
    ObjectPoolLoader(ObjectPool<Obj> &pool) : m_pool(pool) {}
    void load(IdType id) { m_seen.insert(id); }

    ~ObjectPoolLoader() {
      if (!m_seen.size()) return;
      IdType mx = *--m_seen.end();
      m_pool.m_objs.resize(mx + 1);
      m_pool.m_used = m_seen;
      REP (i, mx)
        if (!m_seen.count(i)) m_pool.m_free.push(i);
    }

  private:
    ObjectPool<Obj> &m_pool;
    std::set<IdType> m_seen;
  };

  Obj &get_new2() {
    if (m_free.empty()) {
      int inc = std::max(16, (int)m_objs.size());
      REP (i, inc)
        m_free.push(m_objs.size() + i);
      m_objs.resize(m_objs.size() + inc);
    }
    IdType cur = m_free.front();
    m_free.pop();
    m_objs[cur].id = cur;
    m_used.insert(cur);
    return m_objs[cur];
  }

  IdType get_new() { return get_new2().id; }

  std::vector<IdType> get_new_vec(int n) {
    std::vector<IdType> ids;
    REP (i, n)
      ids.pb(get_new());
    return ids;
  }

  void get_new(int n, ...) {
    auto ids = get_new_vec(n);

    va_list ap;
    va_start(ap, n);
    REP (i, n) {
      Obj **dest = va_arg(ap, Obj **);
      *dest = &get(ids[i]);
    }
    va_end(ap);
  }

  u32 size() const { return m_used.size(); }
  u32 tot_size() const { return m_objs.size(); }
  void remove(IdType id) {
    m_free.push(id);
    m_used.erase(id);
  }
  Obj &get(IdType id) {
    OPA_ASSERT0(has(id));
    return m_objs[id];
  }
  bool has(IdType id) const { return m_used.count(id); }
  const Obj &get(IdType id) const {
    OPA_CHECK0(has(id));
    return m_objs[id];
  }

  const Obj *get_start() const { return m_objs.data(); }
  OPA_ACCESSOR_R(std::set<IdType>, m_used, used);

  void clear() {
    m_objs.clear();
    m_used.clear();
    while (m_free.size() > 0) m_free.pop();
  }

private:
  std::vector<Obj> m_objs;
  std::set<IdType> m_used;
  std::queue<IdType> m_free;
};

template <class Obj> class ContinuousBufferContainer {
public:
  virtual IdType get_new() {
    IdType id = m_gen.get_new();
    get_pos(id) = m_objs.size();
    m_objs.push_back(Obj());
    m_objs.back().id = id;
    return id;
  }

  void move(IdType id, u32 dest_pos) {
    OPA_CHECK0(has(id));
    u32 src_pos = get_pos(id);
    move_pos(src_pos, dest_pos);
  }

  void move_pos(u32 src_pos, u32 dest_pos) {
    if (src_pos == dest_pos) return;

    IdType x = m_objs[dest_pos].id;
    std::swap(m_objs[dest_pos], m_objs[src_pos]);
    get_pos(m_objs[dest_pos].id) = dest_pos;
    get_pos(x) = src_pos;
  }

  virtual void remove(IdType id) {
    move(id, m_objs.size() - 1);
    m_posmap.erase(id);
    m_gen.free(id);
    m_objs.pop_back();
  }

  const Obj *get_start() const { return m_objs.data(); }
  bool has(IdType id) const { return m_posmap.count(id); }

  Obj &get(IdType id) { return m_objs[get_pos(id)]; }
  const Obj &get(IdType id) const { return m_objs[get_pos(id)]; }

  u32 get_size() const { return m_objs.size(); }

protected:
  u32 &get_pos(IdType id) { return m_posmap[id]; }
  u32 get_pos(IdType id) const {
    OPA_CHECK0(m_posmap.count(id));
    return m_posmap.find(id)->ND;
  }

  std::vector<Obj> m_objs;
  IdGenerator m_gen;
  std::map<IdType, u32> m_posmap;
};

template <class Obj>
class BucketContBufferContainer : public ContinuousBufferContainer<Obj> {
public:
  typedef int BucketType;

  BucketContBufferContainer() {}

  virtual IdType get_new(BucketType bucket) {
    IdType id = ContinuousBufferContainer<Obj>::get_new();
    m_buckets.back().second = this->get_size();
    return id;
  }

  void create_bucket(BucketType bucket) {
    m_bucket_pos[bucket] = m_buckets.size();
    m_buckets.push_back(MP(this->get_size(), bucket));
  }

  void set_bucket(IdType id, BucketType dest) {
    u32 pos = get_bucket_pos(id);
    u32 dest_pos = m_bucket_pos[dest];

    if (pos < dest_pos) {
      for (; pos != dest_pos; ++pos) {
        this->move(id, m_buckets[pos].first);
        --m_buckets[pos].first;
      }

    } else {
      for (; pos != dest_pos; --pos) {
        this->move(id, m_buckets[pos - 1].first);
        ++m_buckets[pos - 1].first;
      }
    }
  }

  const Obj *get_start(BucketType bucket) const {
    return ContinuousBufferContainer<Obj>::get_start() +
           get_bucket_start(bucket);
  }

  int get_bucket_size(BucketType bucket) const {
    return m_buckets[get_bucket_id(bucket)].ND - get_bucket_start(bucket);
  }

  virtual void remove(IdType id) override {
    int pos = get_bucket_pos(id);

    for (; pos < m_buckets.size() - 1; ++pos)
      set_bucket(id, m_buckets[pos + 1].second);
    ContinuousBufferContainer<Obj>::remove(id);
    --m_buckets.back().second;
  }

  BucketType get_bucket(IdType id) const {
    return m_buckets[get_bucket_pos(id)].second;
  }

private:
  void check_bucket(BucketType bucket) const {
    OPA_CHECK0(m_bucket_pos.count(bucket));
  }
  u32 get_bucket_id(BucketType bucket) const {
    check_bucket(bucket);
    return m_bucket_pos.find(bucket)->ND;
  }

  u32 get_bucket_start(BucketType bucket) const {
    u32 pos = get_bucket_id(bucket);
    return pos == 0 ? 0 : m_buckets[pos - 1].ND;
  }

  u32 get_bucket_pos(IdType id) {
    auto it = std::lower_bound(ALL(m_buckets), MP(this->get_pos(id), -1));
    return it - m_buckets.begin();
  }

  std::vector<std::pair<u32, BucketType> > m_buckets; // sorted
  std::map<BucketType, u32> m_bucket_pos;
};

template <class T> class DataStream {
public:
  virtual T get() { OPA_CHECK0(false); }
  virtual void push(const T &a) { OPA_CHECK0(false); }
  virtual bool has_more() const { OPA_CHECK0(false); }
};

template <class T> class QueueDataStream {
public:
  virtual T get() {
    T res = m_queue.front();
    m_queue.pop();
    return res;
  }
  virtual void push(const T &a) { m_queue.push(a); }
  virtual bool has_more() const { return !m_queue.empty(); }
  void reset() {
    while (has_more()) get();
  }

private:
  std::queue<T> m_queue;
};

class BitDataStream : public QueueDataStream<u32> {
public:
  template <typename X> void push_multiple(const X &a) {
    REP (i, sizeof(X))
      this->push(a >> i & 1);
  }

  void push_multiple(const u64 &a, int n) {
    REP (i, n)
      this->push(a >> i & 1);
  }
};

template <typename K, typename V> class LRUCache {
public:
  LRUCache() {}
  LRUCache(int n) { init(n); }
  void init(int n) { m_n = n; }

  u64 get_fresh_score(const K &key) const {
    auto it = m_freshness.find(key);
    return it == m_freshness.end() ?  0 : *it;
  }

  u64 new_freshness() const { return step++; }

  V *get(const K &key) const {
    u64 fresh_score = get_fresh_score(key);
    if (fresh_score == 0) return nullptr;
    V *cur = m_entries[fresh_score];

    m_entries.erase(fresh_score);

    fresh_score = new_freshness();
    m_freshness[key] = fresh_score;
    m_entries[fresh_score] = cur;
    return cur;
  }

  void add(const K &key, V *val) const {
    u64 fresh_score = get_fresh_score(key);
    m_entries[fresh_score] = val;
    m_freshness[key] = fresh_score;
    if (m_freshness.size() > m_n) {
      u64 bad = m_freshness.begin()->second;
      auto to_del = m_entries.find(bad);
      delete to_del->second;
      m_entries.erase(to_del);
      m_freshness.erase(m_freshness.begin());
    }
  }

  V *get_or_insert(const K &key, std::function<V *()> func) const {
    V *cur = get(key);
    if (cur == nullptr) {
      cur = func();
      add(key, cur);
    }
    return cur;
  }

private:
  int m_n;
  mutable u64 step = 1;
  mutable std::unordered_map<u64, V *> m_entries;
  mutable std::unordered_map<K, u64> m_freshness;
};

template <class K, class V, class CMP = std::less<K> > class MinFinderPair {
public:
  MinFinderPair(CMP cmp = CMP()) : m_cmp(cmp) {}
  void reset() { m_has = false; }

  void update(const K &k, const V &v) {
    if (!m_has) {
      m_best = MP(k, v);
      m_has = true;
    } else {
      if (m_cmp(k, m_best.first)) m_best = MP(k, v);
    }
  }
  bool has() const { return m_has; }
  const V &get() const { return m_best.second; }
  const K &get_cost() const { return m_best.first; }

private:
  CMP m_cmp;
  std::pair<K, V> m_best;
  bool m_has = false;
};

template <class K, class V>
using MaxFinderPair = MinFinderPair<K, V, std::greater<K> >;

template <class T, class Pred = std::less<T> > class MinFinder {
public:
  typedef MinFinder<T, Pred> SelfType;
  bool update(const T &x) {
    if (!m_has) {
      m_best = x;
      m_has = true;
      return true;
    } else {
      if (Pred()(x, m_best)){
        m_best = x;
        return true;
      }
    }
    return false;
  }

  template <class U, T U::*GetField>
  SelfType &update(const std::vector<U> &tb) {
    for (auto &x : tb) {
      this->update(x.*GetField);
    }
    return *this;
  }

  bool has() const { return m_has; }
  const T &get() const { return m_best; }
  const T &get_or(const T &dv) const { return has() ? get() : dv; }
  void reset() { m_has = false; }

private:
  T m_best;
  bool m_has = false;
};

template <class V> using MaxFinder = MinFinder<V, std::greater<V> >;

template <class T,
          typename MapType = std::unordered_map<
            T, int, std::hash<T> /*::opa::utils::hash_utils::hash<T>i*/> >
class Remapper {
public:
  Remapper(int maxn = -1) { reset(maxn); }
  int get(const T &a) {
    if (m_remap.count(a)) return m_remap[a];
    int nid = m_map.size();
    OPA_CHECK0(maxn == -1 || nid < maxn);
    m_map.push_back(a);
    m_remap[a] = nid;
    return nid;
  }

  bool has(const T &a) const { return m_remap.count(a) != 0; }

  template <typename... Args> int get2(const Args &... args) {
    return get(T(args...));
  }

  int get_or_die(const T &a) const {
    const int *v = glib::gtl::FindOrNull(m_remap, a);
    OPA_CHECK(v != nullptr, a);
    return *v;
  }

  const T &rget(int id) const { return m_map[id]; }
  int size() const { return m_remap.size(); }

  OPA_ACCESSOR_R(std::deque<T>, m_map, mp);
  const MapType &remap() const { return m_remap; }
  void clear() { reset(); }

  void reset(int maxn = -1) {
    m_remap.clear();
    m_map.clear();
    this->maxn = maxn;
  }

private:
  int maxn;
  MapType m_remap;
  std::deque<T> m_map;
};

template <class U, class V, class MapType = std::unordered_map<U, V>,
          class MapType2 = std::unordered_map<V, U> >
class TwoWayRemapper {
public:
  U vu(const V &a) const {
    auto it = m_vu.find(a);
    return it->second;
  }

  V uv(const U &a) const {
    auto it = m_uv.find(a);
    return it->second;
  }

  bool hasu(const U &a) const { return m_uv.count(a); }
  bool hasv(const V &a) const { return m_vu.count(a); }
  void add(const U &a, const V &b) {
    m_uv[a] = b;
    m_vu[b] = a;
  }

  int size() const { return m_uv.size(); }

private:
  MapType m_uv;
  MapType2 m_vu;
};

template <class U> class ConstantChecker {
public:
  bool add(const U &u) {
    if (first) {
      cur = u;
      return true;
    }
    return cur == u;
  }

private:
  U cur;
  bool first = true;
};

template <int K, class Key, class V, class CMP = std::less<Key> >
class KMinFinderPair {
public:
  KMinFinderPair(CMP cmp = CMP()) : m_cmp(cmp) { reset(); }

  void update(const Key &k, const V &v) {
    std::pair<Key, V> cur(k, v);
    REP (i, K) {
      if (i == k_pos) {
        m_best[k_pos++] = cur;
        break;
      }

      if (m_cmp(cur.first, m_best[i].first)) {
        std::swap(m_best[i], cur);
      }
    }
  }

  bool has(int k) const { return k_pos > k; }
  const V &get(int k) const { return m_best[k].second; }
  const Key &get_cost(int k) const { return m_best[k].first; }

  void reset() { k_pos = 0; }

private:
  CMP m_cmp;
  std::array<std::pair<Key, V>, K> m_best;
  int k_pos;
};
template <int K, class Key, class V>
using KMaxFinderPair = KMinFinderPair<K, Key, V, std::greater<Key> >;

OPA_NAMESPACE_DECL2_END

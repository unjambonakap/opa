#pragma once

#include <opa/utils/DataStruct.h>
#include <opa_common.h>
#include <type_traits>

OPA_NAMESPACE_DECL2(opa, algo)

class UnionJoin {
public:
  UnionJoin(int n = 0) { reset(n); }
  void merge(int a, int b) {
    set_repr(a, b);
  }

  void set_repr(int a, int b) {
    // a is repr of b
    a = root(a);
    b = root(b);
    if (a != b) {
      --m_count;
      m_tb[a] += m_tb[b];
      m_tb[b] = a;
    }
  }

  bool same(int a, int b) const { return root(a) == root(b); }
  bool is_repr(int a) const { return m_tb[a] < 0; }

  int size(int a) const { return -m_tb[root(a)]; }
  int root(int a) const {
    int x = m_tb[a];
    return x < 0 ? a : (m_tb[a] = root(x));
  }

  void reset(int n) {
    m_tb = std::vector<int>(n, -1);
    m_count = n;
  }

  static UnionJoin FromVec(const std::vector<int> &lst) {
    UnionJoin res(lst.size());
    std::unordered_map<int, int> rmp;
    REP (i, lst.size()) {
      if (!rmp.count(lst[i])) {
        rmp[lst[i]] = i;
      } else {
        res.merge(i, rmp[lst[i]]);
      }
    }
    return res;
  }

  int n() const { return m_tb.size(); }

  std::vector<int> get_clusters() const {
    std::unordered_map<int, int> rmp;
    std::vector<int> res;
    REP (i, n()) {
      int x = root(i);
      int v = rmp.size();
      if (rmp.count(x)) {
        v = rmp[x];
      } else
        rmp[x] = v;
      res.push_back(v);
    }
    return res;
  }

  std::vector<std::vector<int>> get_groups() const {
    std::unordered_map<int, int> rmp;
    std::vector<std::vector<int>> res;
    REP (i, n()) {
      int x = root(i);
      int v = rmp.size();
      if (rmp.count(x)) {
        v = rmp[x];
      } else{
        rmp[x] = v;
        res.emplace_back();
      }
      res[v].push_back(i);
    }
    return res;
  }


  std::unordered_map<int, std::vector<int>> get_groups_map() const {
    std::unordered_map<int, int> rmp;
    std::unordered_map<int, std::vector<int>> res;
    REP (i, n()) {
      int x = root(i);
      res[root(i)].pb(i);
    }
    return res;
  }

  opa::utils::Remapper<int> get_repr() const {
    opa::utils::Remapper<int> res;
    REP (i, n()) {
      if (is_repr(i)) res.get(i);
    }
    return res;
  }

  bool operator==(const UnionJoin &other) const {
    if (n() != other.n()) return false;
    return get_clusters() == other.get_clusters();
  }
  OPA_ACCESSOR_R(int, m_count, count);

private:
  int m_count;
  mutable std::vector<int> m_tb;
};

OPA_NAMESPACE_DECL2_END

#pragma once

#include <opa_common.h>
#include <opa/or/search.h>

OPA_NAMESPACE(opa, OR)

template <class T, class Comparator = std::less<T> >
class BeamSearch : public Search<T> {
public:
  struct Params {
    typename Search<T>::SearchParams search_params;
    Comparator comparator;
    int node_lim = 100;
    bool debug = false;

    Params(typename Search<T>::SearchParams search_params,
           Comparator comparator, int node_lim) {
      this->search_params = search_params;
      this->comparator = comparator;
      this->node_lim = node_lim;
    }
    Params() {}
  };

  BeamSearch() {}
  BeamSearch(const Params &params) { init(params); }

  virtual void init(const Params &params) {
    Search<T>::init(params.search_params);
    m_params = params;
  }
  virtual void add(const T &a) override { next.pb(a); }
  virtual void start() override { go(); }

private:
  void go() {
    int round = 0;
    while (next.size() && !this->should_stop()) {
      std::swap(next, cur);
      sort(ALL(cur), m_params.comparator);
      reverse(ALL(cur));
      if (cur.size() > m_params.node_lim)
        cur.resize(m_params.node_lim);
      if (m_params.debug) {
        OPA_DISP("BeamSearch >> new step ", round, cur.size());
        REP (i, cur.size())
          OPA_DISP0(i, cur[i]);
        puts("");
      }
      next.clear();
      for (auto &x : cur) {
        if (this->should_stop())
          break;
        this->func()(this, x);
      }
      ++round;
    }
  }

  std::vector<T> cur, next;
  Params m_params;
};

OPA_NAMESPACE_END(opa, OR)

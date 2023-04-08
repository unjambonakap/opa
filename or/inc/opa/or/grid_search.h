#pragma once
#include <absl/types/span.h>
#include <opa/or/or_common.h>

OPA_NAMESPACE_OR

struct SubsetSelectorRes {
  bool stop;
  bool stop_branch;
};
template <class T> class SubsetSelector {
public:
  typedef std::function<SubsetSelectorRes(const std::vector<T> &ids)> CB;

  SubsetSelector(const std::vector<T> &elems, const CB &cb, int num_picks = -1)
      : m_elems(elems), m_cb(cb), num_picks(num_picks) {
    go(0);
  }

private:
  void go(int pos) {
    if (m_stop) return;
    if (m_picks.size() + m_elems.size() - pos < num_picks) return;
    if (pos == 0) {
      if (maybe_call_stop_branch()) return;
    }

    if (pos == m_elems.size()) return;
    go(pos + 1);
    m_picks.push_back(m_elems[pos]);
    if (!maybe_call_stop_branch()) go(pos + 1);
    m_picks.pop_back();
  }

  // Returns whether to keep branching here
  bool maybe_call_stop_branch() {
    if (m_stop) return true;
    if (num_picks != -1 && m_picks.size() != num_picks) return false;
    SubsetSelectorRes res = m_cb(m_picks);
    m_stop = res.stop;
    return res.stop_branch || num_picks != -1;
  }

  bool m_stop = false;
  std::vector<T> m_picks;
  const CB &m_cb;
  int num_picks;
  const std::vector<T> &m_elems;
};

template <class T> class CrossProdGen {
public:
  typedef std::function<void(const std::vector<T> &state)> CB;

  CrossProdGen(const std::vector<std::vector<T> > &space, const CB &cb)
      : m_space(space), m_cb(cb) {
    m_state.resize(space.size());
    go(0);
  }

private:
  void go(int pos) {
    if (pos == m_state.size()) {
      m_cb(m_state);
      return;
    }
    for (const auto &e : m_space[pos]) {
      m_state[pos] = e;
      go(pos + 1);
    }
  }

  std::vector<T> m_state;
  const CB &m_cb;
  const std::vector<std::vector<T> > &m_space;
};

struct GridSearchResult {
  std::vector<double> state;
  double score;
};

struct GridSearchBounds {
  std::vector<std::vector<double> > bounds;
  static GridSearchBounds
  Build(const std::vector<std::tuple<double, double, int> > &tb) {
    GridSearchBounds res;
    for (auto &e : tb) {
      res.bounds.push_back(utils::Range<double>::Build_range(
                             std::get<0>(e), std::get<1>(e), std::get<2>(e))
                             .tb());
    }
    return res;
  }
};

typedef std::function<double(const std::vector<double> &state)>
  GridSearchFunc_t;
// minimize
struct GridSearchFunc {
  GridSearchFunc_t func;
  std::function<void(const GridSearchBounds &bounds)> init = nullptr;
};

GridSearchResult DoGridSearch(const GridSearchBounds &bounds,
                              const GridSearchFunc &func);

GridSearchResult do_1d_grid_search(const std::vector<double> &range,
                                   GridSearchFunc_t func);

class DspBinaryMatcher {
public:
  DspBinaryMatcher(const float *data, int n);
  DspBinaryMatcher(const std::vector<float> &data);

  GridSearchFunc get_grid_search_func() const;

  void grid_search_init_func(const GridSearchBounds &bounds) const;
  double grid_search_func_impl(const std::vector<double> &state) const;

private:
  absl::Span<const float> m_data;
  mutable std::vector<float> m_precomp;
};

OPA_NAMESPACE_OR_END

#include "opa/or/grid_search.h"

OPA_NAMESPACE(opa, OR)

GridSearchResult DoGridSearch(const GridSearchBounds &bounds,
                              const GridSearchFunc &func) {

  double minv;
  bool first = true;
  std::vector<double> best_state;
  if (func.init) {
    func.init(bounds);
  }
  CrossProdGen<double>(bounds.bounds, [&](const std::vector<double> &state) {
    double res = func.func(state);
    if (first || minv > res) {
      best_state = state;
      minv = res;
      first = false;
    }
  });

  GridSearchResult result;
  result.state = best_state;
  result.score = minv;

  return result;
}

DspBinaryMatcher::DspBinaryMatcher(const std::vector<float> &data)
    : m_data(data) {}

DspBinaryMatcher::DspBinaryMatcher(const float *data, int n)
    : m_data(data, n) {}

GridSearchFunc DspBinaryMatcher::get_grid_search_func() const {
  GridSearchFunc func;
  func.func = [&](const std::vector<double> &state) {
    return this->grid_search_func_impl(state);
  };
  func.init = [&](const GridSearchBounds &bounds) {
    return this->grid_search_init_func(bounds);
  };
  return func;
}

void DspBinaryMatcher::grid_search_init_func(
  const GridSearchBounds &bounds) const {
  double max_sps = bounds.bounds[0].back();

  m_precomp.resize(m_data.size() + 1);
  m_precomp[0] = 0;
  for (int i = 0; i < m_data.size(); ++i) {
    m_precomp[i + 1] = m_precomp[i] + m_data[i];
  }
  for (int i = 0; i < max_sps + 1; ++i) m_precomp.push_back(m_precomp.back());
}

double DspBinaryMatcher::grid_search_func_impl(
  const std::vector<double> &state) const {
  OPA_CHECK0(state.size() == 2);
  double sps = state[0];
  double phase = state[1] * sps;
  double score = 0;
  score += std::abs(m_precomp[round(sps)]);
  for (double pos = phase; pos < m_data.size(); pos += sps) {
    double v = m_precomp[pos + sps] - m_precomp[pos];
    score += std::abs(v);
  }
  return -score;
}

GridSearchResult do_1d_grid_search(const std::vector<double> &range,
                                   GridSearchFunc_t func) {
  GridSearchFunc tmp;
  tmp.func = func;
  return DoGridSearch(GridSearchBounds{ { range } }, tmp);
}

OPA_NAMESPACE_END(opa, OR)

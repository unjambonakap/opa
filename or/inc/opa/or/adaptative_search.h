#pragma once
#include <opa/or/manifold.h>
#include <opa/or/or_common.h>
#include <opa/utils/clone.h>

OPA_NAMESPACE_OR

template <typename PointType, typename ScoreType = double>
struct PointAndScore {
  PointType point;
  ScoreType score;
};

template <typename Base> class CloneWrapper {
public:
  template <typename T> CloneWrapper(const T &a) { m_base.reset(a.clone()); }
  CloneWrapper(const CloneWrapper<Base> &a) { m_base.reset(a.m_base->clone()); }

  Base *get() const {
    OPA_CHECK0(m_base);
    return m_base.release();
  }
  mutable UPTR(Base) m_base;
};

struct SearchRoundDataForStopping {
  double new_best;
  double point_dist;
};

class SearchStoppingCriteria : public utils::Clonable<SearchStoppingCriteria> {
public:
  void new_round(const SearchRoundDataForStopping &data) {
    ++m_round;
    this->new_round_impl(data);
  }
  void reset() {
    m_round = 0;
    this->reset_impl();
  }
  // yeah, maybe just use uptr
  CloneWrapper<SearchStoppingCriteria> to_wrapper() {
    return CloneWrapper<SearchStoppingCriteria>(*this);
  }

  virtual bool should_stop() const = 0;
  int round() const { return m_round; }

protected:
  virtual void reset_impl() {}
  virtual void new_round_impl(const SearchRoundDataForStopping &data) {}
  int m_round = 0;
};

class TimeStoppingCriteria : public SearchStoppingCriteria {
public:
  TimeStoppingCriteria(int usec) : m_usec(usec) {}
  virtual void reset_impl() override {
    start = std::chrono::steady_clock::now();
  }
  virtual bool should_stop() const override {
    return (std::chrono::duration_cast<std::chrono::microseconds>(
              std::chrono::steady_clock::now() - start)
              .count() >= m_usec);
  }
  OPA_DECL_CLONE(SearchStoppingCriteria, TimeStoppingCriteria);

private:
  u64 m_usec;
  std::chrono::steady_clock::time_point start;
};

class RoundStoppingCriteria : public SearchStoppingCriteria {
public:
  RoundStoppingCriteria(int n_round) : m_n_round(n_round) {}
  virtual bool should_stop() const override { return m_round >= m_n_round; }
  OPA_DECL_CLONE(SearchStoppingCriteria, RoundStoppingCriteria);

private:
  int m_n_round;
};

class MinimaStoppingCriteria : public SearchStoppingCriteria {
public:
  OPA_DECL_CLONE(SearchStoppingCriteria, MinimaStoppingCriteria);
  MinimaStoppingCriteria(int window, double stop_ratio = 1e-6)
      : m_window(window), m_stop_ratio(stop_ratio) {
    reset_impl();
  }
  // TODO: maybe create RollingAverage class
  virtual void new_round_impl(const SearchRoundDataForStopping &data) override {
    if (m_round >= 0) {
      double diff = data.new_best - m_last;
      diff /= m_window;
      diff = std::abs(diff);
      m_sum += diff;
      m_sum -= m_vals[0];
      m_vals.pop_front();
      m_vals.push_back(diff);
    }
    m_last = data.new_best;
  }
  virtual bool should_stop() const override {
    return m_round > m_window && m_sum / (1e-6 + m_last) < m_stop_ratio - 6;
  }
  virtual void reset_impl() override {
    m_vals = std::deque<double>(m_window, 0);
    m_sum = 0;
  }

private:
  std::deque<double> m_vals;
  double m_sum;
  double m_stop_ratio = 1e-6;
  double m_last;
  int m_window;
};

class AggregateStopCriteria : public SearchStoppingCriteria {
public:
  OPA_DECL_CLONE(SearchStoppingCriteria, AggregateStopCriteria);

  AggregateStopCriteria(const AggregateStopCriteria &peer) {
    for (const auto &e : peer.m_criterias) {
      m_criterias.emplace_back(e->clone());
    }
  }

  AggregateStopCriteria(
    const std::vector<CloneWrapper<SearchStoppingCriteria> > &criterias) {
    for (auto &e : criterias) {
      m_criterias.emplace_back(e.get());
    }
  }

  virtual void new_round_impl(const SearchRoundDataForStopping &data) override {
    for (auto &e : m_criterias) {
      e->new_round(data);
    }
  }
  virtual bool should_stop() const override {
    for (auto &e : m_criterias) {
      if (e->should_stop()) return true;
    }
    return false;
  }
  virtual void reset_impl() override {
    for (auto &e : m_criterias) {
      e->reset();
    }
  }

private:
  std::vector<UPTR(SearchStoppingCriteria)> m_criterias;
};

template <typename PointType> class AdaptativeStrategy {
public:
  virtual ~AdaptativeStrategy() {}

  virtual std::vector<PointType> get_next_vec() const = 0;
  virtual bool done() = 0;
  virtual void push(const PointType &point, double val) = 0;
  virtual void start_round();
  virtual void finish_round() = 0;
  virtual PointAndScore<PointType> get_best() const = 0;

protected:
};

template <typename PointType, typename MetricDistance>
class SimpleAdaptativeStrategy {
public:
  SimpleAdaptativeStrategy(SearchStoppingCriteria &stop_criteria,
                           const MetricDistance &metric_distance)
      : m_stop_criteria(stop_criteria), m_metric_distance(metric_distance) {}

  virtual void push(const PointType &point, double val) {
    m_data.emplace_back(PointAndScore<PointType>{ point, val });
  }
  virtual void start_round() {
    ++m_round;
    m_data.clear();
  }
  virtual void finish_round() {
    m_best_pos = 0;
    OPA_CHECK0(m_data.size() > 0);
    REP (i, m_data.size()) {
      if (m_data[i].score < m_data[m_best_pos].score) m_best_pos = i;
    }
    double new_dist = 0;
    m_prev_best = m_new_best;
    m_new_best = m_data[m_best_pos];
    if (m_round != 1) {
      new_dist =
        m_metric_distance.compute_dist(m_new_best.point, m_prev_best.point);
    }

    m_stop_criteria.new_round(
      SearchRoundDataForStopping{ m_data[m_best_pos].score, new_dist });
  }

  virtual PointAndScore<PointType> get_best() const { return m_new_best; }
  virtual bool done() {
    if (m_round == 0) return false;
    return m_stop_criteria.should_stop();
  }

protected:
  int m_round = 0;
  PointAndScore<PointType> m_new_best;
  PointAndScore<PointType> m_prev_best;
  SearchStoppingCriteria &m_stop_criteria;
  const MetricDistance &m_metric_distance;
  std::vector<PointAndScore<PointType> > m_data;
  int m_best_pos;
};

template <typename PointType, typename ManifoldType> class ConstRManagedDesc {
public:
  ConstRManagedDesc(const ManifoldType &manifold, int init_grid_count,
                    int nb_count)
      : m_manifold(manifold), m_init_grid_count(init_grid_count),
        m_nb_count(nb_count) {}

  std::vector<PointType> get_next_vec(const PointType &prev_best,
                                      const PointType &new_best,
                                      int round) const {
    if (round == 0) {
      m_last_req = Manifold1D::Request{ 0., true, 0., m_init_grid_count };
    } else {
      double cur_precision = m_manifold.average_dist(m_last_req);
      m_last_req =
        Manifold1D::Request{ new_best, false, { cur_precision, m_nb_count } };
    }
    return m_manifold.sample(m_last_req);
  }

  mutable Manifold1D::Request m_last_req;
  const ManifoldType &m_manifold;
  int m_init_grid_count;
  int m_nb_count;
};

template <typename PointType, typename ManifoldType>
class ConstNDimManagedDesc {
public:
  ConstNDimManagedDesc(const ManifoldType &manifold, int init_grid_count,
                       int nb_count, PointType pt = PointType(),
                       std::vector<int> dims = {})
      : m_manifold(manifold), m_init_grid_count(init_grid_count),
        m_nb_count(nb_count) {
    if (dims.size() == 0) {
      REP (i, m_manifold.ndim())
        dims.push_back(i);
    }

    m_last_req.pt = pt;
    m_last_req.data.dims.resize(dims.size());
    REP (i, dims.size())
      m_last_req.data.dims[i].dim_id = dims[i];
  }

  std::vector<PointType> get_next_vec(const PointType &prev_best,
                                      const PointType &new_best,
                                      int round) const {

    if (round == 0) {
      m_last_req.is_grid = true;
      for (auto &e : m_last_req.data.dims) e.req.count = m_init_grid_count;
    } else {
      std::vector<double> cur_precision = m_manifold.average_dist(m_last_req);
      m_last_req.is_grid = false;
      m_last_req.pt = new_best;
      REP (i, m_last_req.data.dims.size()) {
        auto &cur = m_last_req.data.dims[i];
        cur.req.count = m_nb_count;
        cur.req.region = cur_precision[i];
      }
    }
    return m_manifold.sample(m_last_req);
  }

  mutable typename NManifold<PointType>::Request m_last_req;
  const ManifoldType &m_manifold;
  int m_init_grid_count;
  int m_nb_count;
};

template <typename PointType, typename MetricDistance, typename ManagedDesc>
class ManagedSimpleAdaptativeStrategy
  : public SimpleAdaptativeStrategy<PointType, MetricDistance> {
public:
  ManagedSimpleAdaptativeStrategy(const ManagedDesc &managed_desc,
                                  SearchStoppingCriteria &stop_criteria,
                                  const MetricDistance &metric_distance)
      : SimpleAdaptativeStrategy<PointType, MetricDistance>(stop_criteria,
                                                            metric_distance),
        m_managed_desc(managed_desc) {}
  virtual std::vector<PointType> get_next_vec() const {
    return m_managed_desc.get_next_vec(this->m_prev_best.point,
                                       this->m_new_best.point, this->m_round);
  }

private:
  ManagedDesc m_managed_desc;
};

template <class PointType, typename Strategy> class AdaptativeSearcher {
public:
  typedef std::function<double(const PointType &pt)> EvalFunc;

  AdaptativeSearcher(const EvalFunc &eval_func, const Strategy &strategy)
      : m_eval_func(eval_func), m_strategy(strategy) {}

  void process_vec(const std::vector<PointType> &tb) {
    std::vector<double> res;
    res.resize(tb.size());
    m_strategy.start_round();
    REP (i, tb.size()) {
      res[i] = m_eval_func(tb[i]);
      m_strategy.push(tb[i], res[i]);
    }
    m_strategy.finish_round();
  }

  PointAndScore<PointType> find_minimum() {
    while (!m_strategy.done()) {
      this->process_vec(m_strategy.get_next_vec());
    }
    return m_strategy.get_best();
  }

private:
  EvalFunc m_eval_func;
  Strategy m_strategy;
};

template <class PointType> class IterativeSearcher {
  typedef double (*DistFunc)(const PointType &, const PointType &);

public:
  struct IterationParams {
    int round;
  };
  typedef std::function<PointAndScore<PointType>(const PointType &pt,
                                                 const IterationParams &params)>
    EvalFunc;

  IterativeSearcher(const EvalFunc &eval_func,
                    const SearchStoppingCriteria &stop_criteria,
                    const DistFunc &metric_distance)
      : m_eval_func(eval_func), m_stop_criteria(stop_criteria.clone()),
        m_metric_distance(metric_distance) {}

  PointAndScore<PointType> find_minimum(PointType pos) {
    PointAndScore<PointType> cur;
    cur.point = pos;
    IterationParams params;
    params.round = 0;
    while (!m_stop_criteria->should_stop()) {
      PointType prev = cur.point;
      cur = m_eval_func(prev, params);
      double new_dist = m_metric_distance(cur.point, prev);
      m_stop_criteria->new_round(
        SearchRoundDataForStopping{ cur.score, new_dist });
      ++params.round;
    }
    return cur;
  }

private:
  EvalFunc m_eval_func;
  UPTR(SearchStoppingCriteria) m_stop_criteria;
  DistFunc m_metric_distance;
};

OPA_NAMESPACE_OR_END

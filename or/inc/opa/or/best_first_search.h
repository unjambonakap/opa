#pragma once

#include <opa_common.h>
#include <opa/or/search.h>

OPA_NAMESPACE(opa, OR)

template <class T, class Comparator = std::less<T> >
class BestFirstSearch : public Search<T> {
public:
  struct Params {
    typename Search<T>::SearchParams search_params;
    Comparator comparator;
    Params(typename Search<T>::SearchParams search_params,
           Comparator comparator) {
      this->search_params = search_params;
      this->comparator = comparator;
    }
    Params() {}
  };

  BestFirstSearch() {}
  BestFirstSearch(const Params &params) { init(params); }


  void init2() {
    Search<T>::init(m_params.search_params);
    m_q = std::priority_queue<T, std::vector<T>, Comparator>(m_params.comparator);
  }
  virtual void init(const Params &params) {
    m_params = params;
    init2();
  }

  virtual void add(const T &a) override { m_q.push(a); }
  virtual void start() override { go(); }

private:
  void go() {
    while (m_q.size() && !this->should_stop()) {
      T top = m_q.top();
      m_q.pop();
      this->func()(this, top);
    }
  }

  std::priority_queue<T, std::vector<T>, Comparator> m_q;
  Params m_params;
};

OPA_NAMESPACE_END(opa, OR)


#pragma once
#include <opa/utils/misc.h>
#include <opa_common.h>
OPA_NAMESPACE(opa, OR)

template <class T> class Search : public opa::utils::Initable {
public:
  typedef std::function<void(Search<T> *self, const T &state)> Func;

  struct SearchParams {
    Func func;
  };

  Search(const SearchParams &params) {}
  Search() {}

  virtual void init(const SearchParams &params) {
    opa::utils::Initable::init();
    m_params = params;
    m_stop = false;
  }

  virtual void add(const T &a) = 0;
  virtual void start() = 0;
  void stop() { m_stop = true; }
  bool should_stop() const { return m_stop; }

protected:
  OPA_ACCESSOR_R(Func, m_params.func, func);

private:
  SearchParams m_params;
  bool m_stop;
};

template <class T> class Search_bfs : public Search<T> {
public:
  virtual void add(const T &a) override { m_q.push(a); }
  virtual void start() override { go(); }

private:
  void go() {
    while (m_q.size() && !this->should_stop()) {
      T top = m_q.front();
      m_q.pop();
      this->func()(this, top);
    }
  }

  std::queue<T> m_q;
};

template <class T> class Search_dfs : public Search<T> {
public:
  virtual void add(const T &a) override {
    if (!this->should_stop()) this->func()(this, a);
  }
  virtual void start() override {}
};
OPA_NAMESPACE_END(opa, OR)

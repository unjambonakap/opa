#pragma once

#include <opa_common.h>
#include <opa/threading/auto_job.h>
#include <opa/threading/dispatcher.h>
#include <opa/threading/runner.h>
//#include <experimental/tuple>

OPA_NAMESPACE(opa, threading)

template <class Ret, class... Args>
struct FuncObj : public opa::utils::BaseStorable {
public:
  typedef std::function<Ret(Args...)> FuncType;
  FuncObj(const FuncType &func) { m_func = func; }

  FuncType m_func;
};

template <class Ret, class... Args>
struct MapperData : public opa::utils::BaseStorable {
  typedef FuncObj<Ret, Args...> FuncType;
  MapperData() {}

  OPA_TGEN_IMPL(func);
  SPTR(FuncType) func;
  int useless;
};

template <class... Args> struct MapperArgs : public opa::utils::ProtobufParams {
  MapperArgs(std::tuple<Args...> data) { this->data = data; }
  MapperArgs() {}
  std::tuple<Args...> data;
  OPA_TGEN_IMPL(data);
};

template <class Ret, class... Args>
class MapJob : public AutoVectorJob<MapperArgs<Args...>, Ret> {
public:
  virtual void auto_worker_do_work(const MapperArgs<Args...> &data,
                                   Ret &out_res) override {
    // OPA_DISP("CALLING FUN ", uintptr_t(m_mapper_data.func.get()));
    assert(0);
    //out_res = std::experimental::apply(m_mapper_data.func->m_func, data.data);
  }

  virtual void auto_server_get_work() override {
    REP (i, m_args.size()) {
      // OPA_DISP("PUSHING ", i, m_args.size());
      bool out_more;
      this->cb()(m_args[i], out_more);
    }
  }

  OPA_TGEN_IMPL(m_mapper_data);

  std::vector<MapperArgs<Args...> > m_args;
  MapperData<Ret, Args...> m_mapper_data;

  OPA_CLOUDY_JOB_DECL
};

template <class Ret, class... Args> JobId MapJob<Ret, Args...>::StaticJobId;
#define OPA_CLOUDY_DECL_MAPPER(clName, func, ...)                              \
  class clName : public opa::threading::FuncObj<__VA_ARGS__> {                 \
  public:                                                                      \
    clName() : opa::threading::FuncObj<__VA_ARGS__>(func) {}                   \
    static UPTR(clName) getCl() {                                              \
      return UPTR(clName)(OPACS_GET(clName, clName));                          \
    }                                                                          \
  };                                                                           \
  OPA_CLASS_STORE_REGISTER_BY_KEY(clName, clName)

template <class Orig, class Ret, class... Args>
std::vector<Ret> do_map_job(Dispatcher *dispatcher, SPTR(Orig) func,
                            const std::vector<MapperArgs<Args...> > &args) {
  typedef MapJob<Ret, Args...> SelfType;

  UPTR(SelfType) job(Runner::GetJobById<SelfType>(SelfType::StaticJobId));

  job->m_args = args;
  job->m_mapper_data.func =
    std::static_pointer_cast<FuncObj<Ret, Args...>, Orig>(func);
  dispatcher->process_job(*job);
  return job->res_list();
}

OPA_NAMESPACE_END(opa, threading)

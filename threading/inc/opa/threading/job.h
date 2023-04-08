#pragma once

#include <opa/threading/data.h>
#include <opa_common.h>

OPA_NAMESPACE_DECL2(opa, threading)

typedef std::function<DataId(JobMsgPtr work, bool &out_more)> ServerGetWorkFunc;

class Job {
public:
  Job(){

  }
  virtual ~Job() {}
  virtual void reset() {}

  // worker part
  virtual void worker_initialize(const JobMsg &init) = 0;
  virtual void worker_finalize() {}

  virtual JobMsgPtr worker_do_work(const JobMsg &data) = 0;

  // dispatcher part
  virtual JobMsgPtr server_initialize() = 0;
  virtual void server_get_work(const ServerGetWorkFunc &cb) = 0;
  virtual void server_finalize() {}

  virtual void server_set_work_result(const JobMsg &res) = 0;

  virtual bool server_want_more_results() const { return true; }

  JobId get_job_id() const;
  void set_job_id(JobId job_id) { m_job_id = job_id; }

  Job &get_job() { return *this; }

#define OPA_CLOUDY_JOB_DECL                                                    \
  static opa::threading::JobId StaticJobId;                                    \
  static std::string JobName;


#define OPA_CLOUDY_JOB_IMPL(cl)                                                \
  opa::threading::JobId cl::StaticJobId;                                       \
  std::string cl::JobName;
#define OPA_CLOUDY_JOB_IMPL_TMPL(cl, ...)                                      \
  template <> opa::threading::JobId cl<__VA_ARGS__>::StaticJobId;              \
  template <> std::string cl<__VA_ARGS__>::JobName;

private:
  JobId m_job_id = Invalid_JobId;
};
OPA_DECL_SPTR(Job, JobPtr);

template <class TInit, class TData, class TRes> class TJob : public Job {
public:
  virtual void tworker_initialize(const JobMsg &base_msg,
                                  const TInit &init) = 0;
  virtual void tworker_do_work(const JobMsg &base_msg, const TData &data,
                               TRes &out_res) = 0;
  virtual void tserver_initialize(TInit &out_init) = 0;
  virtual void tserver_get_work(
    const std::function<DataId(const TData &data, bool &out_more)> &cb) = 0;
  virtual void tserver_set_work_result(const JobMsg &base_msg, const TRes &res,
                                       DataId data_id) = 0;

  virtual void worker_initialize(const JobMsg &init) {
    TInit init_msg;
    OPA_ASSERT0(init.any().Is<TInit>());
    init.any().UnpackTo(&init_msg);
    tworker_initialize(init, init_msg);
  }

  virtual JobMsgPtr server_initialize() {
    auto msg = new JobMsg;
    TInit init_msg;
    tserver_initialize(init_msg);
    msg->mutable_any()->PackFrom(init_msg);
    return JobMsgPtr(msg);
  }

  virtual JobMsgPtr worker_do_work(const JobMsg &data) {
    auto res_msg = new JobMsg;
    TData data2;
    data.any().UnpackTo(&data2);
    TRes res2;

    tworker_do_work(data, data2, res2);
    res_msg->mutable_any()->PackFrom(res2);
    return JobMsgPtr(res_msg);
  }

  virtual void server_get_work(
    const std::function<DataId(JobMsgPtr work, bool &out_more)> &cb) {
    tserver_get_work([cb](const TData &data, bool &out_more) -> DataId {
      auto res = JobMsgPtr(new JobMsg);
      res->mutable_any()->PackFrom(data);
      return cb(res, out_more);
    });
  }

  virtual void server_set_work_result(const opa::threading::JobMsg &res) {
    TRes res_msg;
    res.any().UnpackTo(&res_msg);
    tserver_set_work_result(res, res_msg, res.data_id());
  }
};

template <class TInit, class TData, class TRes>
class VectorJob : public TJob<TInit, TData, TRes> {
public:
  virtual void reset() { m_res_list.clear(); }

  virtual void tserver_set_work_result(const JobMsg &base_msg, const TRes &res,
                                       DataId data_id) {
    if (m_res_list.size() <= data_id) m_res_list.resize(data_id + 1);
    m_res_list[data_id] = res;
  }

protected:
  std::vector<TRes> m_res_list;
};

template <class TInit, class TData, class TRes>
class FindOneJob : public TJob<TInit, TData, TRes> {
public:
  virtual void reset() { m_found = false; }

  virtual void tserver_set_work_result(const JobMsgPtr &base_msg,
                                       const TRes &res, DataId data_id) {
    if (tserver_handle_res(base_msg, res)) m_found = true;
  }

  virtual bool server_want_more_results() const { return !m_found; }

  virtual bool tserver_handle_res(const JobMsgPtr &base_msg,
                                  const TRes &res) = 0;

protected:
  bool m_found;
};
OPA_NAMESPACE_DECL2_END

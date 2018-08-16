#pragma once

#include <opa_common.h>
#include <opa/threading/data.h>
#include <opa/threading/msg_desc.pb.h>
#include <opa/threading/job.h>
#include <opa/threading/auto_loader.h>

OPA_NAMESPACE_DECL2(opa, threading)

template <class TData, class TRes>
class AutoJob
  : public TJob<opa::utils::AnyMsg, opa::utils::AnyMsg, opa::utils::AnyMsg>,
    public opa::utils::ProtobufParams {

public:
  virtual void tworker_initialize(const JobMsg &base_msg,
                                  const opa::utils::AnyMsg &init) override {
    load(init);
  }

  virtual void tworker_do_work(const JobMsg &base_msg,
                               const opa::utils::AnyMsg &data,
                               opa::utils::AnyMsg &out_res) override {
    TData x;
    TRes res;
    x.load(data);
    auto_worker_do_work(x, res);
    opa::utils::anymsg_store(res, out_res);
  }

  virtual void tserver_initialize(opa::utils::AnyMsg &out_init) override {
    store(out_init);
  }

  virtual void tserver_get_work(
    const std::function<opa::threading::DataId(const opa::utils::AnyMsg &data,
                                               bool &out_more)> &cb) override {

    m_cb = [this, cb](const TData &data, bool &out_more) {
      opa::utils::AnyMsg x;
      data.store(x);
      return cb(x, out_more);
    };
    auto_server_get_work();
  }

  virtual void
  tserver_set_work_result(const JobMsg &base_msg, const opa::utils::AnyMsg &res,
                          opa::threading::DataId data_id) override {
    TRes x;
    opa::utils::anymsg_load(x, res);
    auto_server_set_work_result(x, data_id);
  }

  virtual void auto_worker_initialize() {}
  virtual void auto_worker_do_work(const TData &data, TRes &out_res) = 0;
  virtual void auto_server_get_work() = 0;
  virtual void auto_server_set_work_result(const TRes &res,
                                           opa::threading::DataId data_id) = 0;

  OPA_ACCESSOR(
    std::function<opa::threading::DataId(const TData &data, bool &out_more)>,
    m_cb, cb);

private:
  std::function<opa::threading::DataId(const TData &data, bool &out_more)> m_cb;
};

template <class AData, class ARes>
class AutoFindOneJob : public AutoJob<AData, ARes> {
public:
  virtual void reset() override { m_found = false; }

  virtual void
  auto_server_set_work_result(const ARes &res,
                              opa::threading::DataId data_id) override {
    if (auto_server_handle_res(res)) {
      m_result = res;
      m_found = true;
    }
  }

  virtual bool server_want_more_results() const override { return !m_found; }

  virtual bool auto_server_handle_res(const ARes &res) = 0;
  OPA_ACCESSOR(ARes, m_result, result);

protected:
  ARes m_result;
  bool m_found;
};

template <class AData, class ARes>
class AutoVectorJob : public AutoJob<AData, ARes> {
public:
  virtual void reset() override { m_res_list.clear(); }

  virtual void
  auto_server_set_work_result(const ARes &res,
                              opa::threading::DataId data_id) override {
    if (m_res_list.size() <= data_id)
      m_res_list.resize(data_id + 1);
    m_res_list[data_id] = res;
  }
  OPA_ACCESSOR(std::vector<ARes>, m_res_list, res_list);

protected:
  std::vector<ARes> m_res_list;
};

template <class AData, class ARes>
class AutoCollectJob : public AutoJob<AData, std::vector<ARes> > {
public:
  virtual void reset() override { m_res_list.clear(); }

  virtual void
  auto_server_set_work_result(const std::vector<ARes> &res,
                              opa::threading::DataId data_id) override {
    m_res_list.insert(m_res_list.end(), ALL(res));
  }

  OPA_ACCESSOR(std::vector<ARes>, m_res_list, res_list);

protected:
  std::vector<ARes> m_res_list;
};

OPA_NAMESPACE_DECL2_END

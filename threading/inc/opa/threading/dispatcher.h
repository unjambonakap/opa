#pragma once

#include <opa_common.h>
#include <opa/threading/job.h>
#include <opa/threading/queue.h>
#include <opa/threading/threadable.h>
#include <opa/threading/data.h>

OPA_NAMESPACE_DECL2(opa, threading)

namespace internal {
class Comm;
}
class Job;
class Runner;

class Dispatcher : private Threadable {
public:
  Dispatcher();
  ~Dispatcher();

  void initialize(Runner *runner, const std::string &server_info);
  void process_job(Job &job);
  int nthread() const;

  OPA_ACCESSOR_PTR(Runner, m_runner, runner);

private:
  virtual void run();
  DataId get_work_cb(JobMsgPtr work, bool &out_more);
  void send_msg(JobMsg &msg);

  bool m_has_more;
  Queue<std::pair<DataId, JobMsgPtr> > m_data_queue;
  JobMsgPtr m_init_msg;
  JobId m_job_id;
  DataId m_data_id;
  Job *m_job;
  internal::Comm *m_comm;
  std::mutex m_mutex;
  std::string m_server_info;
  JobNonce m_job_nonce;
  std::string m_sys_id;
  Runner *m_runner;
};

OPA_NAMESPACE_DECL2_END

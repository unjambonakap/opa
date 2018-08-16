#pragma once

#include <opa_common.h>
#include <opa/threading/threadable.h>
#include <opa/threading/job.h>
#include <opa/threading/queue.h>

OPA_NAMESPACE_DECL2(opa, threading)

namespace internal {
class Comm;
}

class ClientDispatcher : public Threadable {
  public:
    ClientDispatcher();
    ~ClientDispatcher();
    void initialize(const std::string &server_info, int nthread);
    virtual void run() override;
    virtual void stop() override;


    DispatcherState proc_next(JobId job_id, JobNonce job_nonce,
                              const std::function<JobMsgPtr(JobMsg *data)> &func);
    DispatcherState proc_init(const std::function<void(JobMsg *init)> &func);

  private:
    JobMsgPtr do_comm(JobMsg &msg);

    int m_nthread;
    std::string m_server_info;
    JobMsgPtr m_init_data;
    Queue<JobMsgPtr> m_data_queue;
    JobId m_job_id;
    internal::Comm *m_comm;
    std::mutex m_mutex;
    JobNonce m_job_nonce;
    std::string m_sys_id;
};

OPA_NAMESPACE_DECL2_END

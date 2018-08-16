#pragma once

#include <opa_common.h>
#include <opa/threading/threadable.h>
#include <opa/threading/job.h>

OPA_NAMESPACE_DECL2(opa, threading)

class ClientDispatcher;

class Worker : public Threadable {
  public:
    virtual void run();
    Worker(ClientDispatcher &dispatcher);
    void do_init(JobMsg *init);
    JobMsgPtr do_next(JobMsg *data);

  private:
    ClientDispatcher &m_dispatcher;
    JobPtr m_job;
    JobNonce m_job_nonce;
};

OPA_NAMESPACE_DECL2_END

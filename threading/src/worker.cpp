#include "worker.h"
#include "runner.h"
#include "client_dispatcher.h"
#include "job.h"
#include <opa/math/common/rng.h>

OPA_NAMESPACE_DECL2(opa, threading)

Worker::Worker(ClientDispatcher &dispatcher) : m_dispatcher(dispatcher) {}

void Worker::run() {
  opa::math::common::reseed_rng();


  while (!want_stop()) {

    DispatcherState state;
    state =
      m_dispatcher.proc_init([this](JobMsg *init) { this->do_init(init); });
    if (state == DispatcherState::Done)
      break;

    if (state != DispatcherState::Ok) {
      usleep(1e5);
      continue;
    }

    while (true) {
      state = m_dispatcher.proc_next(
        m_job->get_job_id(), m_job_nonce,
        [this](JobMsg *data) { return this->do_next(data); });

      if (state == DispatcherState::ChangedJob)
        break;
      else if (state == DispatcherState::Done)
        break;
    }
  }
}
void Worker::do_init(JobMsg *init) {
  m_job.reset(Runner::GetJob(init->job_id()));
  m_job_nonce = init->job_nonce();
  m_job->worker_initialize(*init);
}

JobMsgPtr Worker::do_next(JobMsg *data) { return m_job->worker_do_work(*data); }

OPA_NAMESPACE_DECL2_END

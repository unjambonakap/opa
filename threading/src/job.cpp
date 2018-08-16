#include <opa/threading/job.h>
#include <opa/threading/runner.h>

OPA_NAMESPACE(opa, threading)
JobId Job::get_job_id() const {
  OPA_CHECK0(Runner::CheckJobId(m_job_id));
  return m_job_id;
}

OPA_NAMESPACE_END(opa, threading)

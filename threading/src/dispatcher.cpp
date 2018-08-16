#include "dispatcher.h"
#include "internal/comm.h"
#include <opa/math/common/rng.h>
#include <opa/threading/runner.h>

using namespace std;

OPA_NAMESPACE_DECL2(opa, threading)

Dispatcher::Dispatcher() { m_comm = new internal::Comm(); }
Dispatcher::~Dispatcher() { delete m_comm; }

void Dispatcher::initialize(Runner *runner, const std::string &server_info) {
  m_runner = runner;
  m_server_info = server_info;
  m_comm->bind(m_server_info);
  m_job_nonce = opa::math::common::rng();
  m_sys_id = utils::get_process_fingerprint();
}

void Dispatcher::send_msg(JobMsg &msg) {
  msg.set_host(m_sys_id);
  m_comm->send(msg);
}

int Dispatcher::nthread() const { return m_runner->get_nthread(); }

void Dispatcher::run() {

  m_init_msg->set_job_id(m_job_id);
  m_init_msg->set_job_nonce(m_job_nonce);
  u64 res_count = 0;

  while (m_has_more || res_count < m_data_id) {
    auto msg = m_comm->recv();

    switch (msg->type()) {
    case MessageType::InitRequest:
      m_comm->send(*m_init_msg);
      break;
    case MessageType::DataRequest: {
      StatusCode status;
      bool ok = false;

      if (msg->job_id() != m_job_id || msg->job_nonce() != m_job_nonce)
        status = StatusCode::ChangedJob;
      else {

        std::pair<DataId, JobMsgPtr> data;

        if (!m_data_queue.wait(data)) {
          status = StatusCode::Ok;
        } else {
          ok = true;
          data.ND->set_data_id(data.ST);
          data.ND->set_job_id(m_job_id);
          data.ND->set_job_nonce(m_job_nonce);
          data.ND->set_same_vm(msg->host() == m_sys_id);
          data.ND->set_type(MessageType::DataResponse);
          m_comm->send(*data.ND);
        }
      }

      if (!ok)
        m_comm->send(Build_StatusMsg(status));
    } break;

    case MessageType::ResultOutput:

      // if not this, then it's an answer

      // OPA_ASSERT(res_msg->job_id == m_job_id && res_msg->job_nonce ==
      // m_job_nonce,
      //           "fuckup, msg_job_id=%d, wnat=%d", res_msg->job_id, m_job_id);
      ++res_count;
      m_job->server_set_work_result(*msg);

      m_comm->send(Build_StatusMsg(StatusCode::Ok));
      break;
    default:
      OPA_ABORT0(true);
      break;
    }

    if (!m_job->server_want_more_results())
      break;
  }
  m_data_queue.release();
  fflush(stdout);
}

DataId Dispatcher::get_work_cb(JobMsgPtr work, bool &out_more) {
  m_data_queue.push(MP(m_data_id, work));
  fflush(stdout);
  m_data_queue.wait_size([](int size) -> bool { return size < 100; });
  out_more = m_job->server_want_more_results();
  return m_data_id++;
}

void Dispatcher::process_job(Job &job) {
  m_job = &job;
  m_job_id = job.get_job_id();
  OPA_DISP("DISPATCHING JOb .. ", m_job_id);
  m_job->reset();
  m_has_more = true;
  m_data_queue.reset();
  m_data_id = 0;

  m_init_msg = m_job->server_initialize();
  m_init_msg->set_type(MessageType::InitResponse);

  start();

  using namespace std::placeholders;
  m_job->server_get_work(std::bind(&Dispatcher::get_work_cb, this, _1, _2));
  m_has_more = false;
  m_data_queue.release();
  join();

  m_job->server_finalize();
  ++m_job_nonce;
}

OPA_NAMESPACE_DECL2_END

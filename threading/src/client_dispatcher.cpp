#include "client_dispatcher.h"
#include "internal/comm.h"
#include <opa/threading/msg_desc.pb.h>

using namespace std;
using namespace opa::threading::internal;

OPA_NAMESPACE_DECL2(opa, threading)
const int N_BUF = 1;

ClientDispatcher::ClientDispatcher() {
  m_comm = new Comm();
  m_job_id = Invalid_JobId;
  m_init_data = 0;
}
ClientDispatcher::~ClientDispatcher() { delete m_comm; }

void ClientDispatcher::stop() {
  Threadable::stop();
  m_comm->release();
  m_data_queue.release();
}

void ClientDispatcher::initialize(const std::string &server_info, int nthread) {
  m_server_info = server_info;
  m_nthread = nthread;
  // OPA_DISP0("connect to ",server_info);
  m_comm->connect(server_info);
  // OPA_TRACE0;
  m_sys_id = utils::get_process_fingerprint();
}

void ClientDispatcher::run() {

  while (!want_stop()) {
    auto x = JobMsg();
    x.set_type(MessageType::InitRequest);

    auto msg = do_comm(x);
    if (!msg)
      break;

    if (msg->type() == MessageType::Status) {
      OPA_ABORT0(true);
      break;
    }
    OPA_CHECK(msg->type() == MessageType::InitResponse, "got %d %d\n",
              msg->type(), MessageType::InitResponse);

    {
      lock_guard<mutex> lk(m_mutex);
      m_job_id = msg->job_id();
      m_job_nonce = msg->job_nonce();
      m_init_data = msg;
    }

    while (true) {
      m_data_queue.wait_size([](int n) { return n <= N_BUF; });
      if (m_data_queue.released())
        break;

      JobMsgPtr data_msg;
      JobMsg req_msg;
      req_msg.set_type(MessageType::DataRequest);
      req_msg.set_job_id(m_job_id);
      req_msg.set_job_nonce(m_job_nonce);
      data_msg = do_comm(req_msg);
      if (!data_msg)
        break;

      if (data_msg->type() == MessageType::Status) {
        break;
      }

      OPA_CHECK0(data_msg->type() == MessageType::DataResponse);

      m_data_queue.push(data_msg);
    }

    m_data_queue.wait_size([](int n) { return n == 0; });
  }
  m_comm->close();
}

JobMsgPtr ClientDispatcher::do_comm(JobMsg &msg) {
  lock_guard<mutex> lk(m_mutex);
  msg.set_host(m_sys_id);
  m_comm->send(msg);
  auto res = m_comm->recv();
  return res;
}

DispatcherState ClientDispatcher::proc_next(
  JobId job_id, JobNonce job_nonce,
  const std::function<JobMsgPtr(JobMsg *data)> &func) {

  JobMsgPtr data;
  if (job_id != m_job_id || job_nonce != m_job_nonce)
    return DispatcherState::ChangedJob;

  if (!m_data_queue.wait(data))
    return DispatcherState::Done;

  if (data->job_id() != job_id || data->job_nonce() != job_nonce) {
    m_data_queue.push(data);
    return DispatcherState::ChangedJob;
  }

  auto res = func(data.get());
  {
    res->set_data_id(data->data_id());
    res->set_job_id(data->job_id());
    res->set_job_nonce(data->job_nonce());
    res->set_type(MessageType::ResultOutput);

    JobMsgPtr tmp = do_comm(*res);
    if (!tmp)
      return DispatcherState::Done;
    OPA_CHECK0(tmp->type() == MessageType::Status);
    StatusMsg status_msg;
    tmp->any().UnpackTo(&status_msg);

    OPA_ASSERT(status_msg.status() == StatusCode::Ok,
               "unexpected answer, got=%d, wnat=%d\n", status_msg.status(),
               StatusCode::Ok);
  }

  return DispatcherState::Ok;
}

DispatcherState
ClientDispatcher::proc_init(const std::function<void(JobMsg *init)> &func) {
  // OPA_TRACE0;
  lock_guard<mutex> lk(m_mutex);
  if (want_stop())
    return DispatcherState::Done;
  if (!m_init_data)
    return DispatcherState::Wait;
  // OPA_TRACE0;
  func(m_init_data.get());
  // OPA_TRACE0;
  return DispatcherState::Ok;
}

OPA_NAMESPACE_DECL2_END

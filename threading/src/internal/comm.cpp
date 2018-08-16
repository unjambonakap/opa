#include "comm.h"
#include "runner.h"
#include <google/protobuf/io/coded_stream.h>

#include <zmq.hpp>

using namespace std;

OPA_NAMESPACE_DECL3(opa, threading, internal)
Comm::Comm() { m_ctx = make_shared<zmq::context_t>(1); }
Comm::~Comm() {}

void Comm::bind(const std::string &server_info) {
  m_sock = make_shared<zmq::socket_t>(*m_ctx, ZMQ_REP);
  // OPA_TRACE("Binding to server >> %s", server_info.c_str());

  m_sock->bind(server_info.c_str());
  // OPA_TRACE0;
}

void Comm::connect(const std::string &server_info) {
  m_sock = make_shared<zmq::socket_t>(*m_ctx, ZMQ_REQ);
  // OPA_TRACE0;
  m_sock->connect(server_info.c_str());
  // OPA_TRACE0;
}

void Comm::send(const JobMsg &msg) { send_str(msg.SerializeAsString()); }

std::shared_ptr<JobMsg> Comm::recv() {
  try {
    std::string res = recv_str();
    auto *x = new JobMsg;

    //google::protobuf::io::CodedInputStream input(
    //  reinterpret_cast<const u8 *>(res.c_str()), res.size());
    //input.SetTotalBytesLimit(1e9, 4e8);
    x->ParseFromString(res);
    return JobMsgPtr(x);
  } catch (...) {
    return nullptr;
  }
}

std::string Comm::recv_str() {
  zmq::message_t msg;
  OPA_ASSERTNO0(m_sock->recv(&msg));
  return std::string((const char *)msg.data(), msg.size());
}

void Comm::send_str(const std::string &str) {
  zmq::message_t msg(str.size());
  memcpy(msg.data(), str.c_str(), str.size());
  OPA_ASSERTNO0(m_sock->send(msg));
}

void Comm::close() { m_sock->close(); }
void Comm::release() { m_ctx.reset(); }

OPA_NAMESPACE_DECL3_END

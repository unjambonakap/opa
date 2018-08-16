#pragma once

#include <opa_common.h>
#include <opa/threading/data.h>

namespace zmq {
class context_t;
class socket_t;
};

OPA_NAMESPACE_DECL3(opa, threading, internal)

class Comm {
  public:
    Comm();
    ~Comm();
    void bind(const std::string &server_info);
    void connect(const std::string &server_info);
    void send(const JobMsg &msg);
    void send_str(const std::string &str);
    std::string recv_str();

    JobMsgPtr recv();
    void release();
    void close();

  private:
    std::shared_ptr<zmq::context_t> m_ctx;
    std::shared_ptr<zmq::socket_t> m_sock;
};

OPA_NAMESPACE_DECL3_END

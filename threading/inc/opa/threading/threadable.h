#pragma once

#include <opa_common.h>

OPA_NAMESPACE_DECL2(opa, threading)

class Threadable {
  public:
    Threadable();
    virtual ~Threadable();
    virtual void run() = 0;
    void start();

    void stop_and_join();
    bool want_stop() const { return m_stop; }
    virtual void stop() { m_stop = 1; }
    void detach();
    void join();

  private:
    bool m_stop;
    std::thread m_thread;
};

OPA_NAMESPACE_DECL2_END

#include "threadable.h"

OPA_NAMESPACE_DECL2(opa, threading)

Threadable::Threadable() { m_stop = 0; }

Threadable::~Threadable() { stop_and_join(); }
void Threadable::start() { m_thread = std::thread(&Threadable::run, this); }

void Threadable::detach() { m_thread.detach(); }
void Threadable::join() {
    if (m_thread.joinable())
        m_thread.join();
}
void Threadable::stop_and_join() {
    stop();
    join();
}

OPA_NAMESPACE_DECL2_END

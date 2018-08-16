#pragma once

#include <opa_common_base.h>

OPA_NAMESPACE_DECL2(opa, threading)

template <class T> class Queue {
  public:
    Queue();
    ~Queue() { release(); }
    void push(const T &a);
    bool pop(T &a);
    bool wait(T &a, int timeout_ms = -1);

    void release() {
        std::unique_lock<std::mutex> lk(m_lock);
        m_released = true;
        m_cond.notify_all();
        m_size_cond.notify_all();
    }
    void block() { m_released = false; }

    void reset() {
        m_size = 0;
        block();
        while (m_queue.size())
            m_queue.pop();
    }

    int size() const { return m_size; }
    OPA_ACCESSOR_R(bool, m_released, released);

    bool wait_size(const std::function<bool(int)> &func);

  private:
    void do_pop(T &a);
    std::function<bool()> get_pred() const;

    int m_size;
    bool m_released;
    std::mutex m_lock;
    std::queue<T> m_queue;
    std::condition_variable m_cond;
    std::condition_variable m_size_cond;
};

template <class T> Queue<T>::Queue() { reset(); }

template <class T> void Queue<T>::push(const T &a) {
    std::unique_lock<std::mutex> lk(m_lock);
    m_size += 1;
    m_queue.push(a);
    m_cond.notify_one();
}

template <class T> bool Queue<T>::pop(T &a) {
    std::unique_lock<std::mutex> lk(m_lock);
    if (!m_size)
        return false;
    do_pop(a);
}

template <class T> bool Queue<T>::wait(T &a, int timeout_ms) {
    std::unique_lock<std::mutex> lk(m_lock);

    if (timeout_ms == -1)
        m_cond.wait(lk, get_pred());
    else {
        if (!m_cond.wait_until(lk, std::chrono::system_clock::now() +
                                       std::chrono::milliseconds(timeout_ms),
                               get_pred()))
            return false;
    }

    if (!m_size)
        return false;

    do_pop(a);
    return true;
}

template <class T> void Queue<T>::do_pop(T &a) {
    a = m_queue.front();
    m_queue.pop();
    --m_size;
    m_size_cond.notify_all();
}

template <class T> bool Queue<T>::wait_size(const std::function<bool(int)> &func) {
    std::unique_lock<std::mutex> lk(m_lock);
    m_size_cond.wait(lk, [&]() { return m_released || func(size()); });
    return true;
}

template <class T> std::function<bool()> Queue<T>::get_pred() const {
    return [this]() { return this->m_queue.size() > 0 || m_released; };
}
OPA_NAMESPACE_DECL2_END

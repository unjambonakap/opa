#include <opa/utils/time.h>

OPA_NAMESPACE_DECL2(opa, utils)

Timer Timer::now() {
    Timer res;
    res.update();
    return res;
}

void Timer::update() {
    struct timeval tv;
    gettimeofday(&tv, 0);
    m_time_us = tv.tv_sec * 1000000 + tv.tv_usec;
}

OPA_NAMESPACE_DECL2_END

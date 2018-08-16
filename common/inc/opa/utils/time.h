#pragma once

#include <opa_common.h>

OPA_NAMESPACE_DECL2(opa, utils)

class Timer {
public:
  s64 us() const { return m_time_us; }
  s64 ms() const { return m_time_us / 1000; }
  s64 s() const { return m_time_us / 1000000; }
  Timer(s64 v = 0) : m_time_us(v) {}

  static Timer now();
  static Timer froms(s32 v) { return Timer(v * 1000000); }
  static Timer fromms(s32 v) { return Timer(v * 1000); }
  static Timer fromus(s32 v) { return Timer(v); }
  void update();
  Timer &operator-=(const Timer &a) {
    m_time_us -= a.m_time_us;
    return *this;
  }
  Timer &operator+=(const Timer &a) {
    m_time_us += a.m_time_us;
    return *this;
  }
  Timer operator-(const Timer &a) const {
    Timer res = *this;
    res -= a;
    return res;
  }
  Timer operator+(const Timer &a) const {
    Timer res = *this;
    res += a;
    return res;
  }
  bool operator<(const Timer &a) const { return m_time_us < a.m_time_us; };

private:
  s64 m_time_us;
};

class TimeAlloc {
public:
  TimeAlloc() {}
  TimeAlloc(const Timer &alloc) { reset(alloc); }
  void reset(const Timer &alloc) { m_end = Timer::now() + alloc; }

  bool elapsed() const { return m_end < Timer::now(); }

private:
  Timer m_end;
};

OPA_NAMESPACE_DECL2_END

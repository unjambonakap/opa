#pragma once

#include <opa_common.h>
#include <opa/utils/DataStruct.h>

OPA_NAMESPACE_DECL2(opa, utils)

class Time {
  public:
    Time() { m_us = 0; }
    Time(s64 us) { m_us = us; }
    Time(const Time &x) { m_us = x.m_us; }
    static Time now();

    Time operator+(const Time &a) const { return Time(m_us + a.m_us); }
    Time operator-(const Time &a) const { return Time(m_us - a.m_us); }
    OPA_ACCESSOR_R(s64, m_us, us);

    std::string str() const;

  private:
    s64 m_us;
};

class ProgressBar {
  public:
    typedef IdType JobId;
    class ProgressJob {
      public:
        ProgressJob() {}
        ProgressJob(ProgressBar *bar, JobId id);
        void update(int val, const std::string &desc);
        void increment(const std::string &desc);

      private:
        ProgressBar *m_bar;
        JobId m_id;
    };

    struct EventInfo {
        Time tim;
        int score;

        EventInfo(const Time &tim, int score) {
            this->tim = tim;
            this->score = score;
        }
    };

    struct JobDesc : IdObj {

        JobDesc() {}
        JobDesc(int size, const std::string &name) {
            this->size = size;
            this->cur = 0;
            this->name = name;
        }

        void set_cur(int val, const std::string &desc);
        std::string get_eta() const;

        Time evaluate_remaining_time();

        int size;
        int cur;
        std::string name;
        std::vector<EventInfo> m_times;

        std::string disp_str;
    };

    ProgressJob create(const JobDesc &desc);
    JobDesc &get(JobId id);
    void refresh();

    static ProgressBar *get();

  private:
    ProgressBar() {}

    std::mutex m_lock;
    static ProgressBar *g_instance;
    ObjectPool<JobDesc> m_pool;
};

OPA_NAMESPACE_DECL2_END

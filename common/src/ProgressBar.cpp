#include <opa/utils/ProgressBar.h>

#include <opa/utils/string.h>

OPA_NAMESPACE_DECL2(opa, utils)
#define OPA_PB_LOCK std::unique_lock<std::mutex> lk(m_lock);

const int N_US_MS = 1000;
const int N_MS_S = 1000;
const int N_S_M = 60;
const u64 MS_US = N_US_MS;
const u64 S_US = MS_US * N_MS_S;
const u64 M_US = S_US * N_S_M;

ProgressBar *ProgressBar::g_instance = nullptr;
const int BAR_LEN = 10;

Time Time::now() {
    Time res;
    struct timeval tv;
    gettimeofday(&tv, 0);
    res.m_us = tv.tv_sec * S_US + tv.tv_usec;
    return res;
}

std::string Time::str() const {
    int n_m = m_us / M_US;
    int n_s = m_us / S_US % N_S_M;
    int n_ms = m_us / MS_US % N_MS_S;
    int n_us = m_us % N_US_MS;

    std::string buf;
    buf += format("%dm,%ds,%dms", n_m, n_s, n_ms);
    return buf;
}

ProgressBar::ProgressJob::ProgressJob(ProgressBar *bar, JobId id) {
    m_bar = bar;
    m_id = id;
}

void ProgressBar::ProgressJob::update(int val, const std::string &desc) {
    m_bar->get(m_id).set_cur(val, desc);
}
void ProgressBar::ProgressJob::increment(const std::string &desc) {
    update(m_bar->get(m_id).cur + 1, desc);
}

ProgressBar::ProgressJob ProgressBar::create(const JobDesc &desc) {
    OPA_PB_LOCK;

    JobId id = m_pool.get_new();
    m_pool.get(id) = desc;
    m_pool.get(id).id = id;

    return ProgressJob(this, id);
}

ProgressBar::JobDesc &ProgressBar::get(JobId id) { return m_pool.get(id); }

void ProgressBar::JobDesc::set_cur(int val, const std::string &desc) {
    cur = val;
    Time now = Time::now();
    m_times.pb(EventInfo(now, val));
    disp_str = desc;
}

std::string ProgressBar::JobDesc::get_eta() const {
    if (cur == 0)
        return "NA";
    Time past = Time::now() - m_times[0].tim;
    float rate = 1. * cur / past.us();
    float rem = (size - cur) / rate;
    return Time(rem).str();
}

void ProgressBar::refresh() {
    OPA_PB_LOCK;
    OPA_CHECK0(m_pool.size() <= 1); // atm support only one
    if (m_pool.size() == 0)
        return;
    const JobDesc &desc = m_pool.get(*m_pool.used().begin());
    std::string buf;

    buf += "        \r[";
    float rate = 10. * desc.cur / desc.size;
    REP(i, 10) buf += i < rate ? '#' : '.';

    buf += format("], status=%s:%s, eta=%s", desc.name.c_str(), desc.disp_str.c_str(),
                  desc.get_eta().c_str());
    printf("%s", buf.c_str());
    fflush(stdout);
}
ProgressBar *ProgressBar::get() {
    if (!g_instance)
        g_instance = new ProgressBar;
    return g_instance;
}

OPA_NAMESPACE_DECL2_END

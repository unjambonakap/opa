#include "runner.h"
#include "client_dispatcher.h"
#include "dispatcher.h"
#include "job.h"
#include "worker.h"
#include <gflags/gflags.h>
#include <opa/utils/misc.h>
#include <opa/utils/string.h>

using namespace std;
using namespace opa::utils;

OPA_NAMESPACE_DECL2(opa, threading)

DEFINE_int32(cloudy_nthread, 4, "");
DEFINE_string(cloudy_server, "", "");
DEFINE_string(cloudy_action, "both", "");

std::map<std::string, int> Runner::s_job_map;
std::vector<JobCreator> Runner::s_jobs;

UPTR(Runner) g_runner;

JobId Runner::Register_job(const std::string &name,
                           const JobCreator &job_creator) {
  OPA_ASSERT(s_job_map.count(name) == 0, "Job already registered %s",
             name.c_str());
  JobId cur_id = s_jobs.size();
  OPA_DISPERR("Registering job ", name, cur_id);
  s_job_map[name] = cur_id;
  s_jobs.pb(job_creator);
  return cur_id;
}
bool Runner::CheckJobId(JobId id) { return id >= 0 && id < s_jobs.size(); }

Runner::~Runner() { stop(); }

Runner::Runner() {}

void Runner::EasySetup() {
  Runner::Build();
  g_runner.reset(new Runner);
  g_runner->run_both();
}

Dispatcher *Runner::EasyDispatcher() {
  if (!g_runner) return nullptr;
  return g_runner->dispatcher();
}

void Runner::stop() {
  printf("stopping dispatcher\n");
  m_client_dispatcher.stop();
  for (auto &x : m_workers) x->stop();
  for (auto &x : m_workers) x->join();
  m_client_dispatcher.join();
}

void Runner::Build() { OPA_DISP("Runner with ", FLAGS_cloudy_nthread); }

void Runner::run_cmd() {

  if (FLAGS_cloudy_action == "server")
    run_server_cmd();
  else if (FLAGS_cloudy_action == "client") {
    run_client_cmd();
    m_client_dispatcher.join();
  } else if (FLAGS_cloudy_action == "both")
    run_both();
  else
    OPA_ASSERT(0, "");
}

void Runner::run_client_cmd() {
  std::string server_info = get_server_info(FLAGS_cloudy_server);
  run_client(get_nthread(), server_info);
}
int Runner::get_nthread() const { return FLAGS_cloudy_nthread; }

void Runner::run_client(int nthread, const std::string &server_info) {
  OPA_DISP0("running client ", nthread);

  REP (i, nthread) {
    m_workers.pb(std::make_shared<Worker>(m_client_dispatcher));
    m_workers.back()->start();
  }
  m_client_dispatcher.initialize(server_info, nthread);
  m_client_dispatcher.start();
}
std::string Runner::get_server_info(const std::string &hostname) const {
  OPA_DISP0(hostname);
  std::string server_info =
    utils::stdsprintf("tcp://%s:5555", hostname.c_str());
  return server_info;
}

void Runner::run_server_cmd() { run_server(get_server_info("0.0.0.0")); }

void Runner::run_server(const std::string &server_info) {
  OPA_DISP0(server_info);

  m_dispatcher.initialize(this, server_info);
}

void Runner::run_both() {
  std::string server_info = get_server_info("0.0.0.0");
  run_client(get_nthread(), server_info);
  run_server(server_info);
}

Job *Runner::GetJob(JobId job_id) {
  OPA_ASSERT(job_id < s_jobs.size(), "Bad jobid, v=%d, max=%d", job_id,
             (int)s_jobs.size());
  Job *job = s_jobs[job_id]();
  job->set_job_id(job_id);
  return job;
}

OPA_NAMESPACE_DECL2_END

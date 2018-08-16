#pragma once

#include <opa_common.h>
#include <opa/utils/register.h>
#include <opa/threading/data.h>
#include <opa/threading/dispatcher.h>
#include <opa/threading/client_dispatcher.h>
#include <opa/threading/worker.h>

OPA_NAMESPACE_DECL2(opa, threading)

class Job;

class Runner {
public:
  static JobId Register_job(const std::string &name,
                            const JobCreator &job_creator);
  static void EasySetup();
  static Dispatcher *EasyDispatcher();
  static void Build();
  static JobId GetJobId(const std::string &name) { return s_job_map[name]; }
  static Job *GetJob(JobId job_id);
  static bool CheckJobId(JobId id);
  template <class T> static T *GetJob(const std::string &name) {
    return (T *)GetJob(s_job_map[name]);
  }
  template <class T> static T *GetJobById(JobId id) { return (T *)GetJob(id); }

  void stop();
  ~Runner();
  Runner();

  void run_cmd();
  void run_both();
  OPA_ACCESSOR_PTR(Dispatcher, &m_dispatcher, dispatcher);
  int get_nthread() const;

private:
  void run_client(int nthread, const std::string &server_info);
  void run_client_cmd();

  void run_server_cmd();
  void run_server(const std::string &server_info);

  std::string get_server_info(const std::string &hostname) const;

  static std::map<std::string, JobId> s_job_map;
  static std::vector<JobCreator> s_jobs;

  ClientDispatcher m_client_dispatcher;
  Dispatcher m_dispatcher;
  std::vector<std::shared_ptr<Worker> > m_workers;
};

#define OPA_CLOUDY_REGISTER_BASE_FUNC(cl) opa_cloudy_register_##cl
#define OPA_CLOUDY_REGISTER_BASE(cl)                                           \
  void OPA_CLOUDY_REGISTER_BASE_FUNC(cl)() {                                   \
    std::string name = "CloudyRegisterBase_" #cl;                              \
    opa::threading::JobId id =                                                 \
      opa::threading::Runner::Register_job(name, []() { return new cl; });     \
    cl::StaticJobId = id;                                                      \
    cl::JobName = name;                                                        \
  }                                                                            \
  OPA_REGISTER_INIT(FUNC##cl, opa_cloudy_register_##cl);

#define OPA_CLOUDY_REGISTER_BASE_TMPL(jobname, cl, ...)                        \
  void OPA_CLOUDY_REGISTER_BASE_FUNC(jobname)() {                              \
    std::string name = "CloudyRegisterBase_" #jobname;                         \
    opa::threading::JobId id = opa::threading::Runner::Register_job(           \
      name, []() { return new cl<__VA_ARGS__>; });                             \
    cl<__VA_ARGS__>::StaticJobId = id;                                         \
    cl<__VA_ARGS__>::JobName = name;                                           \
  }                                                                            \
  OPA_REGISTER_INIT(FUNC##jobname, opa_cloudy_register_##jobname);

OPA_NAMESPACE_DECL2_END

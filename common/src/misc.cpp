#include <opa/utils/misc.h>

OPA_NAMESPACE_DECL2(opa, utils)

std::string get_hostname() {
  char buf[HOST_NAME_MAX + 1];
  pii x;
  OPA_DISP0(OPA_FIELDS(x, ST, ND));
  int res = gethostname(buf, HOST_NAME_MAX);
  OPA_CHECK0(res == 0);
  return std::string(buf);
}

std::string get_process_fingerprint() {
  int pid = getpid();
  return stdsprintf("%s#%d", get_hostname().c_str(), pid);
}

Initable::~Initable() { this->fini(); }
void Initable::fini() { m_init = false; }
void Initable::init() {
  if (is_init())
    fini();
  // check_not_init();
  m_init = true;
}

OPA_NAMESPACE_DECL2_END

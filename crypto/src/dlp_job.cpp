#include "dlp_job.h"

#include <opa/threading/runner.h>

using namespace opa::math::common;
using namespace opa::threading;
using namespace std;

OPA_NAMESPACE_DECL2(opa, crypto)
const std::string DlpJob::JOB_NAME = "DLP_JOB";

bignum str2bn(const std::string &buf) { return bignum(buf, 16); }
std::string bn2str(const bignum &a) { return a.str(16); }

void DlpJob::Register() {
  threading::Runner::Register_job(JOB_NAME, []() { return new DlpJob(); });
}

void DlpJob::tworker_initialize(const threading::JobMsg &base_msg,
                                const opa::crypto::protobuf::DlpInit &init) {
  dlp->setup(str2bn(init.n()), str2bn(init.order()), str2bn(init.g()),
             str2bn(init.y()));
}

void DlpJob::worker_finalize() {}

void DlpJob::tworker_do_work(const threading::JobMsg &base_msg,
                             const opa::crypto::protobuf::DlpData &data,
                             opa::crypto::protobuf::DlpRes &out_res) {
  bignum x = dlp->do_one(str2bn(data.suborder()), data.suborder_pw());

  out_res.set_res(bn2str(x));
}

// dispatcher part
void DlpJob::tserver_initialize(opa::crypto::protobuf::DlpInit &out_init) {
  out_init.set_n(bn2str(dlp->n));
  out_init.set_g(bn2str(dlp->g));
  out_init.set_order(bn2str(dlp->order));
  out_init.set_y(bn2str(dlp->y));
}

void DlpJob::tserver_get_work(const std::function<
  DataId(const opa::crypto::protobuf::DlpData &data, bool &out_more)> &cb) {

  bool out_more; // ignored
  REP(i, dlp->factors.size()) {
    protobuf::DlpData data;
    data.set_suborder(bn2str(dlp->factors[i].ST));
    data.set_suborder_pw(dlp->factors[i].ND);
    cb(data, out_more);
  }
}

void DlpJob::server_finalize() {
  vector<bignum> tb;
  for (const auto &x : m_res_list) {
    cout << x.res() << endl;
    tb.pb(str2bn(x.res()));
  }
  dlp->compute_result(tb);
}

OPA_NAMESPACE_DECL2_END

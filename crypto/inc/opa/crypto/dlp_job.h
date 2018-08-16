#pragma once

#include <opa_common.h>
#include <opa/threading/job.h>
#include <opa/threading/runner.h>
#include <opa/crypto/crypto_msg_desc.pb.h>
#include <opa/math/common/bignum.h>
#include <opa/math/common/Utils.h>
#include <opa/crypto/dlp.h>

OPA_NAMESPACE_DECL2(opa, crypto)

class DlpJob : public opa::threading::VectorJob<opa::crypto::protobuf::DlpInit,
                                                opa::crypto::protobuf::DlpData,
                                                opa::crypto::protobuf::DlpRes> {
  public:
    static void Register();
    static const std::string JOB_NAME;
    DlpJob() { dlp = std::make_shared<Dlp>(); }

    virtual void tworker_initialize(const threading::JobMsg &base_msg,
                                    const opa::crypto::protobuf::DlpInit &init);
    virtual void worker_finalize();

    virtual void tworker_do_work(const threading::JobMsg &base_msg,
                                 const opa::crypto::protobuf::DlpData &data,
                                 opa::crypto::protobuf::DlpRes &out_res);
    virtual void tserver_initialize(opa::crypto::protobuf::DlpInit &out_init);

    // dispatcher part
    virtual void server_finalize();
    virtual void tserver_get_work(const std::function<
        opa::threading::DataId(const opa::crypto::protobuf::DlpData &data, bool &out_more)> &cb);

    std::shared_ptr<Dlp> dlp;

  private:
};

OPA_NAMESPACE_DECL2_END

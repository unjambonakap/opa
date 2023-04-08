#pragma once

#include <opa_common.h>
#include <opa/threading/msg_desc.pb.h>
#include <opa/utils/misc.h>

OPA_NAMESPACE_DECL2(opa, threading)

OPA_DECL_SPTR(opa::threading::JobMsg, JobMsgPtr);
typedef s64 JobNonce;
typedef s64 DataId;
typedef s32 JobId;
const JobId Invalid_JobId = -1;
class Job;
typedef std::function<Job *()> JobCreator;

enum class DispatcherState : int { Ok, ChangedJob, Wait, Done };

static opa::threading::JobMsg Build_StatusMsg(StatusCode status) {
    JobMsg res;
    res.set_type(MessageType::Status);
    StatusMsg status_msg;
    status_msg.set_status(status);
    res.mutable_any()->PackFrom(status_msg);
    return res;
}

OPA_NAMESPACE_DECL2_END

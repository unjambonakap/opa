#include <opa_common.h>
#include <opa/threading/runner.h>
#include <opa/threading/dispatcher.h>
#include <opa/crypto/dlp.h>
#include <opa/crypto/dlp_job.h>

using namespace opa::threading;
using namespace opa::crypto;

void init_dlp(Dlp *dlp) { dlp->setup_main(23, 3, 13); }

void run_server(Dispatcher *dispatcher) {
    DlpJob dlp_job;
    init_dlp(dlp_job.dlp.get());
    dispatcher->process_job(dlp_job, Runner::GetJobId(DlpJob::JOB_NAME));
}

void test_multi(int argc, char **argv) {
    Runner runner;

    DlpJob::Register();
    Runner::Build();
    Dispatcher *dispatcher;

    runner.run(argc, (char **)argv, run_server);
}

void test_single() {
    Dlp dlp;
    init_dlp(&dlp);
    OPA_BG res = dlp.solve();
    assert(dlp.check(res));
}

int main(int argc, char **argv) {
    opa::math::common::initMathCommon(0);
    puts("HERE");

    if (1) {
        test_multi(argc, argv);

    } else {
        test_single();
    }

    return 0;
}

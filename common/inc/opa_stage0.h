#ifndef _H_OPA_STAGE0
#define _H_OPA_STAGE0

#include "opa_inc.h"
#include "opa_stage1.h"

#define BUF_LEN 128
#define MAX_SEGS 2
struct stage0_data_t {
    char stage1_elf[BUF_LEN];
    char stage1_main_sym[32];

    u64 dlopen_addr;
    u64 dlsym_addr;
    u64 dl_allocate_tls_addr;

    stage1_data_t stage1_data;
    u64 dummy_sysnum;
};

#endif

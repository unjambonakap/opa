#ifndef _H_INJECTOR
#define _H_INJECTOR

#include "opa_inc.h"

#define STAGE1_DATA_LEN 256

struct stage1_data_t {
    u8 stage1_data[STAGE1_DATA_LEN];
    // u64 dynsym_off;
    // int num_syms;
    // u64 dynstr_off;
    // u64 gotplt_off;

    // u64 dlsym_func;
    // u64 stage1_elf_addr;  // filled by stage0
};

// typedef void (*stage1_entry_t)(stage1_data_t *data);
typedef void (*stage1_main_t)(const void *);

#endif

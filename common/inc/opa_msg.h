#ifndef _H_OPA_MSG
#define _H_OPA_MSG

#include "opa_inc.h"
#define OPA_CONTROLLER_DEV_NAME "opa_ctrl"

BEGIN_DECLS

typedef enum opa_MsgType {
    opa_MsgMapping,
    opa_MsgInjectThread,
    opa_MsgKernSym,
    opa_MsgPeek,
    opa_MsgPoke,
    opa_MsgEnd,
} opa_MsgType;

typedef struct opa_Mapping {
    s32 victim_pid;
    u64 victim_addr;
    s32 n;
    u64 out_mapping_addr;
} opa_Mapping;

typedef union {
    opa_Mapping mapping;
    opa_InjectThread thread;
    opa_KernSymData sym_data;
    opa_PeekPokeData peekpoke;
    opa_None none;
} opa_MsgUnion;

typedef struct opa_Msg {
    opa_MsgType type;
    opa_MsgUnion obj;

} opa_Msg;

END_DECLS

#endif

#ifndef _H_OPA_SYSCALL_DEFS
#define _H_OPA_SYSCALL_DEFS

#include "opa_uapi_inc.h"

#define weak_alias(name, alias)     extern __typeof(name) alias __attribute__((weak, alias(#name)))


#endif


#ifndef _H_OPA_UAPI_INC
#define _H_OPA_UAPI_INC

#include "opa_inc.h"
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/capability.h>
#include <sys/signal.h>
#include <unistd.h>
#include <sys/prctl.h>
#include <asm/prctl.h>

__BEGIN_DECLS
long __arch_prctl(int, unsigned long);
__END_DECLS
#endif

#include <opa_asm.h>
    ENTRY(capget)


    movl    $__NR_capget, %eax
    syscall


    ret
END(capget)


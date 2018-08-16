#include <opa_asm.h>
    ENTRY(capget)


    movl    $__NR_capget, %eax
    int     $0x80


    ret
END(capget)


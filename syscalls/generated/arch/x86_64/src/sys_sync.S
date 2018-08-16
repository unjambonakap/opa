#include <opa_asm.h>
    ENTRY(sync)


    movl    $__NR_sync, %eax
    syscall


    ret
END(sync)


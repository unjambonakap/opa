#include <opa_asm.h>
    ENTRY(setpriority)


    movl    $__NR_setpriority, %eax
    syscall


    ret
END(setpriority)


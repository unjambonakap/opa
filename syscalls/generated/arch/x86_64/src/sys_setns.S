#include <opa_asm.h>
    ENTRY(setns)


    movl    $__NR_setns, %eax
    syscall


    ret
END(setns)


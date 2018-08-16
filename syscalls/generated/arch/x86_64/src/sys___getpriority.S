#include <opa_asm.h>
    ENTRY(__getpriority)


    movl    $__NR_getpriority, %eax
    syscall


    ret
END(__getpriority)


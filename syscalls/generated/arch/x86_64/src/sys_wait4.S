#include <opa_asm.h>
    ENTRY(wait4)

mov %rcx, %r10

    movl    $__NR_wait4, %eax
    syscall


    ret
END(wait4)


#include <opa_asm.h>
    ENTRY(setsockopt)

mov %rcx, %r10

    movl    $__NR_setsockopt, %eax
    syscall


    ret
END(setsockopt)


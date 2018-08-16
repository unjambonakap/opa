#include <opa_asm.h>
    ENTRY(getsockopt)

mov %rcx, %r10

    movl    $__NR_getsockopt, %eax
    syscall


    ret
END(getsockopt)


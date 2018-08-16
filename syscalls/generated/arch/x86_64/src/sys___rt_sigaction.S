#include <opa_asm.h>
    ENTRY(__rt_sigaction)

mov %rcx, %r10

    movl    $__NR_rt_sigaction, %eax
    syscall


    ret
END(__rt_sigaction)


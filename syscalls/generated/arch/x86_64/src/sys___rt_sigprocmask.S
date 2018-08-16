#include <opa_asm.h>
    ENTRY(__rt_sigprocmask)

mov %rcx, %r10

    movl    $__NR_rt_sigprocmask, %eax
    syscall


    ret
END(__rt_sigprocmask)


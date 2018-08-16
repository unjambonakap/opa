#include <opa_asm.h>
    ENTRY(__rt_sigtimedwait)

mov %rcx, %r10

    movl    $__NR_rt_sigtimedwait, %eax
    syscall


    ret
END(__rt_sigtimedwait)


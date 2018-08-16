#include <opa_asm.h>
    ENTRY(prctl)

mov %rcx, %r10

    movl    $__NR_prctl, %eax
    syscall


    ret
END(prctl)


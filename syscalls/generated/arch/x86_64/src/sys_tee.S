#include <opa_asm.h>
    ENTRY(tee)

mov %rcx, %r10

    movl    $__NR_tee, %eax
    syscall


    ret
END(tee)


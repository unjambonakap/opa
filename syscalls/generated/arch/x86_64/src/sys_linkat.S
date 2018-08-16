#include <opa_asm.h>
    ENTRY(linkat)

mov %rcx, %r10

    movl    $__NR_linkat, %eax
    syscall


    ret
END(linkat)


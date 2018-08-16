#include <opa_asm.h>
    ENTRY(pwrite64)

mov %rcx, %r10

    movl    $__NR_pwrite64, %eax
    syscall


    ret
END(pwrite64)


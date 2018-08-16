#include <opa_asm.h>
    ENTRY(fallocate)

mov %rcx, %r10

    movl    $__NR_fallocate, %eax
    syscall


    ret
END(fallocate)


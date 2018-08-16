#include <opa_asm.h>
    ENTRY(utimensat)

mov %rcx, %r10

    movl    $__NR_utimensat, %eax
    syscall


    ret
END(utimensat)


#include <opa_asm.h>
    ENTRY(__openat)

mov %rcx, %r10

    movl    $__NR_openat, %eax
    syscall


    ret
END(__openat)


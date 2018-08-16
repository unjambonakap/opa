#include <opa_asm.h>
    ENTRY(__accept4)

mov %rcx, %r10

    movl    $__NR_accept4, %eax
    syscall


    ret
END(__accept4)


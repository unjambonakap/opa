#include <opa_asm.h>
    ENTRY(socketpair)

mov %rcx, %r10

    movl    $__NR_socketpair, %eax
    syscall


    ret
END(socketpair)


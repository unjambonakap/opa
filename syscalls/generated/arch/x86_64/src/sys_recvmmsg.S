#include <opa_asm.h>
    ENTRY(recvmmsg)

mov %rcx, %r10

    movl    $__NR_recvmmsg, %eax
    syscall


    ret
END(recvmmsg)


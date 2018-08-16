#include <opa_asm.h>
    ENTRY(sendmmsg)

mov %rcx, %r10

    movl    $__NR_sendmmsg, %eax
    syscall


    ret
END(sendmmsg)


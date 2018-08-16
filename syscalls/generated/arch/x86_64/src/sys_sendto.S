#include <opa_asm.h>
    ENTRY(sendto)

mov %rcx, %r10

    movl    $__NR_sendto, %eax
    syscall


    ret
END(sendto)


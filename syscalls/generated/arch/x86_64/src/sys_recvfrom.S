#include <opa_asm.h>
    ENTRY(recvfrom)

mov %rcx, %r10

    movl    $__NR_recvfrom, %eax
    syscall


    ret
END(recvfrom)


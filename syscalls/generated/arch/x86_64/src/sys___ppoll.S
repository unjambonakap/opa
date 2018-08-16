#include <opa_asm.h>
    ENTRY(__ppoll)

mov %rcx, %r10

    movl    $__NR_ppoll, %eax
    syscall


    ret
END(__ppoll)


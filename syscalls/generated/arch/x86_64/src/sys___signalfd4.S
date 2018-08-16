#include <opa_asm.h>
    ENTRY(__signalfd4)

mov %rcx, %r10

    movl    $__NR_signalfd4, %eax
    syscall


    ret
END(__signalfd4)


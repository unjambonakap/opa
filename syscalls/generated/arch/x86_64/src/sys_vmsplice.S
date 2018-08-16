#include <opa_asm.h>
    ENTRY(vmsplice)

mov %rcx, %r10

    movl    $__NR_vmsplice, %eax
    syscall


    ret
END(vmsplice)


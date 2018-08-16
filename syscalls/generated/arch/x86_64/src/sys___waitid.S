#include <opa_asm.h>
    ENTRY(__waitid)

mov %rcx, %r10

    movl    $__NR_waitid, %eax
    syscall


    ret
END(__waitid)


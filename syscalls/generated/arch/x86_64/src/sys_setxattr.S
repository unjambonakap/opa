#include <opa_asm.h>
    ENTRY(setxattr)

mov %rcx, %r10

    movl    $__NR_setxattr, %eax
    syscall


    ret
END(setxattr)


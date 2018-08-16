#include <opa_asm.h>
    ENTRY(getxattr)

mov %rcx, %r10

    movl    $__NR_getxattr, %eax
    syscall


    ret
END(getxattr)


#include <opa_asm.h>
    ENTRY(lgetxattr)

mov %rcx, %r10

    movl    $__NR_lgetxattr, %eax
    syscall


    ret
END(lgetxattr)


#include <opa_asm.h>
    ENTRY(lsetxattr)

mov %rcx, %r10

    movl    $__NR_lsetxattr, %eax
    syscall


    ret
END(lsetxattr)


#include <opa_asm.h>
    ENTRY(fsetxattr)

mov %rcx, %r10

    movl    $__NR_fsetxattr, %eax
    syscall


    ret
END(fsetxattr)


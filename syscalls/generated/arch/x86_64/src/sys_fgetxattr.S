#include <opa_asm.h>
    ENTRY(fgetxattr)

mov %rcx, %r10

    movl    $__NR_fgetxattr, %eax
    syscall


    ret
END(fgetxattr)


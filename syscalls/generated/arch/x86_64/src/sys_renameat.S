#include <opa_asm.h>
    ENTRY(renameat)

mov %rcx, %r10

    movl    $__NR_renameat, %eax
    syscall


    ret
END(renameat)


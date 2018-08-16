#include <opa_asm.h>
    ENTRY(mremap)

mov %rcx, %r10

    movl    $__NR_mremap, %eax
    syscall


    ret
END(mremap)


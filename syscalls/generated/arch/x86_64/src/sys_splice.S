#include <opa_asm.h>
    ENTRY(splice)

mov %rcx, %r10

    movl    $__NR_splice, %eax
    syscall


    ret
END(splice)


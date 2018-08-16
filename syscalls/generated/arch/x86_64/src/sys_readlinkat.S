#include <opa_asm.h>
    ENTRY(readlinkat)

mov %rcx, %r10

    movl    $__NR_readlinkat, %eax
    syscall


    ret
END(readlinkat)


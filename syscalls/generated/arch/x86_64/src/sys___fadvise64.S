#include <opa_asm.h>
    ENTRY(__fadvise64)

mov %rcx, %r10

    movl    $__NR_fadvise64, %eax
    syscall


    ret
END(__fadvise64)


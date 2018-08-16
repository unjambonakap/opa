#include <opa_asm.h>
    ENTRY(mknodat)

mov %rcx, %r10

    movl    $__NR_mknodat, %eax
    syscall


    ret
END(mknodat)


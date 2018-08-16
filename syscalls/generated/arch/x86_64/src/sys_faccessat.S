#include <opa_asm.h>
    ENTRY(faccessat)

mov %rcx, %r10

    movl    $__NR_faccessat, %eax
    syscall


    ret
END(faccessat)


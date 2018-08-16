#include <opa_asm.h>
    ENTRY(prlimit64)

mov %rcx, %r10

    movl    $__NR_prlimit64, %eax
    syscall


    ret
END(prlimit64)


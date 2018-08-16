#include <opa_asm.h>
    ENTRY(pread64)

mov %rcx, %r10

    movl    $__NR_pread64, %eax
    syscall


    ret
END(pread64)


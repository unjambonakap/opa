#include <opa_asm.h>
    ENTRY(fstatat64)

mov %rcx, %r10

    movl    $__NR_newfstatat, %eax
    syscall


    ret
END(fstatat64)


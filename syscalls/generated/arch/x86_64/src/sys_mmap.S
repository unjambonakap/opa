#include <opa_asm.h>
    ENTRY(mmap)

mov %rcx, %r10

    movl    $__NR_mmap, %eax
    syscall


    ret
END(mmap)


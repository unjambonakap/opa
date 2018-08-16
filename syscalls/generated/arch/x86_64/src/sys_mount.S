#include <opa_asm.h>
    ENTRY(mount)

mov %rcx, %r10

    movl    $__NR_mount, %eax
    syscall


    ret
END(mount)


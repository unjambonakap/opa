#include <opa_asm.h>
    ENTRY(fchownat)

mov %rcx, %r10

    movl    $__NR_fchownat, %eax
    syscall


    ret
END(fchownat)


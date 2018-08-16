#include <opa_asm.h>
    ENTRY(__pselect6)

mov %rcx, %r10

    movl    $__NR_pselect6, %eax
    syscall


    ret
END(__pselect6)


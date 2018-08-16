#include <opa_asm.h>
    ENTRY(__clock_nanosleep)

mov %rcx, %r10

    movl    $__NR_clock_nanosleep, %eax
    syscall


    ret
END(__clock_nanosleep)


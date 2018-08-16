#include <opa_asm.h>
    ENTRY(__timer_settime)

mov %rcx, %r10

    movl    $__NR_timer_settime, %eax
    syscall


    ret
END(__timer_settime)


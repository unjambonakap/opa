#include <opa_asm.h>
    ENTRY(timerfd_settime)

mov %rcx, %r10

    movl    $__NR_timerfd_settime, %eax
    syscall


    ret
END(timerfd_settime)


#include <opa_asm.h>
    ENTRY(timerfd_create)


    movl    $__NR_timerfd_create, %eax
    syscall


    ret
END(timerfd_create)


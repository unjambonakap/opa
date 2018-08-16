#include <opa_asm.h>
    ENTRY(timerfd_create)


    movl    $__NR_timerfd_create, %eax
    int     $0x80


    ret
END(timerfd_create)


#include <opa_asm.h>
    ENTRY(timerfd_settime)


    movl    $__NR_timerfd_settime, %eax
    int     $0x80


    ret
END(timerfd_settime)


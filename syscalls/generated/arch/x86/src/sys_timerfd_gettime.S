#include <opa_asm.h>
    ENTRY(timerfd_gettime)


    movl    $__NR_timerfd_gettime, %eax
    int     $0x80


    ret
END(timerfd_gettime)


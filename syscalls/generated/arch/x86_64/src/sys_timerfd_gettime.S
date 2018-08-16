#include <opa_asm.h>
    ENTRY(timerfd_gettime)


    movl    $__NR_timerfd_gettime, %eax
    syscall


    ret
END(timerfd_gettime)


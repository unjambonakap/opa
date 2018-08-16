#include <opa_asm.h>
    ENTRY(sched_getscheduler)


    movl    $__NR_sched_getscheduler, %eax
    syscall


    ret
END(sched_getscheduler)


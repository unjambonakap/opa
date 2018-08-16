#include <opa_asm.h>
    ENTRY(sched_getscheduler)


    movl    $__NR_sched_getscheduler, %eax
    int     $0x80


    ret
END(sched_getscheduler)


#include <opa_asm.h>
    ENTRY(sched_setscheduler)


    movl    $__NR_sched_setscheduler, %eax
    int     $0x80


    ret
END(sched_setscheduler)


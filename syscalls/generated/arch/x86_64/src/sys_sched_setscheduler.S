#include <opa_asm.h>
    ENTRY(sched_setscheduler)


    movl    $__NR_sched_setscheduler, %eax
    syscall


    ret
END(sched_setscheduler)


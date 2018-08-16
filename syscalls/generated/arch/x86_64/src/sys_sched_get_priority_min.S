#include <opa_asm.h>
    ENTRY(sched_get_priority_min)


    movl    $__NR_sched_get_priority_min, %eax
    syscall


    ret
END(sched_get_priority_min)


#include <opa_asm.h>
    ENTRY(sched_rr_get_interval)


    movl    $__NR_sched_rr_get_interval, %eax
    syscall


    ret
END(sched_rr_get_interval)


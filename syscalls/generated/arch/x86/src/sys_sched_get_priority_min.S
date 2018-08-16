#include <opa_asm.h>
    ENTRY(sched_get_priority_min)


    movl    $__NR_sched_get_priority_min, %eax
    int     $0x80


    ret
END(sched_get_priority_min)


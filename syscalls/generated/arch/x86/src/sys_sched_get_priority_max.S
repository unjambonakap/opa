#include <opa_asm.h>
    ENTRY(sched_get_priority_max)


    movl    $__NR_sched_get_priority_max, %eax
    int     $0x80


    ret
END(sched_get_priority_max)


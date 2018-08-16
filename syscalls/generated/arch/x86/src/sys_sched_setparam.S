#include <opa_asm.h>
    ENTRY(sched_setparam)


    movl    $__NR_sched_setparam, %eax
    int     $0x80


    ret
END(sched_setparam)


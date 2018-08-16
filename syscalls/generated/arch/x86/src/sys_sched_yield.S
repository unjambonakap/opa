#include <opa_asm.h>
    ENTRY(sched_yield)


    movl    $__NR_sched_yield, %eax
    int     $0x80


    ret
END(sched_yield)


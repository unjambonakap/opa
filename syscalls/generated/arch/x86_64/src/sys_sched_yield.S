#include <opa_asm.h>
    ENTRY(sched_yield)


    movl    $__NR_sched_yield, %eax
    syscall


    ret
END(sched_yield)


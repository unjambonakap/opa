#include <opa_asm.h>
    ENTRY(sched_setparam)


    movl    $__NR_sched_setparam, %eax
    syscall


    ret
END(sched_setparam)


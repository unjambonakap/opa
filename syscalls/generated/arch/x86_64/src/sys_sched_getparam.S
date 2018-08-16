#include <opa_asm.h>
    ENTRY(sched_getparam)


    movl    $__NR_sched_getparam, %eax
    syscall


    ret
END(sched_getparam)


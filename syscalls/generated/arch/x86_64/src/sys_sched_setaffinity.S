#include <opa_asm.h>
    ENTRY(sched_setaffinity)


    movl    $__NR_sched_setaffinity, %eax
    syscall


    ret
END(sched_setaffinity)


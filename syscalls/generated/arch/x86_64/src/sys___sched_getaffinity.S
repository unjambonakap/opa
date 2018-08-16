#include <opa_asm.h>
    ENTRY(__sched_getaffinity)


    movl    $__NR_sched_getaffinity, %eax
    syscall


    ret
END(__sched_getaffinity)


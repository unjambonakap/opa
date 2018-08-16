#include <opa_asm.h>
    ENTRY(sched_getparam)


    movl    $__NR_sched_getparam, %eax
    int     $0x80


    ret
END(sched_getparam)


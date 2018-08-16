#include <opa_asm.h>
    ENTRY(sched_setaffinity)


    movl    $__NR_sched_setaffinity, %eax
    int     $0x80


    ret
END(sched_setaffinity)


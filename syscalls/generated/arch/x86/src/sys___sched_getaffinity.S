#include <opa_asm.h>
    ENTRY(__sched_getaffinity)


    movl    $__NR_sched_getaffinity, %eax
    int     $0x80


    ret
END(__sched_getaffinity)


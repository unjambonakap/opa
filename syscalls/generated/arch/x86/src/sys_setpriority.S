#include <opa_asm.h>
    ENTRY(setpriority)


    movl    $__NR_setpriority, %eax
    int     $0x80


    ret
END(setpriority)


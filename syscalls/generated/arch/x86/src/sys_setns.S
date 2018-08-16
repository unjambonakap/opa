#include <opa_asm.h>
    ENTRY(setns)


    movl    $__NR_setns, %eax
    int     $0x80


    ret
END(setns)


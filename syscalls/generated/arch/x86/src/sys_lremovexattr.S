#include <opa_asm.h>
    ENTRY(lremovexattr)


    movl    $__NR_lremovexattr, %eax
    int     $0x80


    ret
END(lremovexattr)


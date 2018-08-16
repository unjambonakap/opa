#include <opa_asm.h>
    ENTRY(fremovexattr)


    movl    $__NR_fremovexattr, %eax
    int     $0x80


    ret
END(fremovexattr)


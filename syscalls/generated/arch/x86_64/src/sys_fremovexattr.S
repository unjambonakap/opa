#include <opa_asm.h>
    ENTRY(fremovexattr)


    movl    $__NR_fremovexattr, %eax
    syscall


    ret
END(fremovexattr)


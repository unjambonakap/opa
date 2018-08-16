#include <opa_asm.h>
    ENTRY(lremovexattr)


    movl    $__NR_lremovexattr, %eax
    syscall


    ret
END(lremovexattr)


#include <opa_asm.h>
    ENTRY(removexattr)


    movl    $__NR_removexattr, %eax
    syscall


    ret
END(removexattr)


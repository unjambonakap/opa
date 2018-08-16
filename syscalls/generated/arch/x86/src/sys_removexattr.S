#include <opa_asm.h>
    ENTRY(removexattr)


    movl    $__NR_removexattr, %eax
    int     $0x80


    ret
END(removexattr)


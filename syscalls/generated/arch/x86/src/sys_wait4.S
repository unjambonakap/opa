#include <opa_asm.h>
    ENTRY(wait4)


    movl    $__NR_wait4, %eax
    int     $0x80


    ret
END(wait4)


#include <opa_asm.h>
    ENTRY(swapoff)


    movl    $__NR_swapoff, %eax
    int     $0x80


    ret
END(swapoff)


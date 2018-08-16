#include <opa_asm.h>
    ENTRY(swapoff)


    movl    $__NR_swapoff, %eax
    syscall


    ret
END(swapoff)


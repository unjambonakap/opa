#include <opa_asm.h>
    ENTRY(getpeername)


    movl    $__NR_getpeername, %eax
    syscall


    ret
END(getpeername)


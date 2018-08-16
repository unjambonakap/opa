#include <opa_asm.h>
    ENTRY(writev)


    movl    $__NR_writev, %eax
    syscall


    ret
END(writev)


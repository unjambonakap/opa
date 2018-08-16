#include <opa_asm.h>
    ENTRY(getrusage)


    movl    $__NR_getrusage, %eax
    syscall


    ret
END(getrusage)


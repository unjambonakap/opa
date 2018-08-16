#include <opa_asm.h>
    ENTRY(nanosleep)


    movl    $__NR_nanosleep, %eax
    syscall


    ret
END(nanosleep)


#include <opa_asm.h>
    ENTRY(getitimer)


    movl    $__NR_getitimer, %eax
    syscall


    ret
END(getitimer)


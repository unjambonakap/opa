#include <opa_asm.h>
    ENTRY(setitimer)


    movl    $__NR_setitimer, %eax
    syscall


    ret
END(setitimer)


#include <opa_asm.h>
    ENTRY(kill)


    movl    $__NR_kill, %eax
    syscall


    ret
END(kill)


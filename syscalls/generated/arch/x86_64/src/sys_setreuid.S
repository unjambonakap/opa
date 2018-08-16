#include <opa_asm.h>
    ENTRY(setreuid)


    movl    $__NR_setreuid, %eax
    syscall


    ret
END(setreuid)


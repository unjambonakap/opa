#include <opa_asm.h>
    ENTRY(setuid)


    movl    $__NR_setuid, %eax
    syscall


    ret
END(setuid)


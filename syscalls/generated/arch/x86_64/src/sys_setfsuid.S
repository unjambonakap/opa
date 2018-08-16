#include <opa_asm.h>
    ENTRY(setfsuid)


    movl    $__NR_setfsuid, %eax
    syscall


    ret
END(setfsuid)


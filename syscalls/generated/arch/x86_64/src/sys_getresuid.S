#include <opa_asm.h>
    ENTRY(getresuid)


    movl    $__NR_getresuid, %eax
    syscall


    ret
END(getresuid)


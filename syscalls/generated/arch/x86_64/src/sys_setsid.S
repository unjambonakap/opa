#include <opa_asm.h>
    ENTRY(setsid)


    movl    $__NR_setsid, %eax
    syscall


    ret
END(setsid)


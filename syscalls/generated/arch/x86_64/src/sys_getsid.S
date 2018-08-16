#include <opa_asm.h>
    ENTRY(getsid)


    movl    $__NR_getsid, %eax
    syscall


    ret
END(getsid)


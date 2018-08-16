#include <opa_asm.h>
    ENTRY(getppid)


    movl    $__NR_getppid, %eax
    syscall


    ret
END(getppid)


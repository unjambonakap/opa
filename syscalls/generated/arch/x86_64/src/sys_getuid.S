#include <opa_asm.h>
    ENTRY(getuid)


    movl    $__NR_getuid, %eax
    syscall


    ret
END(getuid)


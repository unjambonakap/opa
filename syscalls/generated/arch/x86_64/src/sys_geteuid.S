#include <opa_asm.h>
    ENTRY(geteuid)


    movl    $__NR_geteuid, %eax
    syscall


    ret
END(geteuid)


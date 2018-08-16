#include <opa_asm.h>
    ENTRY(setresuid)


    movl    $__NR_setresuid, %eax
    syscall


    ret
END(setresuid)


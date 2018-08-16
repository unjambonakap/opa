#include <opa_asm.h>
    ENTRY(setpgid)


    movl    $__NR_setpgid, %eax
    syscall


    ret
END(setpgid)


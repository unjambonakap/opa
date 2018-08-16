#include <opa_asm.h>
    ENTRY(setgid)


    movl    $__NR_setgid, %eax
    syscall


    ret
END(setgid)


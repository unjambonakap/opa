#include <opa_asm.h>
    ENTRY(___fchmodat)


    movl    $__NR_fchmodat, %eax
    syscall


    ret
END(___fchmodat)


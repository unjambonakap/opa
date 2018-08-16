#include <opa_asm.h>
    ENTRY(___fchmod)


    movl    $__NR_fchmod, %eax
    syscall


    ret
END(___fchmod)


#include <opa_asm.h>
    ENTRY(mprotect)


    movl    $__NR_mprotect, %eax
    syscall


    ret
END(mprotect)


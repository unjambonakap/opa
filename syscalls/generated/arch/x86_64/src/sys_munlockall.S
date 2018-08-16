#include <opa_asm.h>
    ENTRY(munlockall)


    movl    $__NR_munlockall, %eax
    syscall


    ret
END(munlockall)


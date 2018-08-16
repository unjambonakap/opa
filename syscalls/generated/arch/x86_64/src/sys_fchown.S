#include <opa_asm.h>
    ENTRY(fchown)


    movl    $__NR_fchown, %eax
    syscall


    ret
END(fchown)


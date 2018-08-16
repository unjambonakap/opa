#include <opa_asm.h>
    ENTRY(bind)


    movl    $__NR_bind, %eax
    syscall


    ret
END(bind)


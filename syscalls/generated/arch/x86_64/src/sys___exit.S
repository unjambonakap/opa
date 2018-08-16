#include <opa_asm.h>
    ENTRY(__exit)


    movl    $__NR_exit, %eax
    syscall


    ret
END(__exit)


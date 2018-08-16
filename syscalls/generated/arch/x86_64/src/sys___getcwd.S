#include <opa_asm.h>
    ENTRY(__getcwd)


    movl    $__NR_getcwd, %eax
    syscall


    ret
END(__getcwd)


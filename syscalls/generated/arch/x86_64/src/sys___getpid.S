#include <opa_asm.h>
    ENTRY(__getpid)


    movl    $__NR_getpid, %eax
    syscall


    ret
END(__getpid)


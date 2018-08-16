#include <opa_asm.h>
    ENTRY(__getdents64)


    movl    $__NR_getdents64, %eax
    syscall


    ret
END(__getdents64)


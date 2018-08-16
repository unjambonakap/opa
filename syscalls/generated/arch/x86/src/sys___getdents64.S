#include <opa_asm.h>
    ENTRY(__getdents64)


    movl    $__NR_getdents64, %eax
    int     $0x80


    ret
END(__getdents64)


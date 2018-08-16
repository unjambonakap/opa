#include <opa_asm.h>
    ENTRY(mkdirat)


    movl    $__NR_mkdirat, %eax
    int     $0x80


    ret
END(mkdirat)


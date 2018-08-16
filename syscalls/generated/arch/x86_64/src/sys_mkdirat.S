#include <opa_asm.h>
    ENTRY(mkdirat)


    movl    $__NR_mkdirat, %eax
    syscall


    ret
END(mkdirat)


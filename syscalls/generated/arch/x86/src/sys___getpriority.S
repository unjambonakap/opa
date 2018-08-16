#include <opa_asm.h>
    ENTRY(__getpriority)


    movl    $__NR_getpriority, %eax
    int     $0x80


    ret
END(__getpriority)


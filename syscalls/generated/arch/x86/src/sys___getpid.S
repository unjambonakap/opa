#include <opa_asm.h>
    ENTRY(__getpid)


    movl    $__NR_getpid, %eax
    int     $0x80


    ret
END(__getpid)


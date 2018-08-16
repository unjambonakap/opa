#include <opa_asm.h>
    ENTRY(__exit)


    movl    $__NR_exit, %eax
    int     $0x80


    ret
END(__exit)


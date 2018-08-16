#include <opa_asm.h>
    ENTRY(__brk)


    movl    $__NR_brk, %eax
    int     $0x80


    ret
END(__brk)


#include <opa_asm.h>
    ENTRY(dup3)


    movl    $__NR_dup3, %eax
    int     $0x80


    ret
END(dup3)


#include <opa_asm.h>
    ENTRY(dup)


    movl    $__NR_dup, %eax
    int     $0x80


    ret
END(dup)


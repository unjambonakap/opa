#include <opa_asm.h>
    ENTRY(mincore)


    movl    $__NR_mincore, %eax
    int     $0x80


    ret
END(mincore)


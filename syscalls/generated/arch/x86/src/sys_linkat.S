#include <opa_asm.h>
    ENTRY(linkat)


    movl    $__NR_linkat, %eax
    int     $0x80


    ret
END(linkat)


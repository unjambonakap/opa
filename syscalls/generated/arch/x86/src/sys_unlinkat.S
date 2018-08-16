#include <opa_asm.h>
    ENTRY(unlinkat)


    movl    $__NR_unlinkat, %eax
    int     $0x80


    ret
END(unlinkat)


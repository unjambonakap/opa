#include <opa_asm.h>
    ENTRY(capset)


    movl    $__NR_capset, %eax
    int     $0x80


    ret
END(capset)


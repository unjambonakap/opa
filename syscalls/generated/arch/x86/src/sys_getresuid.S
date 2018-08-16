#include <opa_asm.h>
    ENTRY(getresuid)


    movl    $__NR_getresuid32, %eax
    int     $0x80


    ret
END(getresuid)


#include <opa_asm.h>
    ENTRY(setreuid)


    movl    $__NR_setreuid32, %eax
    int     $0x80


    ret
END(setreuid)


#include <opa_asm.h>
    ENTRY(setuid)


    movl    $__NR_setuid32, %eax
    int     $0x80


    ret
END(setuid)

